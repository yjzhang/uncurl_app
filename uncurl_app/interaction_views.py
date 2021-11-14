# notes: this backend should be frontend-agnostic.
# I'm thinking of writing the frontend entirely in plotly.js, and not have
# any backend Python rendering components.

import contextlib
import json
import os
import re
import shutil
import traceback
import uuid

import numpy as np
import scipy.io
from flask import request, render_template, redirect, url_for, Blueprint, current_app
from uncurl_analysis import enrichr_api, sc_analysis, custom_cell_selection

from . import generate_analysis
from .cache import cache
from .utils import SimpleEncoder
from .views import state_estimation_preproc_simple

interaction_views = Blueprint('interaction_views', __name__,
        template_folder='templates')

# map of user_id to SCAnalysis objects
# map of gene lists (strings of newline-separated gene names) to enrichr IDs
interaction_views.enrichr_gene_list_ids = {}
# map of tuples (top_genes, gene_set) to enrichr results
interaction_views.enrichr_results = {}

def pmid_to_link(pmid):
    return '<a href="https://www.ncbi.nlm.nih.gov/pubmed/{0}">{0}</a>'.format(pmid)

def split_gene_names(names):
    """
    Splits on comma, newline, or space
    """
    names = [x for x in re.split(r'[\s,]', names.strip()) if len(x)>0]
    return names

@contextlib.contextmanager
def lockfile_context(lockfile_name):
    """
    This is a wrapper around a lockfile.
    Example:
        with lockfile_context('f') as _:
            run_function()
    This will prevent other instances of run_function from running in other threads/processes.
    Raises an exception if the lockfile exists.
    """
    if os.path.exists(lockfile_name):
        raise Exception('Lockfile {0} exists'.format(lockfile_name))
    with open(lockfile_name, 'w') as f:
        f.write(' ')
    yield 1
    os.remove(lockfile_name)

def get_sca(user_id):
    path = user_id_to_path(user_id)
    sca = sc_analysis.SCAnalysis(path)
    sca = sca.load_params_from_folder()
    return sca

@cache.memoize()
def get_sca_dim_red(user_id):
    sca = get_sca(user_id)
    return sca.dim_red

@cache.memoize()
def get_sca_baseline_vis(user_id):
    sca = get_sca(user_id)
    return sca.baseline_vis

@cache.memoize()
def get_sca_top_genes(user_id):
    sca = get_sca(user_id)
    return sca.top_genes

@cache.memoize()
def get_sca_pvals(user_id):
    sca = get_sca(user_id)
    return sca.pvals

@cache.memoize()
def get_sca_pairwise_ratios(user_id):
    sca = get_sca(user_id)
    return sca.t_scores

@cache.memoize()
def get_sca_pairwise_pvals(user_id):
    sca = get_sca(user_id)
    return sca.t_pvals

@cache.memoize()
def get_sca_pval_1vr(user_id):
    sca = get_sca(user_id)
    return sca.pvals_1_vs_rest

@cache.memoize()
def get_sca_top_1vr(user_id):
    sca = get_sca(user_id)
    return sca.top_genes_1_vs_rest

@cache.memoize()
def get_sca_top_genes_custom(user_id, color_track, mode='1_vs_rest'):
    """Output is array of shape [k, genes] for 1_vs_rest or [k, k, genes] for pairwise."""
    sca = get_sca(user_id)
    return sca.calculate_diffexp(color_track, mode=mode)

@cache.memoize()
def get_sca_gene_names(user_id):
    sca = get_sca(user_id)
    return sca.gene_names

@cache.memoize()
def get_sca_color_track(user_id, color_track, return_color=False):
    sca = get_sca(user_id)
    if color_track == 'cluster':
        return sca.labels, True
    if return_color:
        return sca.get_color_track(color_track, return_colors=True)
    else:
        return sca.get_color_track(color_track)

@cache.memoize()
def get_sca_data_sampled_all_genes(user_id):
    sca = get_sca(user_id)
    return sca.data_sampled_all_genes

def color_track_map(color_track):
    """
    Returns a map of labels to ints, and a map of ints to labels.
    """
    colors = sorted(list(set(color_track)))
    return {c: i for i, c in enumerate(colors)}, {i: c for i, c in enumerate(colors)}

def array_to_top_genes(data_array, cluster1, cluster2, is_pvals=False, num_genes=10):
    """
    Given a data_array of shape (k, k, genes), this returns two arrays:
        genes and values.
    """
    data_cluster = data_array[cluster1, cluster2, :]
    if is_pvals:
        order = data_cluster.argsort()
    else:
        order = data_cluster.argsort()[::-1]
    genes = order[:num_genes]
    values = data_cluster[order[:num_genes]]
    return genes, values

def user_id_to_path(user_id, use_secondary=True):
    """
    Given a user id, returns the path to the analysis object's base directory.
    """
    if user_id.startswith('test_'):
        user_id = user_id[5:]
        path = os.path.join(current_app.config['TEST_DATA_DIR'], user_id)
        return path
    else:
        path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
        if not os.path.exists(path) and use_secondary:
            path = os.path.join(current_app.config['SECONDARY_USER_DATA_DIR'], user_id)
        return path

def barplot_data(gene_values, gene_names, cluster_name, x_label,
        title=None):
    """
    Converts data for top genes into a json for building the
    bar plot. Output should be formatted in a way that can be plugged into
    Plotly.

    Args:
        gene_values (list): list of tuples (gene_id, gene_value)
        gene_names (list): list of gene names corresponding to
            the genes in gene_values.
        cluster_name: name of the cluster from which the top genes are drawn.
        x_label: label for the x-axis.
        title: plot title
    """
    if gene_values is None:
        gene_values = [(1,1), (2,2), (3,3)]
    if gene_names is None:
        gene_names = ['placeholder 1', 'placeholder 2', 'placeholder 3']
    if title is None:
        title = 'Top genes for cluster {0}'.format(cluster_name)
    return json.dumps({
            'data': [{
                'x': list(x[1] for x in gene_values),
                'y': gene_names,
                'orientation': 'h',
                'type': 'bar',
            }],
            'layout': {
                'title': title,
                'xaxis': {'title': x_label},
                'margin': {'t': 40},
            },
        }, cls=SimpleEncoder)

def histogram_data(gene_values_cluster, gene_values_all, cluster_name, gene_name, title=None):
    """
    Creates a plotly histogram

    Args:
        gene_values_cluster (array): 1d array for a given gene, within cluster
        gene_values_all (array): 1d array for a given gene, all clusters
    """
    if title is None:
        title = 'Histogram for gene {0}'.format(gene_name)
    return json.dumps({
        'data': [{
            'x': gene_values_cluster.tolist(),
            'type': 'histogram',
            'histnorm': 'probability',
            'opacity': 0.5,
            'name': str(cluster_name),
            'marker': {'color': 'green'},
        },
        {
            'x': gene_values_all.tolist(),
            'type': 'histogram',
            'histnorm': 'probability',
            'opacity': 0.5,
            'name': 'all cells',
            'marker': {'color': 'blue'},
        }],
        'layout': {
            'title': title,
            'barmode': 'overlay',
            'showlegend': True,
            'xaxis': {'title': 'Count'},
            'yaxis': {'title': 'Fraction'},
        },
    }, cls=SimpleEncoder)

def calc_size(labels):
    size = 10
    if len(labels) < 50:
        size = 20
    if len(labels) > 1000:
        size = 5
    if len(labels) > 10000:
        size = 2
    return size

def scatterplot_data(dim_red, labels, colorscale='Portland', mode='cluster',
        color_vals=None, label_text=None, color_dict=None):
    """
    Converts data into a form that will be sent as json for building the
    scatterplot. Output should be formatted in a way that can be used by
    Plotly.

    Args:
        dim_red (array): array of shape (2, n)
        labels (array): 1d array of length n
        colorscale (str)
        mode (str): either 'cluster' or 'entropy'
        color_vals (array): 1d array of length n, of real values that will be used for coloring
        label_text (list or array): labels for each point, with length n
        color_dict (None or dict): mapping of labels to rgb values
    """
    # have size depend on data shape
    size = calc_size(labels)
    data = []
    # label_values is a list 0...number of unique labels - 1
    label_values = list(range(len(set(labels))))
    cluster_names = ['cluster ' + str(c) for c in label_values]
    if isinstance(labels[0], str):
        color_to_index, index_to_color = color_track_map(labels)
        label_values = [index_to_color[c] for c in label_values]
        print('label_values: ', label_values)
        cluster_names = label_values
    # cell_ids indicates the ids of the cells used...
    cell_ids = np.arange(len(labels))
    if label_text is None:
        label_text = np.array([str(x) for x in cell_ids])
    else:
        label_text = np.array(label_text)
    plot_type = 'scattergl' if len(label_text) > 5000 else 'scatter'
    # select color scheme
    if mode == 'cluster':
        if len(label_values) > 10 or color_dict is not None:
            from . import colors
            if len(label_values) <= 25:
                scale0 = colors.CL_25
                color_values = scale0
            else:
                scale0 = colors.CL_25 + colors.CL_25_2 + colors.CL_25_3 + colors.CL_25_4
                if len(label_values) <= 100:
                    color_values = scale0
                else:
                    import colorlover as cl
                    color_values = cl.to_rgb(cl.interp(colors.CL_25, len(label_values)))
        else:
            color_values = label_values
        if color_dict is not None:
            color_values = color_values.copy()
            print('getting color_dict values:', color_dict)
            for i, v in enumerate(label_values):
                if v in color_dict and color_dict[v] is not None:
                    color_values[i] = color_dict[v]
        print('scatterplot color_values:', color_values)
        data =  [
            {
                'x': dim_red[0,labels==c].tolist(),
                'y': dim_red[1,labels==c].tolist(),
                'mode': 'markers',
                'type': plot_type,
                'name': cluster_names[i],
                'marker': {
                    'size': size,
                    'color': color_values[i],
                    'colorscale': colorscale,
                },
                'text': list(label_text[labels==c]),
            }
            for i, c in enumerate(label_values)
        ]
    elif mode == 'entropy':
        if colorscale == 'Portland' or colorscale is None:
            colorscale = 'Reds'
        color_values = [color_vals[labels==c] for c in label_values]
        cmin = min(color_vals)
        cmax = max(color_vals)
        data = [
            {
                'x': dim_red[0,labels==c].tolist(),
                'y': dim_red[1,labels==c].tolist(),
                'mode': 'markers',
                'type': plot_type,
                'name': 'cluster ' + str(c),
                'marker': {
                    'size': size,
                    'color': color_values[c],
                    'colorscale': colorscale,
                    'cmin': cmin,
                    'cmax': cmax,
                    'showscale': True if c==0 else False,
                },
                'text': list(map(str, color_values[c])),
            }
            for c in label_values
        ]
    return json.dumps({
            'data': data,
            'layout': {
                'title': 'Cells',
                'xaxis': {'title': 'dim1'},
                'yaxis': {'title': 'dim2'},
                'margin': {'t':30},
                'showlegend': True if mode =='cluster' else False,
                'hovermode': 'closest',
                'legend': {'x': 1, 'y': 1},
            },
    }, cls=SimpleEncoder)


def volcano_plot_data(user_id, colormap, cluster1, cluster2, selected_genes=None):
    """
    Returns plotly json representation of a volcano plot
    """
    sca = get_sca(user_id)
    gene_names = get_sca_gene_names(user_id)
    # get custom colormap
    if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights', 'read_counts']:
        color_track, is_discrete = get_sca_color_track(user_id, colormap)
        lockfile_name = os.path.join(sca.data_dir, colormap + '_writing_diffexp')
        with lockfile_context(lockfile_name) as _lock:
            selected_diffexp, selected_pvals = get_sca_top_genes_custom(user_id, colormap, 'pairwise')
        color_to_index, index_to_color = color_track_map(color_track)
    # default colormap
    else:
        selected_diffexp = get_sca_pairwise_ratios(user_id)
        selected_pvals = get_sca_pairwise_pvals(user_id)
        index_to_color = list(range(selected_diffexp.shape[0]))
    if selected_genes is not None:
        gene_indices = {g:i for i, g in enumerate(gene_names)}
        selected_gene_indices = np.array([gene_indices[g] for g in selected_genes])
        selected_diffexp = selected_diffexp[:,:,selected_gene_indices]
        selected_pvals = selected_pvals[:,:,selected_gene_indices]
        gene_names = selected_genes
    # get pval data, 
    diffexp_data = selected_diffexp[cluster1, cluster2, :]
    pval_data = selected_pvals[cluster1, cluster2, :]
    pval_2v1_data = selected_pvals[cluster2, cluster1, :]
    pval_combined = np.fmin(pval_data, pval_2v1_data)
    if sca.params['use_fdr']:
        y_desc = '-log10 FDR'
    else:
        y_desc = '-log10 p-value'
    data = [{
                'x': np.log2(diffexp_data + 2e-16),
                'y': -np.log10(pval_combined + 2e-16),
                'mode': 'markers',
                'type': 'scattergl',
                'marker': {
                    'size': 5,
                },
                'text': list(gene_names)
    }]
    return json.dumps({
            'data': data,
            'layout': {
                'title': 'Top genes for label {0} vs label {1}'.format(index_to_color[cluster1], index_to_color[cluster2]),
                'xaxis': {'title': 'log2 fold change', 'autorange': True},
                'yaxis': {'title': y_desc, 'autorange': True},
                'hovermode': 'closest',
                'margin': {'t': 40},
            },
    }, cls=SimpleEncoder)


def violin_plot_data(gene_values_all, cell_labels, gene_name,
        title=None, use_all_clusters=False, selected_clusters=None):
    """
    Returns a violin plot for a single gene, possibly over multiple clusters.

    Args:
        gene_values_cluster: 1d array of gene expression for the selected cluster, for the selected gene
        gene_values_all: 1d array of gene expression for all cells, for the selected gene
    """
    if title is None:
        title = 'Plot for gene {0}'.format(gene_name)
    selected_clusters = set(selected_clusters)
    if use_all_clusters:
        selected_clusters = set(cell_labels)
    selected_cells = np.array([x in selected_clusters for x in cell_labels])
    gene_values_cluster = gene_values_all[selected_cells]
    cell_labels = cell_labels[selected_cells]
    data = [{
        'y': gene_values_cluster.tolist(),
        'x': [str(x) for x in cell_labels],
        'type': 'violin',
        'opacity': 0.5,
        #'marker': {'color': 'green'},
        'box': True,
    },
    {
        'y': gene_values_all.tolist(),
        'type': 'violin',
        'opacity': 0.5,
        'name': 'all cells',
        #'marker': {'color': 'blue'},
        'box': True,
    }]
    return json.dumps({
        'data': data,
        'layout': {
            'title': title,
            'barmode': 'overlay',
            'showlegend': True,
            'xaxis': {'title': 'Cluster', 'type': 'category'},
            'yaxis': {'title': 'Gene level', 'zeroline': False},
            'width': max(400, 150*len(selected_clusters))
        },
    }, cls=SimpleEncoder)



@cache.memoize()
def heatmap_data(user_id, label_name_1, label_name_2, **params):
    """
    Returns a heatmap comparing two color tracks.
    """
    print('heatmap_data', user_id, label_name_1, label_name_2)
    from .advanced_plotting import cluster_heatmap
    color_track_1, is_discrete = get_sca_color_track(user_id, label_name_1)
    if not is_discrete:
        return 'should be a discrete colormap'
    color_track_2, is_discrete = get_sca_color_track(user_id, label_name_2)
    if not is_discrete:
        return 'should be a discrete colormap'
    return cluster_heatmap(color_track_1, color_track_2, label_name_1, label_name_2)

@cache.memoize()
def dendrogram_data(user_id, color_track_name, selected_genes, use_log=False, use_normalize=False):
    """
    Returns a dendrogram json
    """
    if color_track_name in ['entropy', 'gene', 'weights', 'read_count']:
        color_track_name = 'cluster'
    color_track, is_discrete = get_sca_color_track(user_id, color_track_name)
    all_genes = get_sca_gene_names(user_id)
    data = get_sca_data_sampled_all_genes(user_id)
    if len(selected_genes) == 0:
        # get top 5 genes from each cluster
        if color_track_name == 'cluster':
            top_genes = get_sca_top_1vr(user_id)
        else:
            top_genes, pvals = get_sca_top_genes_custom(user_id, color_track_name)
        for i, gene_set in top_genes.items():
            selected_top_genes = gene_set[:5]
            selected_gene_names = [all_genes[int(x[0])] for x in selected_top_genes]
            selected_genes += selected_gene_names
    from .advanced_plotting import dendrogram
    return dendrogram(data, all_genes, selected_genes, color_track_name, color_track, use_log=use_log, use_normalize=use_normalize)

@cache.memoize()
def cluster_correlation_heatmap_data(user_id, color_track_name, method='spearman'):
    """
    Correlation between the mean gene expression profiles of all clusters for the given color track.
    """
    if color_track_name in ['entropy', 'gene', 'weights', 'read_count']:
        color_track_name = 'cluster'
    color_track, is_discrete = get_sca_color_track(user_id, color_track_name)
    data = get_sca_data_sampled_all_genes(user_id)
    from .advanced_plotting import cluster_correlation_heatmap
    return cluster_correlation_heatmap(data, color_track, method)

@cache.memoize()
def gene_heatmap_data(user_id, genes_1, genes_2, color_track_name, cluster_id):
    """
    Returns a gene heatmap
    """
    data_sampled_all_genes = get_sca_data_sampled_all_genes(user_id)
    sca = get_sca(user_id)
    try:
        color_track, is_discrete = get_sca_color_track(user_id, color_track_name)
        if not is_discrete:
            raise Exception()
    except:
        color_track = sca.labels
    from .advanced_plotting import gene_similarity
    all_gene_names = get_sca_gene_names(user_id)
    if cluster_id == 'all':
        return gene_similarity(data_sampled_all_genes, all_gene_names, genes_1, genes_2)
    color_to_index, index_to_color = color_track_map(color_track)
    color_label = index_to_color[int(cluster_id)]
    return gene_similarity(data_sampled_all_genes[:,color_track==color_label], all_gene_names, genes_1, genes_2)

@cache.memoize()
def diff_corr_heatmap_data(user_id, genes_1, genes_2, color_track_name, cluster_id_1, cluster_id_2, value='p'):
    """
    Returns a gene heatmap
    """
    data_sampled_all_genes = get_sca_data_sampled_all_genes(user_id)
    sca = get_sca(user_id)
    try:
        color_track, is_discrete = get_sca_color_track(user_id, color_track_name)
        if not is_discrete:
            raise Exception()
    except:
        color_track = sca.labels
    from .advanced_plotting import differential_correlation
    all_gene_names = get_sca_gene_names(user_id)
    color_to_index, index_to_color = color_track_map(color_track)
    color_label_1 = index_to_color[int(cluster_id_1)]
    color_label_2 = index_to_color[int(cluster_id_2)]
    cells_1 = (color_track == color_label_1)
    cells_2 = (color_track == color_label_2)
    return differential_correlation(data_sampled_all_genes, all_gene_names, genes_1, genes_2, cells_1, cells_2, value=value)

@interaction_views.route('/user/<user_id>/stats')
def data_stats(user_id):
    """
    Returns html view showing data stats
    """
    from flask import Markup
    path = user_id_to_path(user_id)
    try:
        with open(os.path.join(path, 'params.json')) as f:
            params = json.load(f)
    except:
        params = {}
    try:
        with open(os.path.join(path, 'vis_summary.html')) as f:
            v = f.read()
        v = Markup(v)
    except:
        v = ''
    try:
        with open(os.path.join(path, 'read_count_hist_data.json')) as f:
            read_count_hist_data = f.read().strip()
        with open(os.path.join(path, 'gene_count_hist_data.json')) as f:
            gene_count_hist_data = f.read().strip()
        with open(os.path.join(path, 'gene_mean_hist_data.json')) as f:
            gene_mean_hist_data = f.read().strip()
    except:
        from . import data_stats
        sca = get_sca(user_id)
        summary = data_stats.Summary(data_paths=None, gene_paths=None, base_path=path, data=sca.data)
        read_count_hist_data, gene_count_hist_data, gene_mean_hist_data = summary.generate_plotly_jsons()
    # mean read count, median read count, mean gene count, median gene count
    # show genes per cell, switch to using plotly instead of bokeh (remove bokeh as a dependency)
    return render_template('stats.html',
            user_id=user_id,
            read_count_hist_data=read_count_hist_data,
            gene_count_hist_data=gene_count_hist_data,
            gene_mean_hist_data=gene_mean_hist_data,
            params=params)


@interaction_views.route('/user/<user_id>/view')
#@cache.memoize()
def view_plots(user_id):
    """
    Returns main HTML view.
    """
    test_or_user = 'user'
    data_user_id = user_id
    if user_id.startswith('test_'):
        test_or_user = 'test'
        data_user_id = user_id[5:]
    sca = get_sca(user_id)
    import cellmesh
    anatomy_id_names = cellmesh.get_all_cell_id_names(db_dir=cellmesh.ANATOMY_DB_DIR,
            include_cell_lines=True, include_chromosomes=True)
    anatomy_names = [x[1] for x in anatomy_id_names]
    cell_id_names = cellmesh.get_all_cell_id_names(include_cell_components=False)
    cell_names = [x[1] for x in cell_id_names]
    return render_template('state_estimation_static.html', user_id=user_id,
            test_or_user=test_or_user,
            data_user_id=data_user_id,
            gene_names=get_sca_gene_names(user_id),
            gene_sets=enrichr_api.ENRICHR_LIBRARIES,
            color_tracks=sca.get_color_track_names(),
            use_bacillus=True,
            anatomy_names=anatomy_names,
            cell_names=cell_names)


@interaction_views.route('/user/<user_id>/view/update_barplot', methods=['GET', 'POST'])
def update_barplot(user_id):
    """
    Updates barplot data.

    Returns:
        json corresponding to Plotly barplot
    """
    top_or_bulk = str(request.form['top_or_bulk'])
    input_value = int(request.form['input_value'])
    num_genes = int(request.form['num_genes'])
    data_form = request.form.copy()
    try:
        return update_barplot_result(user_id, top_or_bulk, input_value, num_genes,
                data_form)
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)

@cache.memoize()
def update_barplot_result(user_id, top_or_bulk, input_value, num_genes,
        data_form=None):
    """
    Generates a json-encoded string representing a plotly barplot.

    Args:
        user_id
        top_or_bulk (str): barplot option
        input_value (int): cell id
        num_genes (int): number of genes to include
        data_form (dict): copy of request.form
    """
    sca = get_sca(user_id)
    selected_gene = ''
    selected_gene_names = None
    gene_names = get_sca_gene_names(user_id)
    if 'selected_gene' in data_form:
        selected_gene = data_form['selected_gene']
    if len(selected_gene.strip()) > 0:
        gns = set(gene_names)
        selected_gene_names = split_gene_names(selected_gene)
        selected_gene_names = [x for x in selected_gene_names if x in gns]
    if top_or_bulk == 'volcano_pairwise':
        print('creating volcano plot')
        colormap = str(data_form['cell_color'])
        cluster1 = int(data_form['cluster1'])
        cluster2 = int(data_form['cluster2'])
        return volcano_plot_data(user_id, colormap, cluster1, cluster2, selected_genes=selected_gene_names)
    elif top_or_bulk == 'top_gene_expression':
        # get top genes by raw average expression
        colormap = str(data_form['cell_color'])
        cluster_id = int(input_value)
        sca = get_sca(user_id)
        data = get_sca_data_sampled_all_genes(user_id)
        input_label = input_value
        if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights', 'read_counts']:
            color_track, is_discrete = get_sca_color_track(user_id, colormap)
            color_to_index, index_to_color = color_track_map(color_track)
            color_label_1 = index_to_color[int(cluster_id)]
            selected_cells = (color_track == color_label_1)
            input_label = index_to_color[input_value]
        # default colormap
        else:
            labels = sca.labels
            selected_cells = (labels == cluster_id)
        data_subset = data[:, selected_cells]
        # calculate average gene expression for the selected cluster
        data_means = np.array(data_subset.mean(1)).flatten()
        print(np.sort(data_means)[::-1][:10])
        if selected_gene_names is not None:
            gene_indices = {g:i for i, g in enumerate(gene_names)}
            selected_gene_indices = np.array([gene_indices[g] for g in selected_gene_names])
            data_means_subset = data_means[selected_gene_indices]
        else:
            selected_gene_indices = np.argsort(data_means)[::-1]
            print(selected_gene_indices[:20])
            selected_gene_indices = selected_gene_indices[:num_genes]
            data_means_subset = data_means[selected_gene_indices]
            selected_gene_names = np.array([gene_names[i] for i in selected_gene_indices])
            print(selected_gene_names)
        return barplot_data(zip(range(len(selected_gene_names)), data_means_subset), selected_gene_names, input_label,
                    x_label='Mean gene expression', title='Top expressed genes in cluster {0}'.format(input_label))
    elif top_or_bulk == 'sep':
        # show separation score
        sep_scores = sca.separation_scores[int(input_value)]
        cluster_names = map(lambda x: 'cluster ' + str(x),
                list(range(sca.separation_scores.shape[0])))
        cluster_names = list(cluster_names)
        sep_scores = [(c, x) for c, x in zip(cluster_names, sep_scores)]
        return barplot_data(sep_scores,
                cluster_names, input_value,
                x_label='separation score',
                title='Inter-cluster separations for cluster {0}'.format(input_value))
    elif top_or_bulk == 'top_1_vs_rest' or top_or_bulk == 'pval_1_vs_rest':
        colormap = str(data_form['cell_color'])
        input_label = int(input_value)
        # custom colormap
        if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights', 'read_counts']:
            color_track, is_discrete = get_sca_color_track(user_id, colormap)
            lockfile_name = os.path.join(sca.data_dir, colormap + '_writing_diffexp')
            with lockfile_context(lockfile_name) as _lock:
                selected_diffexp, selected_pvals = get_sca_top_genes_custom(user_id, colormap)
            color_to_index, index_to_color = color_track_map(color_track)
            input_label = index_to_color[input_value]
            if top_or_bulk == 'pval_1_vs_rest':
                selected_diffexp = selected_pvals
            selected_diffexp = selected_diffexp[input_label]
        # default cluster colormap
        else:
            if top_or_bulk == 'top_1_vs_rest':
                selected_diffexp = get_sca_top_1vr(user_id)[input_label]
            else:
                selected_diffexp = get_sca_pval_1vr(user_id)[input_label]
        x_label = 'Fold change (1 vs rest)'
        if top_or_bulk == 'selected_color_pval':
            x_label = 'p-value of fold change (1 vs rest)'
        if 'pval' in top_or_bulk:
            # make note of whether p-val is FDR-corrected and
            # note that in the x-value
            is_fdr = sca.params['use_fdr']
            if is_fdr:
                x_label = 'FDR'
            else:
                x_label = 'p-value of fold change'
        # selected genes
        if selected_gene_names:
            selected_top_genes = [x for x in selected_diffexp if gene_names[int(x[0])] in set(selected_gene_names)]
        else:
            selected_top_genes = selected_diffexp[:num_genes]
        print('selected_top_genes:', selected_top_genes)
        selected_gene_names = [gene_names[int(x[0])] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_label,
                title='Top genes for cluster {0}'.format(input_label),
                x_label=x_label)
    elif top_or_bulk == 'top_pairwise' or top_or_bulk == 'pval_pairwise':
        # gets pairwise diffexp for a pair of clusters
        colormap = str(data_form['cell_color'])
        cluster1 = int(data_form['cluster1'])
        cluster2 = int(data_form['cluster2'])
        print('getting pairwise diffexp')
        print('colormap: ', str(colormap), ' clusters: ', cluster1, ' ', cluster2)
        use_baseline_clusters = True
        if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights']:
            try:
                color_track, is_discrete = get_sca_color_track(user_id, colormap)
                if is_discrete:
                    use_baseline_clusters = False
            except:
                pass
        # get the selected clusters
        if use_baseline_clusters:
            print('using default clustering')
            # Data is a numpy array of shape (k, k, genes)
            if top_or_bulk == 'top_pairwise':
                data = get_sca_pairwise_ratios(user_id)
                genes, values = array_to_top_genes(data, cluster1, cluster2, is_pvals=False, num_genes=num_genes)
                desc = 'ratios'
            else:
                data = get_sca_pairwise_pvals(user_id)
                genes, values = array_to_top_genes(data, cluster1, cluster2, is_pvals=True, num_genes=num_genes)
                is_fdr = sca.params['use_fdr']
                if is_fdr:
                    desc = 'FDR of ratios'
                else:
                    desc = 'p-value of ratios'
            # generate barplot
            if len(selected_gene.strip()) > 0:
                genes, values = array_to_top_genes(data, cluster1, cluster2, is_pvals=(top_or_bulk=='top_pairwise'), num_genes=1000000)
                gene_data = list(zip(genes, values))
                gene_data = [x for x in gene_data if gene_names[int(x[0])] in set(selected_gene_names)]
            else:
                gene_data = list(zip(genes, values))
            selected_gene_names = [gene_names[int(x[0])] for x in gene_data]
            return barplot_data(gene_data,
                    selected_gene_names, None,
                    title='Top genes for cluster {0} vs cluster {1}'.format(cluster1, cluster2),
                    x_label='Pairwise {0}'.format(desc))
        else:
            # get barplot for pairwise custom labels
            color_track, is_discrete = get_sca_color_track(user_id, colormap)
            color_to_index, index_to_color = color_track_map(color_track)
            print('using custom clustering')
            lockfile_name = os.path.join(sca.data_dir, colormap + '_writing_diffexp')
            with lockfile_context(lockfile_name) as _lock:
                selected_diffexp, selected_pvals = get_sca_top_genes_custom(user_id, colormap, 'pairwise')
            desc = ''
            if top_or_bulk == 'top_pairwise':
                genes, values = array_to_top_genes(selected_diffexp, cluster1, cluster2, is_pvals=False, num_genes=num_genes)
                desc = 'ratios'
            else:
                genes, values = array_to_top_genes(selected_pvals, cluster1, cluster2, is_pvals=True, num_genes=num_genes)
                is_fdr = sca.params['use_fdr']
                if is_fdr:
                    desc = 'FDR of ratios'
                else:
                    desc = 'p-value of ratios'
            if len(selected_gene.strip()) > 0:
                if top_or_bulk == 'top_pairwise':
                    genes, values = array_to_top_genes(selected_diffexp, cluster1, cluster2, is_pvals=(top_or_bulk=='top_pairwise'), num_genes=1000000)
                else:
                    genes, values = array_to_top_genes(selected_pvals, cluster1, cluster2, is_pvals=(top_or_bulk=='top_pairwise'), num_genes=1000000)
                gene_data = list(zip(genes, values))
                gene_data = [x for x in gene_data if gene_names[int(x[0])] in set(selected_gene_names)]
            else:
                gene_data = list(zip(genes, values))
            selected_gene_names = [gene_names[int(x[0])] for x in gene_data]
            return barplot_data(gene_data,
                    selected_gene_names, None,
                    title='Top genes for label {0} vs label {1}'.format(index_to_color[cluster1], index_to_color[cluster2]),
                    x_label='Pairwise {0}'.format(desc))
    elif top_or_bulk == 'hist':
        # generates a histogram
        selected_gene = data_form['selected_gene']
        use_baseline_clusters = True
        colormap = str(data_form['cell_color'])
        cluster_id = input_value
        if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights']:
            try:
                color_track, is_discrete = get_sca_color_track(user_id, colormap)
                if is_discrete:
                    use_baseline_clusters = False
            except:
                pass
        gene_data = get_gene_data(user_id, selected_gene)
        if use_baseline_clusters:
            gene_cluster_data = gene_data[sca.labels == cluster_id]
            return histogram_data(gene_cluster_data, gene_data, cluster_id, selected_gene)
        else:
            color_track, is_discrete = get_sca_color_track(user_id, colormap)
            color_to_index, index_to_color = color_track_map(color_track)
            gene_cluster_data = gene_data[color_track == index_to_color[cluster_id]]
            return histogram_data(gene_cluster_data, gene_data, index_to_color[cluster_id], selected_gene)
    elif top_or_bulk == 'double_pairs_comparison':
        print('double_pairs_comparison')
        print(data_form)
        colormap = str(data_form['cell_color'])
        c1 = str(data_form['cluster1'])
        c2 = str(data_form['cluster2'])
        c3 = str(data_form['cluster3'])
        c4 = str(data_form['cluster4'])
        return get_double_pairs_comparison_data(user_id, colormap, c1, c2, c3, c4, selected_genes=selected_gene_names)
    elif top_or_bulk == 'violin':
        print('violin_plot')
        print(data_form)
        selected_gene = data_form['selected_gene']
        use_log_transform = False
        if 'violin_use_log' in data_form:
            use_log_transform = (data_form['violin_use_log'] == '1')
        use_all_clusters = False
        if 'violin_all_clusters' in data_form:
            use_all_clusters = (data_form['violin_all_clusters'] == '1')
        use_baseline_clusters = True
        colormap = str(data_form['cell_color'])
        cluster_id = input_value
        all_selected_clusters = [cluster_id]
        if 'all_selected_clusters[]' in data_form:
            all_selected_clusters = data_form.getlist('all_selected_clusters[]', type=int)
            print('selected clusters for violin plot:', all_selected_clusters)
        if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights']:
            try:
                color_track, is_discrete = get_sca_color_track(user_id, colormap)
                if is_discrete:
                    use_baseline_clusters = False
            except:
                pass
        gene_data = get_gene_data(user_id, selected_gene)
        if use_log_transform:
            gene_data = np.log10(gene_data + 1)
        if use_baseline_clusters:
            return violin_plot_data(gene_data, sca.labels, selected_gene, use_all_clusters=use_all_clusters, selected_clusters=all_selected_clusters)
        else:
            color_track, is_discrete = get_sca_color_track(user_id, colormap)
            color_to_index, index_to_color = color_track_map(color_track)
            all_selected_clusters = [index_to_color[c] for c in all_selected_clusters]
            return violin_plot_data(gene_data, color_track, selected_gene, use_all_clusters=use_all_clusters, selected_clusters=all_selected_clusters)
    elif top_or_bulk == 'gene_gene':
        pass
    else:
        return 'Error: '


@cache.memoize()
def get_double_pairs_comparison_data(user_id, colormap, c1, c2, c3, c4, nonzero_threshold=0, selected_genes=None):
    """
    Plot a two-dimensional scatterplot: x-axis shows cluster1-cluster2, y-axis shows cluster3-cluster4
    """
    print('get_double_pairs_comparison_data')
    if colormap in ['cluster', 'gene', 'entropy', 'weights']:
        colormap = 'cluster'
    c1 = int(c1)
    c2 = int(c2)
    c3 = int(c3)
    c4 = int(c4)
    print(c1, c2, c3, c4)
    data_sampled_all_genes = get_sca_data_sampled_all_genes(user_id)
    gene_names = get_sca_gene_names(user_id)
    color_track, is_discrete = get_sca_color_track(user_id, colormap)
    color_to_index, index_to_color = color_track_map(color_track)
    # TODO: use log fold change instead of normalized difference?
    # get selected genes
    if selected_genes is not None:
        gene_indices = {g:i for i, g in enumerate(gene_names)}
        selected_gene_indices = np.array([gene_indices[g] for g in selected_genes])
        #selected_diffexp = selected_diffexp[:,:,selected_gene_indices]
        #selected_pvals = selected_pvals[:,:,selected_gene_indices]
        data_sampled_all_genes = data_sampled_all_genes[selected_gene_indices,:]
        gene_names = selected_genes
    # get means for each of the clusters
    gene_nonzero_counts = data_sampled_all_genes.getnnz(1)
    c1_mean = data_sampled_all_genes[:, color_track == index_to_color[c1]].mean(1)
    c1_mean = np.array(c1_mean).flatten()
    c2_mean = data_sampled_all_genes[:, color_track == index_to_color[c2]].mean(1)
    c2_mean = np.array(c2_mean).flatten()
    c3_mean = data_sampled_all_genes[:, color_track == index_to_color[c3]].mean(1)
    c3_mean = np.array(c3_mean).flatten()
    c4_mean = data_sampled_all_genes[:, color_track == index_to_color[c4]].mean(1)
    c4_mean = np.array(c4_mean).flatten()
    # : by sum of c1 and c2
    if nonzero_threshold == 0:
        nonzero_threshold = float(data_sampled_all_genes.shape[1])/200
    c1_c2 = (c1_mean - c2_mean)/(c1_mean + c2_mean + 1e-8)
    c1_c2[gene_nonzero_counts <= nonzero_threshold] = 0
    c3_c4 = (c3_mean - c4_mean)/(c3_mean + c4_mean + 1e-8)
    c3_c4[gene_nonzero_counts <= nonzero_threshold] = 0
    q1_count = sum((c1_c2 < 0) & (c3_c4 > 0))
    q2_count = sum((c1_c2 > 0) & (c3_c4 > 0))
    q3_count = sum((c1_c2 > 0) & (c3_c4 < 0))
    q4_count = sum((c1_c2 < 0) & (c3_c4 < 0))
    # also print quadrant counts somehow?
    output = {
        'data': [
            {
                'x': c1_c2,
                'y': c3_c4,
                'colorscale': 'Reds',
                'type': 'scattergl',
                'mode': 'markers',
                'text': list(gene_names),
                'marker': {'size': 5},
            },
            {
                'x': [-1, 1, 1, -1],
                'y': [1, 1, -1, -1],
                'mode': 'text',
                'text': [str(x) for x in [q1_count, q2_count, q3_count, q4_count]],
                'type': 'scatter',
                'textfont': {
                    'size': 18,
                    'color': '#ff7f0e',
                }
            }
        ],
        'layout': {
            'xaxis': {'title': 'cluster1 - cluster2', 'automargin': True},
            'yaxis': {'title': 'cluster3 - cluster4', 'automargin': True},
            'margin': {'t': 40},
            'hovermode': 'closest',
            'showlegend': False,
        },
        'misc': {
            'q1_count': q1_count,
            'q2_count': q2_count,
            'q3_count': q3_count,
            'q4_count': q4_count,
        }
    }
    return json.dumps(output, cls=SimpleEncoder)

@interaction_views.route('/user/<user_id>/view/update_scatterplot', methods=['GET', 'POST'])
def update_scatterplot(user_id):
    """
    Updates the scatterplot view.

    Returns: a plotly-formated json
        {data: [{'cluster': c, 'x': x, 'y'; y}...],
         labels: [c1, c2, ...],
         colorscale: 'Viridis'}
    """
    plot_type = request.form['scatter_type']
    cell_color_value = request.form['cell_color']
    try:
        return update_scatterplot_result(user_id, plot_type, cell_color_value,
                request.form.copy())
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)

@cache.memoize()
def update_scatterplot_result(user_id, plot_type, cell_color_value, data_form):
    """
    Returns the plotly JSON representation of the scatterplot.
    """
    sca = get_sca(user_id)
    if plot_type == 'Means':
        labels = np.arange(sca.mds_means.shape[1])
        return scatterplot_data(sca.mds_means,
                labels)
    elif plot_type == 'Cluster_heatmap':
        label_name_1 = data_form['heatmap_cluster_name_1']
        label_name_2 = data_form['heatmap_cluster_name_2']
        return heatmap_data(user_id, label_name_1, label_name_2)
    elif plot_type == 'Dendrogram':
        color_track_name = data_form['cell_color']
        color_track, is_discrete = get_sca_color_track(user_id, color_track_name)
        selected_genes = split_gene_names(data_form['dendrogram_genes'])
        use_log = 'dendrogram_use_log' in data_form and data_form['dendrogram_use_log'] != '0'
        use_normalize = 'dendrogram_normalize' in data_form and data_form['dendrogram_normalize'] != '0'
        return dendrogram_data(user_id, color_track_name, selected_genes, use_log=use_log, use_normalize=use_normalize)
    elif plot_type == 'Genes':
        print('plotting genes')
        dim_red = sca.gene_dim_red
        labels = sca.gene_clusters
        gene_names = get_sca_gene_names(user_id)
        return scatterplot_data(dim_red, labels,
                label_text=gene_names)
    elif plot_type == 'Gene_heatmap':
        print('plotting gene heatmap')
        gene_names_1 = split_gene_names(data_form['heatmap_genes_1'])
        gene_names_2 = split_gene_names(data_form['heatmap_genes_2'])
        color_track_name = data_form['cell_color']
        cluster_id = data_form['gene_heatmap_cluster']
        return gene_heatmap_data(user_id, gene_names_1, gene_names_2,
                color_track_name=color_track_name, cluster_id=cluster_id)
    elif plot_type == 'Diffcorr_heatmap':
        print('plotting differential correlation gene heatmap')
        gene_names_1 = split_gene_names(data_form['diffcorr_genes_1'])
        gene_names_2 = split_gene_names(data_form['diffcorr_genes_2'])
        color_track_name = data_form['cell_color']
        cluster_id_1 = data_form['diffcorr_cluster_1']
        cluster_id_2 = data_form['diffcorr_cluster_2']
        value = data_form['diffcorr_value']
        return diff_corr_heatmap_data(user_id, gene_names_1, gene_names_2,
                color_track_name=color_track_name, cluster_id_1=cluster_id_1,
                cluster_id_2=cluster_id_2, value=value)
    elif plot_type == 'Correlation_heatmap':
        print('plotting cluster correlation heatmap')
        color_track_name = data_form['cell_color']
        return cluster_correlation_heatmap_data(user_id, color_track_name)
    else:
        dim_red = None
        if plot_type == 'Cells':
            dim_red = get_sca_dim_red(user_id)
        elif plot_type == 'Baseline':
            dim_red = get_sca_baseline_vis(user_id)
        if cell_color_value == 'entropy':
            return scatterplot_data(dim_red, sca.labels,
                    colorscale='Viridis',
                    mode='entropy', color_vals=sca.entropy)
        elif cell_color_value == 'gene':
            gene_name = data_form['gene_name']
            use_mw = False
            if 'use_mw' in data_form:
                use_mw = bool(int(data_form['use_mw']))
            gene_data = get_gene_data(user_id, gene_name, use_mw)
            if len(gene_data)==0:
                return 'Error: gene not found'
            return scatterplot_data(dim_red, sca.labels,
                    mode='entropy', color_vals=gene_data)
        elif cell_color_value == 'cluster':
            return scatterplot_data(dim_red, sca.labels)
        elif cell_color_value == 'new':
            # this usually happens by mistake...
            return scatterplot_data(dim_red, sca.labels)
        # if the mode is 'cluster', color based on w
        elif cell_color_value == 'weights':
            cluster = int(data_form['cluster_input'])
            w = sca.w_sampled
            if cluster < 0 or cluster >= w.shape[0]:
                return 'Error: invalid cluster ID'
            return scatterplot_data(dim_red, sca.labels,
                    mode='entropy', color_vals=w[cluster, :])
        elif cell_color_value == 'read_counts':
            read_counts = sca.read_counts
            read_counts = read_counts[sca.cell_subset][sca.cell_sample]
            return scatterplot_data(dim_red, sca.labels,
                    colorscale='Viridis',
                    mode='entropy', color_vals=read_counts)
        elif cell_color_value == 'neural_network_classifier':
            # TODO: get NN results
            color_track, is_discrete = get_sca_color_track(user_id, cell_color_value)
            # set cell class...
            if color_track is None:
                from mouse_cell_query import nn_query
                cell_names, results, class_names = nn_query.predict_using_default_classifier(sca.data.T, sca.genes)
                sca.add_color_track('neural_network_classifier', cell_names, is_discrete=True)
                color_track, is_discrete = get_sca_color_track(user_id, cell_color_value)
                return scatterplot_data(dim_red, color_track)
            else:
                return scatterplot_data(dim_red, color_track)
        else:
            # try to get color track
            # TODO: get color values as well
            color_track, is_discrete, color = get_sca_color_track(user_id, cell_color_value, return_color=True)
            print('scatterplot retrieved color:', color)
            if color_track is None:
                return scatterplot_data(dim_red, sca.labels)
            else:
                if is_discrete:
                    return scatterplot_data(dim_red, color_track, color_dict=color)
                else:
                    return scatterplot_data(dim_red, sca.labels,
                            mode='entropy', color_vals=color_track)

@cache.memoize()
def get_gene_data(user_id, gene_name, use_mw=False):
    """
    Returns an array containing data for a given gene name.
    """
    sca = get_sca(user_id)
    if gene_name is None:
        return None
    gene_data = sca.data_sampled_gene(gene_name, use_mw=use_mw)
    if len(gene_data) == 0:
        return None
    return gene_data

@interaction_views.route('/user/<user_id>/view/cell_info', methods=['GET', 'POST'])
def cell_info(user_id):
    """
    Gets basic statistics about a cell:
        - cell_id (int, 0-indexed)
        - read_count (number)
        - genes_count (int)
        - cluster (int)
        - cluster_median_read_count (float)
        - cluster_median_gene_count (float)
        - cluster_total_gene_count (float) - total genes that are nonzero in any cell in the cluster
        - values for all of the uploaded color maps
    """
    print('cell_info: ', request.form)
    selected_cells = request.form['selected_cells']
    selected_cells = selected_cells.split(',')
    selected_cells = [int(x) for x in selected_cells]
    selected_clusters = request.form['selected_clusters']
    selected_clusters = selected_clusters.split(',')
    selected_clusters = [int(x) for x in selected_clusters]
    color_map = request.form['color_map']
    # TODO: return more results: more cluster info - cluster total count
    return cell_info_result(user_id, selected_cells, selected_clusters, color_map)

@cache.memoize()
def cell_info_result(user_id, selected_cells, selected_clusters, color_map):
    # get read count + gene count for all cells
    # TODO: don't really need cell info; need more cluster info
    sca = get_sca(user_id)
    read_counts = []
    gene_counts = []
    data = get_sca_data_sampled_all_genes(user_id)
    #for cell in selected_cells:
    #    cell_data = data[:, cell]
    #    gene_counts.append(cell_data.count_nonzero())
    #    read_counts.append(cell_data.sum())
    cluster_reads, cluster_genes, total_gene_count = get_cluster_stats(sca, data, selected_clusters[0], color_map)
    return json.dumps({'gene_counts': gene_counts,
        'read_counts': read_counts, 'cluster_reads': cluster_reads,
        'cluster_genes': cluster_genes, 'total_gene_count': total_gene_count}, cls=SimpleEncoder)

def get_cluster_stats(sca, data, cluster_id, colormap='cluster'):
    # get cluster median read count, cluster median gene count
    color_track = sca.labels
    if colormap != 'cluster':
        try:
            color_track, is_discrete = sca.get_color_track(colormap)
            if not is_discrete:
                color_track = sca.labels
            else:
                color_to_index, index_to_color = color_track_map(color_track)
                cluster_id = index_to_color[cluster_id]
        except:
            color_track = sca.labels
    cell_data = data[:, color_track==cluster_id]
    read_counts = np.array(cell_data.sum(0)).flatten()
    gene_counts = np.zeros(cell_data.shape[1])
    # TODO: get total count of nonzero genes in cluster
    total_gene_sum = np.array(cell_data.sum(1)).flatten()
    total_gene_count = (total_gene_sum > 0).sum()
    for i in range(len(gene_counts)):
        gene_counts[i] = cell_data[:,i].count_nonzero()
    return np.median(read_counts), np.median(gene_counts), total_gene_count

@interaction_views.route('/user/<user_id>/view/update_enrichr', methods=['GET', 'POST'])
def update_enrichr(user_id):
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = request.form['top_genes']
    # gene_set is a string.
    gene_set = request.form['gene_set']
    return update_enrichr_result(user_id, top_genes, gene_set)

# TODO: cache this function??? but we don't want to cache timeouts
def update_enrichr_result(user_id, top_genes, gene_set):
    user_list_id = 0
    if top_genes not in interaction_views.enrichr_gene_list_ids:
        gene_list = split_gene_names(top_genes)
        user_list_id = enrichr_api.enrichr_add_list(gene_list)
        if user_list_id == 'timeout':
            return 'Error: Enrichr query timed out'
        interaction_views.enrichr_gene_list_ids[top_genes] = user_list_id
    else:
        user_list_id = interaction_views.enrichr_gene_list_ids[top_genes]
    results = []
    if (top_genes, gene_set) in interaction_views.enrichr_results:
        results = interaction_views.enrichr_results[(top_genes, gene_set)]
    else:
        try:
            results = enrichr_api.enrichr_query(user_list_id, gene_set)
            if results == 'timeout':
                return 'Error: Enrichr query timed out'
            interaction_views.enrichr_results[(top_genes, gene_set)] = results[:10]
        except:
            gene_list = split_gene_names(top_genes)
            user_list_id = enrichr_api.enrichr_add_list(gene_list)
            if user_list_id == 'timeout':
                return 'Error: Enrichr query timed out'
            interaction_views.enrichr_gene_list_ids[top_genes] = user_list_id
            results = enrichr_api.enrichr_query(user_list_id, gene_set)
            if results == 'timeout':
                return 'Error: Enrichr query timed out'
            interaction_views.enrichr_results[(top_genes, gene_set)] = results[:10]
    # only take top 10 results (maybe have this value be variable?)
    results = results[:10]
    interaction_views.last_enrichr_results = [['gene set name',
                                 'p-value',
                                 'z-score',
                                 'combined score']] + \
            [[r[1], r[2], r[3], r[4]] for r in results]
    return json.dumps(interaction_views.last_enrichr_results, cls=SimpleEncoder)

@interaction_views.route('/user/<user_id>/view/update_cellmarker', methods=['GET', 'POST'])
def update_cellmarker(user_id):
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = [x.strip().upper() for x in split_gene_names(request.form['top_genes'])]
    print('update_cellmarker:', top_genes)
    # gene_set is a string.
    test_type = request.form['test_type']
    cells_or_tissues = request.form['cells_or_tissues']
    species = request.form['species']
    return update_cellmarker_result(user_id, top_genes, test_type, cells_or_tissues, species)

@cache.memoize()
def update_cellmarker_result(user_id, top_genes, test, cells_or_tissues, species):
    """
    Gets the CellMarker result for a set of genes.

    Args:
        user_id (str)
        top_genes (list of strings representing genes)
        test (str, currently only 'hypergeom')
        cells_or_tissues (str, either 'cells' or 'tissues')
        species (str, one of 'all', 'Human', 'Mouse')
    """

    import cellmarker
    result = []
    if test == 'hypergeom':
        result = cellmarker.hypergeometric_test(top_genes, cells_or_tissues, return_header=True, return_cl=True, species=species)
    cell_types = [result[0]]
    for i in range(1, min(20, len(result))):
        ri = result[i]
        genes = ri[3]
        gene_pmids = []
        for g in genes:
            gene_pmids.append('{0}: {1}'.format(g, ', '.join(pmid_to_link(x) for x in ri[4][g])))
        cell_types.append((ri[0], ri[1], ri[2], ', '.join(ri[3]), ', '.join(gene_pmids)))
    return json.dumps(cell_types, cls=SimpleEncoder)

@interaction_views.route('/user/<user_id>/view/update_cellmesh', methods=['GET', 'POST'])
def update_cellmesh(user_id):
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = [x.strip().upper() for x in split_gene_names(request.form['top_genes'])]
    print('update_cellmesh:', top_genes)
    # gene_set is a string.
    test_type = request.form['mesh_test_type']
    species = request.form['cellmesh_species']
    return update_cellmesh_result(user_id, top_genes, test_type, species)

@cache.memoize()
def update_cellmesh_result(user_id, top_genes, test, species='human', return_json=True):
    """
    Gets the CellMesh result for a set of genes.

    Args:
        user_id (str)
        top_genes (list of strings representing genes)
        test (str, currently only 'hypergeom')
    """
    import cellmesh
    result = []
    if test == 'hypergeom':
        result = cellmesh.hypergeometric_test(top_genes, species=species, return_header=True)
    elif test == 'norm_hypergeom':
        result = cellmesh.normed_hypergeometric_test(top_genes, species=species, return_header=True)
    elif test == 'prob':
        from cellmesh import prob_method
        result = prob_method.prob_test(top_genes, species=species, return_header=True)
    elif test == 'gsva':
        from cellmesh import gsva_ext_method
        result = gsva_ext_method.gsva_ext_test(top_genes, species=species, return_header=True)
    cell_types = [result[0]]
    for i in range(1, min(20, len(result))):
        ri = result[i]
        gene_pmids = []
        genes = ri[3]
        for g in genes:
            gene_pmids.append('{0}: {1}'.format(g, ', '.join(pmid_to_link(x) for x in ri[4][g])))
        cell_types.append([ri[0], ri[1], ri[2], ', '.join(ri[3]), ', '.join(gene_pmids)])
    if return_json:
        return json.dumps(cell_types, cls=SimpleEncoder)
    else:
        return cell_types

@interaction_views.route('/user/<user_id>/view/update_cellmesh_anatomy', methods=['GET', 'POST'])
def update_cellmesh_anatomy(user_id):
    top_genes = [x.strip().upper() for x in split_gene_names(request.form['top_genes'])]
    mesh_subset = request.form['anatomy_mesh_subset'].strip()
    species = request.form['anatomy_species']
    test = request.form['anatomy_mesh_test_type']
    return update_cellmesh_anatomy_result(top_genes, mesh_subset=mesh_subset, species=species, test=test)

@cache.memoize()
def update_cellmesh_anatomy_result(top_genes, mesh_subset=None, species='human', return_json=True, test='hypergeom'):
    import cellmesh
    if len(mesh_subset) > 1:
        mesh_subset = [cellmesh.get_cell_id_from_name(mesh_subset,
            db_dir=cellmesh.ANATOMY_DB_DIR)]
    else:
        mesh_subset = None
    # TODO: validate mesh_subset
    if test == 'hypergeom':
        result = cellmesh.hypergeometric_test(top_genes, species=species, return_header=True, db_dir=cellmesh.ANATOMY_DB_DIR,
                cell_type_subset=mesh_subset)
    elif test == 'norm_hypergeom':
        result = cellmesh.normed_hypergeometric_test(top_genes, species=species, return_header=True,
                db_dir=cellmesh.ANATOMY_DB_DIR)
    elif test == 'prob':
        from cellmesh import prob_method
        result = prob_method.prob_test(top_genes, species=species, return_header=True, db_dir=cellmesh.ANATOMY_DB_DIR, cell_type_subset=mesh_subset)
    elif test == 'gsva':
        from cellmesh import gsva_ext_method
        result = gsva_ext_method.gsva_ext_test(top_genes, species=species, return_header=True, db_dir=cellmesh.ANATOMY_DB_DIR)
    print('update_cellmesh_anatomy_result', result)
    cell_types = [result[0]]
    for i in range(1, min(20, len(result))):
        ri = result[i]
        gene_pmids = []
        genes = ri[3]
        for g in genes:
            gene_pmids.append('{0}: {1}'.format(g, ', '.join(pmid_to_link(x) for x in ri[4][g])))
        cell_types.append([ri[0], ri[1], ri[2], ', '.join(ri[3]), ', '.join(gene_pmids)])
    if return_json:
        return json.dumps(cell_types, cls=SimpleEncoder)
    else:
        return cell_types

@interaction_views.route('/user/<user_id>/view/update_go', methods=['GET', 'POST'])
def update_go(user_id):
    """
    get gene ontology result
    """
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = [x.strip().upper() for x in split_gene_names(request.form['top_genes'])]
    species = request.form['go_species']
    fdr_threshold = 0.2
    if 'go_fdr_threshold' in request.form:
        fdr_threshold = float(request.form['go_fdr_threshold'])
    print('update_go:', top_genes)
    # gene_set is a string.
    return update_go_result(top_genes, species=species, fdr_threshold=fdr_threshold)

@cache.memoize()
def update_go_result(top_genes, species='mouse', fdr_threshold=0.2, **kwargs):
    from cellmesh import go_query
    top_genes = [x.capitalize() for x in top_genes]
    result = go_query.gene_set_query(top_genes, return_header=True, species=species, fdr_threshold=fdr_threshold)
    for r in result[1:]:
        r[3] = ', '.join(r[3])
    return json.dumps(result, cls=SimpleEncoder)

@interaction_views.route('/user/<user_id>/view/update_subtiwiki', methods=['GET', 'POST'])
def update_subtiwiki(user_id):
    """
    get gene ontology result
    """
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = [x for x in split_gene_names(request.form['top_genes'])]
    mode = request.form['subtiwiki_mode']
    print('update_subtiwiki:', top_genes)
    # gene_set is a string.
    return update_subtiwiki_result(top_genes, mode)

@cache.memoize()
def update_subtiwiki_result(top_genes, mode='all', **kwargs):
    import subtiwiki
    # no capitalization change
    if mode == 'gene_info':
        result = subtiwiki.get_gene_info(top_genes, return_header=True)
    else:
        result = subtiwiki.hypergeometric_test(top_genes, return_header=True, mode=mode)
        result = [list(x) for x in result]
        for r in result[1:]:
            r[3] = ', '.join(r[3])
    print('update_subtiwiki_result:', result)
    return json.dumps(result, cls=SimpleEncoder)

@interaction_views.route('/user/<user_id>/view/update_kegg', methods=['GET', 'POST'])
def update_kegg(user_id):
    """
    get KEGG result
    """
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = [x for x in split_gene_names(request.form['top_genes'])]
    species = request.form['kegg_species']
    print('update_kegg:', top_genes)
    # gene_set is a string.
    return update_kegg_result(top_genes, species=species)

@cache.memoize()
def update_kegg_result(top_genes, species='human', **kwargs):
    import kegg_query
    # no capitalization change
    result = kegg_query.hypergeometric_test(top_genes, return_header=True, species=species)
    result = [list(x) for x in result]
    print('update_kegg_result:', result)
    return json.dumps(result, cls=SimpleEncoder)

@interaction_views.route('/user/<user_id>/view/db_query', methods=['POST'])
def db_query(user_id):
    """
    Queries a single cell db for cell types
    """
    import mouse_cell_query
    sca = get_sca(user_id)
    form_data = request.form.copy()
    db = form_data['cell_search_db']
    cell_color = form_data['cell_color']
    cell_label = form_data['cell_search_cluster']
    means = None
    if cell_color == 'cluster':
        cell_label = cell_label.split()[-1]
        means = sca.cluster_means[:, int(cell_label)]
    else:
        try:
            labels, is_discrete = sca.get_color_track(cell_color)
            color_to_index, index_to_color = color_track_map(labels)
            labels = index_to_color[int(cell_label)]
            data = get_sca_data_sampled_all_genes(user_id)
            data_labels = data[:, labels==cell_label]
            means = np.array(data_labels.mean(1)).flatten()
        except Exception as e:
            text = traceback.format_exc()
            print(text)
            means = sca.cluster_means[:, int(cell_label)]
    # TODO: cache db query
    try:
        results = mouse_cell_query.search_db(means, sca.gene_names, method=form_data['method'], db=db)
        results = [('Cell type', 'Score')] + results
        return json.dumps(results, cls=SimpleEncoder)
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)



@interaction_views.route('/user/<user_id>/view/split_or_merge_cluster', methods=['POST'])
def split_or_merge_cluster(user_id):
    if user_id.startswith('test_'):
        return 'Error: test datasets cannot be modified. Copy the dataset if you wish to modify it.'
    split_or_merge = request.form['split_or_merge']
    selected_clusters = request.form['selected_clusters']
    selected_clusters = selected_clusters.split(',')
    selected_clusters = [int(x) for x in selected_clusters]
    print('split_or_merge:', split_or_merge)
    print('selected_clusters:', selected_clusters)
    sca = get_sca(user_id)
    selected_clusters = list(set(selected_clusters))
    if len(selected_clusters) == 0:
        return 'Error: no selected clusters.'
    # TODO: have a more fine-grained key control. don't just clear the entire
    # cache, but clear some keys selectively from redis.
    if 'DEPLOY' in current_app.config and current_app.config['DEPLOY']:
        print('clearing cache')
        cache.delete_memoized(get_sca_top_genes, user_id)
        cache.delete_memoized(get_sca_top_1vr, user_id)
        cache.delete_memoized(get_sca_pvals, user_id)
        cache.delete_memoized(get_sca_pval_1vr, user_id)
        cache.delete_memoized(update_barplot_result)
        cache.delete_memoized(update_scatterplot_result)
        cache.delete_memoized(heatmap_data)
        cache.delete_memoized(dendrogram_data)
        #print('deleting user_id from cache')
        # TODO: this currently doesn't work, since keys are hashed.
        #clear_cache_user_id(user_id)
        cache.clear()
    else:
        print('clearing cache')
        cache.clear()
    # split clusters
    if split_or_merge == 'split':
        try:
            generate_analysis.generate_analysis_resubmit(sca,
                    'split', selected_clusters)
            return 'Finished splitting selected cluster: ' + str(selected_clusters[0])
        except:
            return 'Error in splitting clusters.'
    # merge clusters
    elif  split_or_merge == 'merge':
        try:
            generate_analysis.generate_analysis_resubmit(sca,
                    'merge', selected_clusters)
            return 'Finished merging selected clusters: ' + ' '.join(map(str, selected_clusters))
        except:
            return 'Error in merging clusters.'
    # create new cluster from selected cells
    elif split_or_merge == 'new':
        try:
            generate_analysis.generate_analysis_resubmit(sca,
                    'new', selected_clusters)
            return 'Finished creating new cluster from selected cells: ' + ' '.join(map(str, selected_clusters))
        except:
            return 'Error in creating new cluster.'
    # delete selected cells
    elif split_or_merge == 'delete':
        try:
            generate_analysis.generate_analysis_resubmit(sca,
                    'delete', selected_clusters)
            return 'Finished deleting selected cells: ' + ' '.join(map(str, selected_clusters))
        except:
            return 'Error in deleting selected cells.'


@interaction_views.route('/user/<user_id>/view/upload_color_track', methods=['POST'])
def upload_color_track(user_id):
    from werkzeug.utils import secure_filename
    sca = get_sca(user_id)
    print('upload_color_track:', request.form)
    print(request.files)
    if 'color_track_file' in request.files:
        f = request.files['color_track_file']
        output_filename = secure_filename(f.filename)
        output_filename = output_filename.split('.')[0]
        color_track_type = request.form['color_track_type']
        if color_track_type == 'continuous':
            data = np.loadtxt(f, dtype=np.float)
            data = data.flatten()
            sca.add_color_track(output_filename, data, False)
        elif color_track_type == 'discrete':
            data = np.loadtxt(f, dtype=str, delimiter='\t')
            data = data.flatten()
            sca.add_color_track(output_filename, data, True)
        elif color_track_type == 'table':
            data = np.loadtxt(f, dtype=str, delimiter='\t')
            for i in range(data.shape[1]):
                column = data[:,i]
                column_data = column[1:]
                column_name = column[0]
                is_discrete = True
                try:
                    column_data = column_data.astype(float)
                    is_discrete = False
                except:
                    pass
                sca.add_color_track(column_name, column[1:], is_discrete)
    # return new color tracks
    return redirect(url_for('interaction_views.view_plots', user_id=user_id))


@interaction_views.route('/user/<user_id>/view/custom_color_map', methods=['POST'])
def custom_color_map(user_id):
    """
    Creates a custom color map, based on user-defined gene selections.
    """
    data_form = request.form.copy()
    try:
        name = data_form['name']
        sca = get_sca(user_id)
        sca.create_custom_selection(name)
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)
    return 'success'


@interaction_views.route('/user/<user_id>/view/get_custom_colormap', methods=['POST'])
def get_custom_colormap(user_id):
    """
    This gets the labels and criteria for a given color map...
    """
    data_form = request.form.copy()
    try:
        name = data_form['name']
        sca = get_sca(user_id)
        results = sca.custom_selections[name];
        return custom_cell_selection.create_json(results)
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)


@interaction_views.route('/user/<user_id>/view/get_colormap_label_criteria', methods=['POST'])
def get_colormap_label_criteria(user_id):
    """
    Returns the labels and criteria for a given colormap.
    """
    data_form = request.form.copy()
    try:
        name = data_form['name']
        sca = get_sca(user_id)
        results = sca.custom_selections[name];
        if 'label' in data_form:
            label_name = data_form['label']
            for label in results.labels:
                if label.name == label_name:
                    return custom_cell_selection.create_json(label)
        if name in sca.custom_selections:
            return custom_cell_selection.create_json(results)
        else:
            return 'Error: name not in custom selections'
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)

@interaction_views.route('/user/<user_id>/view/get_colormap_values', methods=['POST'])
def get_colormap_values(user_id):
    """
    Returns the values for a given colormap.
    """
    data_form = request.form.copy()
    try:
        name = data_form['name']
        sca = get_sca(user_id)
        values = sca.get_color_track_values(name)
        return json.dumps(list(values), cls=SimpleEncoder)
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)

def load_criteria_from_dict(json_dict):
    """
    Returns a list of LabelCriterion objects given a dict loaded from json...
    """
    keys = json_dict.keys()
    ids = set([int(x.split('-')[1]) for x in keys])
    all_criteria = []
    for current_id in ids:
        selection_type = json_dict['selection_type-'+str(current_id)]
        comparison = json_dict['selection_comparison-'+str(current_id)]
        target = json_dict['selection_target-'+str(current_id)]
        and_or = json_dict['selection_and_or-'+str(current_id)]
        criterion = custom_cell_selection.LabelCriterion(selection_type=selection_type, comparison=comparison, target=target, and_or=and_or)
        if 'selection_value-' + str(current_id) in json_dict:
            value = json_dict['selection_value-'+str(current_id)]
            criterion.value = value
        all_criteria.append(criterion)
    return all_criteria

@interaction_views.route('/user/<user_id>/view/update_colormap_label_criteria', methods=['POST'])
def update_colormap_label_criteria(user_id):
    """
    Updates the criteria for a given label in a given colormap, or just returns the existing
    label.

    Returns a json representation of the new label.
    """
    data_form = request.form.copy()
    print(data_form)
    sca = get_sca(user_id)
    colormap_name = data_form['name']
    label_name = data_form['label']
    color = data_form['color']
    print('colormap color:', color)
    if color == '#000000':
        color = None
    if 'criteria' in data_form:
        # load criteria from json
        criteria = load_criteria_from_dict(json.loads(data_form['criteria']))
        sca.update_custom_color_track_label(colormap_name, label_name, criteria, color=color)
        # clear cache for scatterplot results
        print('deleting cached results...')
        cache.delete_memoized(update_barplot_result)
        cache.delete_memoized(update_scatterplot_result)
        cache.delete_memoized(get_sca_color_track)
        cache.delete_memoized(dendrogram_data)
        cache.delete_memoized(heatmap_data)
        cache.delete_memoized(get_sca_top_genes_custom)
    else:
        sca.update_custom_color_track_label(colormap_name, label_name)
    colormap = sca.custom_selections[colormap_name]
    for label in colormap.labels:
        if label.name == label_name:
            return custom_cell_selection.create_json(label)
    return ''

@interaction_views.route('/user/<user_id>/view/copy_dataset', methods=['POST'])
def copy_dataset(user_id):
    """
    Copies the input dataset into a new user id.
    """
    sca = get_sca(user_id)
    path = sca.data_dir
    try:
        new_user_id = str(uuid.uuid4())
        new_user_id = new_user_id + user_id[36:]
        shutil.copytree(path, user_id_to_path(new_user_id, use_secondary=False))
        # change user id in json files (this is a bad hack lol)
        import subprocess
        subprocess.call("sed -i 's/{0}/{1}/g' /tmp/uncurl/{1}/*.json".format(user_id, new_user_id), shell=True)
        return new_user_id
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error - copy failed: ' + str(e)

@interaction_views.route('/user/<user_id>/view/delete_rerun', methods=['GET'])
def delete_rerun(user_id):
    """
    Deletes all stored results and re-starts the uncurl_analysis pipeline
    on this same dataset, possibly with different parameters.
    """
    if user_id.startswith('test_'):
        return 'Error: unable to delete test results'
    sca = get_sca(user_id)
    try:
        # clear cache
        cache.clear()
        sca.delete_uncurl_results()
        return redirect(url_for('views.state_estimation_result', user_id=user_id))
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)

@interaction_views.route('/user/<user_id>/view/recluster', methods=['POST'])
def recluster(user_id):
    """
    Re-clusters - re-runs the labeling method...
    """
    print('reclustering')
    sca = get_sca(user_id)
    data_form = request.form.copy()
    try:
        sca.relabel(data_form['clustering_method'])
    except Exception as e:
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)
    # clear caches?
    cache.delete_memoized(get_sca_top_genes, user_id)
    cache.delete_memoized(get_sca_top_1vr, user_id)
    cache.delete_memoized(get_sca_pvals, user_id)
    cache.delete_memoized(get_sca_pval_1vr, user_id)
    cache.delete_memoized(update_barplot_result)
    cache.delete_memoized(update_scatterplot_result)
    cache.delete_memoized(heatmap_data)
    cache.delete_memoized(dendrogram_data)
    return 'success'

@interaction_views.route('/user/<user_id>/view/subset', methods=['POST'])
def rerun_uncurl(user_id):
    """
    Given a list of cell ids, this creates a new uncurl data
    preview for the data subset.
    """
    print('RUN SUBSET')
    sca = get_sca(user_id)
    # True if the input are cell_ids, False if the input are cluster ids
    print(request.form['is_cells'])
    is_cells = bool(int(request.form['is_cells']))

    cell_ids = request.form['cell_ids']
    cell_ids = [int(x) for x in cell_ids.split(',')]
    print(len(cell_ids))

    # create new user_id
    new_user_id = str(uuid.uuid4())
    new_path = user_id_to_path(new_user_id, use_secondary=False)


    #path = sca.data_dir
    #shutil.copytree(path, user_id_to_path(new_user_id))
    #import subprocess
    #subprocess.call("sed -i 's/{0}/{1}/g' /tmp/uncurl/{1}/*.json".format(user_id, new_user_id), shell=True)

    # get data subset from sca
    if is_cells:
        data_subset = sca.get_data_subset(cell_ids)
    else:
        # TODO: make sure that the colormap matches...
        colormap = request.form['color_map']
        if colormap is not None and colormap not in ['cluster', 'gene', 'entropy', 'weights']:
            pass
        data_subset = sca.get_clusters_subset(cell_ids)

    new_data_path = os.path.join(new_path, 'data.mtx')
    os.makedirs(new_path)
    scipy.io.mmwrite(new_data_path, data_subset)

    # copy gene names?
    shutil.copy(sca.gene_names_f, new_path)
    # run state_estimation_preproc - gets data summary stats 
    state_estimation_preproc_simple(new_user_id, new_path, new_data_path)
    return new_user_id
