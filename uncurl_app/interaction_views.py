# notes: this backend should be frontend-agnostic.
# I'm thinking of writing the frontend entirely in plotly.js, and not have
# any backend Python rendering components.

import contextlib
import json
import os
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
from .views import state_estimation_preproc

interaction_views = Blueprint('interaction_views', __name__,
        template_folder='templates')

# map of user_id to SCAnalysis objects
# map of gene lists (strings of newline-separated gene names) to enrichr IDs
interaction_views.enrichr_gene_list_ids = {}
# map of tuples (top_genes, gene_set) to enrichr results
interaction_views.enrichr_results = {}

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
    sca = get_sca(user_id)
    return sca.calculate_diffexp(color_track, mode=mode)

@cache.memoize()
def get_sca_gene_names(user_id):
    sca = get_sca(user_id)
    return sca.gene_names

@cache.memoize()
def get_sca_color_track(user_id, color_track):
    sca = get_sca(user_id)
    return sca.get_color_track(color_track)

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

def user_id_to_path(user_id):
    """
    Given a user id, returns the path to the analysis object's base directory.
    """
    if user_id.startswith('test_'):
        user_id = user_id[5:]
        path = os.path.join(current_app.config['TEST_DATA_DIR'], user_id)
        return path
    else:
        return os.path.join(current_app.config['USER_DATA_DIR'], user_id)

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
        color_vals=None):
    """
    Converts data into a form that will be sent as json for building the
    scatterplot. Output should be formatted in a way that can be used by
    Plotly.

    Args:
        dim_red (array): array of shape (2, n)
        labels (array): 1d array of length n
        colorscale (str)
        mode (str): either 'cluster' or 'entropy'
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
    if mode == 'cluster':
        color_values = label_values
        data =  [
            {
                'x': dim_red[0,labels==c].tolist(),
                'y': dim_red[1,labels==c].tolist(),
                'mode': 'markers',
                'name': cluster_names[i],
                'marker': {
                    'size': size,
                    'color': color_values[i],
                    'colorscale': colorscale,
                },
                'text': list(map(str, cell_ids[labels==c].tolist())),
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
                'name': 'cluster ' + str(c),
                'marker': {
                    'size': size,
                    'color': color_values[c],
                    'colorscale': colorscale,
                    'cmin': cmin,
                    'cmax': cmax,
                    'showscale': True if c==0 else False,
                },
                #'text': list(map(str, color_values[c])),
                'text': list(map(str, cell_ids[labels==c].tolist())),
            }
            for c in range(len(set(labels)))
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
            },
    }, cls=SimpleEncoder)


@cache.memoize()
def gene_gene_data(user_id, gene_name_1, gene_name_2, labels, mode='cluster', use_mw=False, color_vals=None, jitter=True, **params):
    """
    Returns a plotly-formated json thing representing a gene-gene
    scatterplot. Colorscheme based on labels.
    """
    size = calc_size(labels)
    gene_1_data = get_gene_data(user_id, gene_name_1, use_mw=use_mw)
    gene_2_data = get_gene_data(user_id, gene_name_2, use_mw=use_mw)
    if jitter:
        gene_1_data += np.random.random(gene_1_data.shape)/2 - 0.25
        gene_2_data += np.random.random(gene_2_data.shape)/2 - 0.25
        gene_1_data[gene_1_data < 0] = 0
        gene_2_data[gene_2_data < 0] = 0
    cell_ids = np.arange(len(labels))
    if color_vals is not None:
        label_values = list(sorted(list(set(labels))))
        colorscale = 'Reds'
        color_values = [color_vals[labels==c] for c in label_values]
        cmin = min(color_vals)
        cmax = max(color_vals)
        data = [
            {
                'x': gene_1_data[labels==c].tolist(),
                'y': gene_2_data[labels==c].tolist(),
                'mode': 'markers',
                'name': 'cluster ' + str(c),
                'marker': {
                    'size': size,
                    'color': color_values[c],
                    'colorscale': colorscale,
                    'cmin': cmin,
                    'cmax': cmax,
                    'showscale': True if c==0 else False,
                },
                #'text': list(map(str, color_values[c])),
                'text': list(map(str, cell_ids[labels==c].tolist())),
            }
            for c in range(len(set(labels)))
        ]
    else:
        data = [
            {
                'x': gene_1_data[labels==c].tolist(),
                'y': gene_2_data[labels==c].tolist(),
                'mode': 'markers',
                'name': 'cluster ' + str(c),
                'marker': {
                    'size': size,
                },
                #'text': list(map(str, color_values[c])),
                'text': list(map(str, cell_ids[labels==c].tolist())),
            }
            for c in range(len(set(labels)))
        ]
    return json.dumps({
            'data': data,
            'layout': {
                'title': 'Cells',
                'xaxis': {'title': gene_name_1},
                'yaxis': {'title': gene_name_2},
                'margin': {'t':30},
                'showlegend': True
            },
    }, cls=SimpleEncoder)


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
    # TODO: mean read count, median read count, mean gene count, median gene count
    return render_template('stats.html',
            user_id=user_id,
            visualization=v,
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
    return render_template('state_estimation_static.html', user_id=user_id,
            test_or_user=test_or_user,
            data_user_id=data_user_id,
            gene_names=get_sca_gene_names(user_id),
            gene_sets=enrichr_api.ENRICHR_LIBRARIES,
            color_tracks=sca.get_color_track_names())


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
    gene_names = get_sca_gene_names(user_id)
    if 'selected_gene' in data_form:
        selected_gene = data_form['selected_gene']
    if len(selected_gene.strip()) > 0:
        gns = set(gene_names)
        selected_gene_names = [x.strip(', ') for x in selected_gene.split()]
        selected_gene_names = [x for x in selected_gene_names if x in gns]
    if top_or_bulk == 'top':
        if len(selected_gene.strip()) > 0:
            top_genes = get_sca_top_genes(user_id)[int(input_value)]
            selected_top_genes = [x for x in top_genes if gene_names[x[0]] in set(selected_gene_names)]
        else:
            selected_top_genes = get_sca_top_genes(user_id)[int(input_value)][:num_genes]
        selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_value,
                x_label='c-score',
                title='Top genes for cluster {0}'.format(input_value),
                )
    elif top_or_bulk == 'pval':
        if len(selected_gene.strip()) > 0:
            top_genes = get_sca_pvals(user_id)[int(input_value)]
            selected_top_genes = [x for x in top_genes if gene_names[x[0]] in set(selected_gene_names)]
        else:
            selected_top_genes = get_sca_pvals(user_id)[input_value][:num_genes]
        selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_value,
                title='Top genes for cluster {0}'.format(input_value),
                x_label='p-value of c-score')
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
    elif top_or_bulk == 'top_1_vs_rest':
        if len(selected_gene.strip()) > 0:
            top_genes = get_sca_top_1vr(user_id)[int(input_value)]
            selected_top_genes = [x for x in top_genes if gene_names[x[0]] in set(selected_gene_names)]
        else:
            selected_top_genes = get_sca_top_1vr(user_id)[int(input_value)][:num_genes]
        selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_value,
                title='Top genes for cluster {0}'.format(input_value),
                x_label='Fold change (1 vs rest)')
    elif top_or_bulk == 'pval_1_vs_rest':
        if len(selected_gene.strip()) > 0:
            top_genes = get_sca_pval_1vr(user_id)[input_value]
            selected_top_genes = [x for x in top_genes if gene_names[x[0]] in set(selected_gene_names)]
        else:
            selected_top_genes = get_sca_pval_1vr(user_id)[input_value][:num_genes]
            selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_value,
                title='Top genes for cluster {0}'.format(input_value),
                x_label='p-value of log-fold change (1 vs rest)')
    elif top_or_bulk == 'selected_color' or top_or_bulk == 'selected_color_pval':
        colormap = str(data_form['cell_color'])
        print('getting diffexp for selected color')
        print('barplot cell_color: ', colormap)
        color_track, is_discrete = get_sca_color_track(user_id, colormap)
        # this is some kind of lock to prevent this from running more than once?
        # this is a serious hack... there should be a better way of doing this
        lockfile_name = os.path.join(sca.data_dir, colormap + '_writing_diffexp')
        with lockfile_context(lockfile_name) as _lock:
            selected_diffexp, selected_pvals = get_sca_top_genes_custom(user_id, colormap)
        x_label = 'Fold change (1 vs rest)'
        if top_or_bulk == 'selected_color_pval':
            selected_diffexp = selected_pvals
            x_label = 'p-value of fold change (1 vs rest)'
        _, color_map = color_track_map(color_track)
        input_label = color_map[input_value]
        if len(selected_gene.strip()) > 0:
            top_genes = selected_diffexp[input_label]
            selected_top_genes = [x for x in top_genes if gene_names[int(x[0])] in set(selected_gene_names)]
        else:
            selected_top_genes = selected_diffexp[input_label][:num_genes]
        selected_gene_names = [gene_names[int(x[0])] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_label,
                title='Top genes for label {0}'.format(input_label),
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
                desc = 'p-values of ratios'
            # generate barplot
            if len(selected_gene.strip()) > 0:
                genes, values = array_to_top_genes(data, cluster1, cluster2, is_pvals=(top_or_bulk=='top_pairwise'), num_genes=1000000)
                gene_data = list(zip(genes, values))
                gene_data = [x for x in gene_data if gene_names[x[0]] in set(selected_gene_names)]
            else:
                gene_data = list(zip(genes, values))
            selected_gene_names = [gene_names[x[0]] for x in gene_data]
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
                desc = 'p-values of ratios'
            if len(selected_gene.strip()) > 0:
                if top_or_bulk == 'top_pairwise':
                    genes, values = array_to_top_genes(selected_diffexp, cluster1, cluster2, is_pvals=(top_or_bulk=='top_pairwise'), num_genes=1000000)
                else:
                    genes, values = array_to_top_genes(selected_pvals, cluster1, cluster2, is_pvals=(top_or_bulk=='top_pairwise'), num_genes=1000000)
                gene_data = list(zip(genes, values))
                gene_data = [x for x in gene_data if gene_names[x[0]] in set(selected_gene_names)]
            else:
                gene_data = list(zip(genes, values))
            selected_gene_names = [gene_names[x[0]] for x in gene_data]
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
    elif top_or_bulk == 'custom_selection':
        # TODO: custom selection of cells vs rest
        pass
    else:
        return 'Error: '

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
    elif plot_type == 'Gene-gene':
        if 'gene_name_1' not in data_form or 'gene_name_2' not in data_form:
            return None
        use_mw = False
        if 'use_mw' in data_form:
            use_mw = bool(int(data_form['use_mw']))
        gene_name = ''
        if gene_name in data_form:
            gene_name = data_form['gene_name']
        gene_name_1 = data_form['gene_name_1']
        gene_name_2 = data_form['gene_name_2']
        print('gene-gene plot - gene names:', gene_name_1, gene_name_2)
        color_vals = None
        if cell_color_value == 'gene':
            color_vals = get_gene_data(user_id, gene_name)
        elif cell_color_value == 'entropy':
            color_vals = sca.entropy
        return gene_gene_data(user_id, gene_name_1, gene_name_2, sca.labels, use_mw=use_mw, mode=cell_color_value, color_vals=color_vals)
    else:
        dim_red = None
        if plot_type == 'Cells':
            dim_red = sca.dim_red
        elif plot_type == 'Baseline':
            dim_red = sca.baseline_vis
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
        else:
            # try to get color track
            color_track, is_discrete = get_sca_color_track(user_id, cell_color_value)
            if color_track is None:
                return scatterplot_data(dim_red, sca.labels)
            else:
                if is_discrete:
                    return scatterplot_data(dim_red, color_track)
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
        - expressed_genes (list of strings)
        - expressed_genes_values (list of numbers)
        - cluster (int)
        - values for all of the uploaded color maps
    """
    # TODO: is this really necessary?
    print('cell_info: ', request.form)
    selected_cells = request.form['selected_cells']
    selected_cells = selected_cells.split(',')
    selected_cells = [int(x) for x in selected_cells]
    selected_clusters = request.form['selected_clusters']
    selected_clusters = selected_clusters.split(',')
    selected_clusters = [int(x) for x in selected_clusters]
    return cell_info_result(user_id, selected_cells)

@cache.memoize()
def cell_info_result(user_id, selected_cells):
    # get read count + gene count for all cells
    sca = get_sca(user_id)
    read_counts = []
    gene_counts = []
    for cell in selected_cells:
        cell_data = sca.data_sampled_all_genes[:, cell]
        gene_counts.append(cell_data.count_nonzero())
        read_counts.append(cell_data.sum())
    return json.dumps({'gene_counts': gene_counts,
        'read_counts': read_counts}, cls=SimpleEncoder)


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
        gene_list = top_genes.strip().split()
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
            gene_list = top_genes.strip().split()
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
    top_genes = [x.strip().upper() for x in request.form['top_genes'].split('\n')]
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
    def pmid_to_link(pmid):
        return '<a href="https://www.ncbi.nlm.nih.gov/pubmed/{0}">{0}</a>'.format(pmid)
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
    top_genes = [x.strip().upper() for x in request.form['top_genes'].split('\n')]
    print('update_cellmesh:', top_genes)
    # gene_set is a string.
    test_type = request.form['mesh_test_type']
    return update_cellmesh_result(user_id, top_genes, test_type)

@cache.memoize()
def update_cellmesh_result(user_id, top_genes, test):
    """
    Gets the CellMesh result for a set of genes.

    Args:
        user_id (str)
        top_genes (list of strings representing genes)
        test (str, currently only 'hypergeom')
    """
    def pmid_to_link(pmid):
        return '<a href="https://www.ncbi.nlm.nih.gov/pubmed/{0}">{0}</a>'.format(pmid)
    import cellmesh
    result = []
    if test == 'hypergeom':
        result = cellmesh.hypergeometric_test(top_genes, return_header=True)
    elif test == 'norm_hypergeom':
        result = cellmesh.normed_hypergeometric_test(top_genes, return_header=True)
    cell_types = [result[0]]
    for i in range(1, min(20, len(result))):
        ri = result[i]
        gene_pmids = []
        genes = ri[3]
        for g in genes:
            gene_pmids.append('{0}: {1}'.format(g, ', '.join(pmid_to_link(x) for x in ri[4][g])))
        cell_types.append((ri[0], ri[1], ri[2], ', '.join(ri[3]), ', '.join(gene_pmids)))
    return json.dumps(cell_types, cls=SimpleEncoder)

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
            data = sca.data_sampled_all_genes
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
        #print('deleting user_id from cache')
        # TODO: this currently doesn't work, since keys are hashed.
        #clear_cache_user_id(user_id)
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
    from werkzeug import secure_filename
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
    current_id = 1
    has_id = True
    all_criteria = []
    while has_id:
        selection_type = json_dict['selection_type-'+str(current_id)]
        comparison = json_dict['selection_comparison-'+str(current_id)]
        target = json_dict['selection_target-'+str(current_id)]
        and_or = json_dict['selection_and_or-'+str(current_id)]
        criterion = custom_cell_selection.LabelCriterion(selection_type=selection_type, comparison=comparison, target=target, and_or=and_or)
        all_criteria.append(criterion)
        current_id += 1
        if 'selection_type-'+str(current_id) not in json_dict:
            has_id = False
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
    if 'criteria' in data_form:
        # load criteria from json
        criteria = load_criteria_from_dict(json.loads(data_form['criteria']))
        sca.update_custom_color_track_label(colormap_name, label_name, criteria)
        # clear cache for scatterplot results
        print('deleting cached results...')
        cache.delete_memoized(update_barplot_result)
        cache.delete_memoized(update_scatterplot_result)
        cache.delete_memoized(get_sca_color_track)
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
        shutil.copytree(path, user_id_to_path(user_id))
        return new_user_id
    except:
        return 'Error: copy failed.'

@interaction_views.route('/user/<user_id>/view/delete_rerun', methods=['GET'])
def delete_rerun(user_id):
    """
    Deletes all stored results and re-starts the uncurl_analysis pipeline
    on this same dataset, possibly with different parameters.
    """
    if user_id.startswith('test_'):
        return 'Error: unable to delete test results'
    sca = get_sca(user_id)
    path = sca.data_dir
    try:
        os.remove(os.path.join(path, 'sc_analysis.json'))
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

    # TODO: cluster ids is not currently implemented
    cell_ids = request.form['cell_ids']
    cell_ids = [int(x) for x in cell_ids.split(',')]
    print(len(cell_ids))

    # create new user_id
    new_user_id = str(uuid.uuid4())
    new_path = user_id_to_path(new_user_id)

    # get data subset from sca
    if is_cells:
        data_subset = sca.get_data_subset(cell_ids)
    else:
        data_subset = sca.get_clusters_subset(cell_ids)

    new_data_path = os.path.join(new_path, 'data.mtx')
    os.makedirs(new_path)
    scipy.io.mmwrite(new_data_path, data_subset)

    # run state_estimation_preproc - gets data summary stats 
    state_estimation_preproc(new_user_id, new_path, new_data_path, '')
    return new_user_id
