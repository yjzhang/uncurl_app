# notes: this backend should be frontend-agnostic.
# I'm thinking of writing the frontend entirely in plotly.js, and not have
# any backend Python rendering components.

# TODO: have some sort of caching (redis for deployment, 
# simple cache for testing) - have some mechanism to switch between the two

import json
import os

import numpy as np
from flask import request, render_template
from uncurl_analysis import enrichr_api, sc_analysis

from app import app
from cache import cache

# map of user_id to SCAnalysis objects
app.sc_analysis_dict = {}
# map of gene lists (strings of newline-separated gene names) to enrichr IDs
app.enrichr_gene_list_ids = {}
# map of tuples (top_genes, gene_set) to enrichr results
app.enrichr_results = {}

def get_sca(user_id):
    if user_id in app.sc_analysis_dict:
        return app.sc_analysis_dict[user_id]
    else:
        path = user_id_to_path(user_id)
        sca = sc_analysis.SCAnalysis(path)
        sca = sca.load_params_from_folder()
        app.sc_analysis_dict[user_id] = sca
        return sca

def user_id_to_path(user_id):
    if user_id.startswith('test_'):
        user_id = user_id[5:]
        path = os.path.join(app.config['TEST_DATA_DIR'], user_id)
        return path
    else:
        return os.path.join(app.config['USER_DATA_DIR'], user_id)

def barplot_data(gene_values, gene_names, cluster_name, x_label,
        title=None):
    """
    Converts data for top genes into a json for building the
    bar plot. Output should be formatted in a way that can be plugged into
    Plotly.

    Args:
        selected_top_genes (list): list of tuples (gene_id, gene_value)
        selected_gene_names (list): list of gene names corresponding to
            the genes in selected_top_genes.
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
                'x': [x[1] for x in gene_values],
                'y': gene_names,
                'orientation': 'h',
                'type': 'bar',
            }],
            'layout': {
                'title': title,
                'xaxis': {'title': x_label},
                'margin':{'t':40},
            },
        })

def scatterplot_data(dim_red, labels, colorscale='Portland', mode='cluster',
        gene_expression_list=None, entropy=None):
    """
    Converts data into a form that will be sent as json for building the
    scatterplot. Output should be formatted in a way that can be used by
    Plotly.

    Args:
        dim_red (array): array of shape (2, n)
        labels (array): 1d array of length n

    """
    if mode == 'cluster':
        color_values = list(range(len(set(labels))))
    elif mode == 'entropy':
        color_values = [entropy[labels==c].tolist() for c in set(labels)]
    # TODO: add a colorbar for entropy mode.
    # also, use a different view.
    # have size depend on data shape
    size = 10
    if len(labels) < 50:
        size = 20
    if len(labels) > 2000:
        size = 5
    elif len(labels) > 10000:
        size = 1
    return json.dumps({
            'data': [
                {
                    'x': dim_red[0,labels==c].tolist(),
                    'y': dim_red[1,labels==c].tolist(),
                    'mode': 'markers',
                    'name': 'cluster ' + str(c),
                    'marker': {
                        'size': size,
                        'color': color_values[c],
                        'colorscale': colorscale,
                    },
                }
                for c in range(len(set(labels)))
            ],
            'layout': {
                'title':'Cells',
                'xaxis':{'title': 'dim1'},
                'yaxis':{'title': 'dim2'},
                'margin':{'t':30},
            },
    })



@app.route('/user/<user_id>/view')
def view_plots(user_id):
    test_or_user = 'user'
    data_user_id = user_id
    if user_id.startswith('test_'):
        test_or_user = 'test'
        data_user_id = user_id[5:]
    path = user_id_to_path(user_id)
    sca = get_sca(user_id)
    return render_template('state_estimation_static.html', user_id=user_id,
            test_or_user=test_or_user,
            data_user_id=data_user_id,
            gene_sets=enrichr_api.ENRICHR_LIBRARIES)


@app.route('/user/<user_id>/view/update_barplot', methods=['GET', 'POST'])
def update_barplot(user_id):
    """
    Updates barplot data.

    Returns:
        json corresponding to Plotly barplot
    """
    top_or_bulk = str(request.form['top_or_bulk'])
    input_value = int(request.form['input_value'])
    num_genes = int(request.form['num_genes'])
    return update_barplot_result(user_id, top_or_bulk, input_value, num_genes)

# TODO: cache this function
@cache.memoize()
def update_barplot_result(user_id, top_or_bulk, input_value, num_genes):
    sca = get_sca(user_id)
    if top_or_bulk == 'top':
        selected_top_genes = sca.top_genes[int(input_value)][:num_genes]
        selected_gene_names = [sca.gene_names[x[0]] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_value,
                x_label='c-score',
                title='Top genes for cluster {0}'.format(input_value),
                )
    elif top_or_bulk == 'pval':
        selected_top_genes = sca.pvals[input_value][:num_genes]
        selected_gene_names = [sca.gene_names[x[0]] for x in selected_top_genes]
        return barplot_data(selected_top_genes,
                selected_gene_names, input_value,
                title='Top genes for cluster {0}'.format(input_value),
                x_label='p-value of c-score')
    elif top_or_bulk == 'sep':
        # show separation score
        sep_scores = sca.separation_scores[int(input_value)]
        cluster_names = map(lambda x: 'cluster ' + str(x),
                list(range(sca.separation_scores.shape[0])))
        sep_scores = [(c, x) for c, x in zip(cluster_names, sep_scores)]
        return barplot_data(sep_scores,
                cluster_names, input_value,
                x_label='separation score',
                title='Inter-cluster separations for cluster {0}'.format(input_value))
    else:
        # TODO: show bulk correlations
        pass

@app.route('/user/<user_id>/view/update_scatterplot', methods=['GET', 'POST'])
def update_scatterplot(user_id):
    """
    Updates the scatterplot view.

    Returns:
        {data: [{'cluster': c, 'x': x, 'y'; y}...],
         labels: [c1, c2, ...],
         colorscale: 'Viridis'}
    """
    plot_type = request.form['scatter_type']
    cell_color_value = request.form['cell_color']
    return update_scatterplot_result(user_id, plot_type, cell_color_value)

# TODO: cache this function
@cache.memoize()
def update_scatterplot_result(user_id, plot_type, cell_color_value):
    sca = get_sca(user_id)
    if plot_type == 'Means':
        labels = np.arange(sca.mds_means.shape[1])
        return scatterplot_data(sca.mds_means,
                labels)
    elif plot_type == 'Cells':
        if cell_color_value == 'entropy':
            return scatterplot_data(sca.dim_red, sca.labels,
                    colorscale='Viridis',
                    mode='entropy', entropy=sca.entropy)
        else:
            return scatterplot_data(sca.dim_red, sca.labels)
    elif plot_type == 'Baseline':
        if cell_color_value == 'entropy':
            return scatterplot_data(sca.baseline_vis, sca.labels,
                    colorscale='Viridis',
                    mode='entropy', entropy=sca.entropy)
        else:
            return scatterplot_data(sca.baseline_vis, sca.labels)

@app.route('/user/<user_id>/view/update_enrichr', methods=['GET', 'POST'])
def update_enrichr(user_id):
    # top_genes is a newline-separated string, representing the gene list.
    top_genes = request.form['top_genes']
    # gene_set is a string.
    gene_set = request.form['gene_set']
    return update_enrichr_result(user_id, top_genes, gene_set)

# TODO: cache this function??? but we don't want to cache timeouts
def update_enrichr_result(user_id, top_genes, gene_set):
    user_list_id = 0
    if top_genes not in app.enrichr_gene_list_ids:
        gene_list = top_genes.strip().split()
        user_list_id = enrichr_api.enrichr_add_list(gene_list)
        if user_list_id == 'timeout':
            return 'Error: Enrichr query timed out'
        app.enrichr_gene_list_ids[top_genes] = user_list_id
    else:
        user_list_id = app.enrichr_gene_list_ids[top_genes]
    results = []
    if (top_genes, gene_set) in app.enrichr_results:
        results = app.enrichr_results[(top_genes, gene_set)]
    else:
        try:
            results = enrichr_api.enrichr_query(user_list_id, gene_set)
            if results == 'timeout':
                return 'Error: Enrichr query timed out'
            app.enrichr_results[(top_genes, gene_set)] = results[:10]
        except:
            gene_list = top_genes.strip().split()
            user_list_id = enrichr_api.enrichr_add_list(gene_list)
            if user_list_id == 'timeout':
                return 'Error: Enrichr query timed out'
            app.enrichr_gene_list_ids[top_genes] = user_list_id
            results = enrichr_api.enrichr_query(user_list_id, gene_set)
            if results == 'timeout':
                return 'Error: Enrichr query timed out'
            app.enrichr_results[(top_genes, gene_set)] = results[:10]
    # only take top 10 results (maybe have this value be variable?)
    results = results[:10]
    app.last_enrichr_results = [['gene set name',
                                 'p-value',
                                 'z-score',
                                 'combined score']] + \
            [[r[1], r[2], r[3], r[4]] for r in results]
    return json.dumps(app.last_enrichr_results)

@app.route('/user/<user_id>/view/split_or_merge_cluster', methods=['POST'])
def split_or_merge_cluster(user_id):
    sca = get_sca(user_id)
    selected_clusters = []
    cluster_counts = []
    return
    for point in selected_points['points']:
        cluster = point['curveNumber']
        selected_clusters.append(cluster)
    selected_clusters = list(set(selected_clusters))
    for cluster in selected_clusters:
        cluster_counts.append((sca.labels == cluster).sum())
    # split clusters - TODO: have some kind of progress bar?
    if n_click_split > app.split_clicks:
        if test_or_user == 'test':
            return 'Test datasets cannot be modified.'
        return 'Splitting selected cluster: ' + str(selected_clusters[0]) + '...'
    # merge clusters
    elif n_click_merge > app.merge_clicks:
        if test_or_user == 'test':
            return 'Test datasets cannot be modified.'
        return 'Merging selected clusters: ' + ' '.join(map(str, selected_clusters)) + '...'
    cache.clear()
    return 'Selected clusters: ' + ' '.join(map(lambda x: '{0} ({1} cells)'.format(x[0], x[1]), zip(selected_clusters, cluster_counts)))


def update_all_views(selected_points, n_click_split, n_click_merge):
    """
    """
    if test_or_user == 'test':
        raise Exception('test datasets cannot be changed')
    selected_clusters = []
    for point in selected_points['points']:
        cluster = point['curveNumber']
        selected_clusters.append(cluster)
    selected_clusters = list(set(selected_clusters))
    # split clusters
    if n_click_split > app.split_clicks:
        app.split_clicks = n_click_split
        generate_analysis.generate_analysis_resubmit(sca,
                'split', selected_clusters)
        initialize(app, data_dir, permalink, user_id, test_or_user)
        return generate_cluster_view(sca.mds_means)
    # merge clusters
    elif n_click_merge > app.merge_clicks:
        app.merge_clicks = n_click_merge
        generate_analysis.generate_analysis_resubmit(sca,
                'merge', selected_clusters)
        initialize(app, data_dir, permalink, user_id, test_or_user)
        return generate_cluster_view(sca.dim_red)
    else:
        raise Exception('placeholder')
