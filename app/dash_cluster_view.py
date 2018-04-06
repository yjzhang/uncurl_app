import json
import os, os.path

from flask import url_for

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np
from uncurl_analysis import enrichr_api
import generate_analysis

def create_means_figure(dim_red, colorscale='Portland'):
    """
    create a figure for displaying cluster means.

    Args:
        dim_red (array): array of shape (2, k)
    """
    return {
                'data': [
                    go.Scatter(
                        x=dim_red[0,c:c+1],
                        y=dim_red[1,c:c+1],
                        mode='markers',
                        name='cluster ' + str(c),
                        marker={
                            'size': 20,
                            'color': c,
                            'colorscale': colorscale,
                        },
                    )
                    for c in range(dim_red.shape[1])
                ],
                'layout': go.Layout(
                    title='Cluster means',
                    xaxis={'title': 'dim1'},
                    yaxis={'title': 'dim2'},
                    margin={'t':30},
                ),
        }

def create_cells_figure(dim_red, labels, colorscale='Portland',
        mode='cluster',
        gene_expression_list=None, entropy=None):
    """
    create a figure for displaying cells

    Args:
        dim_red (array): array of shape (2, n)
        labels (array): 1d array of length n
    """
    if mode == 'cluster':
        color_values = list(range(len(set(labels))))
    elif mode == 'entropy':
        color_values = [entropy[labels==c] for c in set(labels)]
    # TODO: add a colorbar for entropy mode.
    # also, use a different view.
    # have size depend on data shape
    size = 10
    if size > 2000:
        size = 5
    elif size > 10000:
        size = 1
    return {
                'data': [
                    go.Scatter(
                        x=dim_red[0,labels==c],
                        y=dim_red[1,labels==c],
                        mode='markers',
                        name='cluster ' + str(c),
                        marker={
                            'size': size,
                            'color': color_values[c],
                            'colorscale': colorscale,
                        },
                    )
                    for c in range(len(set(labels)))
                ],
                'layout': go.Layout(
                    title='Cells',
                    xaxis={'title': 'dim1'},
                    yaxis={'title': 'dim2'},
                    margin={'t':30},
                ),
        }

def create_top_genes_figure(selected_top_genes, selected_gene_names,
        cluster_name, x_label='c-score'):
    """
    Creates a figure for displaying top genes

    Args:
        selected_top_genes (list): list of tuples (gene_id, gene_value)
        selected_gene_names (list): list of gene names corresponding to
                the genes in selected_top_genes.
        cluster_name: name of the cluster from which the top genes are drawn.
    """
    if selected_top_genes is None:
        selected_top_genes = [(1,1),(2,2),(3,3)]
    if selected_gene_names is None:
        selected_gene_names = ['placeholder 1', 'placeholder 2', 'placeholder 3']
    return {
                'data': [
                    go.Bar(
                        y=selected_gene_names,
                        x=[x[1] for x in selected_top_genes],
                        orientation='h',
                    )
                ],
                'layout': go.Layout(
                    title='Top genes for cluster {0}'.format(cluster_name),
                    xaxis={'title': x_label},
                    yaxis={'title': 'genes'},
                ),
        }


def create_bulk_correlation_figure(correlations, bulk_names, cluster_name,
        n_datasets=10):
    """
    Creates a figure for displaying correlations with bulk datasets

    Args:
        correlations (list): list of correlation/similarity values
                (higher = more similarity)
        bulk_names (list): names of the bulk datasets.
        cluster_name (str)
        n_datasets (int): number of bulk datasets.
    """
    # TODO: display top n bulk correlations
    return {
                'data': [
                    go.Bar(
                        y=bulk_names,
                        x=[x[1] for x in correlations],
                        orientation='h',
                    )
                ],
                'layout': go.Layout(
                    title='Top bulk datasets for cluster {0}'.format(cluster_name),
                    xaxis={'title': 'similarity'},
                    yaxis={'title': 'dataset'},
                ),
        }

def generate_cluster_view(dim_red, n_genes=10, gene_names_list=None):
    """
    Generates a cluster view: MDS plot of means on the left, updated bar plot
    of top genes with c-scores on the right.

    Args:
        dim_red (array): 2d array of 2 x cells
    """
    colorscale = 'Portland'
    #if dim_red.shape[1] > 10:
    #    colorscale = 'Jet'
    return html.Div([
        dcc.Location(id='url', refresh=False),
        # view 1: graph plot
        html.Div([
            dcc.RadioItems(
                id='means-or-cell',
                options=[{'label': 'Cluster Means', 'value': 'Means'},
                         {'label': 'Processed Cells', 'value': 'Cells'},
                         {'label': 'Unprocessed Cells', 'value': 'Baseline'}],
                value='Means',
                labelStyle={'display': 'inline-block'},
            ),
            dcc.Graph(id='means',
                figure=create_means_figure(dim_red, colorscale),
                style={'width': 700, 'margin-top': 5}
            ),
            # TODO: dropdown to select color scheme  - color by entropy,
            # color by selected gene expression, etc.
            html.Div([
                dcc.Dropdown(id='cell-color',
                    options=[{'label': 'Color: cluster', 'value': 'cluster'},
                             {'label': 'Color: entropy', 'value': 'entropy'}],
                    value='cluster',)],
                style={'width': 400, 'margin-top': 5}),
            # split/merge clusters
            html.Div([
                html.Button(children='Split selected cluster', id='split'),
                html.Button(children='Merge selected clusters', id='merge'),
                html.Div(id='update-area'),
                ],
            ),
            ],
            style={'display': 'inline-block', 'width': 750, 'float':'left'}),
        # view 2: top genes
        html.Div([
            html.Div([
                dcc.Dropdown(
                id='top-or-bulk',
                options=[{'label': 'Display top genes (c-score)', 'value': 'top'},
                         {'label': 'Display top genes (p-value)', 'value': 'pval'},
                         {'label': 'Display bulk correlations', 'value': 'bulk'}],
                value='top',
                #labelStyle={'display': 'inline-block'},
                )],
                style={'margin-top': -25, 'width': 300}),
            html.Div(['Number of top genes: ', dcc.Input(
                id='num-genes',
                value=10,
                type='number',
            )], style={'margin-top': 4, 'margin-bottom': 0}),
            dcc.Graph(id='top-genes',
                figure=create_top_genes_figure(None, None, '0'),
                style={'width': 550}),
            dcc.Textarea(id='top-genes-view', value='Top genes go here',
                cols='50', rows='10',
                readOnly='true'),
            ],
            style={'display': 'inline-block', 'width': 550}),
        # view 3: Enrichr API
        html.Div([
            'Call Enrichr on current top genes',
            html.Div(
                dcc.Dropdown(id='gene-set-library',
                    options=[{'label': x, 'value': x}
                        for x in enrichr_api.ENRICHR_LIBRARIES]
                ),
                style={'width':500},
            ),
            html.Button(children='Submit',
                id='enrichr-submit'),
            html.Table(id='enrichr-result'),
            ],
            style={
                   'margin-top': 50,
                   'float': 'left',
                   'display': 'inline-block'
                   }),
        ],
        style={'width': '100%', 'display':'inline-block',
            'margin-top': 10,
            'margin-left': 60})


def initialize(app, data_dir=None, permalink='test', user_id='test',
        test_or_user='test', uncurl_args={'max_iters': 20, 'threads':2}):
    """
    This function sets app.layout using a directory containing uncurl results.
    """
    app.uncurl_args = uncurl_args
    # hacking the internals of dash - remove all callbacks...
    app.callback_map = {}
    app.initialized = True
    # app variables: scatter_mode is either 'means' or 'cells'
    app.scatter_mode = 'means'
    # 'bar_node' is either 'top' or 'bulk'
    # top indicates that top genes should be displayed. 'bulk' indicates
    # that bulk correlations should be displayed.
    app.bar_mode = 'top'
    # map from strings of genes (joined by '\n') to enrichr IDs
    app.enrichr_gene_list_ids = {}
    # map of tuple (string, enrichr gene set library) to enrichr results
    app.enrichr_results = {}
    # really terrible hack to get enrichr push to only trigger on click,
    # not on gene list update...
    app.n_clicks_enrichr = 0
    app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

    #M = None
    app.labels = None
    app.mds_means = np.array([[1,2,3,4],[1,2,3,4]])
    app.mds_data = None
    app.top_genes = {'0': [(0,100),(1,50),(2,40)], '1': [(0,50),(1,45),(2,30)]}
    app.p_values = {'0': [(0,0.0),(1,0.1),(2,0.2)], '1': [(0,0.0),(1,0.1),(2,0.2)]}
    app.gene_names = None
    app.entropy = None
    #print('initialize ' + data_dir)
    if data_dir != None:
        app.labels = np.loadtxt(os.path.join(data_dir, 'labels.txt')).astype(int)
        app.mds_means = np.loadtxt(os.path.join(data_dir, 'mds_means.txt'))
        app.mds_data = np.loadtxt(os.path.join(data_dir, 'mds_data.txt'))
        # load genes
        with open(os.path.join(data_dir, 'top_genes.txt')) as f:
            app.top_genes = json.load(f)
        # load pvals
        try:
            with open(os.path.join(data_dir, 'gene_pvals.txt')) as f:
                app.p_values = json.load(f)
        except:
            app.p_values = app.top_genes
        # load gene names
        try:
            app.gene_names = np.loadtxt(os.path.join(data_dir, 'gene_names.txt'), dtype=str)
        except:
            M = np.loadtxt(os.path.join(data_dir, 'm.txt'))
            app.gene_names = np.array(['gene ' + str(x) for x in range(M.shape[0])])
            del M
        # load entropy
        try:
            app.entropy = np.loadtxt(os.path.join(data_dir, 'entropy.txt'))
        except:
            W = np.loadtxt(os.path.join(data_dir, 'w.txt'))
            app.entropy = -(W*np.log2(W)).sum(0)
            del W
        # load baseline visualization
        try:
            app.baseline_vis = np.loadtxt(os.path.join(data_dir, 'baseline_vis.txt'))
        except:
            app.baseline_vis = app.mds_data

    # generate layout
    #app.layout.children[1].children = generate_cluster_view(mds_means)
    app.layout = html.Div([
        html.Div([
            html.Div(html.A('permalink: ' + permalink, href=permalink),
                style={'width': 500}),
            # TODO: add links to data downloads
            html.Div(html.A('Download M', href=url_for('state_estimation_file',
                x=test_or_user,
                user_id=user_id, filename='m.txt')),
                style={'width': 200}),
            html.Div(html.A('Download W', href=url_for('state_estimation_file',
                x=test_or_user,
                user_id=user_id, filename='w.txt')),
                style={'width': 200})
            ],
            style={'display':'inline-block',
                'margin-left': 60}
        ),
        html.Div([generate_cluster_view(app.mds_means)],
            id='cluster-view')
    ])

    # create callback for clicking on clusters
    @app.callback(
            Output(component_id='top-genes', component_property='figure'),
            [Input(component_id='top-or-bulk', component_property='value'),
             Input(component_id='means', component_property='clickData'),
             Input(component_id='num-genes', component_property='value')]
    )
    def update_barplot(top_or_bulk, input_value, num_genes):
        #print('input value:')
        #print(input_value)
        num_genes = int(num_genes)
        if input_value is None:
            input_value = '0'
        else:
            input_value = str(input_value['points'][0]['curveNumber'])
        if top_or_bulk == 'top':
            selected_top_genes = app.top_genes[input_value][:num_genes]
            selected_gene_names = [app.gene_names[x[0]] for x in selected_top_genes]
            return create_top_genes_figure(selected_top_genes,
                    selected_gene_names, input_value)
        elif top_or_bulk == 'pval':
            selected_top_genes = app.p_values[input_value][:num_genes]
            selected_gene_names = [app.gene_names[x[0]] for x in selected_top_genes]
            return create_top_genes_figure(selected_top_genes,
                    selected_gene_names, input_value,
                    x_label='p-value of c-score')
        else:
            # TODO: show bulk correlations
            pass

    # create callback for switching to cell view
    @app.callback(
            Output(component_id='means', component_property='figure'),
            [Input(component_id='means-or-cell', component_property='value'),
             Input(component_id='cell-color', component_property='value')]
    )
    def update_scatterplot(input_value, cell_color_value):
        if input_value == 'Means':
            app.scatter_mode = 'means'
            return create_means_figure(app.mds_means)
        elif input_value == 'Cells':
            app.scatter_mode = 'cells'
            if cell_color_value == 'entropy':
                return create_cells_figure(app.mds_data, app.labels,
                        colorscale='Viridis',
                        mode='entropy', entropy=app.entropy)
            else:
                return create_cells_figure(app.mds_data, app.labels)
        elif input_value == 'Baseline':
            app.scatter_mode = 'baseline'
            if cell_color_value == 'entropy':
                return create_cells_figure(app.baseline_vis, app.labels,
                        colorscale='Viridis',
                        mode='entropy', entropy=app.entropy)
            else:
                return create_cells_figure(app.baseline_vis, app.labels)

    # create callback for top genes list
    @app.callback(
            Output(component_id='top-genes-view', component_property='value'),
            [Input(component_id='top-or-bulk', component_property='value'),
             Input(component_id='means', component_property='clickData'),
             Input(component_id='num-genes', component_property='value')]
    )
    def update_genes_list(top_or_bulk, input_value, num_genes):
        num_genes = int(num_genes)
        if input_value is None:
            input_value = '0'
        else:
            input_value = str(input_value['points'][0]['curveNumber'])
        if top_or_bulk == 'top':
            selected_top_genes = app.top_genes[input_value][:num_genes]
            selected_gene_names = [app.gene_names[x[0]] for x in selected_top_genes]
            return '\n'.join(selected_gene_names)
        elif top_or_bulk == 'pval':
            selected_top_genes = app.p_values[input_value][:num_genes]
            selected_gene_names = [app.gene_names[x[0]] for x in selected_top_genes]
            return '\n'.join(selected_gene_names)
        else:
            # TODO: show bulk correlations
            pass

    app.last_enrichr_results = None
    # create callback for enrichr API
    @app.callback(
            Output(component_id='enrichr-result', component_property='children'),
            [Input(component_id='enrichr-submit', component_property='n_clicks'),
             Input(component_id='top-genes-view', component_property='value'),
             Input(component_id='gene-set-library', component_property='value')]
    )
    def update_enrichr(n_clicks, top_genes_value, gene_set):
        if n_clicks > app.n_clicks_enrichr:
            app.n_clicks_enrichr = n_clicks
            user_list_id = 0
            if top_genes_value not in app.enrichr_gene_list_ids:
                gene_list = top_genes_value.strip().split()
                user_list_id = enrichr_api.enrichr_add_list(gene_list)
                app.enrichr_gene_list_ids[top_genes_value] = user_list_id
            else:
                user_list_id = app.enrichr_gene_list_ids[top_genes_value]
            results = []
            if (top_genes_value, gene_set) in app.enrichr_results:
                results = app.enrichr_results[(top_genes_value, gene_set)]
            else:
                try:
                    results = enrichr_api.enrichr_query(user_list_id, gene_set)
                    app.enrichr_results[(top_genes_value, gene_set)] = results[:10]
                except:
                    gene_list = top_genes_value.strip().split()
                    user_list_id = enrichr_api.enrichr_add_list(gene_list)
                    app.enrichr_gene_list_ids[top_genes_value] = user_list_id
                    results = enrichr_api.enrichr_query(user_list_id, gene_set)
                    app.enrichr_results[(top_genes_value, gene_set)] = results[:10]
            # only take top 10 results (maybe have this value be variable?)
            results = results[:10]
            app.last_enrichr_results = [html.Tr([html.Th('gene set name'),
                             html.Th('p-value'),
                             html.Th('z-score'),
                             html.Th('combined score')])] + \
                    [html.Tr([html.Td(r[1]),
                              html.Td(r[2]),
                              html.Td(r[3]),
                              html.Td(r[4])
                              ])
                    for r in results]
        return app.last_enrichr_results

    app.last_backend_update = None
    app.split_clicks = 0
    app.merge_clicks = 0
    # there are two callbacks dealing with splitting/merging clusters;
    # one that updates a text displaying the selected clusters,
    # and another that actually makes the backend call to update m and w.
    @app.callback(
            Output('update-area', 'children'),
            [Input('means', 'selectedData'),
             Input('split', 'n_clicks'),
             Input('merge', 'n_clicks')])
    def split_or_merge_cluster(selected_points, n_click_split, n_click_merge):
        selected_clusters = []
        cluster_counts = []
        for point in selected_points['points']:
            cluster = point['curveNumber']
            selected_clusters.append(cluster)
        selected_clusters = list(set(selected_clusters))
        for cluster in selected_clusters:
            cluster_counts.append((app.labels == cluster).sum())
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
        return 'Selected clusters: ' + ' '.join(map(lambda x: '{0} ({1} cells)'.format(x[0], x[1]), zip(selected_clusters, cluster_counts)))

    # callback for split/merge operations (actually performing the operations)
    @app.callback(
            Output('cluster-view', 'children'),
            [Input('means', 'selectedData'),
             Input('split', 'n_clicks'),
             Input('merge', 'n_clicks')])
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
            generate_analysis.generate_analysis_resubmit(data_dir,
                    'split',
                    selected_clusters,
                    app.gene_names,
                    **app.uncurl_args)
            initialize(app, data_dir, permalink, user_id, test_or_user)
            return generate_cluster_view(app.mds_means)
        # merge clusters
        elif n_click_merge > app.merge_clicks:
            app.merge_clicks = n_click_merge
            generate_analysis.generate_analysis_resubmit(data_dir,
                    'merge',
                    selected_clusters,
                    app.gene_names,
                    **app.uncurl_args)
            initialize(app, data_dir, permalink, user_id, test_or_user)
            return generate_cluster_view(app.mds_means)
        else:
            raise Exception('placeholder')

def initialize_layout(app):
    app.initialized = False
    app.layout = html.Div(
        id='app-layout',
        children=[
            dcc.Location(id='url', refresh=False),
            html.Div(id='page-content', className='container')
        ]
    )


if __name__ == '__main__':
    app = dash.Dash(name='Cluster view', sharing=True)

    initialize(app, 'test/test1_output')
    app.run_server()
