import json
import os, os.path

from flask import url_for

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np

def create_means_figure(dim_red, colorscale='Portland'):
    """
    create a figure for displaying means
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

def create_cells_figure(dim_red, labels, colorscale='Portland'):
    """
    create a figure for displaying cells
    """
    return {
                'data': [
                    go.Scatter(
                        x=dim_red[0,labels==c],
                        y=dim_red[1,labels==c],
                        mode='markers',
                        name='cluster ' + str(c),
                        marker={
                            'size': 10,
                            'color': c,
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
        cluster_name):
    """
    Creates a figure for displaying top genes
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
                    xaxis={'title': 'c-score'},
                    yaxis={'title': 'genes'},
                ),
        }


def create_bulk_correlation_figure(correlations, bulk_names, n_datasets=10):
    """
    Creates a figure for displaying correlations with bulk datasets
    """
    # TODO

def generate_cluster_view(dim_red, top_genes, n_genes=10):
    """
    Generates a cluster view: MDS plot of means on the left, updated bar plot
    of top genes with c-scores on the right.
    """
    colorscale = 'Portland'
    #if dim_red.shape[1] > 10:
    #    colorscale = 'Jet'
    return html.Div([
        # view 1: graph plot
        html.Div([
            dcc.RadioItems(
                id='means-or-cell',
                options=[{'label': 'Display Means', 'value': 'Means'},
                         {'label': 'Display Cells', 'value': 'Cells'}],
                value='Means',
                labelStyle={'display': 'inline-block'},
            ),
            dcc.Graph(id='means',
                figure=create_means_figure(dim_red, colorscale),
                style={'width': 700})
            ],
            style={'display': 'inline-block', 'width': 750, 'float':'left'}),
        # view 2: top genes
        html.Div([
            dcc.RadioItems(
                id='top-or-bulk',
                options=[{'label': 'Display top genes', 'value': 'top'},
                         {'label': 'Display bulk correlations', 'value': 'bulk'}],
                value='top',
                labelStyle={'display': 'inline-block'},
                style={'margin-top': -25},
            ),
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
        ],
        style={'width': '100%', 'display':'inline-block',
            'margin-top': 10})


# TODO: dynamically generate an app given a path???
def initialize(app, data_dir=None, permalink='test', user_id='test',
        test_or_user='test'):
    """
    This function sets app.layout using a directory containing uncurl results.
    """
    app.initialized = True
    # app variables: scatter_mode is either 'means' or 'cells'
    app.scatter_mode = 'means'
    # 'bar_node' is either 'top' or 'bulk'
    # top indicates that top genes should be displayed. 'bulk' indicates
    # that bulk correlations should be displayed.
    app.bar_mode = 'top'
    app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
    #M = None
    labels = None
    mds_means = np.array([[1,2,3,4],[1,2,3,4]])
    mds_data = None
    top_genes = {'0': [(0,100),(1,50),(2,40)], '1': [(0,50),(1,45),(2,30)]}
    gene_names = None
    #print('initialize ' + data_dir)
    if data_dir != None:
        labels = np.loadtxt(os.path.join(data_dir, 'labels.txt')).astype(int)
        mds_means = np.loadtxt(os.path.join(data_dir, 'mds_means.txt'))
        mds_data = np.loadtxt(os.path.join(data_dir, 'mds_data.txt'))
        with open(os.path.join(data_dir, 'top_genes.txt')) as f:
            top_genes = json.load(f)
        try:
            gene_names = np.loadtxt(os.path.join(data_dir, 'gene_names.txt'), dtype=str)
        except:
            M = np.loadtxt(os.path.join(data_dir, 'm.txt'))
            gene_names = np.array(['gene ' + str(x) for x in range(M.shape[0])])

    # generate layout
    #app.layout.children[1].children = generate_cluster_view(M, mds_means, top_genes)
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
            style={'display':'inline-block'}
        ),
        generate_cluster_view(mds_means, top_genes)
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
            selected_top_genes = top_genes[input_value][:num_genes]
            selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
            return create_top_genes_figure(selected_top_genes,
                    selected_gene_names, input_value)
        else:
            # TODO: show bulk correlations
            pass

    # create callback for switching to cell view
    @app.callback(
            Output(component_id='means', component_property='figure'),
            [Input(component_id='means-or-cell', component_property='value')]
    )
    def update_scatterplot(input_value):
        #print('update_scatterplot')
        #print(input_value)
        if input_value == 'Means':
            app.scatter_mode = 'means'
            return create_means_figure(mds_means)
        elif input_value == 'Cells':
            app.scatter_mode = 'cells'
            return create_cells_figure(mds_data, labels)

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
            selected_top_genes = top_genes[input_value][:num_genes]
            selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
            return '\n'.join(selected_gene_names)
        else:
            # TODO: show bulk correlations
            pass


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
