import json
import os, os.path


import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np

def get_test_dirs(base='test'):
    subdirectories = os.listdir(base)
    return subdirectories

def build_view_from_dirs(test_dirs, user_dirs):
    """
    Given a list of paths, builds a view that contains links to each
    of the paths.
    """
    div = html.Div(
            id='views',
            children=[
                html.Table(
                    [html.Tr(['Test examples'])] +
                    [html.Tr([dcc.Link(x, href='/test/{0}'.format(x)) for x in test_dirs])]
                ),
                html.Table(
                    [html.Tr(['User results'])] +
                    [html.Tr([dcc.Link(x, href='/user/{0}'.format(x)) for x in user_dirs])]
                )
            ]
        )
    return div

def initialize_views():
    try:
        test_dirs = get_test_dirs('test')
    except:
        test_dirs = []
    try:
        user_dirs = get_test_dirs('/tmp/uncurl')
    except:
        user_dirs = []
    links_table = build_view_from_dirs(test_dirs, user_dirs)
    return links_table


# code based on https://github.com/plotly/dash/issues/38#issuecomment-364927379
def router(pathname, **kwargs):
    """
    pathname is of the form <test or user>/<uuid>
    """
    #print(pathname)
    _ = pathname.split('/')
    #print(_)
    if len(_) < 3:
        return initialize_views()
    resource = _[1]
    resource_id = _[2]
    if resource == 'test':
        return initialize(os.path.join('test', resource_id))
    else:
        return initialize(os.path.join('/tmp/uncurl', resource_id))

#@app.callback(Output('page-content', 'children'), [Input('url', 'pathname')])
def display_page(pathname):
    #print('display_page')
    return router(pathname)

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
                    for c in range(max(labels))
                ],
                'layout': go.Layout(
                    title='Cells',
                    xaxis={'title': 'dim1'},
                    yaxis={'title': 'dim2'},
                    margin={'t':30},
                ),
        }

def generate_cluster_view(M, dim_red, top_genes, n_genes=10):
    """
    Generates a cluster view: MDS plot of means on the left, updated bar plot
    of top genes with c-scores on the right.
    """
    colorscale = 'Portland'
    #if dim_red.shape[1] > 10:
    #    colorscale = 'Jet'
    return html.Div([
        html.Div([
            dcc.RadioItems(
                id='means-or-cell',
                options=[{'label': 'Display Means', 'value': 'Means'},
                         {'label': 'Display Cells', 'value': 'Cells'}],
                value='Means',
                labelStyle={'display': 'inline-block'},
            )
        ]),
        # view 1: graph plot
        html.Div(dcc.Graph(id='means',
                figure=create_means_figure(dim_red, colorscale),
                style={'width': 700}),
            style={'display': 'inline-block'}),
        # view 2: top genes
        html.Div(dcc.Graph(id='top-genes',
            figure={
                'data': [
                    go.Bar(
                        x=[1,2,3],
                        y=[1,2,3],
                        orientation='h'
                    )],
                'layout': go.Layout(
                    title='Cluster means',
                    xaxis={'title': 'dim1'},
                    yaxis={'title': 'dim2'},
                ),
                }, style={'width': 400}),
            style={'display': 'inline-block'}),
        ], style={'width': '100%', 'display':'inline-block'})


# TODO: dynamically generate an app given a path???
def initialize(app, data_dir=None):
    """
    This function sets app.layout using a directory containing uncurl results.
    """
    app.initialized = True
    M = None
    W = None
    labels = None
    mds_means = np.array([[1,2,3,4],[1,2,3,4]])
    mds_data = None
    top_genes = {'0': [(0,100),(1,50),(2,40)], '1': [(0,50),(1,45),(2,30)]}
    gene_names = None
    #print('initialize ' + data_dir)
    if data_dir != None:
        M = np.loadtxt(os.path.join(data_dir, 'm.txt'))
        W = np.loadtxt(os.path.join(data_dir, 'w.txt'))
        labels = W.argmax(0)
        mds_means = np.loadtxt(os.path.join(data_dir, 'mds_means.txt'))
        mds_data = np.loadtxt(os.path.join(data_dir, 'mds_data.txt'))
        with open(os.path.join(data_dir, 'top_genes.txt')) as f:
            top_genes = json.load(f)
        try:
            gene_names = np.loadtxt(os.path.join(data_dir, 'gene_names.txt'), dtype=str)
        except:
            gene_names = np.array([str(x) for x in range(len(mds_data.shape[1]))])

    # generate layout
    #app.layout.children[1].children = generate_cluster_view(M, mds_means, top_genes)
    app.layout = generate_cluster_view(M, mds_means, top_genes)

    # create callback for clicking on clusters
    @app.callback(
            Output(component_id='top-genes', component_property='figure'),
            [Input(component_id='means', component_property='clickData')]
    )
    def update_output_div(input_value):
        #print('input value:')
        #print(input_value)
        if input_value is None:
            input_value = '0'
        else:
            input_value = str(input_value['points'][0]['curveNumber'])
        selected_top_genes = top_genes[input_value][:10]
        selected_gene_names = [gene_names[x[0]] for x in selected_top_genes]
        return {
                'data': [
                    go.Bar(
                        y=selected_gene_names,
                        x=[x[1] for x in selected_top_genes],
                        orientation='h',
                    )
                ],
                'layout': go.Layout(
                    title='Top genes for cluster {0}'.format(input_value),
                    xaxis={'title': 'c-score'},
                    yaxis={'title': 'genes'},
                ),
        }

    # create callback for switching to cell view
    @app.callback(
        Output(component_id='means', component_property='figure'),
        [Input(component_id='means-or-cell', component_property='value')]
    )
    def update_scatterplot(input_value):
        #print('update_scatterplot')
        #print(input_value)
        if input_value == 'Means':
            return create_means_figure(mds_means)
        elif input_value == 'Cells':
            return create_cells_figure(mds_data, labels)

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
    app.layout = html.Div(
        id='app-layout',
        children=[
            dcc.Location(id='url', refresh=False),
            html.Div(id='page-content', className='container')
        ]
    )

    initialize(app, 'test/test1_output')
    app.run_server()
