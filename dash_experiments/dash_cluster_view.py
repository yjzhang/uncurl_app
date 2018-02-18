import json
import os, os.path

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np

app = dash.Dash(name='Cluster view', sharing=True)
app.layout = html.Div(
        id='app-layout',
        children=[
            dcc.Location(id='url', refresh=False),
            html.Div(id='page-content', className='container')
        ]
    )

# code based on https://github.com/plotly/dash/issues/38#issuecomment-364927379
def router(pathname, **kwargs):
    """
    pathname is of the form <test or id>/<uuid>
    """
    _ = pathname.split('/')
    resource = _[1]
    resource_id = _[2]
    if resource == 'test':
        return initialize(os.path.join('test', resource_id))
    else:
        return initialize(os.path.join('/tmp/uncurl', resource_id))

@app.callback(Output('page-content', 'children'), [Input('url', 'pathname')])
def display_page(pathname):
    return router(pathname)

def generate_cluster_view(M, dim_red, top_genes, n_genes=10):
    """
    Generates a cluster view: MDS plot of means on the left, updated bar plot
    of top genes with c-scores on the right.
    """
    return html.Div([
        # view 1: graph plot
        html.Div(dcc.Graph(id='means',
            figure={
                'data': [
                    go.Scatter(
                        x=dim_red[0,:],
                        y=dim_red[1,:],
                        mode='markers',
                        text=['cluster ' + str(x) for x in range(dim_red.shape[1])],
                        marker={
                            'size': 20
                        },
                    ),
                ],
                'layout': go.Layout(
                    title='Cluster means',
                    xaxis={'title': 'dim1'},
                    yaxis={'title': 'dim2'},
                ),
                }, style={'width': 700}),
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
def initialize(data_dir=None):
    M = None
    mds_means = np.array([[1,2,3,4],[1,2,3,4]])
    mds_data = None
    top_genes = {'0': [(0,100),(1,50),(2,40)], '1': [(0,50),(1,45),(2,30)]}
    gene_names = None
    if data_dir != None:
        M = np.loadtxt(os.path.join(data_dir, 'm.txt'))
        mds_means = np.loadtxt(os.path.join(data_dir, 'mds_means.txt'))
        mds_data = np.loadtxt(os.path.join(data_dir, 'mds_data.txt'))
        with open(os.path.join(data_dir, 'top_genes.txt')) as f:
            top_genes = json.load(f)
        try:
            gene_names = np.loadtxt(os.path.join(data_dir, 'gene_names.txt'), dtype=str)
        except:
            gene_names = np.array([str(x) for x in range(len(mds_data.shape[1]))])
    app.layout = generate_cluster_view(M, mds_means, top_genes)

    # create callback for clicking on clusters
    @app.callback(
            Output(component_id='top-genes', component_property='figure'),
            [Input(component_id='means', component_property='clickData')]
    )
    def update_output_div(input_value):
        print('input value:')
        print(input_value)
        if input_value is None:
            input_value = '0'
        else:
            input_value = str(input_value['points'][0]['pointIndex'])
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
                    autosize=False,
                ),
        }


if __name__ == '__main__':
    initialize('../data/test1_output')
    app.run_server()
