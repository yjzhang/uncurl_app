import json
import os, os.path

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np

app = dash.Dash()

def generate_cluster_view(M, dim_red, top_genes, n_genes=10):
    """
    Generates a cluster view: MDS plot of means on the left, updated bar plot
    of top genes with c-scores on the right.
    """
    return html.Div([
        # view 1: graph plot
        dcc.Graph(id='means',
            figure={
                'data': [
                    go.Scatter(
                        x=dim_red[0,:],
                        y=dim_red[1,:],
                        mode='markers',
                        text=['cluster ' + str(x) for x in range(dim_red.shape[1])],
                        marker={
                            'size': 50
                        },
                    ),
                ],
                'layout': go.Layout(
                    title='Cluster means',
                    xaxis={'title': 'dim1'},
                    yaxis={'title': 'dim2'},
                ),
            }),
        # view 2: top genes
        dcc.Graph(id='top-genes',
            figure={
                'data': [
                    go.Bar(
                        x=[1,2,3],
                        y=[1,2,3],
                        orientation='h'
                    )],
                'layout': {
                    'title': 'Top genes',
                }
            }),
    ])


def initialize(data_dir=None):
    M = None
    mds_output = np.array([[1,2,3,4],[1,2,3,4]])
    top_genes = {0: [(0,100),(1,50),(2,40)], 1: [(0,50),(1,45),(2,30)]}
    if data_dir != None:
        M = np.loadtxt(os.path.join(data_dir, 'm.txt'))
        mds_output = np.loadtxt(os.path.join(data_dir, 'mds.txt'))
        with open(os.path.join(data_dir, 'top_genes.txt')) as f:
            top_genes = json.load(f)
    app.layout = generate_cluster_view(M, mds_output, top_genes)

    # create callback for clicking on clusters
    @app.callback(
            Output(component_id='top-genes', component_property='figure'),
            [Input(component_id='means', component_property='clickData')]
    )
    def update_output_div(input_value):
        print('input value:')
        print(input_value)
        if input_value is None:
            input_value = 0
        else:
            input_value = input_value['points'][0]['pointIndex']
        selected_top_genes = top_genes[input_value][:10]
        return {
                'data': [
                    go.Bar(
                        y=[str(x[0]) for x in selected_top_genes],
                        x=[x[1] for x in selected_top_genes],
                        orientation='h',
                    )
                ],
                'layout': {
                    'title': 'Top genes',
                }
        }
        #return 'You\'ve entered "{0}"'.format(input_value)


if __name__ == '__main__':
    initialize()
    app.run_server()
