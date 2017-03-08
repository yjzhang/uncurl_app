import os

from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.palettes import Accent8
from bokeh.models import ColumnDataSource

import numpy as np

import uncurl

def vis_clustering(data, assignments, centers):
    """
    Visualizes a hard-assignment clustering. Presents data in 2D.
    """
    X = uncurl.dim_reduce(data)

def vis_state_estimation(data, M, W, user_id):
    """
    2d bokeh visualization of 2d data. Outputs an embeddable html fragment to
    /tmp/user_id/vis_state_estimation.html
    """
    X = uncurl.dim_reduce(data, M, W, 2)
    # reduced_data is of dimensionality 2 x cells
    reduced_data = np.dot(X.T, W)
    # zero-center the data
    reduced_data = reduced_data - reduced_data.mean(1, keepdims=True)
    np.savetxt(os.path.join('/tmp/', user_id, 'reduced_data.txt'), reduced_data)
    clusters = W.argmax(0)
    colors = [Accent8[c] for c in clusters]
    source = ColumnDataSource(dict(
        x=reduced_data[0,:],
        y=reduced_data[1,:],
        color=colors,
        label=clusters))
    f1 = figure(active_scroll='wheel_zoom')
    f1.circle(x='x', y='y',
            color='color', fill_alpha=0.7, size=10,
            legend='label', source=source)
    f1.xaxis.axis_label = 'dim1'
    f1.yaxis.axis_label = 'dim2'
    # TODO: add javascript interactions???
    script, div = components(f1)
    with open(os.path.join('/tmp/', user_id, 'vis_state_estimation.html'), 'w') as f:
        f.write(div)
        f.write('\n\n')
        f.write(script)
    return script, div