import os

from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.palettes import Accent8
from bokeh.models import ColumnDataSource

import numpy as np

import uncurl

def vis_clustering(data, assignments, centers, user_id):
    """
    Visualizes a hard-assignment clustering. Presents data in 2D.
    """
    X = uncurl.dim_reduce_data(data,2)
    colors = [Accent8[c] for c in assignments]
    source = ColumnDataSource(dict(
        x=X[:,0],
        y=X[:,1],
        color=colors,
        label=assignments))
    f1 = figure(active_scroll='wheel_zoom')
    f1.circle(x='x', y='y',
            color='color', fill_alpha=0.7, size=10,
            legend='label', source=source)
    f1.xaxis.axis_label = 'dim1'
    f1.yaxis.axis_label = 'dim2'
    # TODO: add javascript interactions???
    script, div = components(f1)
    with open(os.path.join('/tmp/', user_id, 'vis_clustering.html'), 'w') as f:
        f.write(div)
        f.write('\n\n')
        f.write(script)
    return script, div


def vis_state_estimation(data, M, W, user_id):
    """
    2d bokeh visualization of 2d data. Outputs an embeddable html fragment to
    /tmp/user_id/vis_state_estimation.html
    """
    X = uncurl.dim_reduce(M, W, 2)
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

def vis_lineage(M, W, smoothed_data, edges, clusters, user_id):
    """
    2d bokeh visualization of lineage estimation output. Outputs an embeddable
    html fragment to /tmp/user_id/vis_lineage.html
    """
    X = uncurl.dim_reduce(M, W, 2)
    reduced_data = np.dot(X.T, W)
    colors = [Accent8[c] for c in clusters]
    source = ColumnDataSource(dict(
        x=reduced_data[0,:],
        y=reduced_data[1,:],
        color=colors,
        label=clusters))
    source_smoothed = ColumnDataSource(dict(
        x=smoothed_data[0,:],
        y=smoothed_data[1,:],
        color=colors,
        label=clusters))
    f1 = figure(active_scroll='wheel_zoom')
    f1.circle(x='x', y='y',
            color='color', fill_alpha=0.7, size=10,
            legend='label', source=source)
    # only include the legend once
    f1.circle(x='x', y='y',
            color='color', fill_alpha=1.0, size=5,
            source=source_smoothed)
    # use multi_line to connect the edges
    xs = []
    ys = []
    for edge in edges:
        p1_x = smoothed_data[0, edge[0]]
        p1_y = smoothed_data[1, edge[0]]
        p2_x = smoothed_data[0, edge[1]]
        p2_y = smoothed_data[1, edge[1]]
        xs.append([p1_x, p2_x])
        ys.append([p1_y, p2_y])
    f1.multi_line(xs, ys, line_color='black')
    f1.xaxis.axis_label = 'dim1'
    f1.yaxis.axis_label = 'dim2'
    script, div = components(f1)
    with open(os.path.join('/tmp/', user_id, 'vis_lineage.html'), 'w') as f:
        f.write(div)
        f.write('\n\n')
        f.write(script)
    return script, div
