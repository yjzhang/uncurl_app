# TODO: plot dendrograms and heatmaps
# plot cluster-vs-gene dendrograms
# plot cluster-vs-cluster heatmaps, use SpectralCoclustering to reorder to make them look semi-diagonal

import json

import numpy as np
import plotly.figure_factory as ff

from .utils import SimpleEncoder
from .cache import cache

@cache.memoize()
def cluster_heatmap(cluster1, cluster2, cluster_1_name, cluster_2_name, order='coclustering', normalize_row=True, **params):
    """
    Returns a plotly-formated json that plots the two heatmaps.

    Args:
        cluster1 (array): array of strings
        cluster2 (array): array of strings
        order (str): 'coclustering' or 'alphabetical'. Default: 'coclustering'
        normalize_row (bool): whether or not to normalize by row (so that each row sums to 1). Default: True
    """
    if not isinstance(cluster1[0], str):
        cluster1 = np.array(['c' + str(x) for x in cluster1])
    if not isinstance(cluster2[0], str):
        cluster2 = np.array(['c' + str(x) for x in cluster2])
    cluster1_values = [x for x in sorted(set(cluster1))]
    c1_indices = {c1 : i for i, c1 in enumerate(cluster1_values)}
    cluster2_values = [x for x in sorted(set(cluster2))]
    c2_indices = {c1 : i for i, c1 in enumerate(cluster2_values)}
    data = np.zeros((len(cluster1_values), len(cluster2_values)))
    for i in range(len(cluster1)):
        c1 = cluster1[i]
        c2 = cluster2[i]
        data[c1_indices[c1], c2_indices[c2]] += 1
    # create heatmap json
    if normalize_row:
        data = data/data.sum(1, keepdims=True)
    if order == 'coclustering':
        from sklearn.cluster.bicluster import SpectralCoclustering
        spec = SpectralCoclustering(int(max(len(cluster1_values)/2, len(cluster2_values)/2, 2)))
        spec.fit(data + 1e-8)
        row_labels = spec.row_labels_
        column_labels = spec.column_labels_
        row_order = np.argsort(row_labels)
        col_order = np.argsort(column_labels)
        row_labels = np.array(cluster1_values)[row_order]
        column_labels = np.array(cluster2_values)[col_order]
        data = data[row_order, :]
        data = data[:, col_order]
    else:
        row_labels = np.array(cluster1_values)
        column_labels = np.array(cluster2_values)
    output = {
        'data': [{
            'z': data.tolist(),
            'x': column_labels.tolist(),
            'y': row_labels.tolist(),
            'colorscale': 'Reds',
            'type': 'heatmap',
        }],
        'layout': {
            'xaxis': {'title': cluster_2_name, 'automargin': True},
            'yaxis': {'title': cluster_1_name, 'automargin': True},
        }
    }
    print(output)
    return json.dumps(output, cls=SimpleEncoder)

@cache.memoize()
def dendrogram():
    """
    Returns a json dendrogram from plotly...
    """

