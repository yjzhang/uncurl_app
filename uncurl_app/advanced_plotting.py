# plot dendrograms and heatmaps
# plot cluster-vs-gene dendrograms
# plot cluster-vs-cluster heatmaps, use SpectralCoclustering to reorder to make them look semi-diagonal

import json

import numpy as np

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
        spec = SpectralCoclustering(int(max(len(cluster1_values)/1.5, len(cluster2_values)/1.5, 2)))
        spec.fit(data + 1e-8)
        row_labels = spec.row_labels_
        column_labels = spec.column_labels_
        row_order = np.argsort(row_labels)[::-1]
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
            'font': {'size': 16},
            'height': 550,
            'width': 700,

        }
    }
    return json.dumps(output, cls=SimpleEncoder)

def gene_similarity(data_sampled_all_genes, all_gene_names, gene_names_left, gene_names_top, mode='full', method='pearson'):
    """
    Creates a diagonal gene-gene similarity map using data from all sampled cells.

    mode can be either 'full' or 'reduced'. If 'mode' is full, then this uses the full data matrix for comparison.
    If mode is 'reduced', then this uses the M matrix from uncurl.

    Returns a json dendrogram from plotly
    """
    # TODO: this should be able to use either m_full or the full data matrix
    import scipy.stats
    gene_name_indices = {x: i for i, x in enumerate(all_gene_names)}
    selected_gene_names_left = [x for x in gene_names_left if x in gene_name_indices]
    gene_indices_left = np.array([gene_name_indices[x] for x in selected_gene_names_left])
    selected_gene_names_top = [x for x in gene_names_top if x in gene_name_indices]
    gene_indices_top = np.array([gene_name_indices[x] for x in selected_gene_names_top])
    print('gene heatmap selected gene names:', selected_gene_names_left)
    print('gene heatmap selected gene ids:', gene_indices_left)
    data_subset_1 = data_sampled_all_genes[gene_indices_left, :]
    data_subset_2 = data_sampled_all_genes[gene_indices_top, :]
    from scipy import sparse
    if sparse.issparse(data_subset_1):
        data_subset_1 = data_subset_1.toarray()
        data_subset_2 = data_subset_2.toarray()
    # TODO: have different methods for calculating the correlation matrix
    correlations = np.zeros((len(gene_indices_left), len(gene_indices_top)))
    if method == 'pearson':
        for i in range(len(gene_indices_left)):
            for j in range(len(gene_indices_top)):
                correlations[i, j] = scipy.stats.pearsonr(data_subset_1[i], data_subset_2[j])
    output = {
        'data': [{
            'z': correlations.tolist(),
            'x': selected_gene_names_top.tolist(),
            'y': selected_gene_names_left.tolist(),
            'colorscale': 'Reds',
            'type': 'heatmap',
        }],
        'layout': {
            'font': {'size': 16},
            'height': max(550, 100+len(gene_indices_left)*25),
            'width': max(700, 100+len(gene_indices_top)*25),
        }
    }
    return json.dumps(output, cls=SimpleEncoder)



def dendrogram(data_sampled_all_genes, all_gene_names, selected_gene_names, cluster_name, cluster_data, use_log=False,
        use_normalize=False):
    """
    Returns a json dendrogram from plotly...

    Args:
        data_subset (array): csc matrix, created from data_sampled_all_genes
    """
    import plotly.graph_objects as go
    import plotly.figure_factory as ff
    # TODO
    if not isinstance(cluster_data[0], str):
        cluster_data = np.array(['c' + str(x) for x in cluster_data])
    cluster_values = [x for x in sorted(set(cluster_data))]
    cluster_indices = {c1 : i for i, c1 in enumerate(cluster_values)}
    gene_name_indices = {x: i for i, x in enumerate(all_gene_names)}
    selected_gene_names = [x for x in selected_gene_names if x in gene_name_indices]
    gene_indices = np.array([gene_name_indices[x] for x in selected_gene_names])
    print('dendrogram selected gene names:', selected_gene_names)
    print('dendrogram selected gene ids:', gene_indices)
    selected_gene_indices = {x: i for i, x in enumerate(selected_gene_names)}
    # take mean across all clusters
    data_cluster_means = np.zeros((len(all_gene_names), len(cluster_values)))
    for i, c in enumerate(cluster_values):
        data_cluster_means[:,i] = np.array(data_sampled_all_genes[:, cluster_data==c].mean(1)).flatten()
    data_cluster_means = data_cluster_means[gene_indices, :]
    if use_log:
        data_cluster_means = np.log(1+data_cluster_means)
    if use_normalize:
        # divide data by max value for each gene
        data_cluster_means = data_cluster_means/data_cluster_means.max(1, keepdims=True)
    # code based on https://plot.ly/python/dendrogram/
    # create top dendrogram - plot of cells
    fig = ff.create_dendrogram(data_cluster_means.T, orientation='bottom', labels=cluster_values)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'
    # Create Side Dendrogram - plot of genes
    dendro_side = ff.create_dendrogram(data_cluster_means, orientation='right', labels=selected_gene_names)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    # Create Heatmap
    # TODO: reorder selected_gene_names?
    gene_labels_y = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves_y = [selected_gene_indices[i] for i in gene_labels_y]
    dendro_leaves_x = fig['layout']['xaxis']['ticktext']
    dendro_leaves_x = [cluster_indices[i] for i in dendro_leaves_x]
    heat_data = data_cluster_means[dendro_leaves_y,:]
    heat_data = heat_data[:,dendro_leaves_x]

    heatmap = [
        go.Heatmap(
            x = dendro_leaves_x,
            y = selected_gene_names,
            z = heat_data,
            colorscale = 'Reds',
            showscale = False,
        )
    ]
    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)
    # Edit Layout
    fig.update_layout({'width':700, 'height':100+len(selected_gene_names)*25,
                         'showlegend':False, 'hovermode': 'closest',
                         })
    fig.update_layout(xaxis={'domain': [.15, 1],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'ticks':""})
    fig.update_layout(xaxis2={'domain': [0, .15],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""})
    fig.update_layout(yaxis={'domain': [0, .85],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': True,
                                  'side': 'right',
                                  'tickmode': 'array',
                                  'tickvals': dendro_side['layout']['yaxis']['tickvals'],
                                  'ticktext': gene_labels_y,
                                  'ticks': ""
                        })
    fig.update_layout(yaxis2={'domain':[.825, .975],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""})
    fig.update_layout({'font': {'size': 16}})
    return fig.to_json()

