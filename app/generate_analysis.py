# This generates uncurl analyses as static HTML files.
# inputs: raw data (matrix or file), M (or file), W (or file), reduced_data (2d for vis)
import json
import os

import numpy as np
import scipy.io
from scipy import sparse
from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA, TruncatedSVD

import uncurl
from uncurl.sparse_utils import symmetric_kld

import uncurl_analysis
from uncurl_analysis import gene_extraction, relabeling

from simplex_sampling import simplex_sample

def generate_uncurl_analysis(data, output_dir,
        data_type='dense',
        gene_names=None,
        gene_sub=True,
        dim_red_option='mds',
        bulk_data_dir=None,
        normalize=False,
        min_reads=0,
        max_reads=1e10,
        frac=0.2,
        cell_frac=1.0,
        baseline_dim_red='none',
        **uncurl_kwargs):
    """
    Performs an uncurl analysis of the data, writing the results in the given
    directory.

    Outputs:
        output_dir/data.txt or output_dir/data.mtx
        output_dir/m.txt
        output_dir/w.txt
        output_dir/labels.txt (integer labels)
        output_dir/top_genes.txt (json of a dict mapping cluster ids to a list of (gene_id : c_score) sorted by c_score)
        output_dir/mds_means.txt (mds of the means - 2 x k)
        output_dir/mds_data.txt (mds projection of data - 2 x n)
        output_dir/gene_subset.txt (gene subset selected by uncurl)
        output_dir/gene_names.txt (list of all gene names in data subset)
        output_dir/entropy.txt (entropy of cell labels)

    Args:
        data (array or str): either a data array, or a string containing
            the path to a data array..
        output_dir (str): directory to write output to.
        data_type (str): if data is a path, this indicates whether the data is a dense or sparse array.
        gene_names (list or array): list of all gene names
        gene_sub (bool): whether or not to use gene subset selection (max_variance_genes)
        **uncurl_kwargs: arguments to pass to uncurl.run_state_estimation. has to include clusters=k.
    """
    try:
        os.makedirs(output_dir)
    except:
        print('could not make output dir: {0}'.format(output_dir))
    with open(os.path.join(output_dir, 'submitted'), 'w') as f:
        f.write('')
    if isinstance(data, str) or isinstance(data, unicode):
        if data.endswith('.mtx') or data.endswith('.mtx.gz'):
            data = scipy.io.mmread(data)
            data = sparse.csc_matrix(data)
        else:
            try:
                data = np.loadtxt(data)
            except:
                data = scipy.io.mmread(data)
                data = sparse.csc_matrix(data)
    if isinstance(gene_names, str) or isinstance(gene_names, unicode):
        gene_names = np.loadtxt(gene_names, dtype=str)
    if sparse.issparse(data):
        data = sparse.csc_matrix(data)
    # cell subset
    read_counts = np.array(data.sum(0)).flatten()
    cells_subset = (read_counts >= min_reads) & (read_counts <= max_reads)
    np.savetxt(os.path.join(output_dir, 'cells_subset.txt'), cells_subset,
            fmt='%s')
    data = data[:, cells_subset]
    # normalize data
    if normalize:
        data = uncurl.preprocessing.cell_normalize(data)
        with open(os.path.join(output_dir, 'normalize_data.txt'), 'w') as f:
            f.write('')
    # gene subset
    if gene_sub:
        gene_subset = np.array(uncurl.max_variance_genes(data,
            nbins=5,
            frac=frac))
    else:
        gene_subset = np.array(uncurl.max_variance_genes(data, 1, 1.0))
    print(repr(gene_subset))
    np.savetxt(os.path.join(output_dir, 'gene_subset.txt'), gene_subset, fmt='%d')
    data_subset = data[gene_subset,:]
    print(uncurl_kwargs)
    print(repr(data_subset))
    # run uncurl
    m, w, ll = uncurl.run_state_estimation(data_subset, **uncurl_kwargs)
    np.savetxt(os.path.join(output_dir, 'm.txt'), m)
    np.savetxt(os.path.join(output_dir, 'w.txt'), w)
    print('uncurl done')
    # TODO: cell subset
    if cell_frac < 1:
        sampled_cells = cell_subset(w, cell_frac)
    else:
        sampled_cells = np.arange(data_subset.shape[1])
    np.savetxt(os.path.join(output_dir, 'cell_sample.txt'), sampled_cells, fmt='%d')
    # run baseline dimensionality reduction
    baseline_vis(data_subset[:, sampled_cells], baseline_dim_red)
    print('baseline vis done')
    # run postprocessing
    analysis_postprocessing(data, m, w, output_dir, gene_names,
            dim_red_option, cell_frac)


def cell_subset(w, cell_frac):
    """
    Returns a set of indices.
    """
    k, cells = w.shape
    n_samples = int(cells*cell_frac)
    samples = simplex_sample.sample(k, n_samples)
    indices = simplex_sample.data_sample(w, samples)
    return indices


def baseline_vis(data, baseline_dim_red, output_dir):
    """
    Generates a visualization using the unprocessed data.
    """
    baseline_dim_red = baseline_dim_red.lower()
    if baseline_dim_red == 'none':
        return None
    else:
        tsvd = TruncatedSVD(50)
        data_log_norm = uncurl.preprocessing.log1p(uncurl.preprocessing.cell_normalize(data))
        if baseline_dim_red == 'tsne':
            data_tsvd = tsvd.fit_transform(data_log_norm.T)
            tsne = TSNE(2)
            data_dim_red = tsne.fit_transform(data_tsvd)
        elif baseline_dim_red == 'tsvd':
            tsvd2 = TruncatedSVD(2)
            data_dim_red = tsvd2.fit_transform(data_log_norm.T)
        elif baseline_dim_red == 'mds':
            data_tsvd = tsvd.fit_transform(data_log_norm.T)
            mds = MDS(2)
            data_dim_red = mds.fit_transform(data_tsvd)
        with open(os.path.join(output_dir, 'baseline_vis.txt'), 'w') as f:
            np.savetxt(f, data_dim_red.T)
        return data_dim_red


def analysis_postprocessing(data, m, w, output_dir,
        gene_names=None,
        dim_red_option='mds',
        cell_frac=1.0):
    """
    Runs post-processing steps on the results of uncurl.
    """
    labels = w.argmax(0)
    np.savetxt(os.path.join(output_dir, 'labels.txt'), labels, fmt='%d')
    # find overexpressed genes for clusters
    top_genes = uncurl_analysis.find_overexpressed_genes(data, w.argmax(0))
    with open(os.path.join(output_dir, 'top_genes.txt'), 'w') as f:
        json.dump(top_genes, f)
    print('c-scores done')
    # get p-values for c-scores using permutation test
    permutations = gene_extraction.generate_permutations(data, m.shape[1],
            w.argmax(0),
            n_perms=100)
    p_values = gene_extraction.c_scores_to_pvals(top_genes, permutations)
    with open(os.path.join(output_dir, 'gene_pvals.txt'), 'w') as f:
        json.dump(p_values, f)
    print('p vals done')
    # run mds
    mds_output = uncurl.dim_reduce(m, w, 2)
    print(mds_output.shape)
    np.savetxt(os.path.join(output_dir, 'mds_means.txt'), mds_output.T)
    # dim_red_option
    dim_red_option = dim_red_option.lower()
    with open(os.path.join(output_dir, dim_red_option + '.txt'), 'w') as f:
        f.write('')
    if cell_frac < 1:
        if not os.path.exists(os.path.join(output_dir, 'cell_sample.txt')):
            cells_to_select = cell_subset(w, cell_frac)
            np.savetxt(os.path.join(output_dir, 'cell_sample.txt'), cells_to_select)
        else:
            cells_to_select = np.loadtxt(os.path.join(output_dir, 'cell_sample.txt'), dtype=int)
        w = w[:, cells_to_select]
    if dim_red_option == 'mds':
        mds_data = uncurl.mds(m, w, 2)
        np.savetxt(os.path.join(output_dir, 'mds_data.txt'), mds_data)

    elif dim_red_option == 'tsne':
        tsne = TSNE(2, metric=symmetric_kld)
        tsne_w = tsne.fit_transform(w.T)
        np.savetxt(os.path.join(output_dir, 'mds_data.txt'), tsne_w.T)
    elif dim_red_option == 'pca':
        pca = PCA(2)
        pca_w =  pca.fit_transform(w.T)
        np.savetxt(os.path.join(output_dir, 'mds_data.txt'), pca_w.T)
    if gene_names is not None:
        np.savetxt(os.path.join(output_dir, 'gene_names.txt'), gene_names, fmt='%s')
    print('dimensionality reduction done')
    # entropy
    ent = uncurl_analysis.entropy(w)
    np.savetxt(os.path.join(output_dir, 'entropy.txt'), ent)
    # TODO: implement cluster similarities to bulk

def generate_analysis_resubmit(data_dir,
        split_or_merge='split',
        clusters_to_change=[],
        gene_names=None,
        **uncurl_kwargs):
    """
    Re-runs uncurl by splitting a cluster or merging clusters.

    Args:
        data_dir (str): data directory originally generated by generate_uncurl_analysis
        split_or_merge (str): either 'split' or 'merge'
        clusters_to_change (list): list of cluster numbers. If splitting, only the first cluster will be used. If merging, all the clusters will be merged.
    """
    m_new = np.loadtxt(os.path.join(data_dir, 'm.txt'))
    w_new = np.loadtxt(os.path.join(data_dir, 'w.txt'))
    # load data
    try:
        with open(os.path.join(data_dir, 'params.json')) as f:
            params = json.load(f)
    except:
        params = {}
    try:
        data_path = os.path.join(data_dir, 'data.txt')
        data = np.loadtxt(data_path)
    except:
        data_path = os.path.join(data_dir, 'data.mtx')
        data = scipy.io.mmread(data_path)
        data = sparse.csc_matrix(data)
    try:
        cells_subset = np.loadtxt(os.path.join(data_dir, 'cells_subset.txt'),
                dtype=bool)
        data = data[:, cells_subset]
    except:
        pass
    if 'normalize_data' in params:
        data = uncurl.preprocessing.cell_normalize(data)
    cell_frac = 1.0
    gene_frac = 0.2
    if 'cell_frac' in params:
        cell_frac = params['cell_frac']
    if 'gene_frac' in params:
        gene_frac = params['gene_frac']
    gene_subset = np.array(uncurl.max_variance_genes(data, 5, gene_frac))
    data = data[gene_subset, :]
    if split_or_merge == 'split':
        c = clusters_to_change[0]
        m_new, w_new = relabeling.split_cluster(data, m_new, w_new,
                c, **uncurl_kwargs)
    elif split_or_merge == 'merge':
        m_new, w_new = relabeling.merge_clusters(data, m_new, w_new,
                clusters_to_change, **uncurl_kwargs)
    np.savetxt(os.path.join(data_dir, 'm.txt'), m_new)
    np.savetxt(os.path.join(data_dir, 'w.txt'), w_new)
    dim_red_option = 'mds'
    if os.path.exists(os.path.join(data_dir, 'tsne.txt')):
        dim_red_option = 'tsne'
    elif os.path.exists(os.path.join(data_dir, 'pca.txt')):
        dim_red_option = 'pca'
    analysis_postprocessing(data, m_new, w_new, data_dir, gene_names,
            dim_red_option, cell_frac)
