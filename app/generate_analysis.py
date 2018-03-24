# This generates uncurl analyses as static HTML files.
# inputs: raw data (matrix or file), M (or file), W (or file), reduced_data (2d for vis)
import json
import os

import numpy as np
import scipy.io
from scipy import sparse
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import uncurl
from uncurl.sparse_utils import symmetric_kld
import uncurl_analysis
from uncurl_analysis import gene_extraction, relabeling

def generate_uncurl_analysis(data, output_dir,
        data_type='dense',
        gene_names=None,
        gene_sub=True,
        dim_red_option='mds',
        bulk_data_dir=None,
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
    if isinstance(data, str):
        if data_type == 'dense':
            data = np.loadtxt(data)
        elif data_type == 'sparse':
            data = scipy.io.mmread(data)
            data = sparse.csc_matrix(data)
        # TODO: copy file instead of calling mmwrite
    if sparse.issparse(data):
        data = sparse.csc_matrix(data)
        scipy.io.mmwrite(os.path.join(output_dir, 'data.mtx'), data)
    else:
        np.savetxt(os.path.join(output_dir, 'data.txt'), data)
    if isinstance(gene_names, str):
        gene_names = np.loadtxt(gene_names, dtype=str)
    # run uncurl
    if gene_sub:
        genes_subset = np.array(uncurl.max_variance_genes(data))
    else:
        genes_subset = np.array(uncurl.max_variance_genes(data, 1, 1.0))
    np.savetxt(os.path.join(output_dir, 'gene_subset.txt'), genes_subset, fmt='%d')
    data = data[genes_subset,:]
    if gene_names is not None:
        gene_names = gene_names[genes_subset]
    print(uncurl_kwargs)
    m, w, ll = uncurl.run_state_estimation(data, **uncurl_kwargs)
    np.savetxt(os.path.join(output_dir, 'm.txt'), m)
    np.savetxt(os.path.join(output_dir, 'w.txt'), w)
    analysis_postprocessing(data, m, w, output_dir, gene_names,
            dim_red_option)

def analysis_postprocessing(data, m, w, output_dir,
        gene_names=None,
        dim_red_option='mds'):
    """
    Runs post-processing steps on the results of uncurl.
    """
    labels = w.argmax(0)
    np.savetxt(os.path.join(output_dir, 'labels.txt'), labels, fmt='%d')
    # find overexpressed genes for clusters
    top_genes = uncurl_analysis.find_overexpressed_genes(data, w.argmax(0))
    with open(os.path.join(output_dir, 'top_genes.txt'), 'w') as f:
        json.dump(top_genes, f)
    # get p-values for c-scores using permutation test
    permutations = gene_extraction.generate_permutations(data, m.shape[1],
            w.argmax(0),
            n_perms=100)
    p_values = gene_extraction.c_scores_to_pvals(top_genes, permutations)
    with open(os.path.join(output_dir, 'gene_pvals.txt'), 'w') as f:
        json.dump(p_values, f)
    # run mds
    mds_output = uncurl.dim_reduce(m, w, 2)
    print(mds_output.shape)
    np.savetxt(os.path.join(output_dir, 'mds_means.txt'), mds_output.T)
    # dim_red_option
    dim_red_option = dim_red_option.lower()
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
    # entropy
    ent = uncurl_analysis.entropy(w)
    np.savetxt(os.path.join(output_dir, 'entropy.txt'), ent)
    # TODO: implement cluster similarities to bulk

def generate_analysis_resubmit(data_dir,
        split_or_merge='split',
        clusters_to_change=[],
        gene_names=None,
        dim_red_option='mds',
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
        data_path = os.path.join(data_dir, 'data.txt')
        data = np.loadtxt(data_path)
    except:
        data_path = os.path.join(data_dir, 'data.mtx')
        data = scipy.io.mmread(data_path)
        data = sparse.csc_matrix(data)
    genes_subset1 = np.array(uncurl.max_variance_genes(data))
    genes_subset2 = np.array(uncurl.max_variance_genes(data, 1, 1.0))
    if m_new.shape[0]==len(genes_subset1):
        data = data[genes_subset1, :]
    elif m_new.shape[0]==len(genes_subset2):
        data = data[genes_subset2, :]
    if split_or_merge == 'split':
        c = clusters_to_change[0]
        m_new, w_new = relabeling.split_cluster(data, m_new, w_new,
                c, **uncurl_kwargs)
    elif split_or_merge == 'merge':
        m_new, w_new = relabeling.merge_clusters(data, m_new, w_new,
                clusters_to_change, **uncurl_kwargs)
    np.savetxt(os.path.join(data_dir, 'm.txt'), m_new)
    np.savetxt(os.path.join(data_dir, 'w.txt'), w_new)
    analysis_postprocessing(data, m_new, w_new, data_dir, gene_names,
            dim_red_option)
