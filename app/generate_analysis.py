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
from uncurl_analysis import sc_analysis
from uncurl_analysis import gene_extraction, relabeling

import simplex_sample

def generate_uncurl_analysis(data, output_dir,
        data_type='dense',
        clusters=10,
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
        clusters (int): number of clusters (for uncurl state estimation).
        data_type (str): if data is a path, this indicates whether the data is a dense or sparse array.
        gene_names (list or array): list of all gene names
        gene_sub (bool): whether or not to use gene subset selection (max_variance_genes)
        **uncurl_kwargs: arguments to pass to uncurl.run_state_estimation..
    """
    try:
        os.makedirs(output_dir)
    except:
        print('could not make output dir: {0}'.format(output_dir))
    with open(os.path.join(output_dir, 'submitted'), 'w') as f:
        f.write('')
    data_is_sparse = True
    if isinstance(data, str) or isinstance(data, unicode):
        data_filename = data
        if data.endswith('.mtx') or data.endswith('.mtx.gz'):
            data_is_sparse = True
        else:
            data_is_sparse = False
    else:
        pass
    sca = sc_analysis.SCAnalysis(output_dir,
            data_filename=data_filename,
            data_is_sparse=data_is_sparse,
            normalize=normalize,
            min_reads=min_reads,
            max_reads=max_reads,
            frac=frac,
            cell_frac=cell_frac,
            clusters=clusters,
            dim_red_option=dim_red_option,
            baseline_dim_red=baseline_dim_red,
            **uncurl_kwargs)
    sca.run_full_analysis()
    sca.save_pickle_reset()


def generate_analysis_resubmit(sca,
        split_or_merge='split',
        clusters_to_change=[]):
    """
    Re-runs uncurl by splitting a cluster or merging clusters.

    Args:
        sca (SCAnalysis object)
        split_or_merge (str): either 'split' or 'merge'
        clusters_to_change (list): list of cluster numbers. If splitting, only the first cluster will be used. If merging, all the clusters will be merged.
    """
    sca.recluster(split_or_merge, clusters_to_change)
    sca.run_post_analysis()
    sca.save_pickle_reset()
