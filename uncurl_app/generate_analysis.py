# This generates uncurl analyses as static HTML files.
# inputs: raw data (matrix or file), M (or file), W (or file), reduced_data (2d for vis)
import json
import os

import numpy as np
from scipy import sparse

from uncurl_analysis import sc_analysis

def generate_uncurl_analysis(data, output_dir,
        **uncurl_kwargs):
    """
    Performs an uncurl analysis of the data, writing the results in the given
    directory. Assumes that output_dir contains a file named params.json, with
    all the parameters.

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
            contains params.json, data.mtx/.txt/.gz, and optionally gene_names.txt.
        **uncurl_kwargs: arguments to pass to uncurl.run_state_estimation..
    """
    # TODO: what about init?
    try:
        os.makedirs(output_dir)
    except:
        print('could not make output dir: {0}'.format(output_dir))
    with open(os.path.join(output_dir, 'submitted'), 'w') as f:
        f.write('')
    data_is_sparse = True
    if not isinstance(data, np.ndarray) and not isinstance(data, sparse.spmatrix):
        data_filename = data
        if data.endswith('.mtx') or data.endswith('.mtx.gz'):
            data_is_sparse = True
        else:
            data_is_sparse = False
    else:
        pass
    with open(os.path.join(output_dir, 'uncurl_kwargs.json'), 'w') as f:
        json.dump(uncurl_kwargs, f)
    sca = sc_analysis.SCAnalysis(output_dir,
            data_filename=data_filename,
            data_is_sparse=data_is_sparse)
    sca.load_params_json()
    try:
        sca.run_full_analysis()
    except Exception as e:
        import traceback
        text = traceback.format_exc()
        with open(os.path.join(output_dir, 'error.txt'), 'w') as f:
            f.write(text)
        return
    sca.save_json_reset()
    print('done with generate_analysis')


def generate_analysis_resubmit(sca,
        split_or_merge='split',
        clusters_to_change=[],
        **uncurl_kwargs):
    """
    Re-runs uncurl by splitting a cluster or merging clusters.

    Args:
        sca (SCAnalysis object)
        split_or_merge (str): either 'split' or 'merge'
        clusters_to_change (list): list of cluster numbers. If splitting, only the first cluster will be used. If merging, all the clusters will be merged.
    """
    sca.recluster(split_or_merge, clusters_to_change)
    sca.run_post_analysis()
    sca.save_json_reset()

def get_progress(path):
    """
    Returns the current preprocessing/analysis progress for the given path.

    Returns:
        current_task (str): a description of what step the preprocessing is on.
        time_remaining (str): a description of the time remaining.
    """
    if os.path.exists(os.path.join(path, 'error.txt')):
        with open(os.path.join(path, 'error.txt')) as f:
            text = f.read()
        return text, 'error'
    with open(os.path.join(path, 'params.json')) as f:
        preproc = json.load(f)
    genes = int(preproc['genes'])
    frac = float(preproc['genes_frac'])
    cell_frac = float(preproc['cell_frac'])
    cells = int(preproc['cells'])
    k = int(preproc['clusters'])
    # if k is automatically selected, just guess 15 for now...
    if k == 0:
        k = 15
    # calculate time remaining using genes and cells
    # wow this is really arbitrary but better than nothing???
    uncurl_factor = 120.0/(8000.0*3000.0*8)
    uncurl_total_time = k*genes*frac*cells*uncurl_factor
    vis_factor = 30.0/(1500.0*3000.0)
    visualization_time = genes*frac*cells*cell_frac*vis_factor
    pval_time = 70.0*k**2/8**2
    time_remaining = 500
    if os.path.exists(os.path.join(path, 'sc_analysis.json')):
        current_task = 'DONE'
    elif os.path.exists(os.path.join(path, 'top_genes.txt')):
        current_task = 'p-value calculations'
        time_remaining = 5
    elif os.path.exists(os.path.join(path, 'mds_data.txt')):
        current_task = 'differential expression'
        # time for t-tests is ~70 with 20k genes
        time_remaining = pval_time
    elif os.path.exists(os.path.join(path, 'baseline_vis.txt')):
        current_task = 'data visualization'
        time_remaining = pval_time + visualization_time
    elif os.path.exists(os.path.join(path, 'm.txt')):
        current_task = 'baseline visualization'
        time_remaining = pval_time + 2*visualization_time
    elif os.path.exists(os.path.join(path, 'progress.txt')):
        i = 0
        with open(os.path.join(path, 'progress.txt')) as f:
            i = int(f.read().strip())
        current_task = 'UNCURL progress: {0}/20'.format(i)
        time_remaining = pval_time + 2*visualization_time + uncurl_total_time*(20.0-i)/20.0
    else:
        # TODO: loading data time is nonzero
        current_task = 'loading data'
        time_remaining = pval_time + 2*visualization_time + uncurl_total_time
    time_remaining_minutes = int(time_remaining/60) + 1
    return current_task, '{0} minutes'.format(time_remaining_minutes)
