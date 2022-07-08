import json
import os

import numpy as np
import scipy.io
from scipy import sparse
from uncurl.sparse_utils import sparse_means_var_csc

from .utils import SimpleEncoder

class Summary(object):
    """
    This object contains a summary for a single-cell RNASeq dataset.
    """

    def __init__(self, data_paths, gene_paths, base_path, shapes=['gene_cell'], data=None, dataset_names=None, use_batch_correction=False):
        """
        Args:
            data_paths (list of str): list of paths to data files
            gene_paths (list of str): list of paths to gene files
            base_path (str): path to data dir
            shapes (list of str): list of either gene_cell or cell_gene
        """
        # deal with multiple paths
        if data is None:
            if data_paths is None:
                is_gz = False
                if os.path.exists(os.path.join(base_path, 'data.txt')):
                    data_path = os.path.join(base_path, 'data.txt')
                elif os.path.exists(os.path.join(base_path, 'data.mtx')):
                    data_path = os.path.join(base_path, 'data.mtx')
                elif os.path.exists(os.path.join(base_path, 'data.mtx.gz')):
                    data_path = os.path.join(base_path, 'data.mtx.gz')
                    is_gz = True
                else:
                    raise Exception('data not found')
                self.is_gz = is_gz
                if gene_paths is None:
                    gene_path = os.path.join(base_path, 'gene_names.txt')
                else:
                    gene_path = gene_paths[0]
            else:
                self.is_gz = True
                data_paths_new = []
                for i, data_path in enumerate(data_paths):
                    data_path = str(data_path)
                    is_gz = data_path.endswith('gz')
                    data_path_new = data_path
                    # convert data shape
                    # TODO: try to automatically infer the shape (if gene_names.txt is present)
                    # The problem is inefficiency - we don't want to have to update the dataset twice.
                    if shapes[i] == 'cell_gene':
                        try:
                            data = scipy.io.mmread(data_path)
                        except:
                            data = np.loadtxt(data_path)
                            data_path_new = data_path[:-4] + '.mtx'
                        os.remove(data_path)
                        data = data.T
                        if is_gz:
                            data_path_new = data_path_new[:-3]
                        scipy.io.mmwrite(data_path_new, data)
                        if is_gz:
                            import subprocess
                            subprocess.call(['gzip', data_path_new])
                            data_path_new += '.gz'
                    else:
                        # always convert data to mtx
                        if data_path.endswith('.txt'):
                            data = np.loadtxt(data_path)
                            os.remove(data_path)
                            data_path_new = data_path[:-4] + '.mtx'
                            scipy.io.mmwrite(data_path_new, data)
                            import subprocess
                            subprocess.call(['gzip', data_path_new])
                            data_path_new += '.gz'
                    data_paths_new.append(data_path_new)
                # call merge_datasets
                from uncurl_analysis import merge_datasets
                data_path, gene_path = merge_datasets.merge_files(data_paths_new, gene_paths, dataset_names, base_path, use_batch_correction=use_batch_correction)
            try:
                data = scipy.io.mmread(data_path)
            except:
                data = np.loadtxt(data_path)
        else:
            self.is_gz = True
        self.data = sparse.csc_matrix(data)
        print(gene_path)
        print(gene_paths)
        try:
            self.gene_names = np.loadtxt(gene_path, dtype=str)
        except:
            self.gene_names = np.array([str(x) for x in range(data.shape[0])])
        self.cell_read_counts = np.array(data.sum(0)).flatten()
        if sparse.issparse(data):
            self.cell_gene_counts = data.getnnz(0)
        else:
            self.cell_gene_counts = np.count_nonzero(data, 0)
        self.sorted_read_counts = np.sort(self.cell_read_counts)
        self.sorted_gene_counts = np.sort(self.cell_gene_counts)
        self.cells = data.shape[1]
        self.genes = data.shape[0]
        self.is_sparse = True
        m, v = sparse_means_var_csc(self.data.data,
                self.data.indices, self.data.indptr, self.data.shape[1],
                self.data.shape[0])
        self.gene_means = m
        self.gene_vars = v
        self.sorted_gene_means = np.sort(self.gene_means)
        self.path = base_path

    def summary(self):
        return (self.cells, self.genes)

    def preprocessing_params(self):
        """
        Saves preprocessing parameters as 'preprocess.json'

        params: min_reads (bottom 10th percentile), max_reads
        (top 10th percentile), frac (0.2), nbins (5)
        """
        # cell_frac is set so that there will be a max of 25000 points
        cell_frac = min(1.0, 25000.0/self.cells)
        cell_frac = round(cell_frac, 2)
        top_05 = int(self.cells/20) # 5%
        preproc_params = {'min_reads': int(self.sorted_read_counts[top_05]),
                          'max_reads': int(self.sorted_read_counts[-top_05]),
                          'min_unique_genes': 0,
                          'max_unique_genes': 20000,
                          'max_mt_frac': 1.0,
                          'genes_frac': 0.2,
                          'nbins': 5,
                          'cell_frac': cell_frac,
                          'cells': self.cells,
                          'genes': self.genes,
                          'normalize': False,
                          'is_gz': self.is_gz,
                          'is_sparse': self.is_sparse,
                          'median_reads': np.median(self.cell_read_counts),
                          }
        with open(os.path.join(self.path, 'preprocess.json'), 'w') as f:
            json.dump(preproc_params, f, cls=SimpleEncoder)
        return preproc_params

    def generate_plotly_jsons(self):
        """
        Generate 3 plots: read count histogram, gene count histogram, gene mean expression histogram
        Saves three files: read_count_hist_data.json, gene_count_hist_data.json, gene_mean_hist_data.json
        Returns the results as 3 json-formatted strings.
        """
        top_05 = int(self.cells/20) # 5%
        # TODO: separate bulk data from outliers
        read_counts_max = self.sorted_read_counts[-top_05]
        read_count_hist_data = json.dumps({
             'data': [{
                'x': self.sorted_read_counts[:-top_05].tolist(),
                'type': 'histogram',
                #'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Read counts',
                'marker': {'color': 'blue'},
                'nbinsx': 50,
            }, {
                'x': self.sorted_read_counts[top_05:].tolist(),
                'type': 'histogram',
                #'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Read counts (outliers)',
                'marker': {'color': 'blue'},
                'xbins': {'start': read_counts_max, 'end': self.sorted_read_counts[-1], 'size': int(self.sorted_read_counts[-1]/100)},
            }
            ],
            'layout': {
                'title': 'Read counts per cell',
                'barmode': 'overlay',
                'showlegend': False,
                'xaxis': {'title': 'Read/UMI Count', 'range': [0.0, read_counts_max]},
                'yaxis': {'title': 'Cell Counts'},
            },
        }, cls=SimpleEncoder)
        with open(os.path.join(self.path, 'read_count_hist_data.json'), 'w') as f:
            f.write(read_count_hist_data)
        gene_count_max = self.sorted_gene_counts[-top_05]
        gene_count_hist_data = json.dumps({
             'data': [{
                'x': self.sorted_gene_counts[:-top_05].tolist(),
                'type': 'histogram',
                #'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Gene counts',
                'marker': {'color': 'blue'},
                'nbinsx': 50,
            }, {
                'x': self.sorted_gene_counts[top_05:].tolist(),
                'type': 'histogram',
                #'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Gene counts (outliers)',
                'marker': {'color': 'blue'},
                'xbins': {'start': gene_count_max, 'end': self.sorted_gene_counts[-1], 'size': int(self.sorted_gene_counts[-1]/100)},
            }
            ],
            'layout': {
                'title': 'Unique gene counts per cell',
                'barmode': 'overlay',
                'showlegend': False,
                'xaxis': {'title': 'Gene Count', 'range': [0, gene_count_max]},
                'yaxis': {'title': 'Cell Counts'},
            },
        }, cls=SimpleEncoder)
        with open(os.path.join(self.path, 'gene_count_hist_data.json'), 'w') as f:
            f.write(gene_count_hist_data)
        # plot mtRNA frac as a histogram, no need to plot gene means
        mt_genes = map(lambda x: x.startswith('Mt-') or x.startswith('MT-') or x.startswith('mt-'), self.gene_names)
        mt_genes = np.array(list(mt_genes))
        if len(mt_genes) > 0:
            mt_gene_counts = np.array(self.data[mt_genes, :].sum(0)).flatten()
            mt_gene_frac = mt_gene_counts/self.cell_read_counts
            mt_frac_hist_data = json.dumps({
                 'data': [{
                    'x': mt_gene_frac.tolist(),
                    'type': 'histogram',
                    'opacity': 1.0,
                    'name': 'Gene means',
                    'marker': {'color': 'blue'},
                }],
                'layout': {
                    'title': 'Fraction of mitochondrial genes per cell',
                    'barmode': 'overlay',
                    'showlegend': False,
                    'xaxis': {'title': 'Mt Frac'},
                    'yaxis': {'title': 'Cell Counts'},
                },
            }, cls=SimpleEncoder)
            with open(os.path.join(self.path, 'gene_mean_hist_data.json'), 'w') as f:
                f.write(mt_frac_hist_data)
        else:
            mt_frac_hist_data = None
        return read_count_hist_data, gene_count_hist_data, mt_frac_hist_data

    def load_plotly_json(self):
        """
        if json files already exist, load them. else, generate them.
        """
        if os.path.exists(os.path.join(self.path, 'read_count_hist_data.json')):
            with open(os.path.join(self.path, 'read_count_hist_data.json')) as f:
                read_count_hist_data = f.read()
            with open(os.path.join(self.path, 'gene_count_hist_data.json')) as f:
                gene_count_hist_data = f.read()
            with open(os.path.join(self.path, 'gene_mean_hist_data.json')) as f:
                gene_mean_hist_data = f.read()
            return read_count_hist_data, gene_count_hist_data, gene_mean_hist_data
        else:
            return self.generate_plotly_jsons()

