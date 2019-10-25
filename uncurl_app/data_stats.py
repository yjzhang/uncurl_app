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

    def __init__(self, data_path, base_path, is_gz=False, shape='gene_cell', data=None):
        """
        Args:
            data_path (str): path to data file
            base_path (str): path to data dir
        """
        # a python2 thing for dealing with unicode...
        if data is None:
            data_path = str(data_path)
            data_is_sparse = True
            try:
                data = scipy.io.mmread(data_path)
            except:
                data = np.loadtxt(data_path)
                data_is_sparse = False
            if shape == 'cell_gene':
                os.remove(data_path)
                data = data.T
                if data_is_sparse:
                    if is_gz:
                        data_path = data_path[:-3]
                    scipy.io.mmwrite(data_path, data)
                    if is_gz:
                        import subprocess
                        subprocess.call(['gzip', data_path])
                else:
                    np.savetxt(data_path, data)
                # save data...
        self.cell_read_counts = np.array(data.sum(0)).flatten()
        self.cell_gene_counts = np.array((data>0).sum(0)).flatten()
        self.sorted_read_counts = np.sort(self.cell_read_counts)
        self.sorted_gene_counts = np.sort(self.cell_gene_counts)
        self.cells = data.shape[1]
        self.genes = data.shape[0]
        if sparse.issparse(data):
            self.is_sparse = True
            data_csc = sparse.csc_matrix(data)
            m, v = sparse_means_var_csc(data_csc.data,
                    data_csc.indices, data_csc.indptr, data.shape[1],
                    data.shape[0])
            self.gene_means = m
            self.gene_vars = v
        else:
            self.is_sparse = False
            self.gene_means = data.mean(1).flatten()
            self.gene_vars = data.var(1).flatten()
        self.sorted_gene_means = np.sort(self.gene_means)
        self.path = base_path
        self.is_gz = is_gz

    def summary(self):
        return (self.cells, self.genes)

    def preprocessing_params(self):
        """
        Saves preprocessing parameters as 'preprocess.json'

        params: min_reads (bottom 10th percentile), max_reads
        (top 10th percentile), frac (0.2), nbins (5)
        """
        # cell_frac is set so that there will be a max of 10000 points
        cell_frac = min(1.0, 10000.0/self.cells)
        cell_frac = round(cell_frac, 2)
        top_05 = int(self.cells/20) # 5%
        preproc_params = {'min_reads': int(self.sorted_read_counts[top_05]),
                          'max_reads': int(self.sorted_read_counts[-top_05]),
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
        read_counts_max = self.sorted_read_counts[-top_05]
        read_count_hist_data = json.dumps({
             'data': [{
                'x': self.cell_read_counts.tolist(),
                'type': 'histogram',
                'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Read counts',
                'marker': {'color': 'blue'},
                'nbinsx': max(100, int(self.sorted_read_counts[-1]/100)),
            }],
            'layout': {
                'title': 'Read counts per cell',
                'barmode': 'overlay',
                'xaxis': {'title': 'Read/UMI Count', 'range': [0.0, read_counts_max]},
                'yaxis': {'title': 'Cell Fraction'},
            },
        }, cls=SimpleEncoder)
        with open(os.path.join(self.path, 'read_count_hist_data.json'), 'w') as f:
            f.write(read_count_hist_data)
        gene_count_max = self.sorted_gene_counts[-top_05]
        gene_count_hist_data = json.dumps({
             'data': [{
                'x': self.cell_gene_counts.tolist(),
                'type': 'histogram',
                'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Gene counts',
                'marker': {'color': 'blue'},
                'nbinsx': max(100, int(self.sorted_gene_counts[-1]/100)),
            }],
            'layout': {
                'title': 'Gene counts per cell',
                'barmode': 'overlay',
                'xaxis': {'title': 'Gene Count', 'range': [0, gene_count_max]},
                'yaxis': {'title': 'Cell Fraction'},
            },
        }, cls=SimpleEncoder)
        with open(os.path.join(self.path, 'gene_count_hist_data.json'), 'w') as f:
            f.write(gene_count_hist_data)
        top_10 = int(self.cells/10) # 10%
        gene_means_max = self.sorted_gene_means[-top_10]
        gene_mean_hist_data = json.dumps({
             'data': [{
                'x': self.gene_means.tolist(),
                'type': 'histogram',
                'histnorm': 'probability',
                'opacity': 1.0,
                'name': 'Gene means',
                'marker': {'color': 'blue'},
            }],
            'layout': {
                'title': 'Mean expression per gene',
                'barmode': 'overlay',
                'showlegend': True,
                'xaxis': {'title': 'Gene Mean Expression Level', 'range': [0.0, gene_means_max]},
                'yaxis': {'title': 'Gene Fraction'},
            },
        }, cls=SimpleEncoder)
        with open(os.path.join(self.path, 'gene_mean_hist_data.json'), 'w') as f:
            f.write(gene_mean_hist_data)
        return read_count_hist_data, gene_count_hist_data, gene_mean_hist_data

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

