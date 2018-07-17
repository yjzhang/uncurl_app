import json
import os

from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.layouts import gridplot

import numpy as np
import scipy.io
from scipy import sparse
from uncurl.sparse_utils import sparse_means_var_csc

from .utils import SimpleEncoder

class Summary(object):
    """
    This object contains a summary for a single-cell RNASeq dataset.
    """

    def __init__(self, data, path, is_gz=False, is_data=False):
        if not is_data:
            try:
                data = scipy.io.mmread(data)
            except:
                data = np.loadtxt(data)
        self.cell_read_counts = np.array(data.sum(0)).flatten()
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
        self.path = path
        self.is_gz = is_gz
        # TODO: try some k selection?

    def summary(self):
        return (self.cells, self.genes)

    def preprocessing_params(self):
        """
        Saves preprocessing parameters as 'preprocess.json'

        params: min_reads (bottom 10th percentile), max_reads
        (top 10th percentile), frac (0.2), nbins (5)
        """
        # cell_frac is set so that there will be 500 points?
        cell_frac = min(1.0, 1500.0/self.cells)
        cell_frac = round(cell_frac, 2)
        sorted_read_counts = np.sort(self.cell_read_counts)
        top_05 = int(self.cells/20) #5%
        preproc_params = {'min_reads': int(sorted_read_counts[top_05]),
                          'max_reads': int(sorted_read_counts[-top_05]),
                          'frac': 0.2,
                          'nbins': 5,
                          'cell_frac': cell_frac,
                          'cells': self.cells,
                          'genes': self.genes,
                          'normalize': False,
                          'is_gz': self.is_gz,
                          'is_sparse': self.is_sparse,
                          }
        with open(os.path.join(self.path, 'preprocess.json'), 'w') as f:
            json.dump(preproc_params, f, cls=SimpleEncoder)
        return preproc_params

    def visualize(self):
        # alternative view: use an iframe?
        p1 = figure(title='Histogram of cell read counts', plot_width=400,
                plot_height=400)
        hist, edges = np.histogram(self.cell_read_counts, bins=50)
        p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
        p1.xaxis.axis_label = 'Read count'
        p1.yaxis.axis_label = '# cells'
        p2 = figure(title='Gene variance vs mean', plot_width=400,
                plot_height=400)
        p2.circle(x=self.gene_means, y=self.gene_vars, alpha=0.5, size=5)
        p2.xaxis.axis_label = 'Mean'
        p2.yaxis.axis_label = 'Variance'
        plots = gridplot([[p1, p2]])
        script, div = components(plots)
        with open(os.path.join(self.path, 'vis_summary.html'), 'w') as f:
            f.write(script)
            f.write('\n\n')
            f.write(div)
        return script, div

