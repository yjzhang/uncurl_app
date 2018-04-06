import os
import pickle

import numpy as np
import scipy.io
import uncurl

from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE, MDS

import simplex_sample

class SCAnalysis(object):
    """
    This class represents an ongoing single-cell RNA-Seq analysis.
    """

    def __init__(self, data_dir,
            data_filename='data.mtx',
            data_is_sparse=True,
            normalize=False,
            min_reads=0,
            max_reads=1e10,
            frac=0.2,
            cell_frac=1.0,
            dim_red_option='mds',
            baseline_dim_red='none',
            uncurl_kwargs={}):
        """
        Args:
            data_dir (str): directory where data is stored
        """
        # note: each field contains file names, and whether or not
        # the analysis is complete.
        self.data_dir = data_dir
        self.min_reads = min_reads
        self.max_reads = max_reads
        self.normalize = normalize
        self.is_sparse = data_is_sparse
        self.data_f = os.path.join(data_dir, data_filename)
        self._data = None
        self._data_subset = None

        self.gene_names_f = os.path.join(data_dir, 'gene_names.txt')
        self._gene_names = None

        self.frac = frac
        self.gene_subset_f = os.path.join(data_dir, 'gene_subset.txt')
        self.has_gene_subset = False
        self._gene_subset = None

        self.uncurl_kwargs = uncurl_kwargs

        self.w_f = os.path.join(data_dir, 'w.txt')
        self.has_w = False
        self._w = None

        self.m_f = os.path.join(data_dir, 'm.txt')
        self.has_m = False
        self._m = None

        self.cell_subset_f = os.path.join(data_dir, 'cells_subset.txt')
        self.has_cell_subset = False
        self._cell_subset = None

        self.cell_frac = cell_frac
        self.cell_sample_f = os.path.join(data_dir, 'cell_sample.txt')
        self.has_cell_sample = False
        self._cell_sample = None

        self.baseline_dim_red = baseline_dim_red
        self.baseline_vis_f = os.path.join(data_dir, 'baseline_vis.txt')
        self.has_baseline_vis = False
        self._baseline_vis = None

        self.dim_red_option = dim_red_option
        self.dim_red_f = os.path.join(data_dir, 'mds_data.txt')
        self.has_dim_red = False
        self._dim_red = None

        self.mds_means_f = os.path.join(data_dir, 'mds_means.txt')
        self.has_mds_means = False
        self._mds_means = None

        self.pickle_f = os.path.join(data_dir, 'sc_analysis.pkl')


    @property
    def data(self):
        """
        Data - either a sparse csc matrix or a dense numpy array, of shape
        (genes, cells)
        """
        if self._data is None:
            try:
                if self.is_sparse:
                    self._data = scipy.io.mmread(self.data_f)
                else:
                    self._data = np.loadtxt(self.data_f)
                return self._data
            except:
                return None
        else:
            return self._data

    @property
    def gene_subset(self):
        """
        Gene subset (from max_variance_genes)
        """
        if self._cell_subset is None:
            if not self.has_cell_subset:
                data = self.data
                gene_subset = uncurl.max_variance_genes(data, nbins=5, frac=self.frac)
                np.savetxt(self.gene_subset_f, gene_subset)
                self.has_gene_subset = True
            else:
                gene_subset = np.loadtxt(self.gene_subset_f)
            self._gene_subset = gene_subset
            return gene_subset
        else:
            return self._gene_subset

    @property
    def cell_subset(self):
        """
        Cell subset (array of booleans, based on read count filtering)
        """
        if self._cell_subset is None:
            if not self.has_cell_subset:
                data = self.data
                read_counts = np.array(data.sum(0)).flatten()
                cell_subset = (read_counts >= self.min_reads) & (read_counts <= self.max_reads)
                np.savetxt(self.cell_subset_f, cell_subset)
                self.has_cell_subset = True
            else:
                cell_subset = np.loadtxt(self.cell_subset_f, dtype=bool)
            self._cell_subset = cell_subset
            return cell_subset
        else:
            return self._cell_subset

    @property
    def data_subset(self):
        """
        Data after passed through the gene/cell filters
        """
        if self._data_subset is None:
            data = self.data
            data_subset = data[self.gene_subset, :]
            data_subset = data_subset[:, self.cell_subset]
            self._data_subset = data_subset
        return self._data_subset


    @property
    def gene_names(self):
        """
        Array of gene names
        """
        if self._gene_names is None:
            try:
                self._gene_names = np.loadtxt(self.gene_names_f, dtype=str)
                return self._gene_names
            except:
                return None
        else:
            return self._gene_names

    def run_uncurl(self):
        """
        Runs uncurl on self.data_subset.
        """
        m, w, ll = uncurl.run_state_estimation(self.data_subset,
                **self.uncurl_kwargs)
        np.savetxt(self.w_f, w)
        np.savetxt(self.m_f, m)
        self._m = m
        self._w = w
        self.has_w = True
        self.has_m = True

    @property
    def m(self):
        if self._m is None:
            if self.has_m:
                m = np.loadtxt(self.m_f)
                self._m = m
            else:
                self.run_uncurl()
        return self._m

    @property
    def w(self):
        if self._w is None:
            if self.has_w:
                w = np.loadtxt(self.w_f)
                self._w = w
            else:
                self.run_uncurl()
        return self._w

    @property
    def mds_means(self):
        """
        MDS of the post-uncurl cluster means
        """
        if self._mds_means is None:
            if self.has_mds_means:
                self._mds_means = np.loadtxt(self.mds_means_f)
            else:
                self._mds_means = uncurl.dim_reduce(self.m, self.w, 2)
                np.savetxt(self.mds_means_f, self._mds_means)
                self.has_mds_means = True
        return self._mds_means

    @property
    def cell_sample(self):
        """
        Cell sample (after applying data subset)
        """
        if self.cell_frac == 1:
            self._cell_sample = np.arange(self.data_subset.shape[1])
        else:
            if self._cell_sample is None:
                if self.has_cell_sample:
                    self._cell_sample = np.loadtxt(self.cell_sample_f, dtype=int)
                else:
                    k, cells = self.w.shape
                    n_samples = int(cells*self.cell_frac)
                    samples = simplex_sample.sample(k, n_samples)
                    indices = simplex_sample.data_sample(self.w, samples)
                    np.savetxt(self.cell_sample_f, indices)
                    self.has_cell_sample = True
                    self._cell_sample = indices
        return self._cell_sample

    @property
    def data_sampled(self):
        """
        Data after passed through the gene/cell filters, and sampled.
        """
        data_subset = self.data_subset
        cell_sample = self.cell_sample
        return data_subset[:, cell_sample]

    @property
    def baseline_vis(self):
        # TODO
        baseline_dim_red = self.baseline_dim_red.lower()
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

    def save_pickle_reset(self):
        """
        Removes all cached data, saves to pickle
        """
        # TODO
        self._data = None
        self._data_subset = None
        self._gene_names = None
        self._gene_subset = None
        self._w = None
        self._m = None
        self._cell_subset = None
        self._baseline_vis = None
        self._dim_red = None


    # TODO: uncurl re-initialization
