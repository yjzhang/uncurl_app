# how do we deal with the input of split-seq?
import os
import subprocess
import sys

import pandas as pd
import scipy.io
from scipy import sparse

from . import data_stats

def process_split_seq(path):
    """
    Modifies the data so that it can be used in the uncurl_app pipeline.
    """
    # 1. cp DGE.mtx data.mtx
    subprocess.call(['cp', os.path.join(path, 'DGE.mtx'), os.path.join(path, 'data.mtx')])
    # 2. gzip data.mtx
    subprocess.call(['gzip', os.path.join(path, 'data.mtx')])
    # 3. load genes.csv, write out gene_names.txt
    genes = pd.read_csv(os.path.join(path, 'genes.csv'))
    gene_names = genes.gene_name
    gene_names.to_csv(os.path.join(path, 'gene_names.txt'), header=None, index=None)
    # 4. run the data_stats stuff
    summary = data_stats.Summary(os.path.join(path, 'data.mtx.gz'), path, is_gz=True, shape='cell_gene')
    script, div = summary.visualize()
    summary.preprocessing_params()

def add_cell_metadata_color_track(path):
    # create a color track from the round 1 wells in cell_metadata.csv
    # TODO
    cell_metadata = pd.read_csv('cell_metadata.csv')
    round1_wells = cell_metadata.rnd1_well
