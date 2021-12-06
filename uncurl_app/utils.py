import json

import numpy as np

# to encode numpy stuff...
class SimpleEncoder(json.JSONEncoder):

    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        elif isinstance(o, np.ndarray):
            return o.tolist()
        return json.JSONEncoder.default(self, o)


def get_matrix_header(filename):
    """
    Returns the entries, rows, and cols of a matrix market file.
    """
    with open(filename) as f:
        entries = 0
        rows = 0
        cols = 0
        for line in f.readlines():
            if line.startswith('%'):
                continue
            line = line.split()
            entries = int(line[0])
            rows = int(line[1])
            cols = int(line[2])
        return entries, rows, cols
