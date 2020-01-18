import os
import json

for path, dirs, files in os.walk('.'):
    if 'color_tracks.json' in files:
        with open(os.path.join(path, 'color_tracks.json')) as f:
            data = json.load(f)
        for k, d in data.items():
            try:
                del d['1_vs_rest_scores']
                del d['1_vs_rest_pvals']
                del d['pairwise_pvals']
                del d['pairwise_scores']
            except:
                pass
        with open(os.path.join(path, 'color_tracks.json'), 'w') as f:
            json.dump(data, f)
