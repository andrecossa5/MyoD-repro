"""
Fix cells_meta for Cellula analysis.
"""

import os
import numpy as np
import pandas as pd


## 


# Paths
path_main = '/Users/IEO5505/Desktop/MyoD-repro/'
path_meta = os.path.join(path_main, 'data', 'meta')

# Read and format cells meta
meta = pd.read_csv(os.path.join(path_meta, 'cells_meta_orig.csv'), index_col=0)

meta['sample'].unique()

# Origin and dataset
tests = [meta['sample'].str.contains('PT'),  meta['sample'].str.contains('Lu'),  meta['sample'].str.contains('Li'), meta['sample'].str.contains('LN')]
choices = ['PT', 'Lu', 'Li', 'LN']
meta['origin'] = np.select(tests, choices)
meta['mouse'] = meta['sample'].map(lambda x: f'M{x[-1]}')

# Reorder cols
meta = meta[['GBC', 'sample', 'origin', 'mouse', 'nUMIs', 'mito_perc', 'detected_genes',
      'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes']]

# Save
meta.to_csv(os.path.join(path_meta, 'cells_meta.csv'))


##