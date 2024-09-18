"""
Essential bulk analysis script.
"""

import os
import re
import numpy as np
import pandas as pd
from plotting_utils._plotting import *


## 


# Paths
path_main = '/Users/IEO5505/Desktop/MyoD-repro/'
dataset = 'MM23'
path_data = os.path.join(path_main, 'data', dataset)
path_results = os.path.join(path_main, 'results', 'clonal')


##


# SC
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0).iloc[:,:4]
df = (
    meta.groupby(['mouse', 'origin', 'GBC'])
    .size().to_frame('n')
    .reset_index()
    .assign(freq=lambda x: x['n'] / x.groupby(['mouse', 'origin'])['n'].transform('sum'))
)

# Calculate longitudinal clones statistics
grouped = df.groupby('mouse')
L = []
for mouse, df in grouped:
    origins = df['origin'].unique()
    L.append((
        df[['GBC', 'origin', 'n']]
        .pivot(index='GBC', columns='origin', values='n').fillna(0).astype(int)
        .assign(
            n_cells=lambda x: x.sum(axis=1),
            n_sites=lambda x: (x.loc[:,x.columns.isin(origins)]>0).sum(axis=1),
            median_n_cells_per_site=lambda x: x.iloc[:,:4].median(axis=1),
            mouse=mouse,
            filtered=lambda x: (x['median_n_cells_per_site']>=5) & (x['n_sites']>1)
        )
        .query('n_cells>=50 and n_sites>1')
        .sort_values('median_n_cells_per_site', ascending=False)
    ))
pd.concat(L).to_csv(os.path.join(path_results, f'{dataset}_clone_selection.csv'))


##

