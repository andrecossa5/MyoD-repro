"""
Tables for Daria.
"""

import os
import numpy as np
import pandas as pd


## 


# Paths
path_main = '/Users/IEO5505/Desktop/old/MyoD-repro/'
dataset = 'MM23'
path_data = os.path.join(path_main, 'data', dataset)
path_results = os.path.join(path_main, 'results', 'clonal')


##


# SC
cell_state_df = pd.read_csv(os.path.join(path_data, f'clusters_{dataset}.csv'), index_col=0)
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0).iloc[:,:4]
df_freq = (
    meta.groupby(['mouse', 'origin', 'GBC'])
    .size().to_frame('n')
    .reset_index()
    .pivot_table(values='n', index=['GBC', 'mouse'], columns='origin')
    .fillna(0)
    .reset_index()
    .sort_values('mouse')
)
df_freq.to_csv(os.path.join(path_results, f'clone_frequency_all_{dataset}.csv'))

##

# Enrichments
from utils.tests import compute_enrichment

# Add cell state info
df = meta.join(cell_state_df[['cell_states']])

L = []
for t in df['cell_states'].unique():
    L.append(
        compute_enrichment(df, col1='GBC', col2='cell_states', target=t)
        [['group', 'target', 'enrichment']]
    )
(
    pd.concat(L)
    .pivot(index='group', columns='target', values='enrichment')
    .to_csv(os.path.join(path_results, f'enrichments_cell_states_{dataset}.csv'))
)


##


