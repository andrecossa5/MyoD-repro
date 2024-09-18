"""
Essential bulk analysis script.
"""

import os
import numpy as np
import pandas as pd
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula._utils import *
matplotlib.use('MacOSX')


## 


# Paths
dataset = 'MM23'
path_main = '/Users/IEO5505/Desktop/MyoD-repro/'
path_data = os.path.join(path_main, 'data', dataset)

# Read single-cell data and compute clonal prevalences
sc_meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
freq_sc = (
    sc_meta
    .groupby('sample')
    ['GBC'].value_counts(normalize=True).loc[lambda x:x>0]
    .reset_index(name='freq_sc')
)

# Read bulk data and confront with sc
freq_bulk = pd.read_csv(os.path.join(path_data, 'bulk_GBC_reference.csv'), index_col=0)
freq_bulk = ( 
    freq_bulk
    .groupby('sample')
    .apply(lambda x: x['read_count']/x['read_count'].sum())
    .reset_index(name='freq_bulk')
    .rename(columns={'level_1':'GBC'})
)

# Merge
df = freq_bulk.merge(freq_sc, on=['sample', 'GBC'], how='outer')

# Compute numbers
# freq_bulk.groupby('sample')['GBC'].nunique()
# freq_sc.groupby('sample')['GBC'].nunique()
# df.dropna().groupby('sample')['GBC'].nunique()
# df.dropna().groupby('sample').apply(lambda x: np.corrcoef(x['freq_bulk'], x['freq_sc'])[0,1])


##