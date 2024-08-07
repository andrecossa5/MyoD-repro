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

# Bulk
bulk = pd.read_csv(os.path.join(path_data, 'common.csv'), index_col=0)
# for mice in ['1','2','3']:
#     bulk.loc[bulk.index.str.endswith(mice), bulk.index.str.endswith(mice)].loc[f'13_PT{mice}'].sort_values()


##


# SC
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0).iloc[:,:4]
df_ = meta.groupby(['mouse', 'origin', 'GBC']).size().to_frame('n_cells').reset_index()
df_.query('n_cells>=5').groupby(['mouse', 'origin'])['GBC'].nunique().to_frame('n>5').reset_index()
df_.query('n_cells>=5').groupby(['mouse', 'GBC'])['origin'].nunique().to_frame('n_sites').reset_index().query('n_sites>2').groupby('mouse')['GBC'].nunique()


##


# meta['GBC'].nunique()
# df_.query('GBC=="CAATATGTCCTGCAGGCA" and mouse=="M1"')