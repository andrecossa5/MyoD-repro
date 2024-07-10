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
path_bulk = os.path.join(path_main, 'data', 'MM13')      # Or other bulk output folder
path_sc = os.path.join(path_main, 'data', 'meta')  
path_results = os.path.join(path_main, 'results', 'clonal')

# Bulk
bulk = pd.read_csv(os.path.join(path_bulk, 'common.csv'), index_col=0)
# for mice in ['1','2','3']:
#     bulk.loc[bulk.index.str.endswith(mice), bulk.index.str.endswith(mice)].loc[f'13_PT{mice}'].sort_values()


##


# SC
meta = pd.read_csv(os.path.join(path_sc, 'cells_meta.csv'), index_col=0).iloc[:,:4]
df_ = meta.groupby(['mouse', 'origin', 'GBC']).size().to_frame('n_cells').reset_index()
df_.query('n_cells>=5').groupby(['mouse', 'origin'])['GBC'].nunique().to_frame('n>5').reset_index()
df_.query('n_cells>=5').groupby(['mouse', 'GBC'])['origin'].nunique().to_frame('n_sites').reset_index().query('n_sites>1').groupby('mouse')['GBC'].nunique()


##


# meta['GBC'].nunique()
# df_.query('GBC=="CAATATGTCCTGCAGGCA" and mouse=="M1"')