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
path_PT = os.path.join(path_main, 'data', 'bulk')      # Or other bulk output folder
path_mets = os.path.join(path_main, 'data', 'bulk_metastases_Sept2023')      # Or other bulk output folder
path_results = os.path.join(path_main, 'results', 'clonal')

# Read prevalences, PT and mets
PTs = pd.read_csv(os.path.join(path_PT, 'summary', 'bulk_GBC_reference.csv'), index_col=0)
mets = pd.read_csv(os.path.join(path_mets, 'summary', 'all_prevalences.csv'), index_col=0)

# Reformat
PTs = PTs.query('sample!="Ref"')
mets = mets.loc[lambda x: x['found_wo']]
mets = mets[['read_count', 'sample']]

# Calculate common for three sets
L = []
for pattern in ['1', '2', '3']: #pattern = '1'
    PT_clones = set(PTs.loc[PTs['sample'].str.contains(pattern)].index)
    mets_ = mets.loc[mets['sample'].str.contains(pattern)]
    n = []
    for m in mets_['sample'].unique():
        n.append(len(set(mets_.query('sample==@m').index) & PT_clones))
    L.append(pd.Series(n, index=mets_['sample'].unique()))

# Save
pd.concat(L, axis=0).to_frame('n_clones common with PT').to_csv(os.path.join(path_results, 'n_longitudinal.csv'))


##