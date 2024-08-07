"""
Cirle packed plot.
"""

import os
import random
import pickle
import pandas as pd
from itertools import product
from plotting_utils._plotting import *
from utils.plotting import *
matplotlib.use('macOSX')


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


##


# Colors, for each clone sample
df['GBC'].unique().size

# Random colors for clones
# clones = df['GBC'].unique()
# random.seed(1234)               # 1222 MM13, 1234 MM23
# clones_colors = { 
#     clone : color for clone, color in \
#     zip(
#         clones, 
#         list(
#             ''.join( ['#'] + [random.choice('ABCDEF0123456789') for i in range(6)] )  \
#             for _ in range(clones.size)
#         )
#     )
# }
# with open(os.path.join(path_data, f'clones_colors_{dataset}.pickle'), 'wb') as f:
#     pickle.dump(clones_colors, f)

# Read colors
with open(os.path.join(path_data, f'clones_colors_{dataset}.pickle'), 'rb') as f:
    clones_colors = pickle.load(f)


##


# Fig
combos = list(product(df['mouse'].unique(), df['origin'].value_counts().index))

fig, axs = plt.subplots(df['mouse'].unique().size,df['origin'].unique().size,figsize=(11,8))

for combo, ax in zip(combos, axs.ravel()):
    mouse, origin = combo
    max_n = 10 if origin != 'PT' else 3
    df_ = df.query('n>=@max_n and mouse==@mouse and origin==@origin').set_index('GBC')
    packed_circle_plot(
        df_, covariate='freq', ax=ax, color=clones_colors, annotate=True, t_cov=.05,
        alpha=.65, linewidth=2.5, fontsize=8, fontcolor='k', fontweight='medium'
    )
    df_sample = meta.query('mouse==@mouse and origin==@origin')
    ax.set(title=f'{origin}: {df_sample.shape[0]} cells, {df_sample["GBC"].nunique()} clones')
    ax.axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{dataset}.png'), dpi=1000)


##