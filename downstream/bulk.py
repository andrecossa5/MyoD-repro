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
run = 'MM23'
path_main = '/Users/IEO5505/Desktop/MyoD-repro/'
path_data = os.path.join(path_main, 'data', run)      # Or other bulk output folder
path_results = os.path.join(path_main, 'results', 'clonal')

# Read prevalences
df = pd.read_csv(os.path.join(path_data, 'bulk_GBC_reference.csv'), index_col=0)
common = pd.read_csv(os.path.join(path_data, 'common.csv'), index_col=0)

# Reformat
tests = [ df['sample'].str.contains('PT'), df['sample'].str.contains('LN'), df['sample'].str.contains('Kid') ] 
df['origin'] = np.select(tests, ['PT', 'Li', 'KD'], default='Ref')


##


# n clones
df_ = (
    df.groupby('sample')
    .apply(lambda x: x.index.unique().size)
    .sort_values(ascending=False)
    .to_frame('n')
    .reset_index()
)
fig, ax = plt.subplots(figsize=(8,4.5))
bar(df_, 'n', 'sample', s=.75, c='k', a=.7, ax=ax)
format_ax(ax=ax, title='n clones by sample', ylabel='n', xticks=df_['sample'], rotx=90)
ax.spines[['left', 'top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'n_clones_{run}.png'), dpi=500)


##


# Cumulative clone percentage, all samples
colors = create_palette(df, 'origin', ten_godisnot)

fig, ax = plt.subplots(figsize=(4.5,4.5))
for s in df['sample'].unique():
    df_ = df.query('sample==@s')
    x = (df_['read_count'] / df_['read_count'].sum()).cumsum()
    origin = df.query('sample==@s')['origin'].unique()[0]
    ax.plot(range(len(x)), x, c=colors[origin], linewidth=2.5)

ax.set(title='Clone prevalences', xlabel='Ranked clones', ylabel='Cumulative frequence')
add_legend(ax=ax, colors=colors, bbox_to_anchor=(1,0), loc='lower right', ticks_size=8, label_size=10, artists_size=8)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'cum_percentages_{run}.png'), dpi=300)


##


# Shannon entropies
SH = []
for s in df['sample'].unique():
    df_ = df.query('sample==@s')
    x = df_['read_count'] / df_['read_count'].sum()
    SH.append(-np.sum( np.log10(x) * x ))
df_ = (pd.Series(SH, index=df['sample'].unique())
    .to_frame('SH')
    .sort_values(by='SH', ascending=False)
    .reset_index().rename(columns={'index':'sample'})
    .merge(df[['sample', 'origin']], on='sample')
    .drop_duplicates()
    .set_index('sample')
)

fig, ax = plt.subplots(figsize=(4,4))
box(df_, x='origin', y='SH', ax=ax, with_stats=True, 
    pairs=[['PT', 'KD'], ['KD', 'Li'], ['Li', 'PT']], 
    order=['Ref', 'PT', 'Lu', 'Li']
)
strip(df_, x='origin', y='SH', ax=ax, order=['Ref', 'PT', 'Lu', 'Li'], c='k')
format_ax(ax=ax, title='Shannon Entropy samples', ylabel='SH', rotx=90, reduce_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'SH_{run}.png'), dpi=300)


##


# Commmon clones
fig, ax = plt.subplots(figsize=(5,5))
plot_heatmap(common, ax=ax, annot=True, title='n common clones')
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'common_{run}.png'), dpi=300)


##