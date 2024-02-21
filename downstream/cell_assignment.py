"""
Final cell_assignment first MyoD samples.
"""

import os
import sys
import pickle
import pandas as pd
from plotting_utils._plotting_base import *
sys.path.append('/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc')
from helpers import *


##


# Args
path_main = '/Users/IEO5505/Desktop/MyoD-repro'
path_data = os.path.join(path_main, 'data')
path_clones = os.path.join(path_data, 'clonal_info')


##


# Reformat meta: new metadata
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta_orig.csv'), index_col=0)

# Other params
correction_type = 'reference-free'
coverage_treshold = 15
umi_treshold = 5
p_treshold = 1 
max_ratio_treshold = .7,
normalized_abundance_treshold = .7

# Get all cells, and visualization for all samples
CELL_DFs = []
dist_reads = plt.figure(figsize=(16,8))
dist_combos = plt.figure(figsize=(16,8))

# Read saved sample specific treshold
for i,sample in enumerate(meta['sample'].unique()):

    # Get counts
    with open(os.path.join(path_data, 'clonal_info', f'{sample}.pickle'), 'rb') as p:
        COUNTS = pickle.load(p)

    # Filter only QCed cells
    counts = COUNTS['reference-free']
    counts = counts.loc[
        counts['CBC']
        .isin(meta.query('sample==@sample')
        .index.map(lambda x: x.split('_')[0]))
    ].copy()

    # Filter UMIs
    counts = mark_UMIs(counts, coverage_treshold=coverage_treshold, nbins=50)

    # Add viz dist_reads
    ax = dist_reads.add_subplot(2,4,i+1)
    viz_UMIs(counts, by='status', ax=ax, nbins=50)
    ax.set(title=sample)

    # Get combos
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')

    # Filtering CBC-GBC
    M, _ = filter_and_pivot(
        df_combos, 
        umi_treshold=umi_treshold, 
        p_treshold=p_treshold,  
        max_ratio_treshold=max_ratio_treshold,
        normalized_abundance_treshold=normalized_abundance_treshold
    )

    # Add viz dist_combos
    ax = dist_combos.add_subplot(2,4,i+1)
    sns.kdeplot(data=df_combos, x='normalized_abundance', y='max_ratio', ax=ax)
    ax.axvline(x=normalized_abundance_treshold, color='k')
    ax.axhline(y=max_ratio_treshold, color='k')
    ax.set(title=sample)

    # GBC sets checks
    sets = get_clones(M)
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    ax.text(.6,.1, f'n cells: {len(unique_cells)}', transform=ax.transAxes)
    filtered_M = M.loc[unique_cells]
    clones_df = get_clones(filtered_M)
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC_set')

    )
    cells_df.index = cells_df.index.map(lambda x: x + f'_{sample}')

    # Append
    CELL_DFs.append(cells_df)


##
    

# Save figs and update meta
dist_reads.tight_layout()
dist_reads.savefig(os.path.join(path_main, 'results', 'clonal', 'UMI_tresholds.png'), dpi=300)
dist_combos.tight_layout()
dist_combos.savefig(os.path.join(path_main, 'results', 'clonal', 'filtered_combos.png'), dpi=300)

# New meta
meta_new = pd.concat(CELL_DFs).join(meta).rename(columns={'GBC_set':'GBC'})

# Save final metadata
meta_new['seq_run'] = 'run1'
meta_new.to_csv(os.path.join(path_data, 'cells_meta.csv'))


##