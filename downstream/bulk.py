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
path_main = '/Users/IEO5505/Desktop/MyoD-repro/'
path_data = os.path.join(path_main, 'data')


df = pd.read_excel(os.path.join(path_main, 'data', 'NT_counts.xlsx'))


np.sum(df['Sample_S54609_m1_LN']>0)


fig, ax = plt.subplots()
sns.kdeplot(df['Sample_S54609_m1_LN'].sort_values(ascending=False)[100:], ax=ax, fill=True)
df['Sample_S54609_m1_LN'].sort_values(ascending=False)[100:].describe()
fig.tight_layout()
plt.show()


path_results = os.path.join(path_main, 'results', 'pp_summary')


# Read prevalences
df = pd.read_csv(os.path.join(path_results,'clonal_prevalences.csv'), index_col=0)

df_ = (
    df.groupby('sample')
    .apply(lambda x: x['GBC'].unique().size)
    .sort_values(ascending=False)
    .to_frame('n')
    .reset_index()
)

fig, ax = plt.subplots()
bar(df_, 'n', 'sample', s=.75, c='k', ax=ax)
format_ax(ax=ax, title='n clones by sample', ylabel='n', xticks=df_['sample'])
ax.spines[['left', 'top', 'right']].set_visible(False)
fig.savefig(
    os.path.join(
        path_results, 'n_clones.png'
    ), dpi=300
)


##




