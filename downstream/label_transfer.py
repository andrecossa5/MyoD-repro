"""
Label transfer with data from Pozniack et al., 2024
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import anndata
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics import pairwise_distances
from Cellula.preprocessing._pp import seurat_s
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


## 


# Paths
path_main = '/Users/IEO5505/Desktop/MyoD-repro/'
dataset = 'MM13'
path_data = os.path.join(path_main, 'data', dataset)
path_ref = os.path.join(path_main, 'data', 'misc', 'pozniack.h5ad')


##


# Query
query = sc.read(os.path.join(path_data, 'small.h5ad'))                    
query = anndata.AnnData(X=query.layers['raw'], obs=query.obs[['sample', 'leiden']], var=query.var[[]], layers={'counts':query.layers['raw']}) 
query.obs['study'] = 'Duca et al.,'
query.X.A[:10,:10]
query.layers['counts'].A[:10,:10]
query.obs
query.var

# Ref
# expr = pd.read_csv(os.path.join(path_main, 'data', 'misc', 'hvgs.csv'), index_col=0)
# meta = pd.read_csv(os.path.join(path_main, 'data', 'misc', 'meta.csv'), index_col=0)
# 
# ref = anndata.AnnData( 
#     X=csr_matrix(expr.values),
#     obs=meta[['Malignant_clusters', 'sample_ID']].rename(columns={'sample_ID':'sample'}), 
#     var=pd.DataFrame(index=expr.columns), 
#     layers={'counts':csr_matrix(expr.values)}
# )
# ref.layers['counts'].A[:10,:10]
# ref.obs
# ref.var
# ref = ref[ref.obs['Malignant_clusters'].sample(frac=.4).index].copy()  
# 
# ref.write(os.path.join(path_main, 'data', 'misc', 'pozniack.h5ad'))
ref = sc.read( os.path.join(path_main, 'data', 'misc', 'pozniack.h5ad'))  
ref.obs['study'] = 'Pozniack et al.,'
ref.X.A[:10,:10]
ref.layers['counts'].A[:10,:10]
ref.obs
ref.var


##


# Build joint dataset
adata = anndata.concat([query.copy(), ref.copy()])       
adata.X.A[1:10,1:10]                                                                    
adata.obs['nUMIs'] = (adata.X.A).sum(axis=1)                                            # nUMIs, MT-genes % (nope)
# adata.obs['mito_perc'] = (adata[:,adata.var_names.str.startswith('MT')].X.A).sum(axis=1) / (adata.X.A).sum(axis=1)
adata.obs = adata.obs.join(ref.obs['Malignant_clusters']).join(query.obs['leiden'])     # Put back labels
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata 

# Cell cycle
sc.tl.score_genes(adata, gene_list=seurat_s, score_name='seurat_s')

# HVGs
n_hvgs = 2000
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_hvgs, layer="counts", batch_key="study", subset=True)  # HVGs                    

# Build a latent representation of the joint dataset
# max_epochs = 200
# scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="study", continuous_covariate_keys=['nUMIs', 'seurat_s'])
# scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
# Train SCVI
# scvi_model.train(max_epochs=max_epochs)

# Save
# with open(os.path.join(path_main, 'data', 'misc', 'reference_pozniack.pickle'), 'wb') as f:
#     pickle.dump(scvi_model, f)


##


# Label transfer

# Load pre-trained model
with open(os.path.join(path_main, 'data', 'misc', 'reference_pozniack.pickle'), 'rb') as f:
    scvi_model = pickle.load(f)

# Set keys
SCVI_LATENT_KEY = "X_scVI"
SCANVI_CELLTYPE_KEY = "celltype_scanvi" 
# Get latent space and set celltype keys
adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
adata.obs[SCANVI_CELLTYPE_KEY] = "Unknown"                                                      # Add columns for scanvi prediction 
adata.obs[SCANVI_CELLTYPE_KEY].loc[ref.obs_names] = ref.obs['Malignant_clusters'].values        # Set reference cell types
adata.obs[SCANVI_CELLTYPE_KEY].value_counts()

# Initialize scanvi from pre-trained scvi
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    unlabeled_category="Unknown",
    labels_key=SCANVI_CELLTYPE_KEY,
)
# Train SCANVI
scanvi_model.train(max_epochs=20, n_samples_per_label=100)

# Set keys for prediction
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTION_KEY = "predicted_cell_states"

# Predict celltype of unlabeled query cells
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
adata.obs[SCANVI_PREDICTION_KEY] = scanvi_model.predict(adata)
adata.obs[SCANVI_CELLTYPE_KEY].value_counts()
adata.obs[SCANVI_PREDICTION_KEY].value_counts()

# Add predictions to original data and save
query = sc.read(os.path.join(path_data, 'small.h5ad'))  
query.obs['predicted_cell_states'] = adata[query.obs_names].obs['predicted_cell_states'].values
query.write(os.path.join(path_data, 'small.h5ad'))


##


# UMAP
sc.pp.neighbors(adata, use_rep='X_scVI', n_pcs=30, n_neighbors=15)
sc.tl.umap(adata)

# Viz
df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2']).join(adata.obs)
df.columns

fig, ax = plt.subplots(figsize=(5.3,5))
study_colors = {'Pozniack et al.,' : 'grey', 'Duca et al.,' : 'green'}
celltype_colors = create_palette(df, 'leiden', ten_godisnot)
draw_embeddings(df.query('study == "Duca et al.,"'), 
                cat='leiden', 
                # cont='nUMIs',
                ax=ax, 
                title='Unsupervised clusters: MM13 model', 
                legend_kwargs={'colors':celltype_colors, 'loc':'upper left'},
                #query='study == "Pozniack et al.,"'
                )
ax.axis(False)
fig.tight_layout()
plt.show()


##


# Crosstab
df_duca = df.query('study == "Duca et al.,"')
X = pd.crosstab(df_duca['predicted_cell_states'], df_duca['leiden']).values
order = leaves_list(linkage(pairwise_distances(X)))
order_ = leaves_list(linkage(pairwise_distances(X.T)))

fig, ax = plt.subplots(figsize=(7.2,5))
plot_heatmap(X.loc[X.index[order], X.columns[order_]], ax=ax, x_names_size=10, y_names_size=10, 
            xlabel='Leiden clusters', ylabel='Inferred cell states', 
            title='Relationship leiden clusters vs \n inferred cell states', label='n cells')
fig.tight_layout()
plt.show()


##


# nUMIs, cc
fig, axs = plt.subplots(1,2,figsize=(8,5))
box(df, x='study', y='nUMIs', ax=axs[0], c='grey')
format_ax(ax=axs[0], title='nUMIs', ylabel='nUMIs')
box(df, x='study', y='seurat_s', ax=axs[1], c='grey')
format_ax(ax=axs[1], title='Cell cycle signature (S-phase)', ylabel='nUMIs')
fig.tight_layout()
plt.show()


##