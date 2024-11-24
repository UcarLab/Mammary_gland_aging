# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import skewtest
from matplotlib.pyplot import rc_context


FIG_OUT = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/figures/"
sc.set_figure_params(facecolor="white", figsize=(8,8))
sc.settings.verbosity = 3
sample_list = ["SC2300278_3M-2858","SC2300279_3M-2862","SC2300280_18M-2007","SC2300281_18M-2010"]


def PP(sample_id):
    DIR = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/"
    DIR = DIR+sample_id+"/"
    adata = sc.read_visium(path=DIR,
                           source_image_path=DIR+"spatial/")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    mt_DIR = "/Users/kangh/Desktop/vs_code/data/mouse_scRNA/"
    mt_data = pd.read_excel(mt_DIR+"MitoCarta2.0.xlsx")
    mt_gene = mt_data.sort_values(["MCARTA2.0_score"], 
                                ascending=False).iloc[:250,:]
    adata.var["mt"] = adata.var.index.isin(mt_gene["Symbol"])
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata.obs["sample"] = sample_id
    print(adata)

    ## Filtering, cut off from above plot
    sc.pp.filter_genes(adata, min_cells=3)    # for at least minimum variantion
    sc.pp.filter_cells(adata, min_genes=245)
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    ## Solution for TypeError: can't multiply sequence by non-int of type 'float'
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
    
    adata.layers["raw"] = adata.X.copy()
    adata.layers["sqrt_norm"] = sc.pp.normalize_total(adata, inplace=False, target_sum=1e6)["X"]
    adata.layers["sqrt_norm"] = sc.pp.log1p(adata.layers["sqrt_norm"])
    sc.experimental.pp.normalize_pearson_residuals(adata)

    return adata


merged_df = []
for sample in sample_list:
    merged_df.append(PP(sample))
# %%
## Make merged data
DIR = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/"
con_data = sc.concat(merged_df, uns_merge="unique", fill_value=0)
with plt.rc_context({'axes.facecolor': 'white',
                        'figure.figsize': [4, 5]}):
    sc.pl.violin(con_data, "total_counts", groupby="sample", rotation=90, size=1.5)
results_file = DIR+"merged_YBC.h5ad"
con_data.write(results_file)
# %%
