# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import numpy as np
# import scvi
from scipy.stats import skewtest
from matplotlib.pyplot import rc_context

""" 
    This code is same with CLEAN_01_Main_indi_data.py
    Same process with 01_Main_merge_annot.py only different input file
    This script use individual sample, not merged file
    TO DO: Change font style to Arial
"""

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.settings.verbosity = 3
plt.rcParams["figure.figsize"] = (3,3)
FIG_OUT = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/updates"
DIR = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/"
sc.settings.figdir = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/updates/"
plt.rcParams["font.family"] = "Arial"

sam_pal = {"0": sns.color_palette("tab10")[0],
           "1": sns.color_palette("tab10")[1],
           "2": sns.color_palette("tab10")[2],
           "3": sns.color_palette("tab10")[3],
           "4": sns.color_palette("tab10")[4],
           "5": sns.color_palette("tab10")[5],
           "6": sns.color_palette("tab10")[6],
           "7": sns.color_palette("tab10")[7],
           "8": sns.color_palette("tab10")[8],
           "9": sns.color_palette("tab10")[9],
           "10": "#FF0AB5"}

same_col_code = {"Stroma1": sns.color_palette("tab10")[0],
                "Stroma2": sns.color_palette("tab10")[1],
                "Epithelial": sns.color_palette("tab10")[2],
                "Stroma3": sns.color_palette("tab10")[3],
                "Immune": sns.color_palette("tab10")[4],
                "Stroma4": sns.color_palette("tab10")[5],
                "Stroma5": sns.color_palette("tab10")[6]}

def Load(sam_num):
    adata = sc.read_h5ad(DIR+"merged_YBC.h5ad")    # From 00_new_merge_ST.py
    batch_info = {"SC2300278_3M-2858":"batch1 (Younger)",
                "SC2300279_3M-2862":"batch1 (Younger)",
                "SC2300280_18M-2007":"batch2 (Older)",
                "SC2300281_18M-2010":"batch2 (Older)"}
    adata.obs["batch"] = adata.obs["sample"].map(batch_info)
    sam = adata.obs["sample"].unique()[sam_num]
    adata = adata[adata.obs["sample"]==sam,:]

    return adata, sam


def Do_clustering(adata, hgv, npcs, res):
    def Find_HGVs():
        adata.uns["log1p"]["base"] = None
        sc.pp.highly_variable_genes(adata, 
                                    layer="sqrt_norm",
                                    n_top_genes=2500)
        var_genes_all = adata.var.highly_variable
        print("Highly variable genes: %d"%sum(var_genes_all))
        print(adata.var["highly_variable"].value_counts())

    ## Find HGVs
    if hgv == "y" or hgv =="Y":
        Find_HGVs()
    
    ## RUN PCA
    if "highly_variable" in adata.var.columns:
        sc.pp.pca(adata, use_highly_variable=True)
    else:
        sc.pp.pca(adata, use_highly_variable=False)
    
    sns.lineplot(adata.uns["pca"]["variance_ratio"][1:20])
    sns.scatterplot(adata.uns["pca"]["variance_ratio"][1:20])
    ## Construct neighbor graph and UMAP
    sc.pp.neighbors(adata, n_pcs=npcs)
    plt.rcParams["figure.figsize"] = (6,6)
    sc.tl.umap(adata)
    sc.settings.figdir = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/uncorr/"

    ## Add cluster before correction
    sc.tl.leiden(adata, resolution=res, key_added="Single sample clu")

    return adata


def Check_mark(adata, norm, sam_id, Fig_name):
    """ Check the scRNA markgenes on ST cluster

    Args:
        adata (_type_): output of Do_clustering
        norm (_type_): kinds of count in marker gene analysis
        sam_id (_type_): sample library id
        Fig_name (_type_): output Fig name
    """
    
    plt.rcParams["figure.figsize"] = (6,6)
    ## Check the uMAP results
    sc.pl.umap(adata, color=["annotation"], palette=same_col_code)

    ## Check dot plot on the marker genes from scRNA
    marker_dict = {
                  "Epithelial":["Epcam"],
                  "Lum-AV":["Csn3","Wfdc18","Elf5"],
                  "Lum-HS":["Esr1","Cited1","Prlr"],
                  "Myoepithelial":["Krt17","Acta2","Krt14"],
                  "Immune":["Ptprc"],
                  "T Cell":["Cd3d","Cd4","Cd8b1"],
                  "B Cell":["Blnk","Cd79a","Ms4a1"],
                  "DC/Macro":["Cd14","Itgax","Napsa"],
                  "Fibroblast":["Fn1","Col1a1","Col1a2","Pi16","Col15a1"],
                  "Adipocytes":["Adipoq","Plin1","Cidec"]
                  }
    sc.pl.dotplot(adata, marker_dict, groupby="annotation", 
                  layer=norm, standard_scale="var",
                  dendrogram=False,
                  categories_order=["Epithelial", "Immune", "Stroma1", "Stroma2", "Stroma3", "Stroma4"],
                  save="{}_scRNA_marker_{}".format(Fig_name,sam_id))

    ## Find the marker genes based on found cluster (before correction)
    try:
        adata.uns["log1p"]["base"] = None
    except:
        pass

    sc.tl.rank_genes_groups(adata,
                            groupby="annotation",
                            method="wilcoxon",
                            layer=norm
                            )
    
    sc.pl.rank_genes_groups_dotplot(adata,
                                    groupby="annotation",
                                    groups=["Epithelial", "Immune", "Stroma1", "Stroma2", "Stroma3", "Stroma4"],
                                    min_logfoldchange=1.,
                                    n_genes=7,
                                    standard_scale="var",
                                    layer=norm,
                                    dendrogram=False,
                                    categories_order=["Epithelial", "Immune", "Stroma1", "Stroma2", "Stroma3", "Stroma4"],
                                    save="{}_ST_marker_{}".format(Fig_name,sam_id)
                                    )
    
    if len(adata.obs["sample"].unique()) == 1:  # Only draw figure in single sample
        sc.pl.spatial(adata, 
                      library_id=adata.obs["sample"].unique()[0],
                      color="Single sample clu",
                      palette=sam_pal,
                      alpha=.8,
                      size=1.5,
                      frameon=False)
    
    return adata


# %%
## Data load
age_group = "Old" # It should be set as Old or Young
sample_dict = {"Old":[2,3], "Young":[0,1]}
sample1 = sample_dict[age_group][0]
sample2 = sample_dict[age_group][1]
adata3, lib_id3 = Load(sample1)
adata4, lib_id4 = Load(sample2)
np.random.seed(11)
# %%
## Applying BBKNN
plt.rcParams["figure.figsize"] = (5,5)
old = sc.concat([adata3, adata4])
sc.pp.pca(old)
random_seed = 123123
np.random.seed(random_seed)
sc.pp.neighbors(old)
sc.tl.umap(old)
sc.pl.umap(old, color=["sample"])
sc.external.pp.bbknn(old, batch_key="sample", n_pcs=10)
sc.tl.umap(old)

if age_group == "Young":
    sc.tl.leiden(old, resolution=.1, key_added="Single sample clu")  # Younger samples
else:
    sc.tl.leiden(old, resolution=.15, key_added="Single sample clu")  # Older samples 
# %%
## Draw uMAP
sc.pl.umap(old, color=["sample","Single sample clu"])
# %%
## Label transfer to draw uMAP (From merged one to individual sample)
adata3.obs["Single sample clu"] = old[old.obs["sample"]==lib_id3].obs["Single sample clu"]
adata3.obsm["X_pca"] = old[old.obs["sample"]==lib_id3].obsm["X_pca"]
adata3.obsm["X_umap"] = old[old.obs["sample"]==lib_id3].obsm["X_umap"]
# %%
sc.pl.spatial(adata3,
              library_id=lib_id3,
              color=["Single sample clu"])
# %%
adata4.obs["Single sample clu"] = old[old.obs["sample"]==lib_id4].obs["Single sample clu"]
adata4.obsm["X_pca"] = old[old.obs["sample"]==lib_id4].obsm["X_pca"]
adata4.obsm["X_umap"] = old[old.obs["sample"]==lib_id4].obsm["X_umap"]
# %%
sc.pl.spatial(adata4,
              library_id=lib_id4,
              color=["Single sample clu"],
              size=1.5)
# %%
## Set the same color to same cluster
def Annot_young(df):
    if df[0] == "5":
        return "Stroma4"
    elif df[0] == "1":
        return "Stroma2"
    elif df[0] == "2":
        return "Stroma3"
    elif df[0] == "3":
        return "Epithelial"
    elif df[0] == "4":
        return "Immune"
    elif df[0] == "0":
        return "Stroma1"
    elif df[0] == "6":
        return "Stroma5"

def Annot_old(df):
    if df[0] == "3":
        return "Stroma4"
    elif df[0] == "1":
        return "Stroma3"
    elif df[0] == "2":
        return "Stroma2"
    elif df[0] == "4":
        return "Epithelial"
    elif df[0] == "5":
        return "Immune"
    elif df[0] == "0":
        return "Stroma1"
    elif df[0] == "6":
        return "Stroma5"

age_color = {"Young":Annot_young, "Old":Annot_old}
adata3.obs["annotation"] = adata3.obs["Single sample clu"].apply(age_color[age_group])
adata4.obs["annotation"] = adata4.obs["Single sample clu"].apply(age_color[age_group])
old.obs["annotation"] = old.obs["Single sample clu"].apply(age_color[age_group])

## Draw uMAP
plt.rcParams["figure.figsize"] = (5,5)
sc.pl.umap(old, color=["annotation"], 
           save="Supp_Fig8A", 
           s=12,
           title="Cell type annotation by Visium",
           palette=same_col_code)
# %%
## Check the marker genes
clu_adata3 = Check_mark(adata3, "sqrt_norm", lib_id3, "Supp_Fig8B")
clu_adata4 = Check_mark(adata4, "sqrt_norm", lib_id4, "Supp_Fig8B")
# %%
## Check the marker genes on merged set
clu_old = Check_mark(old, "sqrt_norm", "Merged_old", "Fig7B")
# %%
sc.pl.spatial(clu_adata3,
              library_id=clu_adata3.obs["sample"].unique()[0],
              color=["annotation"],
              size=1.7,
              palette=same_col_code,
              save="Fig7a_ST_{}".format(lib_id3))
# %%
sc.pl.spatial(clu_adata4,
              library_id=clu_adata4.obs["sample"].unique()[0],
              color=["annotation"],
              size=1.7,
              palette=same_col_code,
              save="Fig7a_ST_{}".format(lib_id4))
# %%
clu_adata3.write(f"{DIR}/clu_adata3_input_for_fig7_c-i.h5ad")
clu_adata4.write(f"{DIR}/clu_adata4_input_for_fig7_c-i.h5ad")
# %%
