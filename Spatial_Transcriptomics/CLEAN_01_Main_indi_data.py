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
    adata = sc.read_h5ad(DIR+"merged.h5ad")    # From 00_new_merge_ST.py
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


def Check_mark(adata, norm, sam_id):
    """ Check the scRNA markgenes on ST cluster

    Args:
        adata (_type_): output of Do_clustering
        norm (_type_): kinds of count in marker gene analysis
        sam_id (_type_): sample library id
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
                  save="scRNA_marker_{}".format(sam_id))

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
                                    save="ST_marker_{}".format(sam_id)
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


def Find_co_exp(norm, sub_adata, LN):
    """ Check existence of target co-expression

    Args:
        norm (_type_): norm method
        sub_adata (_type_): each sample anndata
        LN (_type_): the number of lymph node cluster

    Returns:
        _type_: proportion of target co-expression
    """

    sam = sub_adata.obs["sample"].unique()[0]
    total_spots = sub_adata.obs.shape[0]
    ## Extract LN
    # sub_adata = sub_adata[sub_adata.obs["annotation"]!=LN]

    ## Make exp dataframe
    exp = pd.DataFrame()
    target_genes = ["Cd3d","Cd3g","Cd3e","Cd247","Gzmk","Pdcd1","Cd8a","Cd8b1","Trdc","Il17a"]
    for gene in target_genes:
        index = sub_adata.var.index.tolist().index(gene)
        if exp.shape[0] == 0:
            exp = sub_adata.layers[norm].toarray()[:,index]
        else:
            exp = np.vstack([exp, sub_adata.layers[norm].toarray()[:,index]])

    exp = pd.DataFrame(exp).T
    exp.index = sub_adata.obs.index
    exp.columns = target_genes

    ## Make OR combi between
    exp["Cd3"] = exp[["Cd3d","Cd3g","Cd3e","Cd247"]].sum(axis=1)    # CD3 markers more general T cell (CD4, CD8 ...)
    exp["Cd8"] = exp[["Cd8a","Cd8b1"]].sum(axis=1)    # Specific T cell
    exp["Cd3 or Cd8"] = exp[["Cd3","Cd8"]].sum(axis=1)  # Cd3 or Cd8
    exp["gamma-delta"] = exp[["Il17a","Trdc"]].sum(axis=1)    # merged CD8+ marker
    
    ### Definition of co-expression COMBI
    ## Find Cd3 with Gzmk OR Pdcd1
    sub_adata.obs["Cd3_Gzmk"] = pd.Series(((exp[["Cd3", "Gzmk"]]>0).astype(int).sum(axis=1)==2).astype(int))
    sub_adata.obs["Cd3_Pdcd1"] = pd.Series(((exp[["Cd3", "Pdcd1"]]>0).astype(int).sum(axis=1)==2).astype(int))

    ## Find Cd8 with Gzmk OR Pdcd1
    sub_adata.obs["Cd8_Gzmk"] = pd.Series(((exp[["Cd8", "Gzmk"]]>0).astype(int).sum(axis=1)==2).astype(int))
    sub_adata.obs["Cd8_Pdcd1"] = pd.Series(((exp[["Cd8", "Pdcd1"]]>0).astype(int).sum(axis=1)==2).astype(int))

    ## Find Cd3 OR Cd8 with Gzmk OR Pdcd1
    sub_adata.obs["Cd3/8_Gzmk"] = pd.Series(((exp[["Cd3 or Cd8", "Gzmk"]]>0).astype(int).sum(axis=1)==2).astype(int))
    sub_adata.obs["Cd3/8_Pdcd1"] = pd.Series(((exp[["Cd3 or Cd8", "Pdcd1"]]>0).astype(int).sum(axis=1)==2).astype(int))

    ## Find Cd3 with Il17a == Gamma delta T cell
    sub_adata.obs["Gamma-deltaT"] = pd.Series(((exp[["Cd3", "gamma-delta"]]>0).astype(int).sum(axis=1)==2).astype(int))

    ## Check # of co-expressed spot
    ########################################################################
    gz_mark = "Cd3/8_Gzmk"
    pdcd_mark = "Cd3/8_Pdcd1"
    gdt_mark = "Gamma-deltaT"
    ########################################################################

    ## Extract spots which have co-expression
    ## USING this data for overlap
    n_gz = sub_adata[sub_adata.obs[gz_mark]==True,:]
    n_pd = sub_adata[sub_adata.obs[pdcd_mark]==True,:]
    gdt = sub_adata[sub_adata.obs[gdt_mark]==True,:]

    def Draw_fig(data, title, out, cell):
        ## Make a subplot
        figure = plt.figure(figsize=(6,6))
        gs = figure.add_gridspec(1,1)
        ax1 = figure.add_subplot(gs[:,:])
        
        ## For same size
        sc.pl.spatial(sub_adata,
                      library_id=sam,
                      color=["annotation"],
                      size=0.5, 
                      alpha_img=1.0,
                      alpha=0.0,
                      layer=norm,
                      frameon=False,
                      show=False,
                      ax=ax1)
        ax1.legend("",frameon=False)
        
        sc.pl.spatial(data,
                      alpha_img=0,
                      library_id=sam,
                      color=["annotation"],
                      size=3.5,
                      title=title,
                      layer=norm,
                      frameon=False,
                      show=False,
                      palette=same_col_code,
                      ax=ax1)
        ax1.legend("",frameon=False)
            
        plt.savefig(FIG_OUT+sam+"_{}_portion.pdf".format(out), bbox_inches="tight")
        
        ## Zoom in interesting point
        # import matplotlib
        # cvals  = [1, 2] # Value interval 
        # colors = ["#519E3E","#519E3E"]
        # norm_fac=plt.Normalize(min(cvals),max(cvals))
        # tuples = list(zip(map(norm_fac,cvals), colors))
        # cus_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

        # target = data[data.obs["annotation"]=="Epithelial",:].obsm["spatial"]
        # for k,coord in enumerate(target):
        #     figure = plt.figure(figsize=(6,6))
        #     gs = figure.add_gridspec(1,1)
        #     ax1 = figure.add_subplot(gs[:,:])
        #     zoom_int = 1500
        #     zoom_int2 = 5000
            
        #     sc.pl.spatial(data,
        #                 crop_coord=[coord[0]-zoom_int,coord[0]+zoom_int,coord[1]-zoom_int,coord[1]+zoom_int],
        #                 library_id=data.obs["sample"].unique()[0],
        #                 color=[cell],
        #                 alpha=.9,
        #                 cmap=cus_cmap,
        #                 show=False,
        #                 frameon=False,
        #                 ax=ax1,
        #                 colorbar_loc=None)
        
        #     plt.savefig(FIG_OUT+sam+"_{}_portion_zoomed_{}.pdf".format(out,k), bbox_inches="tight")

        #     sc.pl.spatial(data,
        #                 crop_coord=[coord[0]-zoom_int2,coord[0]+zoom_int2,coord[1]-zoom_int2,coord[1]+zoom_int2],
        #                 library_id=data.obs["sample"].unique()[0],
        #                 color=[cell],
        #                 alpha=.45,
        #                 cmap=cus_cmap,
        #                 frameon=False,
        #                 colorbar_loc=None)
        

    try:
        Draw_fig(gdt, "Gamma-delta T (Cd3 withe Il17a)", "Gamma-delta", gdt_mark)
        Draw_fig(n_pd, "Cd3/8 with Pdcd1", "Cd3_pdcd1", pdcd_mark)
        Draw_fig(n_gz, "Cd3/8 with Gzmk", "Cd3_gzmk", gz_mark)
    
    except:
        print("There are not target spots")

    max_value = sub_adata.obs["Cd3/8_Pdcd1"].sum() / sub_adata[sub_adata.obs["annotation"]=="Epithelial"].obs.shape[0]
    dp = sc.pl.dotplot(sub_adata, 
                       var_names=["Cd3/8_Gzmk","Cd3/8_Pdcd1","Gamma-deltaT"],
                       groupby=['annotation'],
                       save=sam+"_tissue_residential",
                       dot_max=float(0.02))  # Adjust dot_max 
                    #    dot_max=float(max_value*1.2))  # Adjust dot_max 

    return sub_adata


def LR(sub_adata, gene1, gene2, gene_list, norm, LN, dot_col):
    """ Find L-R pairs in spatial region
    Args:
        sub_adata (_type_): processed anndata (after clustering)
        gene1 (_type_): Gene1
        gene2 (_type_): Gene2
        gene_list (_type_): target gene_list for complex
        norm (_type_): Method of norm
        LN (_type_): Lymph node cluster
    """

    ## Make a subplot
    figure = plt.figure(figsize=(5,5))
    gs = figure.add_gridspec(4,4)
    ax1 = figure.add_subplot(gs[:,:])

    sam = sub_adata.obs["sample"].unique()[0]
    total_spots = sub_adata.obs.shape[0]
    sub_adata = sub_adata[sub_adata.obs["annotation"]!=LN] # LN

    exp = np.array([])
    col = []
    for gene in gene_list:
        try:
            index = sub_adata.var.index.tolist().index(gene)
            col.append(gene)
            if exp.shape[0] == 0:
                exp = sub_adata.layers[norm].toarray()[:,index]
            else:
                exp = np.vstack([exp, sub_adata.layers[norm].toarray()[:,index]])
        
        except:
            print("{} is not list".format(gene))


    exp = pd.DataFrame(exp).T
    exp.index = sub_adata.obs.index
    exp.columns = col

    ## Make AND combi
    LR_name = gene1 +"-"+ gene2
    if not gene2 in gene_list:
        ## complex case
        col.remove(gene1)
        exp[gene2] = exp[col].sum(axis=1)
        sub_adata.obs[LR_name] = ((exp[[gene1, gene2]]>0).astype(int).sum(axis=1)==2).astype(int)
    else:
        ## individual gene case
        sub_adata.obs[LR_name] = ((exp[[gene1, gene2]]>0).astype(int).sum(axis=1)==2).astype(int)
    
    if (sub_adata.obs[LR_name].sum()) > 0:
        ## Check # OF co-expressed spot
        LR = sub_adata[sub_adata.obs[LR_name]==True,:]
        plt.rcParams["figure.figsize"] = (5,5)
        
        ## For same size
        sc.pl.spatial(sub_adata,
                    library_id=sam,
                    color=["annotation"],
                    size=0.5, img_key='hires', 
                    alpha=0.0,
                    vmax=1, vmin=0,
                    layer=norm,
                    frameon=False,
                    show=False,
                    ax=ax1)
        ax1.legend("",frameon=False)

        sc.pl.spatial(LR,
                    library_id=sam,
                    color=["annotation"],
                    size=4.5, img_key='hires', alpha=0.6,
                    vmax=1, vmin=0,
                    # save=sam+"_LR_{}_{}".format(gene1,gene2),
                    layer=norm,
                    title="{}-{}".format(gene1,gene2),
                    palette=same_col_code,
                    show=False,
                    ax=ax1,
                    frameon=False)
        ax1.legend("",frameon=False)

        dot_col.append(LR_name)
        plt.savefig(FIG_OUT+sam+"_LR_{}_{}_portion.pdf".format(gene1,gene2), 
                    bbox_inches="tight")
        
        return sub_adata, dot_col
    
    else:
        dot_col.append(LR_name)

        return sub_adata, dot_col

# %%
## Data load
file_code = {0:"Y1", 1:"Y2", 2:"O1", 3:"O2"}
file1_value = 2
file2_value = 3

adata3, lib_id3 = Load(file1_value)
adata4, lib_id4 = Load(file2_value)
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
sc.tl.leiden(old, resolution=.1, key_added="Single sample clu")  # Younger samples
# sc.tl.leiden(old, resolution=.15, key_added="Single sample clu")  # Older samples 
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

adata3.obs["annotation"] = adata3.obs["Single sample clu"].apply(Annot_young)
adata4.obs["annotation"] = adata4.obs["Single sample clu"].apply(Annot_young)
old.obs["annotation"] = old.obs["Single sample clu"].apply(Annot_young)
## Draw uMAP
plt.rcParams["figure.figsize"] = (5,5)
sc.pl.umap(old, color=["annotation"], save="_Supp11_sample_uMAP_cluster", s=12,
           title="Cell type annotation by Visium",
           palette=same_col_code)
# %%
adata3.obs[["annotation"]].to_csv(DIR+"Young1_add_clu.txt",
                                  sep="\t")
adata4.obs[["annotation"]].to_csv(DIR+"Young2_add_clu.txt",
                                  sep="\t")
# %%
## Check the marker genes
clu_adata3 = Check_mark(adata3, "sqrt_norm", lib_id3)
clu_adata4 = Check_mark(adata4, "sqrt_norm", lib_id4)
# %%
## Check the marker genes on merged set
clu_old = Check_mark(old, "sqrt_norm", "Merged_old")
# %%
sc.pl.spatial(clu_adata3,
              library_id=clu_adata3.obs["sample"].unique()[0],
              color=["annotation"],
              size=1.7,
              palette=same_col_code,
              title="3M #ST1",
              save="ST_{}".format(lib_id3))
# %%
sc.pl.spatial(clu_adata4,
              library_id=clu_adata4.obs["sample"].unique()[0],
              color=["annotation"],
              size=1.7,
              palette=same_col_code,
              title="3M #ST2",
              save="ST_{}".format(lib_id4))
# %%
## Identification of Existence tissue resident immune cells
## dot size can be changed by line 315
mark_adata3 = Find_co_exp("raw", clu_adata3, "Immune")
mark_adata4 = Find_co_exp("raw", clu_adata4, "Immune")
# %%
## Export for Final_df with for revision
mark_adata3.write(DIR+"{}.h5ad".format(file_code[file1_value]))
mark_adata4.write(DIR+"{}.h5ad".format(file_code[file2_value]))
# %%
## Show the proportion for
# gz_mark = "Cd3/8_Gzmk"
# pdcd_mark = "Cd3/8_Pdcd1"
# gdt_mark = "Gamma-deltaT"
print(mark_adata4.obs["annotation"].value_counts())
print(mark_adata3.obs["annotation"].value_counts())
def Cal_pro(ctype):
    pro1 = mark_adata3.obs[["annotation",ctype]]
    pro2 = mark_adata4.obs[["annotation",ctype]]
    print("{}_{}%".format(ctype,(pro1[ctype].sum()/pro1.shape[0])*100))
    print("{}_{}%".format(ctype,(pro2[ctype].sum()/pro2.shape[0])*100))

Cal_pro("Cd3/8_Gzmk")
Cal_pro("Cd3/8_Pdcd1")
Cal_pro("Gamma-deltaT")
# %%
## Empirical P-value dist for enrichment
def Cal_emp(df, ctype):
    px = pd.DataFrame(df.obs[df.obs[ctype]>0]["annotation"]).copy()  # Total number of Gzmk spots
    px["annotation"] = px["annotation"].apply(lambda x : 1 if x == "Epithelial" else 0)
    dist = []
    k = 0
    while k < 10000:
        trial = []
        for _ in range(df.obs[df.obs[ctype]>0].shape[0]):
            ran = np.random.randint(0,df.obs.shape[0])
            exp = pd.DataFrame(df.obs.iloc[ran,:]).loc["annotation"].values[0]
            trial.append(exp)
        dist.append(trial.count("Epithelial")/df.obs[df.obs[ctype]>0].shape[0])
        k += 1

    ## Calculate empirical P-values
    nor_dist = (np.array(dist) - np.mean(dist)) / np.std(dist)
    nor_px = ((px["annotation"].sum()/df.obs[df.obs[ctype]>0].shape[0]) \
            - np.mean(dist)) / np.std(dist)
    plt.figure(figsize=(3,2))
    plt.axvline(nor_px)

    import scipy.stats as stats
    rv = (stats.norm(0,1))
    x = np.linspace(-5,5,100)
    y1 = rv.pdf(x)
    plt.plot(x,y1)
    print(1-rv.cdf(nor_px))


# for c_type in ["Gamma-deltaT","Cd3/8_Gzmk", "Cd3/8_Pdcd1"]:
#     Cal_emp(mark_adata3, c_type)
#     Cal_emp(mark_adata4, c_type)

# %%
def Cal_(cell):
    if cell == "Cd3/8_Gzmk":
        b = (mark_adata4[(mark_adata4.obs[cell]>0) & (mark_adata4.obs["annotation"]!="Stroma1")].obs["annotation"].value_counts() / mark_adata4[mark_adata4.obs[cell]>0].obs["annotation"].shape[0])
    else:
        b = (mark_adata4[mark_adata4.obs[cell]>0].obs["annotation"].value_counts() / mark_adata4[mark_adata4.obs[cell]>0].obs["annotation"].shape[0])
    a = (mark_adata3[mark_adata3.obs[cell]>0].obs["annotation"].value_counts() / mark_adata3[mark_adata3.obs[cell]>0].obs["annotation"].shape[0])
    a = pd.DataFrame(a).sort_values(by="count")
    b = pd.DataFrame(b).sort_values(by="count")
    
    def Plot(input,sam_id):
        plt.figure(figsize=(1,3))
        for i,cell_type in zip(range(0,input.shape[0]+1),input.index.tolist()):
            if i == 0:
                sns.barplot(input.iloc[0:,:].sum(), 
                            color=same_col_code[cell_type],
                            width=.6)
            else:
                sns.barplot(input.iloc[i:,:].sum(), 
                            color=same_col_code[cell_type],
                            width=.6)
        sns.despine()
        plt.title("O_{}".format(sam_id))
        plt.xticks([])
        plt.xlabel(cell)
        plt.savefig(FIG_OUT+"O_{}_{}.pdf".format(sam_id, cell.replace("/","_")),
                    bbox_inches="tight")
    Plot(a,"#ST1")
    Plot(b,"#ST2")

Cal_("Cd3/8_Gzmk")
Cal_("Cd3/8_Pdcd1")
Cal_("Gamma-deltaT")
# %%
## LR enrichment
def Run_LR(adata, LN, norm):
    # Cd49d=Itga4, Cd29=Itgb1, Cd61=Itgb3, Cd11d=Itgad, and Cd18=Itgb2
    dot_col = []
    lr_adata, dot_col = LR(adata,"Plaur","a4b1",["Plaur","Itga4","Itgb1"], norm, LN, dot_col)
    lr_adata, dot_col = LR(lr_adata,"Vcam1","a4b1",["Vcam1","Itga4","Itgb1"], norm, LN, dot_col)
    lr_adata, dot_col = LR(lr_adata,"Vcam1","a4b7",["Vcam1","Itga4","Itgb7"], norm, LN, dot_col)
    
    dot_col2 = []
    # lr_adata, dot_col2 = LR(lr_adata,"Pdcd1","Fam3c",["Pdcd1","Fam3c"], norm, LN, dot_col2)
    lr_adata, dot_col2 = LR(lr_adata,"Ccl5","Ackr4",["Ccl5","Ackr4"], norm, LN, dot_col2)
    # lr_adata, dot_col2 = LR(lr_adata,"Ccl3","Ide",["Ccl3","Ide"], norm, LN, dot_col2)

    dot_col3 = []
    # lr_adata, dot_col3 = LR(lr_adata,"Jag1","Notch3",["Jag1","Notch3"], norm, LN, dot_col3)
    # lr_adata, dot_col3 = LR(lr_adata,"Notch1","Jag1",["Notch1","Jag1"], norm, LN, dot_col3)
    # lr_adata, dot_col3 = LR(lr_adata,"Jag1","Notch2",["Jag1","Notch2"], norm, LN, dot_col3)
    
    LR_dict = {"Gzmk-Epithelial-left":dot_col,
               "Gzmk-Epithelial-right":dot_col2,
               "Gamma Delta-Epithelial":dot_col3}
    # left, right represent of side of Epithelial cell Ep:Gzmk, Gzmk:Ep
    
    return lr_adata, LR_dict


int_adata3, int_dict3 = Run_LR(mark_adata3, "Immune", "raw")
int_adata4, int_dict4 = Run_LR(mark_adata4, "Immune", "raw")
# %%
def Run_LR2(adata, dict, kws):
    """ Draw dot plot using LR fraction

    Args:
        adata (_type_): From Run_LR
        dict (_type_): From Run_LR
        kws (_type_): "Gzmk-Epithelial" or "Gamma Delta-Epithelial"
    """
    
    epi = adata[adata.obs["annotation"]=="Epithelial"]
    prop = [epi.obs[i].sum()/epi.obs.shape[0] for i in dict[kws]]
    max_prop = max(prop)
    dp = sc.pl.dotplot(adata,
                       var_names={kws:dict[kws]},
                       groupby=['annotation'],
                       standard_scale='var',    # between clus in each pairs
                       dot_max=float(max_prop*1.2),
                       save=adata.obs["sample"].unique()[0]+"_"+kws+"_LR_dotplot_gamma")


# %%
Run_LR2(int_adata3, int_dict3, "Gzmk-Epithelial-left")
Run_LR2(int_adata3, int_dict3, "Gzmk-Epithelial-right")
# Run_LR2(int_adata3, int_dict3, "Gamma Delta-Epithelial")
Run_LR2(int_adata4, int_dict4, "Gzmk-Epithelial-left")
Run_LR2(int_adata4, int_dict4, "Gzmk-Epithelial-right")
# Run_LR2(int_adata4, int_dict4, "Gamma Delta-Epithelial")
# %%
## Find LR with tissue residential
def Vis_co_loc(int_adata, cell, lr, only_epi):
    """ Visualize co-occurence

    Args:
        int_adata (_type_): Anndata
        cell (_type_): Target cell markers
        lr (_type_): Ligand-receptor
        only_epi (_type_): focusing only Epithelial or not
    """

    df = int_adata.obs[[cell,lr,"annotation"]]
    df["{}-{}".format(cell,lr)] = pd.Series(((df[[cell, lr]]>0).astype(int).sum(axis=1)==2).astype(int))
    temp = int_adata.copy()
    temp.obs["{}-{}".format(cell,lr)] = df["{}-{}".format(cell,lr)]
    
    ## Make custom cmap
    import matplotlib
    cvals  = [1, 2] # Value interval 
    colors = ["#FFF200","#FFF200"]
    norm_fac=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm_fac,cvals), colors))
    cus_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

    colors = ["#00BFFF","#00BFFF"]
    norm_fac=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm_fac,cvals), colors))
    cus_cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

    colors = ["#35FF46","#35FF46"]
    norm_fac=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm_fac,cvals), colors))
    cus_cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)
    
    colors = ["#A3A3A3","#A3A3A3"]
    norm_fac=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm_fac,cvals), colors))
    bg = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)
    
    figure = plt.figure(figsize=(6,6))
    gs = figure.add_gridspec(1,1)
    ax1 = figure.add_subplot(gs[:,:])
    sam = temp.obs["sample"].unique()[0]

    ## For same size
    sc.pl.spatial(temp,
                  library_id=sam,
                  color=[cell],
                  size=4, img_key='hires',
                  alpha=.0,
                  frameon=False,
                  show=False,
                  ax=ax1,
                  cmap=bg,
                  colorbar_loc=None)
    ax1.legend("",frameon=False)

    if only_epi != "None":
        temp = temp[temp.obs["annotation"]=="Epithelial"] 
    else:
        pass

    sc.pl.spatial(temp[temp.obs[lr]>0], 
                  alpha_img=0.0,
                  color=[lr], 
                  library_id=sam,
                  size=4,
                  alpha=.85,
                  cmap=cus_cmap2,
                  show=False,
                  frameon=False,
                  ax=ax1,
                  colorbar_loc=None
                  )
    
    sc.pl.spatial(temp[temp.obs[cell]>0],
                  alpha_img=0.0,
                  library_id=sam,
                  color=[cell],
                  size=4, img_key='hires', 
                  alpha=.99,
                  cmap=cus_cmap,
                  frameon=False,
                  show=False,
                  ax=ax1,
                  colorbar_loc=None
                  )
    
    sc.pl.spatial(temp[temp.obs["{}-{}".format(cell,lr)]>0], 
                  alpha_img=0.0,
                  color=["{}-{}".format(cell,lr)], 
                  library_id=sam,
                  size=4,
                  alpha=.99,
                  cmap=cus_cmap3,
                  show=False,
                  frameon=False,
                  ax=ax1,
                  colorbar_loc=None
                  )
                  
    ax1.legend("",frameon=False)
    ax1.set_title("{} on {}".format(lr, cell), fontsize=15)
    
    from matplotlib.lines import Line2D
    custom_legend = [Line2D([0],[0], marker="o", linestyle='None', label=cell,
                            markersize=8,
                            markerfacecolor="#FFF200",
                            markeredgecolor='None'),
                     Line2D([0],[0], marker="o", linestyle='None', label=lr,
                            markersize=8,
                            markerfacecolor="#04B9F6",
                            markeredgecolor='None')]
    ax1.legend(handles=custom_legend, frameon=False)
    sns.move_legend(ax1, "upper left", bbox_to_anchor=(0,1))
    plt.savefig(FIG_OUT+sam+"_LR_{}_on_target.pdf".format(lr), 
                    bbox_inches="tight")
    
    ## Zoom in interesting point
    target = temp[temp.obs[cell]>0].obsm["spatial"]
    for k,coord in enumerate(target):
        figure = plt.figure(figsize=(4,4))
        gs = figure.add_gridspec(1,1)
        ax1 = figure.add_subplot(gs[:,:])
        zoom_int = 1500 # 4500
       
        sc.pl.spatial(temp, 
                      crop_coord=[coord[0]-zoom_int,coord[0]+zoom_int,coord[1]-zoom_int,coord[1]+zoom_int],
                      library_id=temp.obs["sample"].unique()[0],
                      alpha=.0,
                      show=False,
                      frameon=False,
                      ax=ax1)
        
        sc.pl.spatial(temp[temp.obs[lr]>0],
                      crop_coord=[coord[0]-zoom_int,coord[0]+zoom_int,coord[1]-zoom_int,coord[1]+zoom_int],
                      img_key="hires",
                      library_id=temp.obs["sample"].unique()[0],
                      color=[lr],
                      alpha=.85,
                      cmap=cus_cmap2,
                      show=False,
                      frameon=False,
                      ax=ax1,
                      colorbar_loc=None)

        sc.pl.spatial(temp[temp.obs[cell]>0],
                      crop_coord=[coord[0]-zoom_int,coord[0]+zoom_int,coord[1]-zoom_int,coord[1]+zoom_int],
                      library_id=temp.obs["sample"].unique()[0],
                      color=[cell],
                      alpha=.99,
                      cmap=cus_cmap,
                      show=False,
                      frameon=False,
                      ax=ax1,
                      colorbar_loc=None)
        
        sc.pl.spatial(temp[temp.obs["{}-{}".format(cell,lr)]>0],
                      crop_coord=[coord[0]-zoom_int,coord[0]+zoom_int,coord[1]-zoom_int,coord[1]+zoom_int],
                      img_key="hires",
                      library_id=temp.obs["sample"].unique()[0],
                      color=["{}-{}".format(cell,lr)],
                      alpha=.99,
                      cmap=cus_cmap3,
                      show=False,
                      frameon=False,
                      ax=ax1,
                      colorbar_loc=None)
    
        plt.savefig(FIG_OUT+sam+"_LR_{}_on_target_zoom{}.pdf".format(lr,k), 
                    bbox_inches="tight")
    

# %%
## Check Gzmk with various LR
gz_lr = int_dict3["Gzmk-Epithelial-left"] + int_dict3["Gzmk-Epithelial-right"]
for query_lr in gz_lr:
    Vis_co_loc(int_adata3, "Cd3/8_Gzmk", query_lr, "None")
    Vis_co_loc(int_adata4, "Cd3/8_Gzmk", query_lr, "None")

# %%
## Check Gamma-delta T with Jag1-Notch3 LR
Vis_co_loc(int_adata3, "Gamma-deltaT", "Jag1-Notch3", "None")
Vis_co_loc(int_adata4, "Gamma-deltaT", "Jag1-Notch3", "None")
# %%