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
FIG_OUT = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/updates/"
DIR = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/"
sc.settings.figdir = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/updates/"
plt.rcParams["font.family"] = "Arial"


same_col_code = {"Stroma1": sns.color_palette("tab10")[0],
                "Stroma2": sns.color_palette("tab10")[1],
                "Epithelial": sns.color_palette("tab10")[2],
                "Stroma3": sns.color_palette("tab10")[3],
                "Immune": sns.color_palette("tab10")[4],
                "Stroma4": sns.color_palette("tab10")[5],
                "Stroma5": sns.color_palette("tab10")[6]}


def Find_co_exp(norm, sub_adata, LN, Fig_name):
    """ Check existence of target co-expression

    Args:
        norm (_type_): norm method
        sub_adata (_type_): each sample anndata
        LN (_type_): the number of lymph node cluster
        Fig_name (_type): Output fig name
    Returns:
        _type_: proportion of target co-expression
    """

    sam = sub_adata.obs["sample"].unique()[0]
    total_spots = sub_adata.obs.shape[0]
    ## Extract LN
    if LN == "Y":
        pass
    
    else:
        sub_adata = sub_adata[sub_adata.obs["annotation"]!=LN]

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
            
        plt.savefig(FIG_OUT+sam+"_{}_{}_portion.pdf".format(Fig_name, out), bbox_inches="tight")
        
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


def LR(sub_adata, gene1, gene2, gene_list, norm, LN, dot_col, Fig_name):
    """ Find L-R pairs in spatial region
    Args:
        sub_adata (_type_): processed anndata (after clustering)
        gene1 (_type_): Gene1
        gene2 (_type_): Gene2
        gene_list (_type_): target gene_list for complex
        norm (_type_): Method of norm
        LN (_type_): Lymph node cluster
        Fig_name (_type_): Output fig name
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
        plt.savefig(FIG_OUT+sam+"_{}_LR_{}_{}_portion.pdf".format(Fig_name,gene1,gene2), 
                    bbox_inches="tight")
        
        return sub_adata, dot_col
    
    else:
        dot_col.append(LR_name)

        return sub_adata, dot_col

# %%
## Data load
clu_adata3 = sc.read_h5ad(f"{DIR}/clu_adata3_input_for_fig7_c-i.h5ad")
clu_adata4 = sc.read_h5ad(f"{DIR}/clu_adata4_input_for_fig7_c-i.h5ad")

## Identification of Existence tissue resident immune cells
## dot size can be changed by line 315
mark_adata3 = Find_co_exp("raw", clu_adata3, "Y", "Supp_Fig8C")
mark_adata4 = Find_co_exp("raw", clu_adata4, "Y", "Supp_Fig8C")

mark_adata3 = Find_co_exp("raw", clu_adata3, "Immune", "Fig7C")
mark_adata4 = Find_co_exp("raw", clu_adata4, "Immune", "Fig7C")
# %%
## LR enrichment 
## Supp Fig8D
def Run_LR(adata, LN, norm):
    # Cd49d=Itga4, Cd29=Itgb1, Cd61=Itgb3, Cd11d=Itgad, and Cd18=Itgb2
    dot_col = []
    lr_adata, dot_col = LR(adata,"Plaur","a4b1",["Plaur","Itga4","Itgb1"], norm, LN, dot_col, "Fig7E")
    lr_adata, dot_col = LR(lr_adata,"Vcam1","a4b1",["Vcam1","Itga4","Itgb1"], norm, LN, dot_col, "Supp_Fig8D")
    lr_adata, dot_col = LR(lr_adata,"Vcam1","a4b7",["Vcam1","Itga4","Itgb7"], norm, LN, dot_col, "Supp_Fig8D")
    
    dot_col2 = []
    lr_adata, dot_col2 = LR(lr_adata,"Pdcd1","Fam3c",["Pdcd1","Fam3c"], norm, LN, dot_col2, "Fig7H")
    lr_adata, dot_col2 = LR(lr_adata,"Ccl5","Ackr4",["Ccl5","Ackr4"], norm, LN, dot_col2, "Supp_Fig8D")
    lr_adata, dot_col2 = LR(lr_adata,"Ccl3","Ide",["Ccl3","Ide"], norm, LN, dot_col2, "Supp_Fig8D")

    dot_col3 = []
    lr_adata, dot_col3 = LR(lr_adata,"Jag1","Notch3",["Jag1","Notch3"], norm, LN, dot_col3, "Fig7K")
    lr_adata, dot_col3 = LR(lr_adata,"Notch1","Jag1",["Notch1","Jag1"], norm, LN, dot_col3, "None")
    lr_adata, dot_col3 = LR(lr_adata,"Jag1","Notch2",["Jag1","Notch2"], norm, LN, dot_col3, "None")
    
    
    LR_dict = {"Gzmk-Epithelial-left":dot_col,
               "Gzmk-Epithelial-right":dot_col2,
               "Gamma Delta-Epithelial":dot_col3}
    
    # left, right represent of side of Epithelial cell Ep:Gzmk, Gzmk:Ep
    
    return lr_adata, LR_dict


int_adata3, int_dict3 = Run_LR(mark_adata3, "Immune", "raw")
int_adata4, int_dict4 = Run_LR(mark_adata4, "Immune", "raw")

# %%
def Run_LR2(adata, dict, kws, figname):
    """ Draw dot plot using LR fraction

    Args:
        adata (_type_): From Run_LR
        dict (_type_): From Run_LR
        kws (_type_): "Gzmk-Epithelial" or "Gamma Delta-Epithelial"
        figname (_type_): output fig name
    """
    
    epi = adata[adata.obs["annotation"]=="Epithelial"]
    prop = [epi.obs[i].sum()/epi.obs.shape[0] for i in dict[kws]]
    max_prop = max(prop)
    dp = sc.pl.dotplot(adata,
                       var_names={kws:dict[kws]},
                       groupby=['annotation'],
                       standard_scale='var',    # between clus in each pairs
                       dot_max=float(max_prop*1.2),
                       save=figname+"_"+adata.obs["sample"].unique()[0]+"_"+kws+"_LR_dotplot_gamma")


# %%
Run_LR2(int_adata3, int_dict3, "Gzmk-Epithelial-left", "Fig7E")
Run_LR2(int_adata3, int_dict3, "Gzmk-Epithelial-right", "Fig7H")
Run_LR2(int_adata3, int_dict3, "Gamma Delta-Epithelial", "Fig7K")
Run_LR2(int_adata4, int_dict4, "Gzmk-Epithelial-left", "Fig7E")
Run_LR2(int_adata4, int_dict4, "Gzmk-Epithelial-right", "Fig7H")
Run_LR2(int_adata4, int_dict4, "Gamma Delta-Epithelial", "Fig7K")

# %%
## Find LR with tissue residential
## Supp Fig8F
def Vis_co_loc(int_adata, cell, lr, only_epi, Fig_name):
    """ Visualize co-occurence

    Args:
        int_adata (_type_): Anndata
        cell (_type_): Target cell markers
        lr (_type_): Ligand-receptor
        only_epi (_type_): focusing only Epithelial or not
        Fig_name (_type_): Output figname
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
    plt.savefig(FIG_OUT+sam+"_{}_LR_{}_on_target.pdf".format(Fig_name,lr), 
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
    
        plt.savefig(FIG_OUT+sam+"_{}_LR_{}_on_target_zoom{}.pdf".format(Fig_name,lr,k), 
                    bbox_inches="tight")
    

# %%
## Supp Fig8F
## Check Gzmk with various LR
gz_lr = int_dict3["Gzmk-Epithelial-left"] + int_dict3["Gzmk-Epithelial-right"]
Fig_name = {"Plaur-a4b1":"Fig7F", "Vcam1-a4b1":"Supp_Fig8F", "Vcam1-a4b7":"Supp_Fig8F", "Pdcd1-Fam3c":"Fig7I","Ccl5-Ackr4":"Supp_Fig8F"}
for query_lr in gz_lr:
    Vis_co_loc(int_adata3, "Cd3/8_Gzmk", query_lr, "None", Fig_name[query_lr])
    Vis_co_loc(int_adata4, "Cd3/8_Gzmk", query_lr, "None", Fig_name[query_lr])
Vis_co_loc(int_adata3, "Gamma-deltaT", int_dict3["Gamma Delta-Epithelial"][0], "None", "Fig7L")
Vis_co_loc(int_adata4, "Gamma-deltaT", int_dict3["Gamma Delta-Epithelial"][0], "None", "Fig7L")
# %%
