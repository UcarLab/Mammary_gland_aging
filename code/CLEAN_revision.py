# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import numpy as np
# import scvi
from scipy.stats import skewtest, mannwhitneyu
from matplotlib.pyplot import rc_context
""" 
    This code is same with CLEAN_01_Main_indi_data.py
    Same process with 01_Main_merge_annot.py only different input file
    This script use individual sample, not merged file
    TO DO: Change font style to Arial
"""
# %%
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.settings.verbosity = 3
plt.rcParams["figure.figsize"] = (3,3)
FIG_OUT = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/"
DIR = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/"
sc.settings.figdir = "/Users/kangh/Desktop/vs_code/data/our_Visium/Dataset_52_Visium/For_revision/figures/uncorr/"
plt.rcParams["font.family"] = "Arial"


same_col_code = {"Stroma1": sns.color_palette("tab10")[0],
                "Stroma2": sns.color_palette("tab10")[1],
                "Epithelial": sns.color_palette("tab10")[2],
                "Stroma3": sns.color_palette("tab10")[3],
                "Immune": sns.color_palette("tab10")[4],
                "Stroma4": sns.color_palette("tab10")[5],
                "Stroma5": sns.color_palette("tab10")[6]}


def Proportion(sid):
    adata = sc.read_h5ad(DIR+"{}.h5ad".format(sid))    # Marked ST data From 01_Main_indi_.py
    table = pd.DataFrame(adata.obs["annotation"].value_counts()).T
    table.index = [sid]
    
    def Cal_pro(ctype):
        pro1 = adata.obs[["annotation",ctype]]
        print("{}_{}%".format(ctype, round((pro1[ctype].sum()/pro1.shape[0])*100,3)))

        return round((pro1[ctype].sum()/pro1.shape[0])*100,3)

    gzmk = Cal_pro("Cd3/8_Gzmk")
    pdcd = Cal_pro("Cd3/8_Pdcd1")
    delta = Cal_pro("Gamma-deltaT")

    def Cal_emp(df, ctype):
    ## Empirical P-value dist for enrichment
    ## Estimate a significance of enrichment for tissue residential immune cells
        df = df[df.obs["annotation"]!="Immune"]
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
        nor_px = ((px["annotation"].sum()/df.obs[df.obs[ctype]>0].shape[0]) \
                - np.mean(dist)) / np.std(dist)
        plt.figure(figsize=(3,2))
        plt.axvline(nor_px)

        import scipy.stats as stats
        rv = (stats.norm(0,1))
        x = np.linspace(-5,5,100)
        y1 = rv.pdf(x)
        plt.plot(x,y1)
        
        return (nor_px, 1-rv.cdf(nor_px))

    out_p = []
    out_z = []
    for c_type in ["Gamma-deltaT","Cd3/8_Gzmk", "Cd3/8_Pdcd1"]:
        out_z.append(Cal_emp(adata, c_type)[0])
        out_p.append(Cal_emp(adata, c_type)[1])

    def Cal_(df, cell):
        df = df[df.obs["annotation"]!="Immune"]
        a = (df[df.obs[cell]>0].obs["annotation"].value_counts() / df[df.obs[cell]>0].obs["annotation"].shape[0])
        a = pd.DataFrame(a).sort_values(by="count")
    
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
            plt.title("{}".format(sam_id))
            plt.xticks([])
            plt.xlabel(cell)
            plt.savefig(FIG_OUT+"O_{}_{}.pdf".format(sam_id, cell.replace("/","_")),
                        bbox_inches="tight")
        
        Plot(a,sid)

    Cal_(adata, "Cd3/8_Gzmk")
    Cal_(adata, "Cd3/8_Pdcd1")
    Cal_(adata, "Gamma-deltaT")


    return [gzmk, pdcd, delta], out_p, out_z, table


# o1, o1p, o1z, o1t = Proportion("O1")
# o2, o2p, o2z, o2t= Proportion("O2")
y1, y1p, y1z, y1t = Proportion("Y1")
y2, y2p, y2z, y2t = Proportion("Y2")
# %%
merged_table = pd.concat([o1t,o2t,y1t,y2t])
sns.clustermap(merged_table,
               col_cluster=False,
               row_cluster=False,
               figsize=(5,2.5),
               cmap="viridis")
# merged_table.to_csv(DIR+"total_cell_type.txt",
#                     sep="\t")
print(merged_table.sum(axis=1))
# %%
merged = np.vstack([o1,o2,y1,y2])
merged_p = np.vstack([o1p,o2p,y1p,y2p])
merged_z = np.vstack([o1z,o2z,y1z,y2z])

def Draw_merged_bar(merged, *ext):
    age_color = {"Old":"#7A7A7A",
                 "Young":"#CFCFCF"}
    merged_df = pd.DataFrame(merged, 
                             index=["O1","O2","Y1","Y2"],
                             columns=["Gzmk","Pdcd1","Gamma-delta"])
    if ext[0] == "p":
        merged_df = -np.log10(merged_df)

    merged_df["group"] = ["Old","Old","Young","Young"]
    figure, ax = plt.subplots(1,3,figsize=(4,4),
                            sharey=True)
    for axs, ctype in zip(ax.flatten(), merged_df.columns[:-1]):
        print(merged_df[merged_df["group"]=="Old"][ctype])
        print(merged_df[merged_df["group"]!="Old"][ctype])
        log2FC = merged_df[merged_df["group"]=="Old"][ctype].mean() / \
                 merged_df[merged_df["group"]!="Old"][ctype].mean()
        print("{}: {}".format(ctype, log2FC))
        print(mannwhitneyu(merged_df[merged_df["group"]=="Old"][ctype],
                           merged_df[merged_df["group"]!="Old"][ctype],
                           alternative="greater"))
        # sns.barplot(x="group",
        #             y=ctype,
        #             data=merged_df,
        #             ax=axs,
        #             width=.7,
        #             order=["Young","Old"],
        #             palette=age_color,)
        sns.boxplot(x="group",
                    y=ctype,
                    data=merged_df,
                    ax=axs,
                    order=["Young","Old"],
                    palette=age_color,)
        sns.stripplot(x="group",
                      y=ctype,
                      data=merged_df,
                      ax=axs,
                      order=["Young","Old"],
                      color="#000000",
                      s=6)
        axs.set_ylabel("")
        axs.set_xlabel("")
        # axs.set_title(ctype)

    # figure.suptitle(ext[1])
    sns.despine()
    plt.ylim(bottom=0)
    plt.savefig(FIG_OUT+"{}.pdf".format(ext[1]),
                bbox_inches="tight")

# Draw_merged_bar(merged_p, "p")
# Draw_merged_bar(merged_z, "z", "Z-score")
Draw_merged_bar(merged, "np", "Proportion of Epithelial \n (Of whole spots)")
# %%
