# %%
import scanpy as sc
import senepy as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({
'axes.titlesize': 15,     # 제목 글꼴 크기
'axes.labelsize': 14,     # x, y축 라벨 글꼴 크기
'xtick.labelsize': 12,    # x축 틱 라벨 글꼴 크기
'ytick.labelsize': 12,    # y축 틱 라벨 글꼴 크기
'legend.fontsize': 12,    # 범례 글꼴 크기
'figure.titlesize': 15    # figure 제목 글꼴 크기
})
## Load our Fibroblast data
DIR = "/home/kangh/lab-server/script/10X_cellranger/analysis_python/"
adata = sc.read_h5ad(DIR+"Fibro_RNA.h5ad")
print(adata.obs["seurat_clusters"].unique())
def Clu_corr(df):
    if df["seurat_clusters"] == "0":
        return "0"
    elif df["seurat_clusters"] == "6":
        return "1"
    elif df["seurat_clusters"] == "1":
        return "2"
    elif df["seurat_clusters"] == "4":
        return "3"
    elif df["seurat_clusters"] == "3":
        return "4"
    elif df["seurat_clusters"] == "7":
        return "5"
    elif df["seurat_clusters"] == "5":
        return "6"
    elif df["seurat_clusters"] == "10":
        return "7"
    elif df["seurat_clusters"] == "2":
        return "8"
    elif df["seurat_clusters"] == "8":
        return "9"
    elif df["seurat_clusters"] == "9":
        return "10"

adata.obs["new_seurat_clusters"] = adata.obs.apply(Clu_corr, axis=1)
print(adata.obs[["seurat_clusters","new_seurat_clusters"]].head(15))
# %%
## Load Senepy DB (Mouse)
hubs = sp.load_hubs(species = 'Mouse')
hubs.search_hubs_by_genes(gene_set)
# %%
## Make subset of 
filt_meta = hubs.metadata[hubs.metadata.cell == 'fibroblast']
filt_meta = filt_meta[filt_meta["hyp"]<0.05]  # Hyp represents the hypogeometric p-value of the 'known' and novel hub genes
hubs.merge_hubs(filt_meta, new_name = 'Fib')
hubs.hubs['Fib'][:10]
hubs.merge_hubs(filt_meta, new_name = 'Fib_min2', overlap_threshold = 2)
# %%
## Change gene symbol and select fib clusters
translator = sp.translator(hub = hubs.hubs, data = adata)
fib_subset = adata[adata.obs['new_seurat_clusters'].isin(['0','1','2','3','4','5'])].copy()
# %%
## Make gene list for Supp table
pd.DataFrame(set(hubs.get_genes("Fib")) & set(fib_subset.var.index)).to_csv(DIR+"SenePy_genelist.txt",
                                                                            sep="\t", index=None)
print(fib_subset.obs["new_seurat_clusters"].unique())
# %%
## Scaling seneScores
from sklearn.preprocessing import MinMaxScaler
## Using certain hub data to calculate senescnece score
fib_subset.obs['sen_score'] = sp.score_hub(fib_subset, hubs.hubs["Fib"])
sclaer = MinMaxScaler()
print(fib_subset.obs)
# %%
fib_subset.obs
input_df = fib_subset.obs.copy()
input_df["norm_sen_score"] = sclaer.fit_transform(input_df[["sen_score"]])
# %%
def Make_Fig(input_df,score,title):
    from scipy.stats import mannwhitneyu, ttest_ind
    plt.figure(figsize=(5,4))
    for i in input_df["new_seurat_clusters"].unique():
        sub = input_df[input_df["new_seurat_clusters"]==i]
        u,p = mannwhitneyu(sub[sub["orig.ident"]=="mm10_3mths"]["norm_sen_score"],
                           sub[sub["orig.ident"]!="mm10_3mths"]["norm_sen_score"],
                           alternative="less")

        sub = sub.sort_values(by="orig.ident")

        ax = sns.violinplot(x="new_seurat_clusters",
                            y=score,
                            hue="orig.ident",
                            hue_order=["mm10_3mths", "mm10_18mths"],
                            data=sub,
                            inner=None)
        
        ax = sns.stripplot(x="new_seurat_clusters",
                            y=score,
                            hue="orig.ident",
                            hue_order=["mm10_3mths", "mm10_18mths"],
                            data=sub,
                            dodge=True,
                            palette={"mm10_3mths":"#000000",
                                     "mm10_18mths":"#000000"},
                            s=3)
                            

        if p > 0.05:
            star = "ns"
        elif p < 0.05 and p > 0.01:
            star = "*"
        elif p < 0.01 and p > 0.001:
            star = "**"
        elif p < 0.001 and p > 0.0001:
            star = "***"
        elif p < 0.0001:
            star = "****"

        plt.text(0.17*(int(i)+1)-0.12, 0.95,
                "{}".format(star),
                transform=ax.transAxes,
                fontsize=10)
        sns.despine()
        plt.legend(frameon=False,
                labels=["mm10_3mths", "mm10_18mths"],
                bbox_to_anchor=(1.02,1.02),
                loc="upper left")
    plt.tick_params(bottom=False)  
    plt.xticks([0,1,2,3,4,5],
            ["3M 18M\nFib-C0","3M 18M\nFib-C1","3M 18M\nFib-C2","3M 18M\nFib-C3","3M 18M\nFib-C4","3M 18M\nFib-C5"],
            )
    plt.xlabel('')
    plt.title(title)
    plt.savefig(DIR+"{}.pdf".format(title),
                bbox_inches="tight")


# Make_Fig(df,"mean_score",geneset_name)
Make_Fig(input_df,"norm_sen_score","Senepy")
# %%
import numpy as np
import matplotlib.pyplot as plt
from math import pi

# Data
categories = ['Sen', '1/red', 'FDR', 'nrPre', 'Pre', 'PDR']
values = [75, 50, 60, 70, 80, 65]

# Number of variables
N = len(categories)

# What will be the angle of each axis in the plot?
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# Initialise the spider plot
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# Draw one axe per variable and add labels
plt.xticks(angles[:-1], categories)

# Draw ylabels
ax.set_rscale('log')
plt.yticks([10, 30, 60, 90], ["10", "30", "60", "90"], color="grey", size=7)
plt.ylim(0, 100)

# Plot data
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid')

# Fill area
ax.fill(angles, values, 'b', alpha=0.1)

# Customizing the grid to make it hexagonal
ax.set_rgrids([10, 30, 60, 90], angle=angles[0], labels=["10", "30", "60", "90"])
ax.set_thetagrids(np.degrees(angles[:-1]), labels=categories)
ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

# Customizing the gridlines to be hexagonal
for gridline in ax.yaxis.get_gridlines():
    gridline.set_linestyle('-')
    gridline.set_color('grey')
    gridline.set_linewidth(0.5)

# Ensure gridlines are hexagonal
ax.xaxis.grid(True, linestyle='-', color='grey', linewidth=0.5)
ax.yaxis.grid(True, linestyle='-', color='grey', linewidth=0.5)

plt.title('RNA-Bloom')
plt.show()

# %%
import numpy as np
import matplotlib.pyplot as plt
from math import pi

# Data
categories = ['Sen', '1/red', 'FDR', 'nrPre', 'Pre', 'PDR']
values = [75, 50, 60, 70, 80, 65]

# Number of variables
N = len(categories)

# What will be the angle of each axis in the plot?
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# Initialise the spider plot
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# Draw one axe per variable and add labels
plt.xticks(angles[:-1], categories)

# Draw ylabels
ax.set_rscale('log')
plt.yticks([10, 30, 60, 90], ["10", "30", "60", "90"], color="grey", size=7)
plt.ylim(0, 100)

# Plot data
values += values[:1]
ax.plot(angles, values, linewidth=1, linestyle='solid')

# Fill area
ax.fill(angles, values, 'b', alpha=0.1)

plt.title('RNA-Bloom')
plt.show()

# %%
