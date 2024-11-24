# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

var="LuminalA_Luminal-HS"
DIR2 = "/Users/kangh/Downloads/review/{}/".format(var)
test = pd.read_csv(DIR2+"merged.txt",
                   sep="\t", index_col="Term")  # Made by merged observed and simulated data
test = test.fillna(0.0)
test = test.iloc[:10,:] # Using top10

## P-value check
for i in test.index.tolist():
    plt.figure(figsize=(3,3))
    temp = test.loc[i]
    ## Norm
    temp = (temp-temp.mean()) / temp.std()
    target = temp["Observed"]
    p = ((temp>=target).astype(int).sum()/1000)
    sns.histplot(temp, bins=100)
    plt.axvline(target, linewidth=3,
                color="#FF1C1C")
    plt.title("{} P:{}".format(i,p), fontsize=13)
    sns.despine()
    plt.xlabel(i, fontsize=13)
    plt.ylabel("")
    plt.savefig(DIR2+"{}.pdf".format(i.replace("/","_")),
                bbox_inches="tight")

## Make input frame
df = test.T
df["group"] = df.index
df["group"] = df["group"].apply(lambda x: "Random" if x.startswith("Rand") else "Observed")
df_melt = pd.melt(df, value_vars=df.columns[:-1], var_name='melt', id_vars="group")
plt.figure(figsize=(6,4))
ax = sns.stripplot(x=df_melt[df_melt["group"]!="Observed"]["melt"],
                   y=df_melt[df_melt["group"]!="Observed"]["value"],
                   size=2,
                   alpha=.6,
                   hue=df_melt[df_melt["group"]!="Observed"]["group"],
                   palette=["#D1D1D1"],
                   )
ax = sns.stripplot(x=df_melt[df_melt["group"]=="Observed"]["melt"],
                   y=df_melt[df_melt["group"]=="Observed"]["value"],
                   size=6,
                   hue=df_melt[df_melt["group"]=="Observed"]["group"],
                   palette=["#FF5D51"],
                   )
plt.xticks(rotation=90, fontsize=11)
sns.despine()
plt.xlabel("Terms", fontsize=13)
plt.yticks(fontsize=11)
plt.ylabel("-Log10(p)", fontsize=13)
plt.legend(frameon=False)
plt.title(var, fontsize=14)
plt.savefig(DIR2+"Pvalue_comp.pdf",
            bbox_inches="tight")
# %%
