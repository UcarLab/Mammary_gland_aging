# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
plt.rcParams["font.family"] = "Arial"
DIR="/Users/kangh/Downloads/"
data = pd.read_csv(DIR+"luminalav_luminala.tsv",
                 sep="\t")
excel_gene = data["Mouse Gene Name"].str.upper().tolist()

def GO(target_list, cutoff, term):
    key_set = {"mf":"GO_Molecular_Function_2021",
               "bp":"GO_Biological_Process_2018", 
               "react":"Reactome_2022",
               "Msig":"MSigDB_Hallmark_2020",
               "cc":"GO_Cellular_Component_2021"}

    import gseapy as gp
    enr = gp.enrichr(gene_list=target_list,
                    gene_sets=key_set[term],  #GO_Biological_Process_2023, MSigDB_Hallmark_2020, GO_Molecular_Function_2023, GO_Cellular_Component_2023
                    organism="human", # don't forget to set organism to the one you desired! e.g. Yeast
                    # outdir=DIR+"GO_{}".format(title), # don't write to disk
                    cutoff=cutoff
                    )
    
    input = enr.results[enr.results["Adjusted P-value"]<cutoff][["Term", "Adjusted P-value","Combined Score","Genes"]]
    input["-log(q)"] = -np.log10(input["Adjusted P-value"])
    input["Term"] = input["Term"].str.split("(").str[0]

    return input


def Data(cancer, single):
    DIR="/Users/kangh/Downloads/"
    OUT=DIR+"Brittany_review/{}_{}/".format(cancer,single)
    if not os.path.exists(OUT):
        os.mkdir(OUT)

    cancer_deg = pd.read_csv(DIR+"TCGA_{}vsNormal_DE.csv".format(cancer),
                            sep=",")
    sig = cancer_deg[(cancer_deg["padj"]<=0.05) &
                     (np.abs(cancer_deg["log2FoldChange"])>=0.5)]
    up_sig = sig[sig["direction"]=="up"][["external_gene_name"]]
    down_sig = sig[sig["direction"]=="down"][["external_gene_name"]]
    frac = (sig.shape[0]/cancer_deg.shape[0])   # Number of cancer DEGs

    ## Make sc Sig genes
    sc_deg = pd.read_csv(DIR+"scdeg_whole.txt",
                        sep="\t")  # pseudo bulk
    sc_deg = sc_deg[(sc_deg["Significant_in"]=="Pseudobulk") |
                    (sc_deg["Significant_in"]=="Single_cell") |
                    (sc_deg["Significant_in"]=="Both")] # Union gene set
    
    sig_sc = sc_deg[sc_deg["Cluster"]==single]
    sig_sc["Gene Name"] = sig_sc["Gene Name"].str.upper()
    up_sig_sc = sig_sc[sig_sc["Log2 Fold Change"]>0][["Gene Name"]]
    down_sig_sc = sig_sc[sig_sc["Log2 Fold Change"]<0][["Gene Name"]]
    
    ## Real related genes (Observed)
    sim = []
    shared_list = list(set(list(set(up_sig_sc["Gene Name"]) &
                                set(up_sig["external_gene_name"]))))
    shared_list2 = list(set(list(set(down_sig_sc["Gene Name"]) &
                                 set(down_sig["external_gene_name"]))))
    obs = shared_list + shared_list2
    obs_GO = GO(obs, 0.01, "Msig")
    obs_GO = obs_GO[obs_GO["Adjusted P-value"]<0.01][["Term","-log(q)"]].set_index("Term")
    obs_GO.columns = ["Observed"]
    obs_GO.to_csv(OUT+"Obs.txt",
                       sep="\t")
    obs = [len(obs)]

    ## Make random cancer gene with same size of original (Permutation)
    for i in np.random.randint(1,1000000,1000):
        rand = cancer_deg.sample(frac=frac, random_state=i)
        up_rand = rand[rand["direction"]=="up"][["external_gene_name"]]
        down_rand = rand[rand["direction"]=="down"][["external_gene_name"]]
        f_shared_up = list(set(up_sig_sc["Gene Name"]) &
                           set(up_rand["external_gene_name"]))
        f_shared_down = list(set(down_sig_sc["Gene Name"]) &
                           set(down_rand["external_gene_name"]))
        f_shared_list = list(set(f_shared_up + f_shared_down))
        rand_GO = GO(f_shared_list, 1.0, "Msig")
        rand_GO = rand_GO[rand_GO["Term"].isin(obs_GO.index.tolist())][["Term","-log(q)"]].set_index("Term")
        rand_GO.columns = ["Rand{}".format(i)]
        
        ## Do something to increase speed
        rand_GO = pd.merge(obs_GO, rand_GO, left_index=True, right_index=True, how="left").drop("Observed",axis=1)
        rand_GO = rand_GO.reindex(obs_GO.index.tolist())    # Sort by observed data index
        rand_GO.to_csv(OUT+"Rand{}.txt".format(i),
                       sep="\t", index=None)


# Data("LuminalA","Luminal-AV")
# Data("LuminalA","Luminal-HS")
# Data("LuminalB","Luminal-AV")
Data("LuminalB","Luminal-HS")
# %%
DIR2 = "/Users/kangh/Downloads/Brittany_review/LuminalA_Luminal-AV/"
test = pd.read_csv(DIR2+"merged.txt",
                   sep="\t")
# print(test.fillna(0.0))
sns.clustermap(test,
               row_cluster=False,
               col_cluster=False,
               standard_scale=0,
               cmap="viridis",
               figsize=(10,4))
# %%
