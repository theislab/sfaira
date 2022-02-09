import os
import pandas as pd
import scanpy as sc

# This is provided in plain text on GEO
sample_map = {"Sample1": "nCoV 1 scRNA-seq",
              "Sample2": "nCoV 2 scRNA-seq",
              "Sample3": "Flu 1 scRNA-seq",
              "Sample4": "Flu 2 scRNA-seq",
              "Sample5": "Normal 1 scRNA-seq",
              "Sample6": "Flu 3 scRNA-seq",
              "Sample7": "Flu 4 scRNA-seq",
              "Sample8": "Flu 5 scRNA-seq",
              "Sample9": "nCoV 3 scRNA-seq",
              "Sample10": "nCoV 4 scRNA-seq",
              "Sample11": "nCoV 5 scRNA-seq",
              "Sample12": "nCoV 6 scRNA-seq",
              "Sample13": "Normal 2 scRNA-seq",
              "Sample14": "Normal 3 scRNA-seq",
              "Sample15": "nCoV 7 scRNA-seq",
              "Sample16": "nCoV 8 scRNA-seq",
              "Sample17": "nCoV 9 scRNA-seq",
              "Sample18": "nCoV 10 scRNA-seq",
              "Sample19": "Normal 4 scRNA-seq",
              "Sample20": "nCoV 11 scRNA-seq"
              }


def load(data_dir, sample_fn, **kwargs):
    adata = sc.read_10x_mtx(data_dir, prefix="GSE149689_")
    adata.obs["sample"] = "Sample" + adata.obs.index.str.split("-").str[1]
    adata.obs["GEO upload info"] = adata.obs["sample"].map(sample_map)

    fn_meta = os.path.join(data_dir, "Table_S1.xlsx")
    metadata = pd.read_excel(fn_meta, engine="openpyxl")
    metadata.fillna(method="ffill", inplace=True)
    metadata.replace(to_replace="(\\n)", value=" ", regex=True, inplace=True)
    metadata.rename(columns={"Experimental\nbatch": "Experimental batch"}, inplace=True)
    metadata.drop_duplicates(inplace=True)

    adata.obs = adata.obs.reset_index().merge(metadata, how="left").set_index("index")
    adata.obs["Age"] = [str(x) for x in adata.obs["Age"].values]
    return adata
