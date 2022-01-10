import anndata
import os
import pandas as pd
import scipy.io
import gzip

sex_map = {"M": "male", "F": "female"}

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

disease_map = {"severe COVID-19": "COVID-19",
               "mild COVID-19": "COVID-19",
               "severe influenza": "influenza",
               "Healthy donor": "healthy",
               "mild COVID-19 (asymptomatic)": "COVID-19"
               }

severness_map = {"severe COVID-19": "severe",
                 "mild COVID-19": "mild",
                 "severe influenza": "severe",
                 "Healthy donor": "healthy",
                 "mild COVID-19 (asymptomatic)": "asymptomatic"
                 }


def load(data_dir, sample_fn, **kwargs):
    fn = [os.path.join(data_dir, "GSE149689_matrix.mtx.gz"),
          os.path.join(data_dir, "GSE149689_barcodes.tsv.gz"),
          os.path.join(data_dir, "GSE149689_features.tsv.gz"),
          os.path.join(data_dir, "Table_S1.xlsx")]

    with gzip.open(fn[0]) as file:
        x = scipy.io.mmread(file).T.tocsr()

    obs = pd.read_csv(fn[1], header=None, sep="\t", index_col=0)
    obs.index.name = None
    var = pd.read_csv(fn[2], header=None, sep="\t")
    var.columns = ["ensembl", "names", "type"]
    var.index = var["ensembl"].values
    adata = anndata.AnnData(X=x, obs=obs, var=var.drop("type", axis=1))

    adata.obs["sample"] = "Sample" + adata.obs.index.str.split("-").str[1]
    adata.obs["GEO upload info"] = adata.obs["sample"].map(sample_map)

    metadata = pd.read_excel(fn[3], engine="openpyxl")
    metadata.fillna(method="ffill", inplace=True)
    metadata.replace(to_replace="(\\n)", value=" ", regex=True, inplace=True)
    metadata.rename(columns={"Experimental\nbatch": "Experimental batch"}, inplace=True)
    metadata.drop_duplicates(inplace=True)

    adata.obs = adata.obs.reset_index().merge(metadata, how="left").set_index("index")

    adata.obs["Age"] = [f"{int(x)}-year-old human stage" for x in adata.obs["Age"].values]
    adata.obs["Sex"] = adata.obs["Sex"].map(sex_map)
    adata.obs["severeness"] = adata.obs["Disease group"].map(severness_map)
    adata.obs["Disease group"] = adata.obs["Disease group"].map(disease_map)
    return adata
