import os
import gzip
import anndata
import shutil


def load(data_dir, **kwargs):

    sex_dict = {
        "F": "female",
        "M": "male"
    }

    dev_stage_dict = {
        72: "11th week post-fertilization homosapiens stage",
        74: "11th week post-fertilization homosapiens stage",
        85: "13th week post-fertilization homosapiens stage",
        89: "13th week post-fertilization homosapiens stage",
        90: "13th week post-fertilization homosapiens stage",
        94: "14th week post-fertilization homosapiens stage",
        96: "14th week post-fertilization homosapiens stage",
        100: "15th week post-fertilization homosapiens stage",
        110: "16th week post-fertilization homosapiens stage",
        112: "17th week post-fertilization homosapiens stage",
        113: "17th week post-fertilization homosapiens stage",
        115: "17th week post-fertilization homosapiens stage",
        117: "17th week post-fertilization homosapiens stage",
        119: "18th week post-fertilization homosapiens stage",
        120: "18th week post-fertilization homosapiens stage",
        122: "18th week post-fertilization homosapiens stage",
        125: "18th week post-fertilization homosapiens stage",
        129: "19th week post-fertilization homosapiens stage",
    }

    organ_dict = {
        "Adrenal": "adrenal gland",
        "Cerebellum": "cerebellum",
        "Cerebrum": "telencephalon",
        "Eye": "eye",
        "Heart": "heart",
        "Intestine": "intestine",
        "Kidney": "kidney",
        "Liver": "liver",
        "Lung": "lung",
        "Muscle": "muscle organ",
        "Pancreas": "pancreas",
        "Placenta": "placenta",
        "Spleen": "spleen",
        "Stomach": "stomach",
        "Thymus": "thymus",
    }

    fn = os.path.join(data_dir, "GSE156793_S3_gene_count.loom.gz")
    fn_tmp = os.path.join(os.path.expanduser("~"), "tmp.loom")
    with gzip.open(fn, 'rb') as f_in:
        with open(fn_tmp, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    adata = anndata.read_loom(fn_tmp)
    os.remove(fn_tmp)

    adata.obs["Sex"] = [sex_dict[x] for x in adata.obs["Sex"]]
    adata.obs["Organ"] = [organ_dict[x] for x in adata.obs["Organ"]]
    adata.obs["Developmental_stage"] = [dev_stage_dict[x] for x in adata.obs["Development_day"]]

    return adata
