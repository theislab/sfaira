import anndata
import os
from typing import Union

from sfaira.data import DatasetBaseGroupLoadingManyFiles

SAMPLE_FNS = [
    "tabula-muris-senis-droplet-processed-official-annotations-Fat.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-BAT.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-GAT.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-MAT.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-SCAT.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Bladder.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Bladder.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Marrow.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Brain_Myeloid.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Large_Intestine.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Large_Intestine.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Diaphragm.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Heart_and_Aorta.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Heart.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Aorta.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Kidney.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Kidney.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Liver.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Liver.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Lung.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Lung.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Mammary_Gland.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Mammary_Gland.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Pancreas.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Pancreas.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Skin.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Skin.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Spleen.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Spleen.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Thymus.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Thymus.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Tongue.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Tongue.h5ad",
    "tabula-muris-senis-droplet-processed-official-annotations-Trachea.h5ad",
    "tabula-muris-senis-facs-processed-official-annotations-Trachea.h5ad",
]


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def __init__(
            self,
            sample_fn: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, data_path=data_path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        organ = "-".join(sample_fn.split("-")[7:]).split(".")[0].lower()
        organ = "adipose tissue" if organ in ["fat", "bat", "gat", "mat", "scat"] else \
            "aorta" if organ in ["aorta"] else \
            "urinary bladder" if organ in ["bladder"] else \
            "bone marrow" if organ in ["marrow"] else \
            "brain" if organ in ["brain_non-myeloid", "brain_myeloid"] else \
            "colon" if organ in ["large_intestine"] else \
            "diaphragm" if organ in ["diaphragm"] else \
            "heart" if organ in ["heart_and_aorta", "heart"] else \
            "kidney" if organ in ["kidney"] else \
            "liver" if organ in ["liver"] else \
            "lung" if organ in ["lung"] else \
            "mammary gland" if organ in ["mammary_gland"] else \
            "muscle organ" if organ in ["limb_muscle"] else \
            "pancreas" if organ in ["pancreas"] else \
            "skin of body" if organ in ["skin"] else \
            "spleen" if organ in ["spleen"] else \
            "thymus" if organ in ["thymus"] else \
            "tongue" if organ in ["tongue"] else \
            "trachea" if organ in ["trachea"] else organ
        # ToDo: heart_and_aorta could be a distinct UBERON term, e.g. cardiovascular system?

        self.download_url_data = f"https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/{sample_fn}"
        self.download_url_meta = None

        self.obs_key_cellontology_original = "cell_ontology_class"
        self.obs_key_age = "age"
        self.obs_key_dev_stage = "development_stage"  # not given in all data sets
        self.obs_key_sex = "sex"
        # ToDo: further anatomical information for subtissue in "subtissue"?

        self.author = "Pisco"
        self.doi = "10.1101/661728"
        self.healthy = True
        self.normalization = "norm"
        self.organism = "mouse"
        self.organ = organ
        self.protocol = "10X sequencing" if sample_fn.split("-")[3] == "droplet" else "Smart-seq2"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_ensembl_col = None
        self.var_symbol_col = "index"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read_h5ad(fn)
    adata.X = adata.raw.X
    adata.var = adata.raw.var
    del adata.raw
    adata.obsm = {}
    adata.varm = {}
    adata.uns = {}

    return adata
