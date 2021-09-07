import anndata
import os

from sfaira.data import DatasetBase

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


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        organ = "-".join(self.sample_fn.split("-")[7:]).split(".")[0].lower()
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

        self.download_url_data = f"https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/" \
                                 f"{self.sample_fn}"
        self.download_url_meta = None

        self.cell_type_obs_key = "cell_ontology_class"
        self.development_stage_obs_key = "development_stage"
        self.sex_obs_key = "sex"
        # ToDo: further anatomical information for subtissue in "subtissue"?

        self.author = "Pisco"
        self.disease = "healthy"
        self.doi_journal = "10.1038/s41586-020-2496-1"
        self.doi_preprint = "10.1101/661728"
        self.normalization = "norm"
        self.organism = "mouse"
        self.organ = organ
        self.assay_sc = "10x 3' v2" if self.sample_fn.split("-")[3] == "droplet" else "Smart-seq2"
        self.year = 2019
        self.sample_source = "primary_tissue"

        self.gene_id_ensembl_var_key = None
        self.gene_id_symbols_var_key = "index"

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    dev_stage_dict = {
        "18m": "18 month-old stage",
        "1m": "4 weeks",
        "21m": "20 month-old stage and over",
        "24m": "20 month-old stage and over",
        "30m": "20 month-old stage and over",
        "3m": "3 month-old stage",
    }
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read_h5ad(fn)
    adata.X = adata.raw.X
    adata.var = adata.raw.var
    del adata.raw
    adata.obsm = {}
    adata.varm = {}
    adata.uns = {}
    adata.obs['development_stage'] = [dev_stage_dict[i] for i in adata.obs['age']]

    return adata
