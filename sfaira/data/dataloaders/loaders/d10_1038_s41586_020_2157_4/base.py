import anndata
import numpy as np
import os
import pandas as pd
import scipy.sparse
from typing import Union
import zipfile

from sfaira.data import DatasetBaseGroupLoadingOneFile

SAMPLE_IDS = [
    'AdultAdipose_1',
    'AdultAdrenalGland_2',
    'AdultAdrenalGland_3',
    'AdultArtery_1',
    'AdultAscendingColon_1',
    'AdultBladder_1',
    'AdultBladder_2',
    'AdultCerebellum_1',
    'AdultCervix_1',
    'AdultColon_1',
    'AdultDuodenum_1',
    'AdultEpityphlon_1',
    'AdultEsophagus_1',
    'AdultEsophagus_2',
    'AdultFallopiantube_1',
    'AdultGallbladder_1',
    'AdultGallbladder_2',
    'AdultHeart_1',
    'AdultHeart_2',
    'AdultIleum_2',
    'AdultJejunum_2',
    'AdultKidney_2',
    'AdultKidney_3',
    'AdultKidney_4',
    'AdultLiver_1',
    'AdultLiver_2',
    'AdultLiver_4',
    'AdultLung_1',
    'AdultLung_2',
    'AdultLung_3',
    'AdultMuscle_1',
    'AdultOmentum_1',
    'AdultOmentum_2',
    'AdultOmentum_3',
    'AdultPancreas_1',
    'AdultPeripheralBlood_3',
    'AdultPeripheralBlood_4',
    'AdultPleura_1',
    'AdultProstate_1',
    'AdultRectum_1',
    'AdultSigmoidColon_1',
    'AdultSpleenParenchyma_1',
    'AdultSpleen_1',
    'AdultStomach_1',
    'AdultStomach_2',
    'AdultStomach_3',
    'AdultTemporalLobe_1',
    'AdultThyroid_1',
    'AdultThyroid_2',
    'AdultTrachea_2',
    'AdultTransverseColon_2',
    'AdultUreter_1',
    'AdultUterus_1',
    'BoneMarrow_1',
    'BoneMarrow_2',
    'ChorionicVillus_1',
    'CordBloodCD34P_1',
    'CordBloodCD34P_2',
    'CordBlood_1',
    'CordBlood_2',
    'FetalAdrenalGland_2',
    'FetalAdrenalGland_3',
    'FetalAdrenalGland_4',
    'FetalBrain_3',
    'FetalBrain_4',
    'FetalBrain_5',
    'FetalBrain_6',
    'FetalCalvaria_1',
    'FetalEyes_1',
    'FetalFemaleGonad_1',
    'FetalFemaleGonad_2',
    'FetalHeart_1',
    'FetalHeart_2',
    'FetalIntestine_1',
    'FetalIntestine_2',
    'FetalIntestine_3',
    'FetalIntestine_4',
    'FetalIntestine_5',
    'FetalKidney_3',
    'FetalKidney_4',
    'FetalKidney_5',
    'FetalKidney_6',
    'FetalLung_1',
    'FetalLung_2',
    'FetalMaleGonad_1',
    'FetalMaleGonad_2',
    'FetalMuscle_1',
    'FetalPancreas_1',
    'FetalPancreas_2',
    'FetalPancreas_3',
    'FetalRib_2',
    'FetalRib_3',
    'FetalSkin_2',
    'FetalSkin_3',
    'FetalSpinalCord_1',
    'FetalStomach_1',
    'FetalStomach_2',
    'FetalThymus_1',
    'FetalThymus_2',
    'HESC_1',
    'Liver_1',
    'Liver_2',
    'NeonatalAdrenalGland_1',
    'PeripheralBlood_1',
    'Placenta_1'
]


class Dataset(DatasetBaseGroupLoadingOneFile):

    def __init__(
            self,
            sample_id: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):

        super().__init__(
            sample_id=sample_id,
            data_path=data_path,
            meta_path=meta_path,
            cache_path=cache_path,
            **kwargs
        )

        sample_organ_dict = {
            'AdultAdipose_1': 'adipose tissue of abdominal region',
            'AdultAdrenalGland_2': 'adrenal gland',
            'AdultAdrenalGland_3': 'adrenal gland',
            'AdultArtery_1': 'artery',
            'AdultAscendingColon_1': 'ascending colon',
            'AdultBladder_1': 'urinary bladder',
            'AdultBladder_2': 'urinary bladder',
            'AdultCerebellum_1': 'cerebellum',
            'AdultCervix_1': 'uterine cervix',
            'AdultColon_1': 'colon',
            'AdultDuodenum_1': 'duodenum',
            'AdultEpityphlon_1': 'caecum',
            'AdultEsophagus_1': 'esophagus',
            'AdultEsophagus_2': 'esophagus',
            'AdultFallopiantube_1': 'fallopian tube',
            'AdultGallbladder_1': 'gall bladder',
            'AdultGallbladder_2': 'gall bladder',
            'AdultHeart_1': 'heart',
            'AdultHeart_2': 'heart',
            'AdultIleum_2': 'ileum',
            'AdultJejunum_2': 'jejunum',
            'AdultKidney_2': 'kidney',
            'AdultKidney_3': 'kidney',
            'AdultKidney_4': 'kidney',
            'AdultLiver_1': 'liver',
            'AdultLiver_2': 'liver',
            'AdultLiver_4': 'liver',
            'AdultLung_1': 'lung',
            'AdultLung_2': 'lung',
            'AdultLung_3': 'lung',
            'AdultMuscle_1': 'skeletal muscle organ',
            'AdultOmentum_1': 'omentum',
            'AdultOmentum_2': 'omentum',
            'AdultOmentum_3': 'omentum',
            'AdultPancreas_1': 'pancreas',
            'AdultPeripheralBlood_3': 'blood',
            'AdultPeripheralBlood_4': 'blood',
            'AdultPleura_1': 'pleura',
            'AdultProstate_1': 'prostate gland',
            'AdultRectum_1': 'rectum',
            'AdultSigmoidColon_1': 'sigmoid colon',
            'AdultSpleenParenchyma_1': 'parenchyma of spleen',
            'AdultSpleen_1': 'spleen',
            'AdultStomach_1': 'stomach',
            'AdultStomach_2': 'stomach',
            'AdultStomach_3': 'stomach',
            'AdultTemporalLobe_1': 'temporal lobe',
            'AdultThyroid_1': 'thyroid gland',
            'AdultThyroid_2': 'thyroid gland',
            'AdultTrachea_2': 'trachea',
            'AdultTransverseColon_2': 'transverse colon',
            'AdultUreter_1': 'ureter',
            'AdultUterus_1': 'uterus',
            'BoneMarrow_1': 'bone marrow',
            'BoneMarrow_2': 'bone marrow',
            'ChorionicVillus_1': 'chorionic villus',
            'CordBloodCD34P_1': 'umbilical cord blood',
            'CordBloodCD34P_2': 'umbilical cord blood',
            'CordBlood_1': 'umbilical cord blood',
            'CordBlood_2': 'umbilical cord blood',
            'FetalAdrenalGland_2': 'adrenal gland',
            'FetalAdrenalGland_3': 'adrenal gland',
            'FetalAdrenalGland_4': 'adrenal gland',
            'FetalBrain_3': 'brain',
            'FetalBrain_4': 'brain',
            'FetalBrain_5': 'brain',
            'FetalBrain_6': 'brain',
            'FetalCalvaria_1': 'vault of skull',
            'FetalEyes_1': 'eye',
            'FetalFemaleGonad_1': 'ovary',
            'FetalFemaleGonad_2': 'ovary',
            'FetalHeart_1': 'heart',
            'FetalHeart_2': 'heart',
            'FetalIntestine_1': 'intestine',
            'FetalIntestine_2': 'intestine',
            'FetalIntestine_3': 'intestine',
            'FetalIntestine_4': 'intestine',
            'FetalIntestine_5': 'intestine',
            'FetalKidney_3': 'kidney',
            'FetalKidney_4': 'kidney',
            'FetalKidney_5': 'kidney',
            'FetalKidney_6': 'kidney',
            'FetalLung_1': 'lung',
            'FetalLung_2': 'lung',
            'FetalMaleGonad_1': 'testis',
            'FetalMaleGonad_2': 'testis',
            'FetalMuscle_1': 'skeletal muscle organ',
            'FetalPancreas_1': 'pancreas',
            'FetalPancreas_2': 'pancreas',
            'FetalPancreas_3': 'pancreas',
            'FetalRib_2': 'rib',
            'FetalRib_3': 'rib',
            'FetalSkin_2': 'skin of body',
            'FetalSkin_3': 'skin of body',
            'FetalSpinalCord_1': 'spinal cord',
            'FetalStomach_1': 'stomach',
            'FetalStomach_2': 'stomach',
            'FetalThymus_1': 'thymus',
            'FetalThymus_2': 'thymus',
            'HESC_1': '',
            'Liver_1': 'liver',
            'Liver_2': 'liver',
            'NeonatalAdrenalGland_1': 'adrenal gland',
            'PeripheralBlood_1': 'blood',
            'Placenta_1': 'placenta',
        }

        self.download_url_data = "https://ndownloader.figshare.com/files/17727365"
        self.download_url_meta = [
            "https://ndownloader.figshare.com/files/21758835",
            "https://ndownloader.figshare.com/files/22447898",
        ]

        self.obs_key_sample = "sample"

        self.organ = sample_organ_dict[self.sample_id]
        self.id = f"human_{self.organ}_2020_microwell_han_{str(SAMPLE_IDS.index(self.sample_id)).zfill(3)}" \
                  f"_10.1038/s41586-020-2157-4"

        self.author = "Guo"
        self.doi = "10.1038/s41586-020-2157-4"
        self.healthy = True
        self.normalization = "raw"
        self.organism = "human"
        self.protocol = "microwell-seq"
        self.state_exact = "healthy"
        self.year = 2020

        self.obs_key_cellontology_original = "cell_ontology_class"
        self.obs_key_dev_stage = "dev_stage"
        self.obs_key_sex = "gender"
        self.obs_key_age = "age"

        self.var_symbol_col = "index"

    def _load_full(self):
        self.adata = anndata.read(os.path.join(self.path, "human", self.directory_formatted_doi, "HCL_Fig1_self.adata.h5ad"))
        # convert to sparse matrix
        self.adata.X = scipy.sparse.csr_matrix(self.adata.X).copy()

        # harmonise annotations
        for col in ["batch", "tissue"]:
            self.adata.obs[col] = self.adata.obs[col].astype("str")
        self.adata.obs.index = self.adata.obs.index.str.replace("AdultJeJunum", "AdultJejunum", regex=True).str.replace(
            "AdultGallBladder", "AdultGallbladder", regex=True).str.replace(
            "FetalFemaleGonald", "FetalFemaleGonad", regex=True)
        self.adata.obs.replace({"AdultJeJunum": "AdultJejunum", "AdultGallBladder": "AdultGallbladder",
                           "FetalFemaleGonald": "FetalFemaleGonad"}, regex=True, inplace=True)
        self.adata.obs.index = ["-".join(i.split("-")[:-1]) for i in self.adata.obs.index]

        # load celltype labels and harmonise them
        # This pandas code should work with pandas 1.2 but it does not and yields an empty data frame:
        fig1_anno = pd.read_excel(
            os.path.join(self.data_dir_base, "human", self.directory_formatted_doi, "HCL_Fig1_cell_Info.xlsx"),
            index_col="cellnames",
            engine="xlrd",  # ToDo: Update when pandas xlsx reading with openpyxl is fixed: yields empty tables
        )
        fig1_anno.index = fig1_anno.index.str.replace("AdultJeJunum", "AdultJejunum", regex=True).str.replace(
            "AdultGallBladder", "AdultGallbladder", regex=True).str.replace(
            "FetalFemaleGonald", "FetalFemaleGonad", regex=True)

        # check that the order of cells and cell labels is the same
        assert np.all(fig1_anno.index == self.adata.obs.index)

        # add annotations to self.adata object and rename columns
        self.adata.obs = pd.concat([self.adata.obs, fig1_anno[["cluster", "stage", "donor", "celltype"]]], axis=1)
        self.adata.obs.columns = ["sample", "tissue", "n_genes", "n_counts", "cluster_global", "stage", "donor",
                             "celltype_global"]

        # add sample-wise annotations to the full self.adata object
        df = pd.DataFrame(
            columns=["Cell_barcode", "Sample", "Batch", "Cell_id", "Cluster_id", "Ages", "Development_stage", "Method",
                     "Gender", "Source", "Biomaterial", "Name", "ident", "Celltype"])
        archive = zipfile.ZipFile(
            os.path.join(self.data_dir_base, "human", self.directory_formatted_doi, "annotation_rmbatch_data_revised417.zip")
        )
        for f in archive.namelist():
            df1 = pd.read_csv(archive.open(f), encoding="unicode_escape")
            df = pd.concat([df, df1], sort=True)
        df = df.set_index("Cell_id")
        self.adata = self.adata[[i in df.index for i in self.adata.obs.index]].copy()
        a_idx = self.adata.obs.index.copy()
        self.adata.obs = pd.concat([self.adata.obs, df[["Ages", "Celltype", "Cluster_id", "Gender", "Method", "Source"]]], axis=1)
        assert np.all(a_idx == self.adata.obs.index)

        # remove mouse cells from the object  # ToDo: add this back in as mouse data sets?
        self.adata = self.adata[self.adata.obs["Source"] != "MCA2.0"].copy()

        # tidy up the column names of the obs annotations
        self.adata.obs.columns = ["sample", "sub_tissue", "n_genes", "n_counts", "cluster_global", "dev_stage",
                             "donor", "celltype_global", "age", "celltype_specific", "cluster_specific", "gender",
                             "protocol", "source"]
