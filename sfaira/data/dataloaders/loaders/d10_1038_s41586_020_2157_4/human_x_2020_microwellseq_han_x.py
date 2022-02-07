import anndata
import numpy as np
import os
import pandas as pd
import scipy.sparse
import zipfile

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://ndownloader.figshare.com/files/17727365"
        self.download_url_meta = [
            "https://ndownloader.figshare.com/files/21758835",
            "https://ndownloader.figshare.com/files/22447898",
        ]

        self.author = "Han"
        self.doi_journal = "10.1038/s41586-020-2157-4"
        self.healthy = True
        self.layer_counts = "X"
        self.organism = "Homo sapiens"
        self.assay_sc = "microwell-seq"
        self.state_exact = "healthy"
        self.year = 2020
        self.sample_source = "primary_tissue"

        self.bio_sample_obs_key = "sample"
        self.cell_type_obs_key = "celltype_specific"
        self.development_stage_obs_key = "dev_stage"
        self.organ_obs_key = "organ"
        self.primary_data = True
        self.sex_obs_key = "sex"
        self.age_obs_key = "age"

        self.feature_symbol_var_key = "index"
        self.feature_type = "rna"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
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
        'FetalIntetsine_3': 'intestine',
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
        'HESC_1': 'blastocyst',
        'Liver_1': 'liver',
        'Liver_2': 'liver',
        'NeonatalAdrenalGland_1': 'adrenal gland',
        'PeripheralBlood_1': 'blood',
        'Placenta_1': 'placenta',
    }
    sample_dev_stage_dict = {
        'AdultAdipose_1': 'homosapiens adult stage',
        'AdultAdrenalGland_2': 'homosapiens adult stage',
        'AdultAdrenalGland_3': 'homosapiens adult stage',
        'AdultArtery_1': 'homosapiens adult stage',
        'AdultAscendingColon_1': 'homosapiens adult stage',
        'AdultBladder_1': 'homosapiens adult stage',
        'AdultBladder_2': 'homosapiens adult stage',
        'AdultCerebellum_1': 'homosapiens adult stage',
        'AdultCervix_1': 'homosapiens adult stage',
        'AdultColon_1': 'homosapiens adult stage',
        'AdultDuodenum_1': 'homosapiens adult stage',
        'AdultEpityphlon_1': 'homosapiens adult stage',
        'AdultEsophagus_1': 'homosapiens adult stage',
        'AdultEsophagus_2': 'homosapiens adult stage',
        'AdultFallopiantube_1': 'homosapiens adult stage',
        'AdultGallbladder_1': 'homosapiens adult stage',
        'AdultGallbladder_2': 'homosapiens adult stage',
        'AdultHeart_1': 'homosapiens adult stage',
        'AdultHeart_2': 'homosapiens adult stage',
        'AdultIleum_2': 'homosapiens adult stage',
        'AdultJejunum_2': 'homosapiens adult stage',
        'AdultKidney_2': 'homosapiens adult stage',
        'AdultKidney_3': 'homosapiens adult stage',
        'AdultKidney_4': 'homosapiens adult stage',
        'AdultLiver_1': 'homosapiens adult stage',
        'AdultLiver_2': 'homosapiens adult stage',
        'AdultLiver_4': 'homosapiens adult stage',
        'AdultLung_1': 'homosapiens adult stage',
        'AdultLung_2': 'homosapiens adult stage',
        'AdultLung_3': 'homosapiens adult stage',
        'AdultMuscle_1': 'homosapiens adult stage',
        'AdultOmentum_1': 'homosapiens adult stage',
        'AdultOmentum_2': 'homosapiens adult stage',
        'AdultOmentum_3': 'homosapiens adult stage',
        'AdultPancreas_1': 'homosapiens adult stage',
        'AdultPeripheralBlood_3': 'homosapiens adult stage',
        'AdultPeripheralBlood_4': 'homosapiens adult stage',
        'AdultPleura_1': 'homosapiens adult stage',
        'AdultProstate_1': 'homosapiens adult stage',
        'AdultRectum_1': 'homosapiens adult stage',
        'AdultSigmoidColon_1': 'homosapiens adult stage',
        'AdultSpleenParenchyma_1': 'homosapiens adult stage',
        'AdultSpleen_1': 'homosapiens adult stage',
        'AdultStomach_1': 'homosapiens adult stage',
        'AdultStomach_2': 'homosapiens adult stage',
        'AdultStomach_3': 'homosapiens adult stage',
        'AdultTemporalLobe_1': 'homosapiens adult stage',
        'AdultThyroid_1': 'homosapiens adult stage',
        'AdultThyroid_2': 'homosapiens adult stage',
        'AdultTrachea_2': 'homosapiens adult stage',
        'AdultTransverseColon_2': 'homosapiens adult stage',
        'AdultUreter_1': 'homosapiens adult stage',
        'AdultUterus_1': 'homosapiens adult stage',
        'BoneMarrow_1': 'homosapiens adult stage',
        'BoneMarrow_2': 'homosapiens adult stage',
        'ChorionicVillus_1': 'homosapiens adult stage',
        'CordBloodCD34P_1': 'homosapiens adult stage',
        'CordBloodCD34P_2': 'homosapiens adult stage',
        'CordBlood_1': 'homosapiens adult stage',
        'CordBlood_2': 'homosapiens adult stage',
        'FetalAdrenalGland_2': 'fetal stage',
        'FetalAdrenalGland_3': 'fetal stage',
        'FetalAdrenalGland_4': 'fetal stage',
        'FetalBrain_3': 'fetal stage',
        'FetalBrain_4': 'fetal stage',
        'FetalBrain_5': 'fetal stage',
        'FetalBrain_6': 'fetal stage',
        'FetalCalvaria_1': 'fetal stage',
        'FetalEyes_1': 'fetal stage',
        'FetalFemaleGonad_1': 'fetal stage',
        'FetalFemaleGonad_2': 'fetal stage',
        'FetalHeart_1': 'fetal stage',
        'FetalHeart_2': 'fetal stage',
        'FetalIntestine_1': 'fetal stage',
        'FetalIntestine_2': 'fetal stage',
        'FetalIntetsine_3': 'fetal stage',
        'FetalIntestine_4': 'fetal stage',
        'FetalIntestine_5': 'fetal stage',
        'FetalKidney_3': 'fetal stage',
        'FetalKidney_4': 'fetal stage',
        'FetalKidney_5': 'fetal stage',
        'FetalKidney_6': 'fetal stage',
        'FetalLung_1': 'fetal stage',
        'FetalLung_2': 'fetal stage',
        'FetalMaleGonad_1': 'fetal stage',
        'FetalMaleGonad_2': 'fetal stage',
        'FetalMuscle_1': 'fetal stage',
        'FetalPancreas_1': 'fetal stage',
        'FetalPancreas_2': 'fetal stage',
        'FetalPancreas_3': 'fetal stage',
        'FetalRib_2': 'fetal stage',
        'FetalRib_3': 'fetal stage',
        'FetalSkin_2': 'fetal stage',
        'FetalSkin_3': 'fetal stage',
        'FetalSpinalCord_1': 'fetal stage',
        'FetalStomach_1': 'fetal stage',
        'FetalStomach_2': 'fetal stage',
        'FetalThymus_1': 'fetal stage',
        'FetalThymus_2': 'fetal stage',
        'HESC_1': 'blastula stage',
        'Liver_1': 'homosapiens adult stage',
        'Liver_2': 'homosapiens adult stage',
        'NeonatalAdrenalGland_1': 'newborn homosapiens stage',
        'PeripheralBlood_1': 'homosapiens adult stage',
        'Placenta_1': 'homosapiens adult stage',
    }
    sex_dict = {
        'Male': "male",
        'Female': "female",
        'nan': "unknown",
        'FeM=male': "unknown",
    }

    adata = anndata.read(os.path.join(data_dir, "HCL_Fig1_adata.h5ad"))
    # convert to sparse matrix
    adata.X = scipy.sparse.csr_matrix(adata.X).copy()

    # harmonise annotations
    for col in ["batch", "tissue"]:
        adata.obs[col] = adata.obs[col].astype("str")
    adata.obs.index = adata.obs.index.str.replace("AdultJeJunum", "AdultJejunum", regex=True).str.replace(
        "AdultGallBladder", "AdultGallbladder", regex=True).str.replace(
        "FetalFemaleGonald", "FetalFemaleGonad", regex=True)
    adata.obs.replace({"AdultJeJunum": "AdultJejunum", "AdultGallBladder": "AdultGallbladder",
                       "FetalFemaleGonald": "FetalFemaleGonad"}, regex=True, inplace=True)
    adata.obs.index = ["-".join(i.split("-")[:-1]) for i in adata.obs.index]

    # load celltype labels and harmonise them
    # This pandas code should work with pandas 1.2 but it does not and yields an empty data frame:
    fig1_anno = pd.read_excel(
        os.path.join(data_dir, "HCL_Fig1_cell_Info.xlsx"),
        index_col="cellnames",
        engine="xlrd",  # ToDo: Update when pandas xlsx reading with openpyxl is fixed: yields empty tables
    )
    fig1_anno.index = fig1_anno.index.str.replace("AdultJeJunum", "AdultJejunum", regex=True).str.replace(
        "AdultGallBladder", "AdultGallbladder", regex=True).str.replace(
        "FetalFemaleGonald", "FetalFemaleGonad", regex=True)

    # check that the order of cells and cell labels is the same
    assert np.all(fig1_anno.index == adata.obs.index)

    # add annotations to adata object and rename columns
    adata.obs = pd.concat([adata.obs, fig1_anno[["cluster", "stage", "donor", "celltype"]]], axis=1)
    adata.obs.columns = ["sample", "tissue", "n_genes", "n_counts", "cluster_global", "stage", "donor",
                         "celltype_global"]

    # add sample-wise annotations to the full adata object
    df = pd.DataFrame(
        columns=["Cell_barcode", "Sample", "Batch", "Cell_id", "Cluster_id", "Ages", "Development_stage", "Method",
                 "Gender", "Source", "Biomaterial", "Name", "ident", "Celltype"])
    archive = zipfile.ZipFile(os.path.join(data_dir, "annotation_rmbatch_data_revised417.zip"))
    for f in archive.namelist():
        df1 = pd.read_csv(archive.open(f), encoding="unicode_escape")
        df = pd.concat([df, df1], sort=True)
    df = df.set_index("Cell_id")
    adata = adata[[i in df.index for i in adata.obs.index]].copy()
    a_idx = adata.obs.index.copy()
    adata.obs = pd.concat([adata.obs, df[
        ["Ages", "Celltype", "Cluster_id", "Gender", "Method", "Source"]
    ]], axis=1)
    assert np.all(a_idx == adata.obs.index)

    # remove musmusculus cells from the object  # ToDo: add this back in as musmusculus data sets?
    adata = adata[adata.obs["Source"] != "MCA2.0"].copy()

    # tidy up the column names of the obs annotations
    adata.obs.columns = [
        "sample", "sub_tissue", "n_genes", "n_counts", "cluster_global", "dev_stage", "donor", "celltype_global",
        "age", "celltype_specific", "cluster_specific", "sex", "assay_sc", "source"]
    # Remove new line characters from cell type:
    adata.obs["celltype_specific"] = [
        x.replace("\n", "").rstrip()
        for x in adata.obs["celltype_specific"].values
    ]
    adata.obs["organ"] = [sample_organ_dict[x] for x in adata.obs["sample"].values]
    adata.obs["sex"] = [sex_dict[str(x)] for x in adata.obs["sex"].values]
    # TODO are the more exact developmental stages in dev_stage?
    adata.obs["dev_stage"] = [sample_dev_stage_dict[str(x)] for x in adata.obs["sample"].values]

    return adata
