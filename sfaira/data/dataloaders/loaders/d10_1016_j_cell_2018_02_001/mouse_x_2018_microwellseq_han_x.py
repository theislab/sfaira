import anndata
import numpy as np
import pandas
import zipfile
import tarfile
import os

from sfaira.data import DatasetBase

SAMPLE_FNS = [
    "Bladder_dge.txt.gz",
    "BoneMarrow1_dge.txt.gz",
    # "BoneMarrow2_dge.txt.gz",  # not annotated, potentially bad quality
    # "BoneMarrow3_dge.txt.gz",  # not annotated, potentially bad quality
    "BoneMarrowcKit1_dge.txt.gz",
    "BoneMarrowcKit2_dge.txt.gz",
    "BoneMarrowcKit3_dge.txt.gz",
    "Brain1_dge.txt.gz",
    "Brain2_dge.txt.gz",
    # "CJ7.EB14.Ezh2.1_dge.txt.gz",  # TODO: sort out meta data for these
    # "CJ7.EB14.WT.1_dge.txt.gz",  # TODO: sort out meta data for these
    # "CJ7.EB14.WT.2_dge.txt.gz",  # TODO: sort out meta data for these
    # "EB.Ezh2_dge.txt.gz",  # TODO: sort out meta data for these
    # "EB.WT_dge.txt.gz",  # TODO: sort out meta data for these
    # "EmbryonicMesenchymeE14.5_dge.txt.gz",  # TODO: sort out meta data for these
    # "EmbryonicStemCell.CJ7_Deep_dge.txt.gz",  # TODO: sort out meta data for these
    # "EmbryonicStemCells_dge.txt.gz",  # TODO: sort out meta data for these
    "FetalBrain_dge.txt.gz",
    "FetalFemaleGonad_dge.txt.gz",
    "FetalIntestine_dge.txt.gz",
    "FetalKidney1_dge.txt.gz",
    "FetalKidney2_dge.txt.gz",
    "FetalLiverE14.1_dge.txt.gz",
    "FetalLung_dge.txt.gz",
    "FetalMaleGonad_dge.txt.gz",
    "FetalPancreas_dge.txt.gz",
    "FetalStomach_dge.txt.gz",
    # "human-293T_dge.txt.gz",  # ToDo: sort out meta data for these
    "Kidney1_dge.txt.gz",
    "Kidney2_dge.txt.gz",
    "Liver1_dge.txt.gz",
    "Liver2_dge.txt.gz",
    "Lung1_dge.txt.gz",
    "Lung2_dge.txt.gz",
    "Lung3_dge.txt.gz",
    # "MammaryGland.Involution.CD45.1_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Involution.CD45.2_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Involution1_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Involution2_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Lactation1_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Lactation2_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Pregnancy_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Virgin.CD45.1_dge.txt.gz",  # TODO not annotated?
    # "MammaryGland.Virgin.CD45.2_dge.txt.gz",  # TODO not annotated?
    "MammaryGland.Virgin1_dge.txt.gz",
    "MammaryGland.Virgin2_dge.txt.gz",
    "MammaryGland.Virgin3_dge.txt.gz",
    "MammaryGland.Virgin4_dge.txt.gz",
    # "mES.CJ7_dge.txt.gz",  # TODO: sort out meta data for these
    "MesenchymalStemCells_dge.txt.gz",
    "MesenchymalStemCellsPrimary_dge.txt.gz",
    # "mouse-3T3_dge.txt.gz",  # TODO: sort out meta data for these
    "Muscle_dge.txt.gz",
    "NeonatalCalvaria1_dge.txt.gz",
    "NeonatalCalvaria2_dge.txt.gz",
    "NeonatalHeart_dge.txt.gz",
    "NeonatalMuscle1_dge.txt.gz",
    "NeonatalMuscle2_dge.txt.gz",
    # "NeonatalPancreas_dge.txt.zip",  # TODO enable zip file here
    "NeonatalRib1_dge.txt.gz",
    "NeonatalRib2_dge.txt.gz",
    "NeonatalRib3_dge.txt.gz",
    "NeonatalSkin_dge.txt.gz",
    "NeontalBrain1_dge.txt.gz",
    "NeontalBrain2_dge.txt.gz",
    "Ovary1_dge.txt.gz",
    "Ovary2_dge.txt.gz",
    "Pancreas_dge.txt.gz",
    "PeripheralBlood1_dge.txt.gz",
    "PeripheralBlood2_dge.txt.gz",
    "PeripheralBlood3_dge.txt.gz",
    "PeripheralBlood4_dge.txt.gz",
    "PeripheralBlood5_dge.txt.gz",
    "PeripheralBlood6_dge.txt.gz",
    "PlacentaE14.1_dge.txt.gz",
    "PlacentaE14.2_dge.txt.gz",
    "Prostate1_dge.txt.gz",
    "Prostate2_dge.txt.gz",
    "SmallIntestine.CD45_dge.txt.gz",
    "SmallIntestine1_dge.txt.gz",
    "SmallIntestine2_dge.txt.gz",
    "SmallIntestine3_dge.txt.gz",
    "Spleen_dge.txt.gz",
    "Stomach_dge.txt.gz",
    "Testis1_dge.txt.gz",
    "Testis2_dge.txt.gz",
    "Thymus1_dge.txt.gz",
    "Thymus2_dge.txt.gz",
    "TrophoblastStemCells_dge.txt.gz",
    "Uterus1_dge.txt.gz",
    "Uterus2_dge.txt.gz",
]


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        sample_organ_dict = {
            "Bladder_dge.txt.gz": "urinary bladder",
            "BoneMarrow1_dge.txt.gz": "bone marrow",
            "BoneMarrow2_dge.txt.gz": "bone marrow",
            "BoneMarrow3_dge.txt.gz": "bone marrow",
            "BoneMarrowcKit1_dge.txt.gz": "bone marrow",
            "BoneMarrowcKit2_dge.txt.gz": "bone marrow",
            "BoneMarrowcKit3_dge.txt.gz": "bone marrow",
            "Brain1_dge.txt.gz": "brain",
            "Brain2_dge.txt.gz": "brain",
            "CJ7.EB14.Ezh2.1_dge.txt.gz": None,
            "CJ7.EB14.WT.1_dge.txt.gz": None,
            "CJ7.EB14.WT.2_dge.txt.gz": None,
            "EB.Ezh2_dge.txt.gz": None,
            "EB.WT_dge.txt.gz": None,
            "EmbryonicMesenchymeE14.5_dge.txt.gz": "mesenchyme",
            "EmbryonicStemCell.CJ7_Deep_dge.txt.gz": "blastocyst",
            "EmbryonicStemCells_dge.txt.gz": "blastocyst",
            "FetalBrain_dge.txt.gz": "brain",
            "FetalFemaleGonad_dge.txt.gz": "ovary",
            "FetalIntestine_dge.txt.gz": "intestine",
            "FetalKidney1_dge.txt.gz": "kidney",
            "FetalKidney2_dge.txt.gz": "kidney",
            "FetalLiverE14.1_dge.txt.gz": "liver",
            "FetalLung_dge.txt.gz": "lung",
            "FetalMaleGonad_dge.txt.gz": "testis",
            "FetalPancreas_dge.txt.gz": "pancreas",
            "FetalStomach_dge.txt.gz": "stomach",
            "human-293T_dge.txt.gz": None,
            "Kidney1_dge.txt.gz": "kidney",
            "Kidney2_dge.txt.gz": "kidney",
            "Liver1_dge.txt.gz": "liver",
            "Liver2_dge.txt.gz": "liver",
            "Lung1_dge.txt.gz": "lung",
            "Lung2_dge.txt.gz": "lung",
            "Lung3_dge.txt.gz": "lung",
            "MammaryGland.Involution.CD45.1_dge.txt.gz": "mammary gland",
            "MammaryGland.Involution.CD45.2_dge.txt.gz": "mammary gland",
            "MammaryGland.Involution1_dge.txt.gz": "mammary gland",
            "MammaryGland.Involution2_dge.txt.gz": "mammary gland",
            "MammaryGland.Lactation1_dge.txt.gz": "mammary gland",
            "MammaryGland.Lactation2_dge.txt.gz": "mammary gland",
            "MammaryGland.Pregnancy_dge.txt.gz": "mammary gland",
            "MammaryGland.Virgin.CD45.1_dge.txt.gz": "mammary gland",
            "MammaryGland.Virgin.CD45.2_dge.txt.gz": "mammary gland",
            "MammaryGland.Virgin1_dge.txt.gz": "mammary gland",
            "MammaryGland.Virgin2_dge.txt.gz": "mammary gland",
            "MammaryGland.Virgin3_dge.txt.gz": "mammary gland",
            "MammaryGland.Virgin4_dge.txt.gz": "mammary gland",
            "mES.CJ7_dge.txt.gz": "blastocyst",
            "MesenchymalStemCells_dge.txt.gz": "mesenchyme",
            "MesenchymalStemCellsPrimary_dge.txt.gz": "mesenchyme",
            "mouse-3T3_dge.txt.gz": None,
            "Muscle_dge.txt.gz": "skeletal muscle organ",
            "NeonatalCalvaria1_dge.txt.gz": "vault of skull",
            "NeonatalCalvaria2_dge.txt.gz": "vault of skull",
            "NeonatalHeart_dge.txt.gz": "heart",
            "NeonatalMuscle1_dge.txt.gz": "skeletal muscle organ",
            "NeonatalMuscle2_dge.txt.gz": "skeletal muscle organ",
            "NeonatalPancreas_dge.txt.zip": "pancreas",
            "NeonatalRib1_dge.txt.gz": "rib",
            "NeonatalRib2_dge.txt.gz": "rib",
            "NeonatalRib3_dge.txt.gz": "rib",
            "NeonatalSkin_dge.txt.gz": "skin of body",
            "NeontalBrain1_dge.txt.gz": "brain",
            "NeontalBrain2_dge.txt.gz": "brain",
            "Ovary1_dge.txt.gz": "ovary",
            "Ovary2_dge.txt.gz": "ovary",
            "Pancreas_dge.txt.gz": "pancreas",
            "PeripheralBlood1_dge.txt.gz": "blood",
            "PeripheralBlood2_dge.txt.gz": "blood",
            "PeripheralBlood3_dge.txt.gz": "blood",
            "PeripheralBlood4_dge.txt.gz": "blood",
            "PeripheralBlood5_dge.txt.gz": "blood",
            "PeripheralBlood6_dge.txt.gz": "blood",
            "PlacentaE14.1_dge.txt.gz": "placenta",
            "PlacentaE14.2_dge.txt.gz": "placenta",
            "Prostate1_dge.txt.gz": "prostate gland",
            "Prostate2_dge.txt.gz": "prostate gland",
            "SmallIntestine.CD45_dge.txt.gz": "small intestine",
            "SmallIntestine1_dge.txt.gz": "small intestine",
            "SmallIntestine2_dge.txt.gz": "small intestine",
            "SmallIntestine3_dge.txt.gz": "small intestine",
            "Spleen_dge.txt.gz": "spleen",
            "Stomach_dge.txt.gz": "stomach",
            "Testis1_dge.txt.gz": "testis",
            "Testis2_dge.txt.gz": "testis",
            "Thymus1_dge.txt.gz": "testis",
            "Thymus2_dge.txt.gz": "testis",
            "TrophoblastStemCells_dge.txt.gz": "trophoblast",
            "Uterus1_dge.txt.gz": "uterus",
            "Uterus2_dge.txt.gz": "uterus",
        }
        sample_dev_stage_dict = {
            "Bladder_dge.txt.gz": "adult",
            "BoneMarrow1_dge.txt.gz": "adult",
            "BoneMarrow2_dge.txt.gz": "adult",
            "BoneMarrow3_dge.txt.gz": "adult",
            "BoneMarrowcKit1_dge.txt.gz": "adult",
            "BoneMarrowcKit2_dge.txt.gz": "adult",
            "BoneMarrowcKit3_dge.txt.gz": "adult",
            "Brain1_dge.txt.gz": "adult",
            "Brain2_dge.txt.gz": "adult",
            "CJ7.EB14.Ezh2.1_dge.txt.gz": None,
            "CJ7.EB14.WT.1_dge.txt.gz": None,
            "CJ7.EB14.WT.2_dge.txt.gz": None,
            "EB.Ezh2_dge.txt.gz": None,
            "EB.WT_dge.txt.gz": None,
            "EmbryonicMesenchymeE14.5_dge.txt.gz": "embryonic",
            "EmbryonicStemCell.CJ7_Deep_dge.txt.gz": "embryonic",
            "EmbryonicStemCells_dge.txt.gz": "embryonic",
            "FetalBrain_dge.txt.gz": "fetal",
            "FetalFemaleGonad_dge.txt.gz": "fetal",
            "FetalIntestine_dge.txt.gz": "fetal",
            "FetalKidney1_dge.txt.gz": "fetal",
            "FetalKidney2_dge.txt.gz": "fetal",
            "FetalLiverE14.1_dge.txt.gz": "fetal",
            "FetalLung_dge.txt.gz": "fetal",
            "FetalMaleGonad_dge.txt.gz": "fetal",
            "FetalPancreas_dge.txt.gz": "fetal",
            "FetalStomach_dge.txt.gz": "fetal",
            "human-293T_dge.txt.gz": None,
            "Kidney1_dge.txt.gz": "adult",
            "Kidney2_dge.txt.gz": "adult",
            "Liver1_dge.txt.gz": "adult",
            "Liver2_dge.txt.gz": "adult",
            "Lung1_dge.txt.gz": "adult",
            "Lung2_dge.txt.gz": "adult",
            "Lung3_dge.txt.gz": "adult",
            "MammaryGland.Involution.CD45.1_dge.txt.gz": "adult",
            "MammaryGland.Involution.CD45.2_dge.txt.gz": "adult",
            "MammaryGland.Involution1_dge.txt.gz": "adult",
            "MammaryGland.Involution2_dge.txt.gz": "adult",
            "MammaryGland.Lactation1_dge.txt.gz": "adult",
            "MammaryGland.Lactation2_dge.txt.gz": "adult",
            "MammaryGland.Pregnancy_dge.txt.gz": "adult",
            "MammaryGland.Virgin.CD45.1_dge.txt.gz": "adult",
            "MammaryGland.Virgin.CD45.2_dge.txt.gz": "adult",
            "MammaryGland.Virgin1_dge.txt.gz": "adult",
            "MammaryGland.Virgin2_dge.txt.gz": "adult",
            "MammaryGland.Virgin3_dge.txt.gz": "adult",
            "MammaryGland.Virgin4_dge.txt.gz": "adult",
            "mES.CJ7_dge.txt.gz": "embryonic",
            "MesenchymalStemCells_dge.txt.gz": "embryonic",
            "MesenchymalStemCellsPrimary_dge.txt.gz": "embryonic",
            "mouse-3T3_dge.txt.gz": None,
            "Muscle_dge.txt.gz": "adult",
            "NeonatalCalvaria1_dge.txt.gz": "neonatal",
            "NeonatalCalvaria2_dge.txt.gz": "neonatal",
            "NeonatalHeart_dge.txt.gz": "neonatal",
            "NeonatalMuscle1_dge.txt.gz": "neonatal",
            "NeonatalMuscle2_dge.txt.gz": "neonatal",
            "NeonatalPancreas_dge.txt.zip": "neonatal",
            "NeonatalRib1_dge.txt.gz": "neonatal",
            "NeonatalRib2_dge.txt.gz": "neonatal",
            "NeonatalRib3_dge.txt.gz": "neonatal",
            "NeonatalSkin_dge.txt.gz": "neonatal",
            "NeontalBrain1_dge.txt.gz": "neonatal",
            "NeontalBrain2_dge.txt.gz": "neonatal",
            "Ovary1_dge.txt.gz": "adult",
            "Ovary2_dge.txt.gz": "adult",
            "Pancreas_dge.txt.gz": "adult",
            "PeripheralBlood1_dge.txt.gz": "adult",
            "PeripheralBlood2_dge.txt.gz": "adult",
            "PeripheralBlood3_dge.txt.gz": "adult",
            "PeripheralBlood4_dge.txt.gz": "adult",
            "PeripheralBlood5_dge.txt.gz": "adult",
            "PeripheralBlood6_dge.txt.gz": "adult",
            "PlacentaE14.1_dge.txt.gz": "fetal",
            "PlacentaE14.2_dge.txt.gz": "fetal",
            "Prostate1_dge.txt.gz": "adult",
            "Prostate2_dge.txt.gz": "adult",
            "SmallIntestine.CD45_dge.txt.gz": "adult",
            "SmallIntestine1_dge.txt.gz": "adult",
            "SmallIntestine2_dge.txt.gz": "adult",
            "SmallIntestine3_dge.txt.gz": "adult",
            "Spleen_dge.txt.gz": "adult",
            "Stomach_dge.txt.gz": "adult",
            "Testis1_dge.txt.gz": "adult",
            "Testis2_dge.txt.gz": "adult",
            "Thymus1_dge.txt.gz": "adult",
            "Thymus2_dge.txt.gz": "adult",
            "TrophoblastStemCells_dge.txt.gz": "embryonic",
            "Uterus1_dge.txt.gz": "adult",
            "Uterus2_dge.txt.gz": "adult",
        }

        self.organ = sample_organ_dict[self.sample_fn]

        self.download_url_data = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.download_url_meta = None

        self.author = "Han"
        self.dev_stage = sample_dev_stage_dict[self.sample_fn]
        self.disease = "healthy"
        self.doi_journal = "10.1016/j.cell.2018.02.001"
        self.normalization = "raw"
        self.organism = "mouse"
        self.assay_sc = "microwell-seq"
        self.year = 2018
        self.sample_source = "primary_tissue"

        self.gene_id_symbols_var_key = "index"

        # Only adult and neonatal samples are annotated:
        self.cell_types_obs_key = "Annotation" \
            if sample_dev_stage_dict[self.sample_fn] in ["adult", "neonatal"] and \
            self.sample_fn not in [
                "NeontalBrain1_dge.txt.gz",
                "NeontalBrain2_dge.txt.gz",
                "SmallIntestine.CD45_dge.txt.gz",
                "Thymus2_dge.txt.gz",
            ] else None

        self.set_dataset_id(idx=1)


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, '5435866.zip')
    with zipfile.ZipFile(fn) as archive:
        celltypes = pandas.read_csv(archive.open('MCA_CellAssignments.csv'), index_col=1)
        celltypes = celltypes.drop(["Unnamed: 0"], axis=1)

        with tarfile.open(fileobj=archive.open('MCA_500more_dge.tar.gz')) as tar:
            data = pandas.read_csv(tar.extractfile(f'500more_dge/{sample_fn}'),
                                   compression="gzip",
                                   sep=" ",
                                   header=0
                                   )

    adata = anndata.AnnData(data.T)
    annotated_cells = np.array([x in celltypes.index for x in adata.obs_names])
    # Add annotation if available for this data set:
    if np.sum(annotated_cells) > 0:
        # Subset to annotated cells if any are annotated:
        adata = adata[annotated_cells].copy()
        # Clean nans in data frame to avoid issues with cell type annotation:
        celltypes["Annotation"] = [
            x if x not in [np.nan, "nan"] else "unknown"
            for x in celltypes["Annotation"].values
        ]
        adata.obs = celltypes.loc[adata.obs_names, :]

    return adata
