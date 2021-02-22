import os
import anndata as ad
from typing import Union

from sfaira.data import DatasetBaseGroupLoadingOneFile


# SFAIRA TODO: Add correct sample IDs here.
SAMPLE_IDS = [
    # "your_sample_id_1",
    # "your_sample_id_2"
]


class Dataset(DatasetBaseGroupLoadingOneFile):

    def __init__(
            self,
            sample_fn: str,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(sample_fn=sample_fn, sample_ids=SAMPLE_IDS, meta_path=meta_path, cache_path=cache_path,
                         path=path, **kwargs)

        # SFAIRA TODO: Add you meta data here.
        self.id = '{{ cookiecutter.id }}'  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).

        self.author = {{ cookiecutter.author }}  # author (list) who sampled / created the data set
        self.doi = '{{ cookiecutter.doi }}'  # doi of data set accompanying manuscript

        self.download_url_data = '{{ cookiecutter.download_url_data }}'  # download website(s) of data files
        # self.download_url_meta = 'x'  # download website(s) of meta data files

        self.organ = '{{ cookiecutter.organ }}'  # organ (anatomical structure)
        self.organism = '{{ cookiecutter.organism }}'  # (*) species / organism
        self.protocol = '{{ cookiecutter.protocol }}'  # (*, optional) protocol used to sample data (e.g. smart-seq2)
        self.year = {{ cookiecutter.year }}  # year in which sample was acquired
        # self.age = 'x'  # (*, optional) age of sample
        # self.dev_stage = x  # (*, optional) developmental stage of organism
        # self.ethnicity = x  # (*, optional) ethnicity of sample
        # self.healthy = x  # (*, optional) whether sample represents a healthy organism
        # self.normalisation = x  # (optional) normalisation applied to raw data loaded (ideally counts, "raw")
        # self.sex = x  # (*, optional) sex
        # self.state_exact = x  # (*, optional) exact disease, treatment or perturbation state of sample

        # SFAIRA: The following meta data may instead also be supplied on a cell level if an appropriate column
        # SFAIRA: is present in the anndata instance (specifically in .obs) after loading. You need to make sure this is loaded in the loading script)!
        # SFAIRA: See above for a description what these meta data attributes mean. If these attributes are note available, you can simply leave this out.
        # self.obs_key_age = x  # (optional, see above, do not provide if .age is provided)
        # self.obs_key_dev_stage = x  # (optional, see above, do not provide if .dev_stage is provided)
        # self.obs_key_ethnicity = x  # (optional, see above, do not provide if .ethnicity is provided)
        # self.obs_key_healthy = x  # (optional, see above, do not provide if .healthy is provided)
        # self.obs_key_organ = x  # (optional, see above, do not provide if .organ is provided)
        # self.obs_key_organism = x  # (optional, see above, do not provide if .organism is provided)
        # self.obs_key_protocol = x  # (optional, see above, do not provide if .protocol is provided)
        # self.obs_key_sex = x  # (optional, see above, do not provide if .sex is provided)
        # self.obs_key_state_exact = x  # (optional, see above, do not provide if .state_exact is provided)
        # SFAIRA: Additionally, cell type annotation is ALWAYS provided per cell in .obs, this annotation is optional though.
        # SFAIRA: name of column which contain streamlined cell ontology cell type classes:
        # self.obs_key_cellontology_original = x  # (optional)

        # SFAIRA TODO: Make sure to include this attribute which indicates the column in self.adata in which you saved the sample IDs.
        self.obs_key_sample = 'x'

    def _load_full(self) -> ad.AnnData:
        # fn = os.path.join(self.data_dir,)  # SFAIRA TODO: add the name of the raw file
        # SFAIRA TODO: load full data into AnnData object (no subsetting!) and return it
        pass
