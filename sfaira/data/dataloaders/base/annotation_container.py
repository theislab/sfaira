import os
from typing import List, Tuple, Union

from sfaira.data.dataloaders.utils import read_yaml
from sfaira.data.dataloaders.obs_utils import get_ontology, value_protection


class AnnotationContainer:

    _dataset_index: Union[None, int] = None

    layer_counts: Union[None, int] = None
    layer_processed: Union[None, str] = None
    layer_spliced_counts: Union[None, str] = None
    layer_spliced_processed: Union[None, str] = None
    layer_unspliced_counts: Union[None, str] = None
    layer_unspliced_processed: Union[None, str] = None
    layer_velocity: Union[None, str] = None

    _assay_sc: Union[None, str] = None
    _assay_differentiation: Union[None, str] = None
    _assay_type_differentiation: Union[None, str] = None
    author: Union[None, str] = None
    _bio_sample: Union[None, str] = None
    _cell_line: Union[None, str] = None
    _cell_type: Union[None, str] = None
    _default_embedding: Union[None, str] = None
    _development_stage: Union[None, str] = None
    _disease: Union[None, str] = None
    doi_journal: Union[None, str] = None
    doi_preprint: Union[None, str] = None
    _download_url_data: Union[Tuple[List[None]], Tuple[List[str]], None] = None
    _download_url_meta: Union[Tuple[List[None]], Tuple[List[str]], None] = None
    _ethnicity: Union[None, str] = None
    gm: Union[None, str] = None
    individual: Union[None, str] = None
    _organ: Union[None, str] = None
    _organism: Union[None, str] = None
    _primary_data: Union[None, bool] = None
    _sex: Union[None, str] = None
    source: Union[None, str] = None
    source_doi: Union[None, str] = None
    _sample_source: Union[None, str] = None
    state_exact: Union[None, str] = None
    _tech_sample: Union[None, str] = None
    treatment: Union[None, str] = None
    _year: Union[None, int] = None

    _bio_sample_obs_key: Union[None, str] = None
    _tech_sample_obs_key: Union[None, str] = None
    assay_sc_obs_key: Union[None, str] = None
    assay_differentiation_obs_key: Union[None, str] = None
    assay_type_differentiation_obs_key: Union[None, str] = None
    cell_line_obs_key: Union[None, str] = None
    cell_type_obs_key: Union[None, str] = None
    development_stage_obs_key: Union[None, str] = None
    disease_obs_key: Union[None, str] = None
    ethnicity_obs_key: Union[None, str] = None
    gm_obs_key: Union[None, str] = None
    individual_obs_key: Union[None, str] = None
    organ_obs_key: Union[None, str] = None
    organism_obs_key: Union[None, str] = None
    sample_source_obs_key: Union[None, str] = None
    sex_obs_key: Union[None, str] = None
    source_doi_obs_key: Union[None, str] = None
    state_exact_obs_key: Union[None, str] = None
    treatment_obs_key: Union[None, str] = None

    _feature_reference: Union[None, str] = None
    _feature_type: Union[None, str] = None

    feature_id_var_key: Union[None, str] = None
    feature_reference_var_key: Union[None, str] = None
    feature_symbol_var_key: Union[None, str] = None
    feature_type_var_key: Union[None, str] = None

    spatial_x_coord_obs_key: Union[None, str] = None
    spatial_y_coord_obs_key: Union[None, str] = None
    spatial_z_coord_obs_key: Union[None, str] = None
    vdj_vj_1_obs_key_prefix: Union[None, str] = None
    vdj_vj_2_obs_key_prefix: Union[None, str] = None
    vdj_vdj_1_obs_key_prefix: Union[None, str] = None
    vdj_vdj_2_obs_key_prefix: Union[None, str] = None
    vdj_c_call_obs_key_suffix: Union[None, str] = None
    vdj_consensus_count_obs_key_suffix: Union[None, str] = None
    vdj_d_call_obs_key_suffix: Union[None, str] = None
    vdj_duplicate_count_obs_key_suffix: Union[None, str] = None
    vdj_j_call_obs_key_suffix: Union[None, str] = None
    vdj_junction_obs_key_suffix: Union[None, str] = None
    vdj_junction_aa_obs_key_suffix: Union[None, str] = None
    vdj_locus_obs_key_suffix: Union[None, str] = None
    vdj_productive_obs_key_suffix: Union[None, str] = None
    vdj_v_call_obs_key_suffix: Union[None, str] = None

    def read_from_yaml(self, yaml_path: Union[str, None], sample_fn: Union[None, str] = None):
        # Check if YAML files exists, read meta data from there if available:
        if yaml_path is not None:
            assert os.path.exists(yaml_path), f"did not find yaml {yaml_path}"
            yaml_vals = read_yaml(fn=yaml_path)
            # Set organism first as this is required to disambiguate valid entries for other meta data.
            self.organism = yaml_vals["attr"]["organism"]
            for k, v in yaml_vals["attr"].items():
                if v is not None and k not in ["organism", "sample_fns"]:
                    if isinstance(v, dict):  # v is a dictionary over file-wise meta-data items
                        if sample_fn in v.keys():
                            # only set value if field exists
                            v = v[sample_fn]
                    # Catches spelling errors in meta data definition (yaml keys).
                    if not hasattr(self, k) and not hasattr(self, "_" + k):
                        raise ValueError(f"Tried setting unavailable property {k}.")
                    try:
                        setattr(self, k, v)
                    except AttributeError as e:
                        raise ValueError(f"ValueError when setting {k} as {v}: {e}")

    @property
    def assay_sc(self) -> Union[None, str]:
        return self._assay_sc

    @assay_sc.setter
    def assay_sc(self, x: str):
        x = value_protection(attr="assay_sc", allowed=self.get_ontology(k="assay_sc"), attempted=x)
        self._assay_sc = x

    @property
    def assay_differentiation(self) -> Union[None, str]:
        return self._assay_differentiation

    @assay_differentiation.setter
    def assay_differentiation(self, x: str):
        x = value_protection(attr="assay_differentiation", allowed=self.get_ontology(k="assay_differentiation"),
                             attempted=x)
        self._assay_differentiation = x

    @property
    def assay_type_differentiation(self) -> Union[None, str]:
        return self._assay_type_differentiation

    @assay_type_differentiation.setter
    def assay_type_differentiation(self, x: str):
        x = value_protection(attr="assay_type_differentiation",
                             allowed=self.get_ontology(k="assay_type_differentiation"), attempted=x)
        self._assay_type_differentiation = x

    @property
    def bio_sample(self) -> Union[None, str]:
        if self._bio_sample is not None:
            return self._bio_sample
        else:
            # Define biological sample automatically.
            bio_sample = "*".join([x for x in [
                self.assay_sc,
                self.assay_differentiation,
                self.assay_type_differentiation,
                self.development_stage,
                self.disease,
                self.ethnicity,
                self.individual,
                self.organ,
                self.organism,
                self.sex,
            ] if x is not None])
            return bio_sample if bio_sample != "" else None

    @bio_sample.setter
    def bio_sample(self, x: str):
        self._bio_sample = x

    @property
    def bio_sample_obs_key(self) -> Union[None, str]:
        if self._bio_sample_obs_key is not None:
            return self._bio_sample_obs_key
        else:
            # Define biological sample automatically.
            bio_sample_obs_key = "*".join([x for x in [
                self.assay_sc_obs_key,
                self.assay_differentiation_obs_key,
                self.assay_type_differentiation_obs_key,
                self.development_stage_obs_key,
                self.disease_obs_key,
                self.ethnicity_obs_key,
                self.individual_obs_key,
                self.organ_obs_key,
                self.organism_obs_key,
                self.sex_obs_key,
                self.state_exact_obs_key,
            ] if x is not None])
            return bio_sample_obs_key if bio_sample_obs_key != "" else None

    @bio_sample_obs_key.setter
    def bio_sample_obs_key(self, x: str):
        self._bio_sample_obs_key = x

    @property
    def cell_line(self) -> Union[None, str]:
        return self._cell_line

    @cell_line.setter
    def cell_line(self, x: str):
        # TODO add value protection here.
        self._cell_line = x

    @property
    def cell_type(self) -> Union[None, str]:
        return self._cell_type

    @cell_type.setter
    def cell_type(self, x: str):
        x = value_protection(attr="cell_type", allowed=self.get_ontology(k="cell_type"), attempted=x)
        self._cell_type = x

    @property
    def dataset_index(self) -> Union[None, int]:
        x = self._dataset_index
        # Catch if index was omitted:
        if x is None:
            x = 1
        return x

    @dataset_index.setter
    def dataset_index(self, x: int):
        self._dataset_index = x

    @property
    def default_embedding(self) -> Union[None, str]:
        return self._default_embedding

    @default_embedding.setter
    def default_embedding(self, x: str):
        x = value_protection(attr="default_embedding", allowed=self.get_ontology(k="default_embedding"), attempted=x)
        self._default_embedding = x

    @property
    def development_stage(self) -> Union[None, str]:
        return self._development_stage

    @development_stage.setter
    def development_stage(self, x: str):
        x = value_protection(attr="development_stage", allowed=self.get_ontology(k="development_stage"), attempted=x)
        self._development_stage = x

    @property
    def disease(self) -> Union[None, str]:
        return self._disease

    @disease.setter
    def disease(self, x: str):
        x = value_protection(attr="disease", allowed=self.get_ontology(k="disease"), attempted=x)
        self._disease = x

    @property
    def download_url_data(self) -> Union[Tuple[List[str]], Tuple[List[None]]]:
        """
        Data download website(s).

        Save as tuple with single element, which is a list of all download websites relevant to dataset.
        :return:
        """
        x = self._download_url_data
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        return x

    @download_url_data.setter
    def download_url_data(self, x: Union[str, None, List[str], Tuple[str], List[None], Tuple[None]]):
        # Formats to tuple with single element, which is a list of all download websites relevant to dataset,
        # which can be used as a single element column in a pandas data frame.
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        self._download_url_data = x

    @property
    def download_url_meta(self) -> Union[Tuple[List[str]], Tuple[List[None]]]:
        """
        Meta data download website(s).

        Save as tuple with single element, which is a list of all download websites relevant to dataset.
        :return:
        """
        x = self._download_url_meta
        # if self._download_url_meta is not None:  # TODO add this back in once download_meta is set in all datasets
        #    x = self._download_url_meta
        # else:
        #    if self.meta is None:
        #        self.load_meta(fn=None)
        #    x = self.meta[self._adata_ids.download_url_meta]
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        return x

    @download_url_meta.setter
    def download_url_meta(self, x: Union[str, None, List[str], Tuple[str], List[None], Tuple[None]]):
        # Formats to tuple with single element, which is a list of all download websites relevant to dataset,
        # which can be used as a single element column in a pandas data frame.
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        self._download_url_meta = x

    @property
    def ethnicity(self) -> Union[None, str]:
        return self._ethnicity

    @ethnicity.setter
    def ethnicity(self, x: str):
        x = value_protection(attr="ethnicity", allowed=self.get_ontology(k="ethnicity"), attempted=x)
        self._ethnicity = x

    @property
    def feature_reference(self) -> str:
        return self._feature_reference

    @feature_reference.setter
    def feature_reference(self, x: str):
        # TODO add value protection here.
        self._feature_reference = x

    @property
    def feature_type(self) -> str:
        return self._feature_type

    @feature_type.setter
    def feature_type(self, x: str):
        x = value_protection(attr="feature_type", allowed=self.get_ontology(k="feature_type"), attempted=x)
        self._feature_type = x

    @property
    def primary_data(self) -> Union[None, bool]:
        return self._primary_data

    @primary_data.setter
    def primary_data(self, x: bool):
        x = value_protection(attr="primary_data", allowed=self.get_ontology(k="primary_data"), attempted=x)
        self._primary_data = x

    @property
    def organ(self) -> Union[None, str]:
        return self._organ

    @organ.setter
    def organ(self, x: str):
        x = value_protection(attr="organ", allowed=self.get_ontology(k="organ"), attempted=x)
        self._organ = x

    @property
    def organism(self) -> Union[None, str]:
        return self._organism

    @organism.setter
    def organism(self, x: str):
        x = value_protection(attr="organism", allowed=self.get_ontology(k="organism"), attempted=x)
        self._organism = x

    @property
    def sample_source(self) -> Union[None, str]:
        return self._sample_source

    @sample_source.setter
    def sample_source(self, x: str):
        x = value_protection(attr="sample_source", allowed=self.get_ontology(k="sample_source"), attempted=x)
        self._sample_source = x

    @property
    def sex(self) -> Union[None, str]:
        return self._sex

    @sex.setter
    def sex(self, x: str):
        x = value_protection(attr="sex", allowed=self.get_ontology(k="sex"), attempted=x)
        self._sex = x

    @property
    def tech_sample(self) -> Union[None, str]:
        if self._tech_sample is not None:
            return self._tech_sample
        else:
            # Define technical batch automatically as biological sample.
            return self.bio_sample

    @tech_sample.setter
    def tech_sample(self, x: str):
        self._tech_sample = x

    @property
    def tech_sample_obs_key(self) -> Union[None, str]:
        if self._tech_sample_obs_key is not None:
            return self._tech_sample_obs_key
        else:
            # Define technical batch automatically as biological sample.
            return self.bio_sample_obs_key

    @tech_sample_obs_key.setter
    def tech_sample_obs_key(self, x: str):
        self._tech_sample_obs_key = x

    @property
    def year(self) -> Union[None, int]:
        return self._year

    @year.setter
    def year(self, x: int):
        x = value_protection(attr="year", allowed=self.get_ontology(k="year"), attempted=x)
        self._year = x

    def get_ontology(self, **kwargs):
        return get_ontology(organism=self.organism, **kwargs)
