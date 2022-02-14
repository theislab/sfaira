import logging
import numpy as np
import os

import sys
from rich import print

from sfaira.commands.utils import get_ds
from sfaira.consts.utils import clean_doi

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None

log = logging.getLogger(__name__)


class DataloaderTester:

    def __init__(self, path_loader, path_data, doi):
        self.WD = os.path.dirname(__file__)
        self.path_loader = path_loader
        self.path_data = path_data
        self.cwd = os.getcwd()
        self.doi = doi
        self.doi_sfaira_repr = ''

    def test_dataloader(self, clean_tsvs: bool, in_phase_3: bool):
        """
        Runs a predefined unit test on a given dataloader.
        """
        self.doi_sfaira_repr = clean_doi(self.doi)
        self._test_dataloader(clean_tsvs=clean_tsvs, in_phase_3=in_phase_3)

    def _get_ds(self):
        return get_ds(doi_sfaira_repr=self.doi_sfaira_repr, path_data=self.path_data, path_loader=self.path_loader)

    def _test_dataloader(self, clean_tsvs: bool, in_phase_3: bool):
        """
        Tests the dataloader.
        """
        print('[bold blue]Conflicts are not automatically resolved.')
        print('[bold blue]In case of coflicts, please go back to [bold]https://www.ebi.ac.uk/ols/ontologies/cl[blue] '
              'and add the correct cell ontology class name into the .tsv "target" column.')

        ds, cache_path = self._get_ds()
        if clean_tsvs:
            ds.clean_ontology_class_maps()

        ds, cache_path = self._get_ds()
        ds.load(load_raw=True, allow_caching=False)
        # Test that not too much of the count matrix is lost during feature streamlining:
        # This would indicate an issue with the feature assignments.
        for k, v in ds.datasets.items():
            signal_raw = np.asarray(v.adata.X.sum()).sum()
            # Check that fields are there:
            # .obs:
            for x in [
                "assay_sc_obs_key",
                "assay_differentiation_obs_key",
                "assay_type_differentiation_obs_key",
                "bio_sample_obs_key",
                "cell_line_obs_key",
                "cell_type_obs_key",
                "development_stage_obs_key",
                "disease_obs_key",
                "ethnicity_obs_key",
                "gm_obs_key",
                "individual_obs_key",
                "organ_obs_key",
                "organism_obs_key",
                "sample_source_obs_key",
                "sex_obs_key",
                "source_doi_obs_key",
                "state_exact_obs_key",
                "tech_sample_obs_key",
                "treatment_obs_key",
                "spatial_x_coord_obs_key",
                "spatial_y_coord_obs_key",
                "spatial_z_coord_obs_key",
                "vdj_vj_1_obs_key_prefix",
                "vdj_vj_2_obs_key_prefix",
                "vdj_vdj_1_obs_key_prefix",
                "vdj_vdj_2_obs_key_prefix",
                "vdj_c_call_obs_key_suffix",
                "vdj_consensus_count_obs_key_suffix",
                "vdj_d_call_obs_key_suffix",
                "vdj_duplicate_count_obs_key_suffix",
                "vdj_j_call_obs_key_suffix",
                "vdj_junction_obs_key_suffix",
                "vdj_junction_aa_obs_key_suffix",
                "vdj_locus_obs_key_suffix",
                "vdj_productive_obs_key_suffix",
                "vdj_v_call_obs_key_suffix",
            ]:
                val = getattr(v, x)
                if val is not None:
                    if x in ["bio_sample_obs_key", "individual_obs_key", "tech_sample_obs_key"]:
                        for vali in val.split("*"):
                            if vali not in v.adata.obs.columns:
                                print(f'[bold red]Did not find column {vali} for {x} in data set {k}, found: '
                                      '{v.adata.obs.columns}.')
                                sys.exit()
                    else:
                        if val not in v.adata.obs.columns:
                            print(f'[bold red]Did not find column {val} for {x} in data set {k}, found: '
                                  f'{v.adata.obs.columns}.')
                            sys.exit()
            # .var:
            for x in [
                "feature_id_var_key",
                "feature_symbol_var_key",
                "feature_reference_var_key",
                "feature_type_var_key",
            ]:
                val = getattr(v, x)
                if val is not None:
                    if val not in v.adata.var.columns and val != "index":
                        print(f'[bold red]Did not find column {val} for {x} in data set {k}, found: '
                              f'{v.adata.var.columns}.')
                        sys.exit()
            v.streamline_features(schema="cellxgene:" + "2.0.0")
            signal_proc = np.asarray(v.adata.X.sum()).sum()
            if signal_proc < 0.01 * signal_raw and v.feature_type != "peak":
                print('[bold red]Mapping your feature space to a reference annotation resulted in a heavy loss of '
                      f'counts in dataset {k}.')
                print(f'[bold red]From {signal_raw} total counts before streamlining, {signal_proc} are left after.')
                print('[bold red]Consider revising feature meta data.')
                sys.exit()
            v.streamline_metadata(schema="cellxgene")
        print("[bold blue]Completed testing of data loader, the data loader is now ready for use.")
        if in_phase_3:
            print('[bold orange]Sfaira butler: "You data loader is finished!"')
            print('[bold orange]               "Proceed to phase 4 (publish) or use data loader."')
            print('[bold orange]               "Copy the following lines as a post into the pull request:"')
            print('[bold blue]=========================')
            print('[bold blue]Data loader test passed:')
            for x in ds.datasets.keys():
                print(f'[bold blue]    - {x}')
            print('[bold blue]=========================')
        else:
            print('[bold orange]Sfaira butler: "You data loader works!"')
