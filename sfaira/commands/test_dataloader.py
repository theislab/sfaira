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
