import logging
import os
from rich import print
import shutil

from sfaira.commands.utils import get_pydoc
from sfaira.consts.utils import clean_doi
from sfaira.data import DatasetGroupDirectoryOriented

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
        file_path, _ = get_pydoc(path_loader=self.path_loader, doi_sfaira_repr=self.doi_sfaira_repr)
        cache_path = None
        # Clear dataset cache
        shutil.rmtree(cache_path, ignore_errors=True)

        ds = DatasetGroupDirectoryOriented(
            file_base=file_path,
            data_path=self.path_data,
            meta_path=None,
            cache_path=None
        )

        return ds, cache_path

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
