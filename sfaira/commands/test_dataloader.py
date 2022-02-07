import logging
import os
import shutil
import pydoc

from rich import print
from sfaira.commands.questionary import sfaira_questionary
from sfaira.consts.utils import clean_doi
from sfaira.data import DatasetGroupDirectoryOriented

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None

log = logging.getLogger(__name__)


class DataloaderTester:

    def __init__(self, path, test_data, doi):
        self.WD = os.path.dirname(__file__)
        self.path = path
        self.test_data = test_data
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
        dir_loader_sfaira = "sfaira.data.dataloaders.loaders."
        file_path_sfaira = os.path.dirname(str(pydoc.locate(dir_loader_sfaira + "FILE_PATH")))

        dir_loader_sfairae = "sfaira_extension.data.dataloaders.loaders." if sfairae else None
        file_path_sfairae = os.path.dirname(str(pydoc.locate(dir_loader_sfairae + "FILE_PATH"))) if sfairae else None

        # Check if loader name is a directory either in sfaira or sfaira_extension loader collections:
        if self.doi_sfaira_repr in os.listdir(file_path_sfaira):
            dir_loader = dir_loader_sfaira + "." + self.doi_sfaira_repr
        elif file_path_sfairae and self.doi_sfaira_repr in os.listdir(file_path_sfairae):
            dir_loader = dir_loader_sfairae + "." + self.doi_sfaira_repr
        else:
            raise ValueError("data loader not found in sfaira and also not in sfaira_extension")
        file_path = str(pydoc.locate(dir_loader + ".FILE_PATH"))
        cache_path = None
        # Clear dataset cache
        shutil.rmtree(cache_path, ignore_errors=True)

        ds = DatasetGroupDirectoryOriented(
            file_base=file_path,
            data_path=self.test_data,
            meta_path=None,
            cache_path=None
        )

        return ds, cache_path

    def _test_dataloader(self, clean_tsvs: bool, in_phase_3: bool):
        """
        Tests the dataloader.
        """
        print('[bold blue]Conflicts are not automatically resolved.')
        print('[bold blue]Please go back to [bold]https://www.ebi.ac.uk/ols/ontologies/cl[blue] for every mismatch or '
              'conflicts and add the correct cell ontology class name into the .tsv "target" column.')

        ds, cache_path = self._get_ds()
        if clean_tsvs:
            ds.clean_ontology_class_maps()

        # TODO try-except with good error description saying that the data loader is broken here:
        ds.load(load_raw=True, allow_caching=True)
        # Try loading from cache:
        ds, cache_path = self._get_ds()
        # TODO try-except with good error description saying that the data loader is broken here:
        ds.load(load_raw=False, allow_caching=True)
        shutil.rmtree(cache_path, ignore_errors=True)
        print("[bold blue]Completed testing of data loader, the data loader is now ready for use.")
        if in_phase_3:
            print('[bold orange]Sfaira butler: "You data loader is finished!"')
            print('[bold orange]"Proceed to phase 4 (publish) or use data loader."')
        else:
            print('[bold orange]Sfaira butler: "You data loader works!"')
