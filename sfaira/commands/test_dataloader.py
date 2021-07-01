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

    def test_dataloader(self):
        """
        Runs a predefined unit test on a given dataloader.
        """
        if not self.doi:
            self._prompt_doi()
        self.doi_sfaira_repr = clean_doi(self.doi)
        print(f'[bold blue]Please ensure that your dataloader is in sfaira/dataloaders/loaders/{self.doi_sfaira_repr}.')
        self._test_dataloader()

    def _prompt_doi(self):
        self.doi = sfaira_questionary(function='text',
                                      question='Enter your DOI',
                                      default='10.1000/j.journal.2021.01.001')

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

    def _test_dataloader(self):
        """
        Tests the dataloader.
        """
        print('[bold blue]Conflicts are not automatically resolved.')
        print('[bold blue]Please go back to [bold]https://www.ebi.ac.uk/ols/ontologies/cl[blue] for every mismatch or '
              'conflicts and add the correct cell ontology class name into the .tsv "target" column.')

        ds, cache_path = self._get_ds()
        ds.clean_ontology_class_map()

        # TODO try-except with good error description saying that the data loader is broken here:
        ds.load(
            remove_gene_version=True,
            # match_to_reference=TODO get organism here,
            load_raw=True,
            allow_caching=True
        )
        # Try loading from cache:
        ds, cache_path = self._get_ds()
        # TODO try-except with good error description saying that the data loader is broken here:
        ds.load(
            remove_gene_version=True,
            # match_to_reference=TODO get organism here,
            load_raw=False,
            allow_caching=True
        )
        shutil.rmtree(cache_path, ignore_errors=True)
