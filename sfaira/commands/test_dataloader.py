import logging
import os
from subprocess import Popen

from rich import print
from sfaira.commands.questionary import sfaira_questionary

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
        print('[bold blue]Please ensure that your dataloader is in sfaira/dataloaders/loaders/<doi_flattened>.')
        if not self.doi:
            self._prompt_doi()
        self.doi_sfaira_repr = f'd{self.doi.translate({ord(c): "_" for c in r"!@#$%^&*()[]/{};:,.<>?|`~-=_+"})}'
        self._run_unittest()

    def _prompt_doi(self):
        self.doi = sfaira_questionary(function='text',
                                      question='Enter your DOI',
                                      default='10.1000/j.journal.2021.01.001')

    def _run_unittest(self):
        """
        Runs the actual integration test by invoking pytest on it.
        """
        print('[bold blue]Conflicts are not automatically resolved.')
        print('[bold blue]Please go back to [bold]https://www.ebi.ac.uk/ols/ontologies/cl[blue] for every mismatch or conflicts '
              'and add the correct cell ontology class name into the .tsv "target" column.')

        os.chdir(f'{self.path}/sfaira/unit_tests/data_contribution')

        pytest = Popen(['pytest', 'test_data_template.py', '--doi_sfaira_repr', self.doi_sfaira_repr, '--test_data', self.test_data],
                       universal_newlines=True, shell=False, close_fds=True)
        (pytest_stdout, pytest_stderr) = pytest.communicate()
        if pytest_stdout:
            print(pytest_stdout)
        if pytest_stderr:
            print(pytest_stderr)

        os.chdir(self.cwd)
