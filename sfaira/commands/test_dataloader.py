import logging
import os
from subprocess import Popen

from rich import print
from sfaira.commands.questionary import sfaira_questionary

log = logging.getLogger(__name__)


class DataloaderTester:

    def __init__(self, path):
        self.WD = os.path.dirname(__file__)
        self.path = path
        self.doi = ''
        self.doi_sfaira_repr = ''

    def test_dataloader(self):
        """
        Runs a predefined unit test on a given dataloader.
        """
        print('[bold red]This command is currently disabled.')
        # print('[bold blue]Please ensure that your dataloader is in sfaira/dataloaders/loaders/<doi_flattened>.')
        # print('[bold blue]Please ensure that your test data is in sfaira/unit_tests/template_data/<doi_flattened>.')
        # self._prompt_doi()
        # self._run_unittest()

    def _prompt_doi(self):
        self.doi = sfaira_questionary(function='text',
                                      question='Enter your DOI',
                                      default='10.1000/j.journal.2021.01.001')
        self.doi_sfaira_repr = f'd{self.doi.translate({ord(c): "_" for c in r"!@#$%^&*()[]/{};:,.<>?|`~-=_+"})}'

    def _run_unittest(self):
        print('[bold blue]Conflicts are not automatically resolved.')
        print('[bold blue]Please go back to [bold]https://www.ebi.ac.uk/ols/ontologies/cl[blue] for every mismatch or conflicts '
              'and add the correct cell ontology class name into the .csv "target" column.')
        pytest = Popen(['pytest', '-s', self.path, '--doi_sfaira_repr', self.doi_sfaira_repr],
                       universal_newlines=True, shell=False, close_fds=True)
        (pytest_stdout, pytest_stderr) = pytest.communicate()
        print(pytest_stderr)
