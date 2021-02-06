import logging

from rich import print

log = logging.getLogger(__name__)


class DataloaderTester:

    @classmethod
    def test_dataloader(cls):
        """
        Runs a predefined unit test on a given dataloader.
        """
        doi = cls._prompt_doi()
        cls._run_unittest(doi)

    @classmethod
    def _prompt_doi(cls):
        pass

    @classmethod
    def _run_unittest(cls, doi: str):
        pass


