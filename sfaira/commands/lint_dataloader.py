import logging

from rich import print

log = logging.getLogger(__name__)


class DataloaderLinter:

    def __init__(self, path='.'):
        self.path: str = path
        self.passed: list = []
        self.warned: list = []
        self.failed: list = []
        self.linting_functions: list = [
            'lint_required_attributes',
            'run_flake8'
        ]

    def lint_dataloader(self) -> None:
        """
        Statically verifies a dataloader against a predefined set of rules.
        Every rule is a function defined in this class, which must be part of this class' linting_functions.
        """
        pass

    def lint_required_attributes(self):
        """
        Verifies that all required attributes for every dataloader are present.
        These are
        1.
        2.
        3.
        """
        pass

    def run_flake8(self):
        """
        Runs flake8 on the provided dataloader.
        Check
            passes if flake8 returns exitcode 0
            fails if flake8 returns an exitcode != 0
        """
        pass
