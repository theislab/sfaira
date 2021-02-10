import logging

import rich
from rich import print
from rich.panel import Panel

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

    def lint(self) -> None:
        """
        Statically verifies a dataloader against a predefined set of rules.
        Every rule is a function defined in this class, which must be part of this class' linting_functions.
        """
        pass

    def _lint_required_attributes(self):
        """
        Verifies that all required attributes for every dataloader are present.
        These are
        1.
        2.
        3.
        """
        pass

    def _lint_header(self):
        """

        :return:
        """

    def _print_results(validator: Validator):
        console = rich.console.Console()
        console.print()
        console.rule("[bold green] LINT RESULTS")
        console.print()
        console.print(
            f'     [bold green][[\u2714]] {len(validator.passed):>4} tests passed\n     [bold yellow][[!]] {len(validator.warnings):>4} tests had warnings\n'
            f'     [bold red][[\u2717]] {len(validator.errors):>4} tests failed',
            overflow="ellipsis",
            highlight=False,
        )

        def format_negative_result(linting_results, color):
            results = []
            for code, result in linting_results.items():
                results.append(f'[bold blue]Line: [not bold {color}]{result.lines} [bold blue]'
                               f'Error: [not bold {color}]{result.error} [bold blue]')
            return "\n".join(results)

        def format_positive_result(linting_results):
            results = []
            for result in linting_results:
                results.append(f'[bold]{result}')
            return "\n".join(results)

        if len(validator.passed) > 0:
            console.print()
            console.rule("[bold green][[\u2714]] Tests Passed", style="green")
            console.print(Panel(format_positive_result(validator.passed), style="green"), no_wrap=False, overflow="ellipsis")
        if len(validator.warnings) > 0:
            console.print()
            console.rule("[bold yellow][[!]] Test Warnings", style="yellow")
            console.print(Panel(format_negative_result(validator.warnings, 'yellow'), style="yellow"), no_wrap=False, overflow="ellipsis")
        if len(validator.errors) > 0:
            console.print()
            console.rule("[bold red][[\u2717]] Test Failures", style="red")
            console.print(Panel(format_negative_result(validator.errors, 'red'), style="red"), no_wrap=False, overflow="ellipsis")
