import logging

import rich
import yaml
from rich.panel import Panel
from flatten_dict import flatten
from flatten_dict.reducer import make_reducer
from rich.progress import Progress, BarColumn

log = logging.getLogger(__name__)


class DataloaderValidator:

    def __init__(self, path='.'):
        self.path: str = path
        self.content: dict = {}
        self.passed: dict = {}
        self.warned: dict = {}
        self.failed: dict = {}
        self.validation_functions: list = [
            '_validate_required_attributes',
        ]

    def validate(self) -> None:
        """
        Statically verifies a yaml dataloader file against a predefined set of rules.
        Every rule is a function defined in this class, which must be part of this class' linting_functions.
        """
        with open(self.path) as yaml_file:
            self.content = yaml.load(yaml_file, Loader=yaml.FullLoader)

        progress = Progress("[bold green]{task.description}", BarColumn(bar_width=None),
                            "[bold yellow]{task.completed} of {task.total}[reset] [bold green]{task.fields[func_name]}")
        with progress:
            lint_progress = progress.add_task("Running lint checks",
                                              total=len(self.validation_functions),
                                              func_name=self.validation_functions)
            for fun_name in self.validation_functions:
                progress.update(lint_progress, advance=1, func_name=fun_name)
                getattr(self, fun_name)()

        self._print_results()

    def _validate_required_attributes(self):
        """
        Verifies that all required attributes for every dataloader are present.
        """
        passed_required_attributes = True

        attributes = ['dataset_structure:sample_fns',
                      'dataset_wise:author',
                      'dataset_wise:doi',
                      'dataset_wise:download_url_data',
                      'dataset_wise:download_url_meta',
                      'dataset_wise:normalization',
                      'dataset_wise:year',
                      'dataset_or_observation_wise:assay',
                      'dataset_or_observation_wise:organ',
                      'dataset_or_observation_wise:organism']

        flattened_dict = flatten(self.content, reducer=make_reducer(delimiter=':'))
        for attribute in attributes:
            try:
                detected = False
                for key in flattened_dict.keys():
                    if key.startswith(attribute):
                        detected = True
                if not detected:
                    passed_required_attributes = False
                    self.failed['-1'] = f'Missing attribute: {attribute}'
            except KeyError:
                passed_required_attributes = False
                self.failed['-1'] = f'Missing attribute: {attribute}'

        if passed_required_attributes:
            self.passed[0] = 'Passed required dataloader attributes checks.'

    def _print_results(self):
        console = rich.console.Console()
        console.print()
        console.rule("[bold green] LINT RESULTS")
        console.print()
        console.print(
            f'     [bold green][[\u2714]] {len(self.passed):>4} tests passed\n     [bold yellow][[!]] {len(self.warned):>4} tests had warnings\n'
            f'     [bold red][[\u2717]] {len(self.failed):>4} tests failed',
            overflow="ellipsis",
            highlight=False,
        )

        def format_result(linting_results: dict, color):
            results = []
            for line, result in linting_results.items():
                results.append(f'[bold {color}]Result: {result}')
            return "\n".join(results)

        if len(self.passed) > 0:
            console.print()
            console.rule("[bold green][[\u2714]] Tests Passed", style='green')
            console.print(Panel(format_result(self.passed, 'green'), style='green'), no_wrap=False, overflow='ellipsis')
        if len(self.warned) > 0:
            console.print()
            console.rule("[bold yellow][[!]] Test Warnings", style='yellow')
            console.print(Panel(format_result(self.warned, 'yellow'), style="yellow"), no_wrap=False, overflow='ellipsis')
        if len(self.failed) > 0:
            console.print()
            console.rule("[bold red][[\u2717]] Test Failures", style='red')
            console.print(Panel(format_result(self.failed, 'red'), style='red'), no_wrap=False, overflow='ellipsis')
