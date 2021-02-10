import logging

import rich
from rich.panel import Panel
from rich.progress import Progress, BarColumn

log = logging.getLogger(__name__)


class DataloaderLinter:

    def __init__(self, path='.'):
        self.path: str = path
        self.content: list = []
        self.passed: dict = {}
        self.warned: dict = {}
        self.failed: dict = {}
        self.linting_functions: list = [
            '_lint_header',
            '_lint_dataloader_object',
            '_lint_required_attributes',
            '_lint_sfaira_todos',
            '_lint_load'
        ]

    def lint(self, path) -> None:
        """
        Statically verifies a dataloader against a predefined set of rules.
        Every rule is a function defined in this class, which must be part of this class' linting_functions.
        :param path: Path to an existing dataloader
        """
        with open(path, 'r') as f:
            self.content = list(map(lambda line: line.strip(), f.readlines()))

        progress = Progress("[bold green]{task.description}", BarColumn(bar_width=None),
                            "[bold yellow]{task.completed} of {task.total}[reset] [bold green]{task.fields[func_name]}")
        with progress:
            lint_progress = progress.add_task("Running lint checks",
                                              total=len(self.linting_functions),
                                              func_name=self.linting_functions)
            for fun_name in self.linting_functions:
                progress.update(lint_progress, advance=1, func_name=fun_name)
                getattr(self, fun_name)()

        self._print_results()

    def _lint_header(self):
        """
        Verifies the docstring header against all requirements.
        Expected are
            Author: str
            Email: valid mail
            Version: semantic ver
        """
        passed_lint_header = True

        try:
            line, author = list(filter(lambda line_author: line_author[1].startswith('Author:'), enumerate(self.content)))[0]
            line, email = list(filter(lambda line_author: line_author[1].startswith('Email:'), enumerate(self.content)))[0]
            line, version = list(filter(lambda line_author: line_author[1].startswith('Version:'), enumerate(self.content)))[0]
        except IndexError:
            passed_lint_header = False
            self.failed['-1'] = 'Docstring for Author, Email and Version is invalid.'

        if passed_lint_header:
            self.passed[line] = 'Passed header checks.'

    def _lint_dataloader_object(self):
        """
        Verifies that the Dataloader Object itself (no the attributes) is valid
        """
        # TODO Could be more strict by checking also whether the constructor is valid, but too much of a hazzle with Black formatting.
        passed_lint_dataloader_object = True

        try:
            line, dl_object = list(filter(lambda line_dl_object: line_dl_object[1].startswith(('class Dataset(DatasetBaseGroupLoadingManyFiles):',
                                                                                               'class Dataset(DatasetBase):')), enumerate(self.content)))[0]
        except IndexError:
            passed_lint_dataloader_object = False
            self.failed['-1'] = 'Missing one of class Dataset(DatasetBase) or class Dataset(DatasetBaseGroupLoadingManyFiles)'

        if passed_lint_dataloader_object:
            self.passed[line] = 'Passed dataloader object checks.'

    def _lint_load(self):
        """
        Verifies that the method _load_any_object(self, fn=None) is present.
        """
        passed_load = True

        try:
            line, dl_object = list(filter(lambda line_dl_object: line_dl_object[1].startswith(('def _load_any_object(self, fn=None):', 'def _load(self, fn):')),
                                          enumerate(self.content)))[0]
        except IndexError:
            passed_load = False
            self.failed['-1'] = 'Missing one of methods    _load_any_object(self, fn=None)  or    def _load(self, fn)'

        if passed_load:
            self.passed[line] = 'Passed dataloader object checks.'

    def _lint_required_attributes(self):
        """
        Verifies that all required attributes for every dataloader are present.
        """
        passed_required_attributes = True

        attributes = ['self.id',
                      'self.author',
                      'self.doi',
                      'self.download_url_data',
                      'self.organ',
                      'self.organism',
                      'self.protocol',
                      'self.year']

        for attribute in attributes:
            try:
                line, attribute = list(filter(lambda line_attribute: line_attribute[1].startswith(attribute), enumerate(self.content)))[0]
            except IndexError:
                passed_required_attributes = False
                self.failed['-1'] = 'One of required attributes  id, author, doi, download_url_data, organ, organism, protocol, year   is missing.'

        if passed_required_attributes:
            self.passed[0] = 'Passed required dataloader attributes checks.'

    def _lint_sfaira_todos(self):
        """
        Warns if any SFAIRA TODO: statements were found
        """
        passed_sfaira_todos = True

        for index, line in enumerate(self.content):
            if 'SFAIRA TODO' in line:
                passed_sfaira_todos = False
                self.warned[f'{index}'] = f'Line {index}: {line[2:]}'

        if passed_sfaira_todos:
            self.passed['0'] = 'Passed sfaira TODOs checks.'

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
