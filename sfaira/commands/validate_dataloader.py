import sys
import logging
import os
import pydoc
import rich
from rich.panel import Panel
from rich.progress import Progress, BarColumn
from typing import Dict, List
import yaml

from sfaira.commands.utils import get_pydoc
from sfaira.consts.utils import clean_doi

log = logging.getLogger(__name__)


class DataloaderValidator:
    fns_loaders: List[str]
    paths_yamls: List[str]
    yaml_dicts: Dict[str, dict]
    passed: List[str]
    warned: List[str]
    failed: List[str]

    def __init__(self, path_loader, doi, schema):
        self.doi = doi
        self._clean_doi = clean_doi(doi)

        loader_filenames = [str(x) for x in os.listdir(os.path.join(path_loader, self._clean_doi))
                            if str(x).endswith(".py") and str(x) != "__init__.py"]
        yaml_filenames = [x for x in os.listdir(os.path.join(path_loader, self._clean_doi))
                          if str(x).endswith(".yaml")]
        self.path_loader = path_loader
        self.fns_loaders = loader_filenames
        self.paths_yamls = [os.path.join(path_loader, self._clean_doi, x) for x in yaml_filenames]
        self.yaml_dicts = {}
        self.passed = []
        self.warned = []
        self.failed = []
        self.loader_validation_functions: list = [
            '_validate_namespace',
        ]
        self.yaml_validation_functions: list = [
            '_validate_yaml_sections',
            '_validate_required_attributes',
        ]
        self.schema = schema

    def validate(self) -> None:
        """
        Statically verifies data loader files against a predefined set of rules.
        """
        progress = Progress("[bold green]{task.description}", BarColumn(bar_width=None),
                            "[bold yellow]{task.completed} of {task.total}[reset] [bold green]{task.fields[func_name]}")
        self._validate_loader(progress=progress)
        self._validate_yaml(progress=progress)
        self._print_results()

    def _validate_loader(self, progress) -> None:
        """
        Statically verifies a loader (.py) file against a predefined set of rules.
        Every rule is a function defined in this class, which must be part of this class' linting_functions.
        """
        with progress:
            lint_progress = progress.add_task("Running loader lint checks",
                                              total=len(self.yaml_validation_functions),
                                              func_name=self.yaml_validation_functions)
            for fun_name in self.loader_validation_functions:
                progress.update(lint_progress, advance=1, func_name=fun_name)
                getattr(self, fun_name)()

    def _validate_namespace(self):
        """
        Verifies that namespace of .py files contain all necessary elements.

            - .py files contain load() function.
        """
        passed_required_elements = True

        elements = ['load', ]

        for fn in self.fns_loaders:
            file_module = ".".join(fn.split(".")[:-1])
            _, pydoc_handle = get_pydoc(path_loader=self.path_loader, doi_sfaira_repr=self._clean_doi)
            pydoc_handle = pydoc_handle + "." + file_module
            for x in elements:
                query_element = pydoc.locate(pydoc_handle + "." + x)
                if query_element is None:
                    passed_required_elements = False
                    self.failed.append(f'Missing element in namespace of {fn}.py file: {x}')

        if passed_required_elements:
            self.passed.append('Passed required data loader load() namespace checks.')

    def _validate_yaml(self, progress) -> None:
        """
        Statically verifies the yaml data loader files against a predefined set of rules.
        Every rule is a function defined in this class, which must be part of this class' linting_functions.
        """
        for x in self.paths_yamls:
            with open(x) as yaml_file:
                self.yaml_dicts[x] = yaml.load(yaml_file, Loader=yaml.FullLoader)

        with progress:
            lint_progress = progress.add_task("Running yaml lint checks",
                                              total=len(self.yaml_validation_functions),
                                              func_name=self.yaml_validation_functions)
            for fun_name in self.yaml_validation_functions:
                progress.update(lint_progress, advance=1, func_name=fun_name)
                getattr(self, fun_name)()

    def _validate_yaml_sections(self):
        """
        Verifies that .yaml files contain correct sections and keys within each section.

        This can prevent copy or delete mistakes, for example.
        """
        passed_required_sections = True

        # Note: Each key in 'dataset_or_observation_wise' appears as written below and as f"{x}_obs_key" in the .yaml.
        # Each key is accordingly listed twice and the list of lists is flattened in this dictionary:
        section_keys = {
            'dataset_structure': [
                'dataset_index',
                'sample_fns'],
            'dataset_wise': [
                'author',
                'default_embedding',
                'doi_journal',
                'doi_preprint',
                'download_url_data',
                'download_url_meta',
                'primary_data',
                'year'],
            'layers': [
                'layer_counts',
                'layer_processed',
                'layer_spliced_counts',
                'layer_spliced_processed',
                'layer_unspliced_counts',
                'layer_unspliced_processed',
                'layer_velocity'],
            'dataset_or_feature_wise': [
                'feature_reference',
                'feature_reference_var_key',
                'feature_type',
                'feature_type_var_key'],
            'dataset_or_observation_wise': [z for y in [[x, x + "_obs_key"] for x in [
                'assay_sc',
                'assay_differentiation',
                'assay_type_differentiation',
                'bio_sample',
                'cell_line',
                'cell_type',
                'development_stage',
                'disease',
                'ethnicity',
                'gm',
                'individual',
                'organ',
                'organism',
                'sample_source',
                'sex',
                'source_doi',
                'state_exact',
                'tech_sample',
                'treatment']] for z in y],
            'feature_wise': [
                'feature_id_var_key',
                'feature_symbol_var_key'],
            'observation_wise': [
                'spatial_x_coord_obs_key',
                'spatial_y_coord_obs_key',
                'spatial_z_coord_obs_key',
                'vdj_vj_1_obs_key_prefix',
                'vdj_vj_2_obs_key_prefix',
                'vdj_vdj_1_obs_key_prefix',
                'vdj_vdj_2_obs_key_prefix',
                'vdj_c_call_obs_key_suffix',
                'vdj_consensus_count_obs_key_suffix',
                'vdj_d_call_obs_key_suffix',
                'vdj_duplicate_count_obs_key_suffix',
                'vdj_j_call_obs_key_suffix',
                'vdj_junction_obs_key_suffix',
                'vdj_junction_aa_obs_key_suffix',
                'vdj_locus_obs_key_suffix',
                'vdj_productive_obs_key_suffix',
                'vdj_v_call_obs_key_suffix'],
            'meta': [
                'version'],
        }

        for fn, content in self.yaml_dicts.items():
            for k, v in section_keys.items():
                if k not in content.keys():
                    passed_required_sections = False
                    self.failed.append(f'Missing section {k} in file {fn}')
                else:
                    for kk in v:
                        if kk not in content[k].keys():
                            passed_required_sections = False
                            self.failed.append(f'Missing key {kk} in section {k} in file {fn}')

        if passed_required_sections:
            self.passed.append('Passed required data loader section checks.')

    def _validate_required_attributes(self):
        """
        Verifies that all required attributes for every data loader are present.
        """
        passed_required_attributes = True

        attributes = {
            'dataset_structure': [
                'dataset_index'],
            'dataset_wise': [
                'author',
                ['doi_journal', 'doi_preprint'],
                'download_url_data',
                'primary_data',
                'year'],
            'layers': [
                ['layer_counts', 'layer_processed']],
            'dataset_or_feature_wise': [
                ['feature_type', 'feature_type_var_key']],
            'feature_wise': [
                ['feature_id_var_key', 'feature_symbol_var_key']],
            'meta': [
                'version'],
        }
        if self.schema == "sfaira":
            attributes['dataset_or_observation_wise'] = [[x, x + "_obs_key"] for x in [
                'assay_sc',
                'organism']]
        elif self.schema == "cellxgene":
            attributes['dataset_wise'].append("default_embedding")
            attributes['dataset_or_observation_wise'] = [[x, x + "_obs_key"] for x in [
                'assay_sc',
                'cell_type',
                'development_stage',
                'disease',
                'ethnicity',
                'organ',
                'organism',
                'sex']]
        else:
            print(f"[bold red]Did not recognize schema {self.schema} in validate.")
            sys.exit()

        for fn, yaml_dict in self.yaml_dicts.items():
            for section, attrs in attributes.items():
                for attr in attrs:
                    detected = False
                    # Lists of attributes are handled in the following way:
                    # One of the two keys has to have a value.
                    if isinstance(attr, list):
                        for sub_attr in attr:
                            val = yaml_dict[section][sub_attr]
                            if val is not None:
                                detected = True
                        if not detected:
                            passed_required_attributes = False
                            self.failed.append('Missing any of the following keys to describe the required meta ' +
                                               f'data {section}:{attr} in file {fn}')
                    else:
                        val = yaml_dict[section][attr]
                        if val is not None:
                            detected = True
                        if not detected:
                            passed_required_attributes = False
                            self.failed.append(f'Missing required meta data {section}:{attr} in file {fn}')

        if passed_required_attributes:
            self.passed.append('Passed required meta data checks.')

    def _print_results(self):
        console = rich.console.Console()
        console.print()
        console.rule("[bold green] LINT RESULTS")
        console.print()
        console.print(
            f'     [bold green][[\u2714]] {len(self.passed):>4} tests passed\n'
            f'     [bold yellow][[!]] {len(self.warned):>4} tests had warnings\n'
            f'     [bold red][[\u2717]] {len(self.failed):>4} tests failed',
            overflow="ellipsis",
            highlight=False,
        )

        def format_result(linting_results: dict, color):
            results = []
            for result in linting_results:
                results.append(f'[bold {color}]Result: {result}')
            return "\n".join(results)

        if len(self.passed) > 0:
            console.print()
            console.rule("[bold green][[\u2714]] Tests Passed", style='green')
            console.print(Panel(format_result(self.passed, 'green'), style='green'), no_wrap=False, overflow='ellipsis')
        if len(self.warned) > 0:
            console.print()
            console.rule("[bold yellow][[!]] Test Warnings", style='yellow')
            console.print(Panel(format_result(self.warned, 'yellow'), style="yellow"), no_wrap=False,
                          overflow='ellipsis')
        if len(self.failed) > 0:
            console.print()
            console.rule("[bold red][[\u2717]] Test Failures", style='red')
            console.print(Panel(format_result(self.failed, 'red'), style='red'), no_wrap=False, overflow='ellipsis')
