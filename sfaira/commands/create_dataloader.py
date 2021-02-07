import logging
import os
import re
from dataclasses import dataclass
from typing import Union

from sfaira.commands.questionary import sfaira_questionary
from rich import print
from switchlang import switch

log = logging.getLogger(__name__)


@dataclass
class TemplateAttributes:
    dataloader_type: str = ''
    id: str = ''  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).

    author: Union[str, list] = ''  # author (list) who sampled / created the data set
    doi: str = ''  # doi of data set accompanying manuscript

    download_url_data: str = ''  # download website(s) of data files
    download_url_meta: str = ''  # download website(s) of meta data files

    age: str = '0'  # (*, optional) age of sample
    dev_stage: str = ''  # (*, optional) developmental stage of organism
    ethnicity: str = ''  # (*, optional) ethnicity of sample
    healthy: bool = True  # (*, optional) whether sample represents a healthy organism
    normalisation: str = ''  # (optional) normalisation applied to raw data loaded (ideally counts, "raw")
    organ: str = ''  # (*, optional) organ (anatomical structure)
    organism: str = ''  # (*) species / organism
    protocol: str = ''  # (*, optional) protocol used to sample data (e.g. smart-seq2)
    sex: str = ''  # (*, optional) sex
    state_exact: str = ''  # (*, optional) exact disease, treatment or perturbation state of sample
    year: int = 2021  # year in which sample was acquired

    # The following meta data may instead also be supplied on a cell level if an appropriate column is present in the
    # anndata instance (specifically in .obs) after loading. You need to make sure this is loaded in the loading script)!
    obs_key_age: int = 0  # (optional, see above, do not provide if .age is provided)
    obs_key_dev_stage: str = ''  # (optional, see above, do not provide if .dev_stage is provided)
    obs_key_ethnicity: str = ''  # (optional, see above, do not provide if .ethnicity is provided)
    obs_key_healthy: str = ''  # (optional, see above, do not provide if .healthy is provided)
    obs_key_organ: str = ''  # (optional, see above, do not provide if .organ is provided)
    obs_key_organism: str = ''  # (optional, see above, do not provide if .organism is provided)
    obs_key_protocol: str = ''  # (optional, see above, do not provide if .protocol is provided)
    obs_key_sex: str = ''  # (optional, see above, do not provide if .sex is provided)
    obs_key_state_exact: str = ''  # (optional, see above, do not provide if .state_exact is provided)
    # Additionally, cell type annotation is ALWAYS provided per cell in .obs, this annotation is optional though.
    # name of column which contain streamlined cell ontology cell type classes:
    obs_key_cellontology_original: str = ''  # (optional)


class DataloaderCreator:

    def __init__(self):
        self.WD = os.path.dirname(__file__)
        self.TEMPLATES_PATH = f'{self.WD}/templates'
        self.template_attributes = TemplateAttributes()

    def create_dataloader(self):
        """
        Prompts and guides the user through a number of possible dataloader choices.
        Prompts the user for required attributes which must be present in the dataloader.
        Finally creates the specific cookiecutter dataloader template.
        """
        self._prompt_dataloader_template()
        self._prompt_dataloader_configuration()
        self._create_dataloader_template()

    def _prompt_dataloader_template(self) -> None:
        """
        Guides the user to select the appropriate dataloader template for his dataset.
        Sets the dataloader_type
        """
        number_datasets = sfaira_questionary(function='select',
                                             question='How many datasets does your project have?',
                                             choices=['One', 'More than one'])
        # One dataset
        if number_datasets == 'One':
            self.template_attributes.dataloader_type = 'single_dataset'

        # More than one dataset
        dataset_counts = sfaira_questionary(function='select',
                                            question='Are your datasets in a single file or is there one file per dataset?',
                                            choices=['Single dataset file', 'Multiple dataset files'])
        if dataset_counts == 'Single dataset file':
            self.template_attributes.dataloader_type = 'multiple_datasets_single_file'

        # streamlined?
        streamlined_datasets = sfaira_questionary(function='select',
                                                  question='Are your datasets in a similar format?',
                                                  choices=['Same format', 'Different formats'])
        if streamlined_datasets == 'Same format':
            self.template_attributes.dataloader_type = 'multiple_datasets_streamlined'
        else:
            self.template_attributes.dataloader_type = 'multiple_datasets_not_streamlined'

    def _prompt_dataloader_configuration(self):
        """
        Prompts the user for all required attributes for a dataloader such as DOI, author, etc.
        """
        # Prompts shared by all templates
        # self.author = "Last name, first name"  # or ["Last name, first name", "Last name, first name"]

        doi = sfaira_questionary(function='text',
                                 question='DOI:',
                                 default='10.1000/j.journal.2021.01.001')
        while not re.match(r'\b10\.\d+/[\w.]+\b', doi):
            print('[bold red]The entered DOI is malformed!')  # noqa: W605
            doi = sfaira_questionary(function='text',
                                     question='DOI:',
                                     default='10.1000/j.journal.2021.01.001')
        self.template_attributes.doi = doi

        id: str = ''  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).

        author: str = ''  # author (list) who sampled / created the data set
        doi: str = ''  # doi of data set accompanying manuscript

        download_url_data: str = ''  # download website(s) of data files
        download_url_meta: str = ''  # download website(s) of meta data files

        age: str = '0'  # (*, optional) age of sample
        dev_stage: str = ''  # (*, optional) developmental stage of organism
        ethnicity: str = ''  # (*, optional) ethnicity of sample
        healthy: bool = True  # (*, optional) whether sample represents a healthy organism
        normalisation: str = ''  # (optional) normalisation applied to raw data loaded (ideally counts, "raw")
        organ: str = ''  # (*, optional) organ (anatomical structure)
        organism: str = ''  # (*) species / organism
        protocol: str = ''  # (*, optional) protocol used to sample data (e.g. smart-seq2)
        sex: str = ''  # (*, optional) sex
        state_exact: str = ''  # (*, optional) exact disease, treatment or perturbation state of sample
        year: int = 2021  # year in which sample was acquired

        # Prompts which are dataloader type specific
        def _single_dataset_prompts():
            pass

        def _multiple_datasets_single_file_prompts():
            pass

        def _multiple_datasets_streamlined():
            pass

        def _multiple_datasets_not_streamlined():
            pass

        with switch(self.template_attributes.dataloader_type) as s:
            s.case('single_dataset', _single_dataset_prompts)
            s.case('multiple_datasets_single_file', _multiple_datasets_single_file_prompts)
            s.case('multiple_datasets_streamlined', _multiple_datasets_streamlined)
            s.case('multiple_datasets_not_streamlined', _multiple_datasets_not_streamlined)
            s.default(lambda error: print('[bold red] Invalid dataloader type select. Internal error!'))

    @classmethod
    def _create_dataloader_template(cls):
        pass
        # ensure that it handles multiple dataloaders on the same DOI

        # mb add a clean command which gets rids of outcommented self.whatever -> provide everything by default and then just clean it up
