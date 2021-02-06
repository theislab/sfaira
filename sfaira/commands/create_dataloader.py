import logging
import os

from sfaira.commands.questionary import sfaira_questionary
from rich import print
from switchlang import switch

log = logging.getLogger(__name__)


class DataloaderCreator:
    WD = os.path.dirname(__file__)
    TEMPLATES_PATH = f'{WD}/templates'

    @classmethod
    def create_dataloader(cls):
        """
        Prompts and guides the user through a number of possible dataloader choices.
        Prompts the user for required attributes which must be present in the dataloader.
        Finally creates the specific cookiecutter dataloader template.
        """
        dataloader_template_type = cls._prompt_dataloader_template()
        cls._prompt_dataloader_configuration(dataloader_template_type)
        cls._create_dataloader_template()

    @classmethod
    def _prompt_dataloader_template(cls) -> str:
        """
        Guides the user to select the appropriate dataloader template for his dataset.

        :return: Type of desired dataloader. One of
            single_dataset
            multiple_datasets_single_file
            multiple_datasets_streamlined
            multiple_datasets_not_streamlined
        """
        number_datasets = sfaira_questionary(function='select',
                                             question='How many datasets does your project have?',
                                             choices=['One', 'More than one'])
        # One dataset
        if number_datasets == 'One':
            return 'single_dataset'

        # More than one dataset
        dataset_counts = sfaira_questionary(function='select',
                                            question='Are your datasets in a single file or is there one file per dataset?',
                                            choices=['Single dataset file', 'Multiple dataset files'])
        if dataset_counts == 'Single dataset file':
            return 'multiple_datasets_single_file'

        # streamlined?
        streamlined_datasets = sfaira_questionary(function='select',
                                                  question='Are your datasets in a similar format?',
                                                  choices=['Same format', 'Different formats'])
        if streamlined_datasets == 'Same format':
            return 'multiple_datasets_streamlined'
        else:
            return 'multiple_datasets_not_streamlined'

    @classmethod
    def _prompt_dataloader_configuration(cls, dataloader_template_type: str):
        """
        Prompts the user for all required attributes for a dataloader such as DOI, author, etc.

        :param dataloader_template_type: One of
            single_dataset
            multiple_datasets_single_file
            multiple_datasets_streamlined
            multiple_datasets_not_streamlined
        """
        with switch(dataloader_template_type) as s:
            s.case('single_dataset', lambda: print('selected single dataset'))
            s.case('multiple_datasets_single_file', lambda: print('selected multiple datasets single file'))
            s.case('multiple_datasets_streamlined', lambda: print('selected streamlined'))
            s.case('multiple_datasets_not_streamlined', lambda: print('selected not streamlined'))
            s.default(lambda error: print('[bold red] Invalid dataloader type select. Internal error!'))

    @classmethod
    def _create_dataloader_template(cls):
        pass
        # ensure that it handles multiple dataloaders on the same DOI

        # mb add a clean command which gets rids of outcommented self.whatever -> provide everything by default and then just clean it up
