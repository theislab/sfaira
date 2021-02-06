import logging
import os
from dataclasses import dataclass

from sfaira.commands.questionary import sfaira_questionary
from rich import print
from switchlang import switch

log = logging.getLogger(__name__)


@dataclass
class TemplateAttributes:
    dataloader_type: str = ''


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
        # TODO

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
