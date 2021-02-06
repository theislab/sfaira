import logging
import os

from sfaira.commands.questionary import sfaira_questionary
from rich import print

log = logging.getLogger(__name__)


class CreateDataloader:
    WD = os.path.dirname(__file__)
    TEMPLATES_PATH = f'{WD}/templates'

    @classmethod
    def create_dataloader(cls):
        cls._prompt_dataloader_template()
        cls._prompt_dataloader_configuration()
        cls._create_dataloader_template()

    @classmethod
    def _prompt_dataloader_template(cls):
        answer = sfaira_questionary(function='select',
                                    question='what do you want',
                                    choices=['Apple', 'Banana'],
                                    default='Apple')
        print(answer)
        with open(f'{cls.TEMPLATES_PATH}/test.txt') as f:
            print(f.readlines())

    @classmethod
    def _prompt_dataloader_configuration(cls):
        pass

    @classmethod
    def _create_dataloader_template(cls):
        pass
