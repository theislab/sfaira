import logging
import os

from sfaira.commands.questionary import sfaira_questionary
from rich import print

log = logging.getLogger(__name__)


class CreateDataloader:
    def __init__(self):
        self.WD = os.path.dirname(__file__)
        self.TEMPLATES_PATH = f'{self.WD}/templates'

    def test(self):
        answer = sfaira_questionary(function='select',
                                    question='what do you want',
                                    choices=['Apple', 'Banana'],
                                    default='Apple')
        print(answer)
        with open(f'{self.WD}/templates/test.txt') as f:
            print(f.readlines())
