import logging

import yaml
from boltons.iterutils import remap

log = logging.getLogger(__name__)


class DataloaderCleaner:

    def __init__(self, path):
        self.path = path

    def clean_dataloader(self) -> None:
        """
        Removes unused keys from the yaml file
        """
        with open(self.path) as yaml_file:
            content = yaml.load(yaml_file, Loader=yaml.FullLoader)
            drop_falsey = lambda path, key, value: bool(value)
            clean = remap(content, visit=drop_falsey)

            with open(self.path, 'w') as file:
                yaml.dump(clean, file)
