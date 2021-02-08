import logging

log = logging.getLogger(__name__)


class DataloaderCleaner:

    def __init__(self, path):
        self.path = path

    def clean_dataloader(self):
        cleaned_content = []

        with open(self.path, 'r') as data_loader_file:
            content = data_loader_file.readlines()
            for line in content:
                line_stripped = line.strip()
                if line_stripped.startswith('# self.'):
                    continue
                elif line_stripped.startswith('# SFAIRA:'):
                    continue
                else:
                    if '#' in line:
                        cleaned_content += f'{line.split("#")[0]}\n'
                    else:
                        cleaned_content += line

        with open(self.path, 'w') as data_loader_file:
            for line in cleaned_content:
                data_loader_file.write(line)
            data_loader_file.write('\n')
