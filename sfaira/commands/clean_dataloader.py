import logging
from subprocess import Popen

log = logging.getLogger(__name__)


class DataloaderCleaner:

    def __init__(self, path):
        self.path = path

    def clean_dataloader(self) -> None:
        """
        Removes any unwanted artifacts from a dataloader Python script and formats the code.
        1. Any line that starts with # self.      <- outcommented attribute
        2. Any line that starts with # SFARA:     <- explicitly marked full comments
        3. Any line with             # something  <- comments after attributes
        4. Runs black
        """
        # Remove all unwanted artifacts
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
                        if len(line.split('#')) > 1:
                            try:
                                cleaned_content += f'{line.split("#")[0]}\n'
                            except KeyError:
                                cleaned_content += line
                    else:
                        cleaned_content += line

        with open(self.path, 'w') as data_loader_file:
            for line in cleaned_content:
                data_loader_file.write(line)
            data_loader_file.write('\n')

        # run black
        black = Popen(['black', self.path],
                      universal_newlines=True, shell=False, close_fds=True)
        (black_stdout, black_stderr) = black.communicate()
