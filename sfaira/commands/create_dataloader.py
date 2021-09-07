import logging
import os
import re
from dataclasses import dataclass, asdict
from typing import Union, Dict

from sfaira.consts.utils import clean_doi, clean_id_str
from sfaira.commands.questionary import sfaira_questionary
from rich import print
from cookiecutter.main import cookiecutter

log = logging.getLogger(__name__)


@dataclass
class TemplateAttributes:
    dataloader_type: str = ''  # One of single_dataset, multiple_datasets_single_file, multiple_datasets_streamlined, multiple_datasets_not_streamlined
    id: str = ''  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).
    id_without_doi: str = ''  # complete id without the doi -> usually used to name the python scripts
    create_extra_description: str = ''  # Whether to create an optional extra description file or not

    author: Union[str, list] = ''  # author (list) who sampled / created the data set
    doi: str = ''  # doi of data set accompanying manuscript
    doi_sfaira_repr: str = ''  # internal representation with any special characters replaced with underscores

    sample_fns: Union[str, Dict[str, list]] = ''  # file name of the first *.h5ad file
    download_url_data: str = ''  # download website(s) of data files
    download_url_meta: str = ''  # download website(s) of meta data files
    organ: str = ''  # (*) organ (anatomical structure)
    organism: str = ''  # (*) species / organism
    assay_sc: str = ''  # (*, optional) protocol used to sample data (e.g. smart-seq2)
    normalization: str = ''  # raw or the used normalization technique
    default_embedding: str = ''  # Default embedding of the data
    primary_data: str = ''  # Is this a primary dataset?
    disease: str = ''  # name of the disease of the condition
    ethnicity: str = ''  # ethnicity of the sample
    sample_source: str = ''  # source of the sample
    state_exact: str = ''  # state of the sample
    year: str = 2021  # year in which sample was acquired
    number_of_datasets: str = 1  # Required to determine the file names

    cell_types_original_obs_key: str = ''  # Original cell type key in obs


class DataloaderCreator:

    def __init__(self, path_loader, doi):
        self.WD = os.path.dirname(__file__)
        self.TEMPLATES_PATH = f'{self.WD}/templates'
        self.template_attributes = TemplateAttributes()
        self.out_path = path_loader
        self.doi = doi

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
        else:
            self.template_attributes.dataloader_type = 'multiple_datasets'

    def _prompt_dataloader_configuration(self):
        """
        Prompts the user for all required attributes for a dataloader such as DOI, author, etc.
        """
        author = sfaira_questionary(function='text',
                                    question='Author(s):',
                                    default='Einstein, Albert; Hawking, Stephen')
        self.template_attributes.author = author.split(';') if ';' in author else author
        if self.doi:
            doi = self.doi
        else:
            doi = sfaira_questionary(function='text',
                                     question='DOI:',
                                     default='10.1000/j.journal.2021.01.001')
            while not re.match(r'\b10\.\d+/[\w.]+\b', doi):
                print('[bold red]The entered DOI is malformed!')
                doi = sfaira_questionary(function='text',
                                         question='DOI:',
                                         default='10.1000/j.journal.2021.01.001')
        self.template_attributes.doi = doi
        self.template_attributes.doi_sfaira_repr = clean_doi(doi)

        self.template_attributes.number_of_datasets = sfaira_questionary(function='text',
                                                                         question='Number of datasets:',
                                                                         default='1')

        # Differentiate between a single dataset or multiple datasets to get sample file names
        if self.template_attributes.dataloader_type == 'multiple_datasets':
            self.template_attributes.sample_fns = {'fns': []}
            for ds in range(int(self.template_attributes.number_of_datasets)):
                fn = sfaira_questionary(function='text',
                                        question='Sample file name:',
                                        default=f'data_{ds}.h5ad')
                self.template_attributes.sample_fns['fns'].append(fn)
        else:
            self.template_attributes.sample_fns = sfaira_questionary(function='text',
                                                                     question='Sample file name of the first dataset:',
                                                                     default='data.h5ad')

        self.template_attributes.primary_data = str(sfaira_questionary(function='confirm',
                                                                       question='Primary data:',
                                                                       default='Yes'))
        self.template_attributes.default_embedding = sfaira_questionary(function='text',
                                                                        question='Default embedding:',
                                                                        default='NA')
        self.template_attributes.organism = sfaira_questionary(function='text',
                                                               question='Organism:',
                                                               default='NA')
        self.template_attributes.organ = sfaira_questionary(function='text',
                                                            question='Organ:',
                                                            default='NA')
        self.template_attributes.assay_sc = sfaira_questionary(function='text',
                                                               question='Assay:',
                                                               default='NA')
        self.template_attributes.normalization = sfaira_questionary(function='text',
                                                                    question='Normalization:',
                                                                    default='raw')
        self.template_attributes.disease = sfaira_questionary(function='text',
                                                              question='Disease:',
                                                              default='healthy')
        self.template_attributes.state_exact = sfaira_questionary(function='text',
                                                                  question='Sample state:',
                                                                  default='healthy')
        self.template_attributes.sample_source = sfaira_questionary(function='text',
                                                                    question='Sample source:',
                                                                    default='NA')
        is_cell_type_annotation = sfaira_questionary(function='confirm',
                                                     question='Does your dataset have a cell type annotation?',
                                                     default='No')
        if is_cell_type_annotation:
            self.template_attributes.cell_types_original_obs_key = sfaira_questionary(function='text',
                                                                                      question='Cell type annotation obs key:',
                                                                                      default='')
        self.template_attributes.year = sfaira_questionary(function='text',
                                                           question='Year:',
                                                           default='2021')
        first_author = author[0] if isinstance(author, list) else author
        try:
            first_author_lastname = first_author.split(',')[0]
        except KeyError:
            print('[bold yellow] First author was not in the expected format. Using full first author for the id.')
            first_author_lastname = first_author
        self.template_attributes.id_without_doi = f'{clean_id_str(self.template_attributes.organism)}_' \
                                                  f'{clean_id_str(self.template_attributes.organ)}_' \
                                                  f'{clean_id_str(self.template_attributes.year)}_' \
                                                  f'{clean_id_str(self.template_attributes.assay_sc)}_' \
                                                  f'{clean_id_str(first_author_lastname)}_001'
        self.template_attributes.id = f'{self.template_attributes.id_without_doi}_' \
                                      f'{self.template_attributes.doi_sfaira_repr}'
        if self.template_attributes.dataloader_type == 'single_dataset':
            self.template_attributes.download_url_data = sfaira_questionary(function='text',
                                                                            question='URL to download the data',
                                                                            default='https://ftp.ncbi.nlm.nih.gov/geo/')
            self.template_attributes.download_url_meta = sfaira_questionary(function='text',
                                                                            question='URL to download the meta data',
                                                                            default='https://ftp.ncbi.nlm.nih.gov/geo/')
        self.template_attributes.create_extra_description = sfaira_questionary(function='confirm',
                                                                               question='Do you want to add additional custom metadata?',
                                                                               default='Yes')
        if is_cell_type_annotation:
            print('[bold blue]You will have to run \'sfaira annotate-dataloader\' after the template has been created and filled.')
        else:
            print('[bold blue]You can skip \'sfaira annotate-dataloader\'.')

    def _template_attributes_to_dict(self) -> dict:
        """
        Create a dict from the our Template Structure dataclass
        :return: The dict containing all key-value pairs with non empty values
        """
        return {key: val for key, val in asdict(self.template_attributes).items() if val != ''}

    def _create_dataloader_template(self):
        template_path = f'{self.TEMPLATES_PATH}/{self.template_attributes.dataloader_type}'
        cookiecutter(f'{template_path}',
                     output_dir=self.out_path,
                     no_input=True,
                     overwrite_if_exists=True,
                     extra_context=self._template_attributes_to_dict())

    def create_datadir(self, path_data):
        os.makedirs(os.path.join(path_data, self.template_attributes.doi_sfaira_repr))
