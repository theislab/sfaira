import logging
import numpy as np
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
    doi_preprint: str = ''  # doi of data set accompanying preprint manuscript
    doi_journal: str = ''  # doi of data set accompanying journal manuscript
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
    year: int = 2021  # year in which sample was acquired
    number_of_datasets: str = 1  # Required to determine the file names

    assay_sc_obs_key: str = ''  # Assay key in obs
    cell_type_obs_key: str = ''  # Original cell type key in obs
    development_stage_obs_key: str = ''  # Development stage key in obs
    disease_obs_key: str = ''  # Disease key in obs
    organ_obs_key: str = ''  # Organ key in obs
    sex_obs_key: str = ''  # Sex key in obs

    gene_id_ensembl_var_key: str = ''  # Gene id ensembl key in var
    gene_id_symbols_var_key: str = ''  # Gene id symbols key in var


class DataloaderCreator:

    def __init__(self, path_loader):
        self.WD = os.path.dirname(__file__)
        self.TEMPLATES_PATH = f'{self.WD}/templates'
        self.template_attributes = TemplateAttributes()
        self.out_path = path_loader

    def create_dataloader(self):
        """
        Prompts and guides the user through a number of possible dataloader choices.
        Prompts the user for required attributes which must be present in the dataloader.
        Finally creates the specific cookiecutter dataloader template.
        """
        self._prompt_dataloader_template()
        self._prompt_dataloader_configuration(path_data)
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

    def _prompt_dataloader_configuration(self, path_data):
        """
        Prompts the user for all required attributes for a dataloader such as DOI, author, etc.
        """
        author = sfaira_questionary(function='text',
                                    question='Author(s):',
                                    default='Einstein, Albert; Hawking, Stephen')
        self.template_attributes.author = author.split(';') if ';' in author else author
        doi = ""
        counter = 0
        while doi == "":
            if counter > 0:
                print('[bold red]You need to supply either a preprint or a jounral DOI!'
                      'Use no_doi_AUTHOR_SOME-NAME to name data sets that do not have a corresponding publication')
            doi_preprint = sfaira_questionary(
                function='text',
                question='DOI of the preprint publication: [10.1000/j.journal.2021.01.001]',
                default='')
            if doi_preprint != "":
                while not re.match(r'\b10\.\d+/[\w.]+\b', doi_preprint) and not doi_preprint.startswith("no_doi"):
                    print(f'[bold red]The entered DOI {doi_preprint} is malformed and needs to be an exact DOI! '
                          'You may leave this field empty if you supply doi_journal.')
                    doi_preprint = sfaira_questionary(
                        function='text',
                        question='DOI of the preprint publication: [10.1000/j.journal.2021.01.001]',
                        default='')
            self.template_attributes.doi_preprint = doi_preprint
            doi_journal = sfaira_questionary(
                function='text',
                question='DOI of the journal publication: [10.1000/j.journal.2021.01.001]',
                default='')
            if doi_journal != "":
                while not re.match(r'\b10\.\d+/[\w.]+\b', doi_journal) and not doi_journal.startswith("no_doi"):
                    print(f'[bold red]The entered DOI {doi_preprint} is malformed and needs to be an exact DOI! '
                          'You may leave this field empty if you supply doi_preprint.')
                    doi_journal = sfaira_questionary(
                        function='text',
                        question='DOI of the journal publication: [10.1000/j.journal.2021.01.001]',
                        default='')
            self.template_attributes.doi_journal = doi_journal
            doi = doi_journal if doi_journal != "" else doi_preprint
            counter += 1
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
            self.template_attributes.sample_fns = sfaira_questionary(
                function='text',
                question='Sample file name of the first dataset (these will be made accessible to load() via '
                         'sample_fn):',
                default='data.h5ad')

        print('[bold blue]Meta data dataset-wise (encouraged).')
        print('[bold blue]The following meta data queries are on dataset-resolution. '
              'These can also later be modified manually in the .yaml which has the same effect as setting them here.')
        first_author = author[0] if isinstance(author, list) else author
        try:
            first_author_lastname = first_author.split(',')[0]
        except KeyError:
            print('[bold yellow] First author was not in the expected format. Using full first author for the id.')
            first_author_lastname = first_author
        self.template_attributes.id_without_doi = f'{clean_id_str(self.template_attributes.organism)}_' \
                                                  f'{clean_id_str(self.template_attributes.organ)}_' \
                                                  f'{clean_id_str(str(self.template_attributes.year))}_' \
                                                  f'{clean_id_str(self.template_attributes.assay_sc)}_' \
                                                  f'{clean_id_str(first_author_lastname)}_001'
        self.template_attributes.id = f'{self.template_attributes.id_without_doi}_' \
                                      f'{self.template_attributes.doi_sfaira_repr}'
        self.template_attributes.download_url_data = sfaira_questionary(function='text',
                                                                        question='URL to download the data',
                                                                        default='')
        self.template_attributes.download_url_meta = sfaira_questionary(function='text',
                                                                        question='URL to download the meta data',
                                                                        default='')
        self.template_attributes.primary_data = str(sfaira_questionary(function='confirm',
                                                                       question='Primary data:',
                                                                       default='Yes'))
        self.template_attributes.default_embedding = str(sfaira_questionary(
            function='text',
            question='Key of default embedding in obsm (if available):',
            default='X_umap'))
        self.template_attributes.year = sfaira_questionary(function='text',
                                                           question='Year:',
                                                           default="2021")
        self.template_attributes.organism = sfaira_questionary(function='text',
                                                               question='Organism (from NCBItaxon):',
                                                               default='')
        self.template_attributes.organ = sfaira_questionary(function='text',
                                                            question='Organ [from UBERON, can also be set per cell later]:',
                                                            default='')
        self.template_attributes.assay_sc = sfaira_questionary(function='text',
                                                               question='Assay (from EFO, can also be set per cell later):',
                                                               default='')
        self.template_attributes.normalization = sfaira_questionary(function='text',
                                                                    question='Normalization:',
                                                                    default='raw')
        self.template_attributes.disease = sfaira_questionary(function='text',
                                                              question='Disease (from MONDO, can also be set per cell later):',
                                                              default='healthy')
        self.template_attributes.sample_source = sfaira_questionary(
            function='text',
            question='Sample source ["primary_tissue", "2d_culture", "3d_culture", "tumor"]:',
            default='')

        # Encouraged meta data that tend to be in .obs and would require mapping .tsv:
        print('[bold blue]Meta data cell-wise (encouraged).')
        print('[bold blue]The following meta data queries are on cell-resolution and do not need to be set if they can '
              'be annotated on a dataset-level. Skip cell- and dataset-level entries if the corresponding item is not '
              'available all-together. These items can also later be modified manually in the .yaml which has the same '
              'effect as setting them here, we do additionally check here if you need to run '
              '\'sfaira annotate-dataloader\' though.')
        self.template_attributes.assay_sc_obs_key = sfaira_questionary(
            function='text',
            question='Key of assay annotation field in .obs (from EFO, not required if "Assay" was set above):',
            default='')
        if self.template_attributes.assay_sc != "" and self.template_attributes.assay_sc_obs_key != "":
            print("[bold yellow]WARNING: you set 'Assay' before already, this new entry on assay is ignored."
                  "Define this meta data item either data set wide (as before) or cell-specific (here).")
        self.template_attributes.cell_type_obs_key = sfaira_questionary(
            function='text',
            question='Key of cell type annotation field in .obs (from CL):',
            default='')
        self.template_attributes.development_stage_obs_key = sfaira_questionary(
            function='text',
            question='Key of developmental stage annotation field in .obs (from hsapdv for human or mmusdv for mouse):',
            default='')
        self.template_attributes.disease_obs_key = sfaira_questionary(
            function='text',
            question='Key of disease annotation field in .obs (from MONDO, not required if "Disease" was set above):',
            default='')
        if self.template_attributes.disease != "" and self.template_attributes.disease_obs_key != "":
            print("[bold yellow]WARNING: you set 'Disease' before already, this new entry on disease is ignored."
                  "Define this meta data item either data set wide (as before) or cell-specific (here).")
        self.template_attributes.organ_obs_key = sfaira_questionary(
            function='text',
            question='Key of organ annotation field in .obs (from UBERON, not required if "Organ" was set above]:',
            default='')
        if self.template_attributes.organ != "" and self.template_attributes.organ_obs_key != "":
            print("[bold yellow]WARNING: you set 'Organ' before already, this new entry on organ is ignored."
                  "Define this meta data item either data set wide (as before) or cell-specific (here).")
        self.template_attributes.sex_obs_key = sfaira_questionary(
            function='text',
            question='Key of sex annotation field in .obs (from PATO, not required if "Sex" will be set manually):',
            default='')
        requires_annotate = np.any([x != "" for x in [
            self.template_attributes.assay_sc_obs_key,
            self.template_attributes.cell_type_obs_key,
            self.template_attributes.development_stage_obs_key,
            self.template_attributes.disease_obs_key,
            self.template_attributes.organ_obs_key,
            self.template_attributes.sex_obs_key,
        ]])

        print('[bold blue]Meta data gene-wise (encouraged).')
        print('[bold blue]The following meta data queries are on gene-resolution. '
              'These can also later be modified manually in the .yaml which has the same effect as setting them here.')
        self.template_attributes.gene_id_symbols_var_key = sfaira_questionary(
            function='text',
            question='Key of gene symbol field in .var:',
            default='index')
        self.template_attributes.gene_id_ensembl_var_key = sfaira_questionary(
            function='text',
            question='Key of gene id ensembl field in .var:',
            default='')

        print('[bold orange]Sfaira butler: "Up next:"')
        self.action_counter = 1
        path_loader = os.path.join(self.out_path, self.template_attributes.doi_sfaira_repr)
        print(f'[bold orange]{self.action_counter}) Proceed to modify the .yaml and .py files in {path_loader}')
        self.action_counter += 1
        self.create_datadir(path_data=path_data)
        if requires_annotate:
            print(f'[bold orange]{self.action_counter}) You will have to run \'sfaira annotate-dataloader\'.')
            self.action_counter += 1
        else:
            print(f'[bold orange]{self.action_counter}) Proceed to run \'sfaira test-dataloader\', you do not need to '
                  f'run \'sfaira annotate-dataloader\'.')
            self.action_counter += 1

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
        path_data = os.path.join(path_data, self.template_attributes.doi_sfaira_repr)
        os.makedirs(path_data)
        print(f'[bold orange]{self.action_counter}) Proceed to copy unmodified downloaded data files into {path_data} '
              f'for annotation and/or testing.')
        self.action_counter += 1
