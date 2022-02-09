from cookiecutter.main import cookiecutter
from dataclasses import dataclass, asdict
import logging
import numpy as np
import os
import re
from rich import print
from typing import Union, Dict

from sfaira.consts.utils import clean_doi, clean_id_str
from sfaira.commands.questionary import sfaira_questionary

log = logging.getLogger(__name__)


@dataclass
class TemplateAttributes:
    create_extra_description: str = ''  # Whether to create an optional extra description file or not
    dataloader_type: str = ''  # One of single_dataset, multiple_datasets_single_file, multiple_datasets_streamlined, multiple_datasets_not_streamlined
    doi_sfaira_repr: str = ''
    id: str = ''  # unique identifier of data set (Organism_Organ_Year_Protocol_NumberOfDataset_FirstAuthorLastname_doi).
    id_without_doi: str = ''  # complete id without the doi -> usually used to name the python scripts
    number_of_datasets: str = 1

    sample_fns: Union[str, Dict[str, list]] = ''

    author: Union[str, list] = ''
    default_embedding: str = ''
    doi_preprint: str = ''
    doi_journal: str = ''
    download_url_data: str = ''
    download_url_meta: str = ''
    organism: str = ''
    primary_data: str = ''
    year: int = 0

    layer_counts: str = ''
    layer_processed: str = ''
    layer_spliced_counts: str = ''
    layer_spliced_processed: str = ''
    layer_unspliced_counts: str = ''
    layer_unspliced_processed: str = ''
    layer_velocity: str = ''

    assay_sc: str = ''
    cell_type: str = ''
    development_stage: str = ''
    disease: str = ''
    ethnicity: str = ''
    organ: str = ''
    sample_source: str = ''
    sex: str = ''

    assay_sc_obs_key: str = ''
    cell_type_obs_key: str = ''
    development_stage_obs_key: str = ''
    disease_obs_key: str = ''
    ethnicity_obs_key: str = ''
    organ_obs_key: str = ''
    sample_source_obs_key: str = ''
    sex_obs_key: str = ''

    spatial_x_coord_obs_key: str = ''
    spatial_y_coord_obs_key: str = ''
    spatial_z_coord_obs_key: str = ''
    vdj_vj_1_obs_key_prefix: str = ''
    vdj_vj_2_obs_key_prefix: str = ''
    vdj_vdj_1_obs_key_prefix: str = ''
    vdj_vdj_2_obs_key_prefix: str = ''
    vdj_c_call_obs_key_suffix: str = ''
    vdj_consensus_count_obs_key_suffix: str = ''
    vdj_d_call_obs_key_suffix: str = ''
    vdj_duplicate_count_obs_key_suffix: str = ''
    vdj_j_call_obs_key_suffix: str = ''
    vdj_junction_obs_key_suffix: str = ''
    vdj_junction_aa_obs_key_suffix: str = ''
    vdj_locus_obs_key_suffix: str = ''
    vdj_productive_obs_key_suffix: str = ''
    vdj_v_call_obs_key_suffix: str = ''

    feature_reference: str = ''
    feature_type: str = ''

    feature_id_var_key: str = ''
    feature_symbol_var_key: str = ''


class DataloaderCreator:

    def __init__(self, path_loader):
        self.WD = os.path.dirname(__file__)
        self.TEMPLATES_PATH = os.path.join(self.WD, "templates")
        self.template_attributes = TemplateAttributes()
        self.out_path = path_loader

    def create_dataloader(self, path_data):
        """
        Prompts and guides the user through a number of possible data loader choices.
        Prompts the user for required attributes which must be present in the data loader.
        Finally creates the specific cookiecutter data loader template.
        """
        self._prompt_dataloader_template()
        self._prompt_dataloader_configuration(path_data)
        self._create_dataloader_template()

    def _prompt_dataloader_template(self) -> None:
        """
        Guides the user to select the appropriate data loader template for his dataset.
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
        Prompts the user for all required attributes for a data loader such as DOI, author, etc.
        """
        author = sfaira_questionary(function='text',
                                    question='Leading author "Last name, first name" (e.g.: Einstein, Albert)',
                                    default='')
        self.template_attributes.author = author.split(';') if ';' in author else author
        doi = ""
        counter = 0
        while doi == "":
            if counter > 0:
                print('[bold red]You need to supply either a preprint or a journal DOI!'
                      'Use "no_doi_AUTHOR_SOME-NAME" to name data sets that do not have a corresponding publication.')
            doi_preprint_try = 0
            doi_preprint = ""
            while doi_preprint_try == 0 or (
                    doi_preprint != "" and
                    not re.match(r'\b10\.\d+/[\w.]+\b', doi_preprint) and
                    not doi_preprint.startswith("no_doi")
            ):
                if doi_preprint_try > 0:
                    print(f'[bold red]The entered DOI {doi_preprint} is malformed and needs to be an exact DOI! '
                          'You may leave this field empty if you supply doi_journal.')
                doi_preprint = sfaira_questionary(
                    function='text',
                    question='DOI of the preprint publication (e.g.: 10.1000/j.journal.2022.01.001)',
                    default='')
                doi_preprint_try += 1
            self.template_attributes.doi_preprint = doi_preprint
            doi_journal_try = 0
            doi_journal = ""
            while doi_journal_try == 0 or (
                    doi_journal != "" and
                    not re.match(r'\b10\.\d+/[\w.]+\b', doi_journal) and
                    not doi_journal.startswith("no_doi")
            ):
                if doi_journal_try > 0:
                    print(f'[bold red]The entered DOI {doi_preprint} is malformed and needs to be an exact DOI! '
                          'You may leave this field empty if you supply doi_preprint.')
                doi_journal = sfaira_questionary(
                    function='text',
                    question='DOI of the journal publication (e.g.: 10.1000/j.journal.2022.01.001)',
                    default='')
                doi_journal_try += 1
            self.template_attributes.doi_journal = doi_journal
            doi = doi_journal if doi_journal != "" else doi_preprint
            counter += 1
        self.template_attributes.doi_sfaira_repr = clean_doi(doi)
        if os.path.exists(self.path_loader):
            print(f"[bold orange]A data loader for this DOI already exists: {self.path_loader}. "
                  f"If you want to write a second data loader for this publication, "
                  f"proceed with this curation workflow, take care that the old one is not overwritten, though. "
                  f"If you were not aware of the other data loader, consider aborting this phase 1a and "
                  f"decide if want delete the old version of this loader "
                  f"or consider completing the previous curation attempt in phase 1b. "
                  )

        self.template_attributes.number_of_datasets = sfaira_questionary(function='text',
                                                                         question='Number of datasets',
                                                                         default='1')

        # Differentiate between a single dataset or multiple datasets to get sample file names
        if self.template_attributes.dataloader_type == 'multiple_datasets':
            self.template_attributes.sample_fns = {'fns': []}
            for ds in range(int(self.template_attributes.number_of_datasets)):
                fn = sfaira_questionary(function='text',
                                        question=f'Filename for dataset {ds}',
                                        default=f'data_{ds}.h5ad')
                self.template_attributes.sample_fns['fns'].append(fn)

        print('[bold blue]Cell-wise meta data that are shared across the dataset.')
        print('[bold blue]The following meta data queries are on dataset-resolution. '
              'These can also later be modified manually in the .yaml which has the same effect as setting them here.')
        first_author = author[0] if isinstance(author, list) else author
        try:
            first_author_lastname = first_author.split(',')[0]
        except KeyError:
            print('[bold yellow] First author was not in the expected format. Using full first author for the id.')
            first_author_lastname = first_author

        self.template_attributes.default_embedding = str(sfaira_questionary(
            function='text',
            question='Key of default embedding in .obsm:',
            default=''))
        download_url_data_try = 0
        download_url_data = ''
        while download_url_data_try == 0 or download_url_data == '':
            if download_url_data_try > 0:
                print('[bold red]You need to supply a download url of the data.')
            download_url_data = sfaira_questionary(
                function='text',
                question='URL to download the data',
                default='')
            download_url_data_try += 1
        self.template_attributes.download_url_data = download_url_data
        self.template_attributes.download_url_meta = sfaira_questionary(
            function='text',
            question='URL to download the meta data (only necessary if different from download_url_data)',
            default='')
        organism_try = 0
        organism = ''
        while organism_try == 0 or organism == '':
            if organism_try > 0:
                print('[bold red]You need to supply the main organism.')
            organism = sfaira_questionary(
                function='text',
                question='Main organism (from NCBItaxon):',
                default='')
            organism_try += 1
        self.template_attributes.organism = organism
        self.template_attributes.primary_data = str(sfaira_questionary(
            function='confirm',
            question='Is this primary data?',
            default='Yes'))
        year_try = 0
        year = ''
        while year_try == 0 or year == '':
            if year_try > 0:
                print('[bold red]You need to supply the year of first publication.')
            year = sfaira_questionary(
                function='text',
                question='Year of first publication (e.g. 2022)',
                default="")
            year_try += 1
        self.template_attributes.year = int(year)

        # Data matrices in the dataset:
        print('[bold blue]Data matrices.')
        print('[bold blue]An AnnData object may contain multiple data matrices: raw and processed gene expression '
              'counts, or spliced and unspliced count data and velocity estimates, for example. '
              'Minimally, you need to supply either of the matrices "counts" or "processed". '
              'In the following, "*counts" refers to the INTEGER count of alignment events (e.g. transcripts for RNA). '
              'In the following, "*processed" refers to any processing that modifies these counts, for example: '
              'normalization, batch correction, ambient RNA correction. '
              'These items can also later be modified manually in the .yaml which has the same '
              'effect as setting them here.')

        def format_q_mat_key(attr) -> str:
            return f"Layer that contains {attr} (either 'X', 'raw', or a .layers key)"

        main_layer = ""
        counter = 0
        while main_layer == "":
            if counter > 0:
                print('[bold red]You need to supply either a matrix for "counts" or for "processed"!')
            self.template_attributes.layer_counts = sfaira_questionary(
                function='text',
                question=format_q_mat_key("counts"),
                default='')
            self.template_attributes.layer_processed = sfaira_questionary(
                function='text',
                question=format_q_mat_key("processed counts"),
                default='')
            main_layer = self.template_attributes.layer_counts if self.template_attributes.layer_counts != "" else \
                self.template_attributes.layer_processed
            counter += 1
        self.template_attributes.layer_spliced_counts = sfaira_questionary(
            function='text',
            question=format_q_mat_key("spliced counts"),
            default='')
        self.template_attributes.layer_spliced_processed = sfaira_questionary(
            function='text',
            question=format_q_mat_key("processed spliced counts"),
            default='')
        self.template_attributes.layer_unspliced_counts = sfaira_questionary(
            function='text',
            question=format_q_mat_key("unspliced counts"),
            default='')
        self.template_attributes.layer_unspliced_processed = sfaira_questionary(
            function='text',
            question=format_q_mat_key("processed unspliced counts"),
            default='')
        self.template_attributes.layer_velocity = sfaira_questionary(
            function='text',
            question=format_q_mat_key("gene-wise velocities"),
            default='')

        # Meta data that may also be observation-wise:
        print('[bold blue]Cell-wise meta data that are sometimes shared across the dataset.')
        print('[bold blue]The following meta data queries are on a dataset-resolution and should be set '
              'if they do not vary across the cells in that dataset. '
              'Skip the corresponding cell-wise below annotation if you annotate metadata here. '
              'If these meta data vary across cells in a data set, skip them here and annotate them in the next '
              'section. '
              'These items can also later be modified manually in the .yaml which has the same '
              'effect as setting them here. '
              'A lot of these meta data are ontology constrained.'
              'You should input symbols, ie. readable words and not IDs here.'
              'You can look up term symbols here https://www.ebi.ac.uk/ols/index.')

        def format_q_uns_key(attr, onto) -> str:
            return f"Dataset-wide {attr} annotation (from {onto})"

        self.template_attributes.assay_sc = sfaira_questionary(
            function='text',
            question=format_q_uns_key("assay", "EFO"),
            default='')
        self.template_attributes.cell_type = sfaira_questionary(
            function='text',
            question=format_q_uns_key("cell type", "CL (Cell ontology)"),
            default='')
        self.template_attributes.development_stage = sfaira_questionary(
            function='text',
            question=format_q_uns_key("developmental stage", "hsapdv for human, mmusdv for mouse"),
            default='')
        self.template_attributes.disease = sfaira_questionary(
            function='text',
            question=format_q_uns_key("disease", "MONDO"),
            default='')
        self.template_attributes.ethnicity = sfaira_questionary(
            function='text',
            question=format_q_uns_key("ethnicity", "HANCESTRO for human, skip for non-human"),
            default='')
        self.template_attributes.organ = sfaira_questionary(
            function='text',
            question=format_q_uns_key("organ/tissue", "UBERON"),
            default='')
        self.template_attributes.sample_source = sfaira_questionary(
            function='text',
            question='Dataset-wide sample source annotation (from '
                     '["primary_tissue", "2d_culture", "3d_culture", "tumor"]):',
            default='')
        self.template_attributes.sex = sfaira_questionary(
            function='text',
            question=format_q_uns_key("sex", "PATO"),
            default='')

        # Encouraged meta data that tend to be in .obs and would require mapping .tsv:
        print('[bold blue]Cell-wise meta data.')
        print('[bold blue]The following meta data queries are on cell-resolution '
              'and do not need to be set if they were annotated on a dataset-level above. '
              'Skip cell- and dataset-level entries if the corresponding item is not available at all. '
              'These items can also later be modified manually in the .yaml which has the same effect as setting them '
              'here except of that we do additionally check here if you need to run \'sfaira annotate-dataloader\'.')

        def format_q_obs_key(attr, onto) -> str:
            return f"Key of {attr} annotation field in .obs (from {onto})"

        def format_warning_double_curation(attr) -> str:
            return f"[bold yellow]WARNING: you set '{attr}' before already, this new entry is ignored. " \
                   "Define this meta data item either data set wide (as before) or cell-specific (here)."

        # assay_sc:
        self.template_attributes.assay_sc_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("assay", "EFO"),
            default='')
        if self.template_attributes.assay_sc != "" and self.template_attributes.assay_sc_obs_key != "":
            print(format_warning_double_curation('assay_sc'))
        # cell_type:
        self.template_attributes.cell_type_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("cell type", "CL"),
            default='')
        if self.template_attributes.cell_type != "" and self.template_attributes.cell_type_obs_key != "":
            print(format_warning_double_curation('cell type'))
        # development_stage:
        self.template_attributes.development_stage_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("development stage", "hsapdv for human, mmusdv for mouse"),
            default='')
        if self.template_attributes.development_stage != "" and \
                self.template_attributes.development_stage_obs_key != "":
            print(format_warning_double_curation('development_stage'))
        # disease:
        self.template_attributes.disease_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("disease", "MONDO"),
            default='')
        if self.template_attributes.disease != "" and self.template_attributes.disease_obs_key != "":
            print(format_warning_double_curation('disease'))
        # ethnicity:
        self.template_attributes.ethnicity_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("ethnicity", "HANCESTRO for human, skip for non-human"),
            default='')
        if self.template_attributes.ethnicity != "" and self.template_attributes.ethnicity_obs_key != "":
            print(format_warning_double_curation('ethnicity'))
        # organ:
        self.template_attributes.organ_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("organ/tissue", "UBERON"),
            default='')
        if self.template_attributes.organ != "" and self.template_attributes.organ_obs_key != "":
            print(format_warning_double_curation('organ'))
        # sample_source:
        self.template_attributes.sample_source_obs_key = sfaira_questionary(
            function='text',
            question='Key of sample source annotation field in .obs (from '
                     '["primary_tissue", "2d_culture", "3d_culture", "tumor"]):',
            default='')
        if self.template_attributes.sample_source != "" and self.template_attributes.sample_source_obs_key != "":
            print(format_warning_double_curation('sex'))
        # sex:
        self.template_attributes.sex_obs_key = sfaira_questionary(
            function='text',
            question=format_q_obs_key("sex", "PATO"),
            default='')
        if self.template_attributes.sex != "" and self.template_attributes.sex_obs_key != "":
            print(format_warning_double_curation('sex'))
        requires_annotate = np.any([x != "" for x in [
            self.template_attributes.assay_sc_obs_key,
            self.template_attributes.cell_type_obs_key,
            self.template_attributes.development_stage_obs_key,
            self.template_attributes.disease_obs_key,
            self.template_attributes.ethnicity_obs_key,
            self.template_attributes.organ_obs_key,
            self.template_attributes.sex_obs_key,
        ]])

        # Modality-specific data:

        def format_q_modality_obs_key(modality, attr) -> str:
            return f"{modality}: Key of {attr} annotation field in .obs"

        has_spatial = sfaira_questionary(
            function='confirm',
            question='Modality spatial: Does this data set contain spatial coordinates of observations?',
            default='No')
        if has_spatial:
            self.template_attributes.spatial_x_coord_obs_key = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("Spatial", "x coordinate"),
                default='')
            self.template_attributes.spatial_y_coord_obs_key = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("Spatial", "y coordinate"),
                default='')
            self.template_attributes.spatial_z_coord_obs_key = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("Spatial", "z coordinate"),
                default='')

        has_vdj = sfaira_questionary(
            function='confirm',
            question='Modality V(D)J: Does this data set contain V(D)J gene reconstructions by observation?',
            default='No')
        if has_vdj:
            print('[bold blue]V(D)J annotation: Meta data definitions are documented here '
                  'https://docs.airr-community.org/en/latest/datarep/rearrangements.html.')
            print('[bold blue]Below, default column names in scirpy are chose, use deafults if you read the V(D)J data '
                  'with scirpy, for example.')
            self.template_attributes.vdj_vj_1_obs_key_prefix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "prefix of first VJ locus/chain"),
                default='IR_VJ_1_')
            self.template_attributes.vdj_vj_2_obs_key_prefix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "prefix of second VJ locus/chain"),
                default='IR_VJ_2_')
            self.template_attributes.vdj_vdj_1_obs_key_prefix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "prefix of first VDJ locus/chain"),
                default='IR_VDJ_1_')
            self.template_attributes.vdj_vdj_2_obs_key_prefix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "prefix of second VDJ locus/chain"),
                default='IR_VDJ_2_')
            self.template_attributes.vdj_c_call_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "C gene"),
                default='c_call')
            self.template_attributes.vdj_d_call_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "D gene"),
                default='d_call')
            self.template_attributes.vdj_j_call_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "J gene"),
                default='j_call')
            self.template_attributes.vdj_v_call_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("VDJ", "V gene"),
                default='v_call')
            self.template_attributes.vdj_duplicate_count_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "number of duplicate UMIs (duplicate count)"),
                default='duplicate_count')
            self.template_attributes.vdj_junction_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "junction nt sequence"),
                default='junction')
            self.template_attributes.vdj_junction_aa_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "junction aa sequence"),
                default='junction_aa')
            self.template_attributes.vdj_locus_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "gene locus (e.g elements from IGH, IGK, IGL, TRA, TRB, "
                                                            "TRD, or TRG)"),
                default='locus')
            self.template_attributes.vdj_productive_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "'Is the gene productive?'"),
                default='productive')
            self.template_attributes.vdj_consensus_count_obs_key_suffix = sfaira_questionary(
                function='text',
                question=format_q_modality_obs_key("V(D)J", "number of reads contributing to consensus "
                                                            "(consensus count)"),
                default='consensus_count')

        print('[bold blue]Feature-wise meta data that are shared across the dataset.')
        print('[bold blue]The following meta data queries are on feature(gene)-wise. '
              'In some cases, these meta data need to be set per feature, '
              'skip here and annotate in the .yaml if that is the case. '
              'These can also later be modified manually in the .yaml which has the same effect as setting them here.')

        self.template_attributes.feature_reference = str(sfaira_questionary(
            function='text',
            question='Reference genome annotation release (e.g. Homo_sapiens.GRCh38.104)',
            default=''))
        self.template_attributes.feature_type = sfaira_questionary(
            function='text',
            question='Feature type of features annotated (either rna, protein, or peak)',
            default='rna')

        print('[bold blue]Feature-wise meta data.')
        print('[bold blue]The following meta data queries are on gene-resolution. '
              'These can also later be modified manually in the .yaml which has the same effect as setting them here.')
        feature_id = ""
        counter = 0
        while feature_id == "":
            if counter > 0:
                print('[bold red]You need to supply either a .var key for feature IDs or for symbols!')
            self.template_attributes.feature_id_var_key = sfaira_questionary(
                function='text',
                question='Key of feature ID (e.g. ENSEMBL ID) field in .var',
                default='')
            self.template_attributes.feature_symbol_var_key = sfaira_questionary(
                function='text',
                question='Key of feature symbol (e.g. gene names) field in .var',
                default='index')
            feature_id = self.template_attributes.feature_id_var_key \
                if self.template_attributes.feature_id_var_key != "" else \
                self.template_attributes.feature_symbol_var_key

        def insert_placeholder(id_str):
            return id_str if id_str != "" else "x"

        self.template_attributes.id_without_doi = \
            f'{insert_placeholder(clean_id_str(self.template_attributes.organism))}_' \
            f'{insert_placeholder(clean_id_str(self.template_attributes.organ))}_' \
            f'{insert_placeholder(clean_id_str(str(self.template_attributes.year)))}_' \
            f'{insert_placeholder(clean_id_str(self.template_attributes.assay_sc))}_' \
            f'{insert_placeholder(clean_id_str(first_author_lastname))}_001'
        self.template_attributes.id = f'{self.template_attributes.id_without_doi}_' \
                                      f'{self.template_attributes.doi_sfaira_repr}'

        print('[bold orange]Sfaira butler: "Up next:"')
        self.action_counter = 1
        print(f'[bold orange]               "{self.action_counter}) Proceed to modify the .yaml and .py files in '
              f'{self.path_loader}"')
        self.action_counter += 1
        self.check_datadir(path_data=path_data)
        if requires_annotate:
            print(f'[bold orange]               "{self.action_counter}) Proceed to phase 2: '
                  f'\'sfaira annotate-dataloader\'."')
            self.action_counter += 1
        else:
            print(f'[bold orange]               "{self.action_counter}) Proceed to phase 3: \'sfaira finalize-dataloader\', '
                  f'you can skip phase 2 (\'sfaira annotate-dataloader\')."')
            self.action_counter += 1

    @property
    def path_loader(self):
        return os.path.join(self.out_path, self.template_attributes.doi_sfaira_repr)

    def _template_attributes_to_dict(self) -> dict:
        """
        Create a dict from the our Template Structure dataclass
        :return: The dict containing all key-value pairs with non empty values
        """
        return {key: val for key, val in asdict(self.template_attributes).items() if val != ''}

    def _create_dataloader_template(self):
        template_path = os.path.join(self.TEMPLATES_PATH, self.template_attributes.dataloader_type)
        cookiecutter(f'{template_path}',
                     output_dir=self.out_path,
                     no_input=True,
                     overwrite_if_exists=True,
                     extra_context=self._template_attributes_to_dict())

    def check_datadir(self, path_data):
        path_data = os.path.join(path_data, self.template_attributes.doi_sfaira_repr)
        if not os.path.exists(path_data):
            print(f"[bold red]The unmodified downloaded data files were anticipated to lie in {path_data} "
                  f"but this path was not found. "
                  f"Create this directory and move the raw data files there.")
        else:
            print(f'[bold orange]{self.action_counter}) Make sure that the unmodified downloaded data files are in  '
                  f'{path_data}.')
            self.action_counter += 1
