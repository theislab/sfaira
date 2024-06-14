import logging
import os
import sys

import click
import rich
import rich.logging
from rich import traceback
from rich import print

from sfaira.commands.annotate_dataloader import DataloaderAnnotater
from sfaira.commands.cache_control import CacheControl
from sfaira.commands.export_h5ad import H5adExport
from sfaira.commands.submit_pullrequest import PullRequestHandler
from sfaira.commands.test_dataloader import DataloaderTester
from sfaira.commands.utils import doi_lint
from sfaira.commands.validate_dataloader import DataloaderValidator
from sfaira.commands.validate_h5ad import H5adValidator
from sfaira.commands.test_dataloader import DataloaderTester
from sfaira.commands.submit_pullrequest import PullRequestHandler

from sfaira import __version__
from sfaira.commands.create_dataloader import DataloaderCreator
from sfaira.commands.upgrade import UpgradeCommand

WD = os.path.dirname(__file__)
log = logging.getLogger()

DEFAULT_DATA_PATH = "./"
PACKAGE_LOADER_PATH = "sfaira/data/dataloaders/loaders/"


def check_paths(pathlist):
    for p in pathlist:
        if p is None or not os.path.isdir(p):
            print(f"Error: Path {p} does not exist.")
            sys.exit()


def set_paths(loader=None, data=None, cache=None):
    env_loader_path = os.getenv('SFAIRA_LOADER_PATH')
    env_data_path = os.getenv('SFAIRA_DATA_PATH')
    env_cache_path = os.getenv('SFAIRA_CACHE_PATH')
    if env_loader_path:
        if not os.getenv('SFAIRA_DOCKER'):  # Skip print if running in sfaira docker where env variables are always set
            print('[bold blue]SFAIRA_LOADER_PATH environment variable detected. Ignoring --path-loader if supplied.')
        loader = env_loader_path
    if env_data_path:
        if not os.getenv('SFAIRA_DOCKER'):  # Skip print if running in sfaira docker where env variables are always set
            print('[bold blue]SFAIRA_DATA_PATH environment variable detected. Ignoring --path-data if supplied.')
        data = env_data_path
    if env_cache_path:
        if not os.getenv('SFAIRA_DOCKER'):  # Skip print if running in sfaira docker where env variables are always set
            print('[bold blue]SFAIRA_CACHE_PATH environment variable detected. Ignoring --path-cache if supplied.')
        cache = env_cache_path
    check_paths([loader, data])
    return loader, data, cache


def main():
    traceback.install(width=200, word_wrap=True)
    print(r"""[bold blue]
███████ ███████  █████  ██ ██████   █████ 
██      ██      ██   ██ ██ ██   ██ ██   ██ 
███████ █████   ███████ ██ ██████  ███████ 
     ██ ██      ██   ██ ██ ██   ██ ██   ██ 
███████ ██      ██   ██ ██ ██   ██ ██   ██ 
                                           
        """)

    print('[bold blue]Run [green]sfaira --help [blue]for an overview of all commands\n')

    # Is the latest sfaira version installed? Upgrade if not!
    if not UpgradeCommand.check_sfaira_latest():
        print('[bold blue]Run [green]sfaira upgrade [blue]to get the latest version.')
    sfaira_cli()


@click.group()
@click.version_option(__version__, message=click.style(f'sfaira Version: {__version__}', fg='blue'))
@click.pass_context
def sfaira_cli(ctx):
    """
    Create and manage sfaira dataloaders.
    """
    pass


@sfaira_cli.command()
@click.option('--reload',
              is_flag=True,
              help='Downloads new ontology versions into cache after clearing it.')
def clear_cache(reload) -> None:
    """Clears sfaira cache, including ontology and genome cache."""
    ctrl = CacheControl()
    if reload:
        ctrl.clear(countdown=False)
        ctrl.reload()
    else:
        ctrl.clear(countdown=True)


@sfaira_cli.command()
@click.option('--path-data',
              default=DEFAULT_DATA_PATH,
              type=click.Path(exists=False),
              help='Absolute path of the desired location of the raw data directory.')
@click.option('--path-loader',
              default=PACKAGE_LOADER_PATH,
              type=click.Path(exists=False),
              help='Relative path from the current directory to the desired location of the data loader.')
def create_dataloader(path_data, path_loader) -> None:
    """
    Interactively create a new sfaira dataloader.
    """
    path_loader, path_data, _ = set_paths(loader=path_loader, data=path_data)
    dataloader_creator = DataloaderCreator(path_loader)
    dataloader_creator.create_dataloader(path_data)


@sfaira_cli.command()
@click.option('--doi', type=str, required=True, help="The doi of the paper that the data loader refers to.")
@click.option('--path-data',
              default=DEFAULT_DATA_PATH,
              type=click.Path(exists=False),
              help='Absolute path of the location of the raw data directory.')
@click.option('--path-loader',
              default=PACKAGE_LOADER_PATH,
              type=click.Path(exists=False),
              help='Relative path from the current directory to the location of the data loader.')
@click.option('--schema', type=str, default="sfaira", help="The curation schema to check meta data availability for.")
def annotate_dataloader(doi, path_data, path_loader, schema) -> None:
    """
    Annotates a dataloader.
    """
    path_loader, path_data, _ = set_paths(loader=path_loader, data=path_data)
    if doi_lint(doi):
        dataloader_validator = DataloaderValidator(path_loader, doi, schema=schema)
        dataloader_validator.validate()
        dataloader_annotater = DataloaderAnnotater()
        dataloader_annotater.annotate(path_loader, path_data, dataloader_validator.doi)
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


@sfaira_cli.command()
@click.option('--doi', type=str, required=True, help="The doi of the paper that the data loader refers to.")
@click.option('--path-loader',
              default=PACKAGE_LOADER_PATH,
              type=click.Path(exists=False),
              help='Relative path from the current directory to the desired location of the data loader.')
@click.option('--schema', type=str, default="sfaira", help="The curation schema to check meta data availability for.")
def validate_dataloader(doi, path_loader, schema) -> None:
    """
    Verifies the dataloader against sfaira's requirements.
    """
    path_loader, _, _ = set_paths(loader=path_loader)
    if doi_lint(doi):
        dataloader_validator = DataloaderValidator(path_loader=path_loader, doi=doi, schema=schema)
        dataloader_validator.validate()
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


def _full_test(path_loader, path_data, doi, schema, clean_tsvs, in_phase_3):
    if doi is None or doi_lint(doi):
        dataloader_validator = DataloaderValidator(path_loader=path_loader, doi=doi, schema=schema)
        dataloader_validator.validate()
        dataloader_tester = DataloaderTester(path_loader, path_data, doi)
        dataloader_tester.test_dataloader(clean_tsvs=clean_tsvs, in_phase_3=in_phase_3)
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


@sfaira_cli.command()
@click.option('--doi', type=str, required=True, help="The doi of the paper that the data loader refers to.")
@click.option('--path-data',
              default=DEFAULT_DATA_PATH,
              type=click.Path(exists=False),
              help='Absolute path of the location of the raw data directory.')
@click.option('--path-loader',
              default=PACKAGE_LOADER_PATH,
              type=click.Path(exists=False),
              help='Relative path from the current directory to the location of the data loader.')
@click.option('--schema', type=str, default="sfaira", help="The curation schema to check meta data availability for.")
def test_dataloader(doi, path_data, path_loader, schema) -> None:
    """
    Runs a full data loader test.
    """
    path_loader, path_data, _ = set_paths(loader=path_loader, data=path_data)
    _full_test(path_loader=path_loader, path_data=path_data, doi=doi, schema=schema, clean_tsvs=False, in_phase_3=False)


@sfaira_cli.command()
@click.option('--doi', type=str, required=True, help="The doi of the paper that the data loader refers to.")
@click.option('--path-data',
              default=DEFAULT_DATA_PATH,
              type=click.Path(exists=False),
              help='Absolute path of the location of the raw data directory.')
@click.option('--path-loader',
              default=PACKAGE_LOADER_PATH,
              type=click.Path(exists=False),
              help='Relative path from the current directory to the location of the data loader.')
@click.option('--schema', type=str, default="sfaira", help="The curation schema to check meta data availability for.")
def finalize_dataloader(doi, path_data, path_loader, schema) -> None:
    """
    Formats .tsvs and runs a full data loader test.
    """
    path_loader, path_data, _ = set_paths(loader=path_loader, data=path_data)
    _full_test(path_loader=path_loader, path_data=path_data, doi=doi, schema=schema, clean_tsvs=True, in_phase_3=True)


@sfaira_cli.command()
def publish_dataloader() -> None:
    """
    Interactively create a GitHub pull request for a newly created data loader.
    This only works when called in the sfaira CLI docker container.
    Runs a full data loader test before starting the pull request.
    """
    path_loader, _, _ = set_paths()
    pullrequest_handler = PullRequestHandler(path_loader)
    pullrequest_handler.submit_pr()


@sfaira_cli.command()
@click.option('--doi', type=str, required=True, help="The doi of the paper that the data loader refers to.")
@click.option('--schema', type=str, required=True, help="Schema to streamline to, e.g. 'cellxgene'")
@click.option('--path-out',
              type=click.Path(exists=True),
              help='Absolute path of the location of the streamlined output h5ads.')
@click.option('--path-data',
              default=DEFAULT_DATA_PATH,
              type=click.Path(exists=True),
              help='Absolute path of the location of the raw data directory.')
@click.option('--contributors',
              type=str,
              default=None,
              help='Contributors for this curated object (listed in cellxgene for example).')
@click.option('--del-obs',
              type=str,
              default=None,
              help='# separated string of columns to remove from .obs, do not add spaces.')
@click.option('--keep-obs',
              type=str,
              default=None,
              help='# separated string of columns to remove from .obs, do not add spaces.')
@click.option('--path-cache',
              type=click.Path(exists=True),
              default=None,
              help='The optional absolute path to cached data library maintained by sfaira. Using such a cache speeds '
                   'up loading in sequential runs but is not necessary.')
@click.option('--path-loader',
              default=PACKAGE_LOADER_PATH,
              type=click.Path(exists=False),
              help='Relative path from the current directory to the location of the data loader.')
@click.option('--title',
              type=str,
              default=None,
              help='Title of dataloader / publication')
def export_h5ad(del_obs, contributors, doi, schema, keep_obs, path_out, path_data, path_loader, path_cache, title) -> \
        None:
    """Creates a collection of streamlined h5ad object for a given DOI."""
    del_obs = del_obs if del_obs is None else del_obs.split("#")
    keep_obs = keep_obs if keep_obs is None else keep_obs.split("#")
    path_loader, path_data, _ = set_paths(loader=path_loader, data=path_data, cache=path_cache)
    h5ad_export = H5adExport(clean_obs=del_obs, contributors=contributors, doi=doi, keep_obs=keep_obs,
                             path_cache=path_cache, path_data=path_data, path_loader=path_loader, path_out=path_out,
                             schema=schema, title=title)
    h5ad_export.write()


@sfaira_cli.command()
@click.option('--h5ad', type=click.Path(exists=True))
@click.option('--schema', type=str, default=None)
def validate_h5ad(h5ad, schema) -> None:
    """Runs a component test on a streamlined h5ad object.

    h5ad is the absolute path of the .h5ad file to test.
    schema is the schema type ("cellxgene",) to test.
    """
    h5ad_validator = H5adValidator(h5ad, schema)
    h5ad_validator.validate()


if __name__ == "__main__":
    traceback.install()
    sys.exit(main())  # pragma: no cover
