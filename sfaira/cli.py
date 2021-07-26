import logging
import os
import sys
import re

import click
import rich
import rich.logging
from rich import traceback
from rich import print

from sfaira.commands.annotate_dataloader import DataloaderAnnotater
from sfaira.commands.test_dataloader import DataloaderTester

from sfaira.commands.validate_dataloader import DataloaderValidator

import sfaira
from sfaira.commands.create_dataloader import DataloaderCreator
from sfaira.commands.upgrade import UpgradeCommand

WD = os.path.dirname(__file__)
log = logging.getLogger()


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
@click.version_option(sfaira.__version__, message=click.style(f'sfaira Version: {sfaira.__version__}', fg='blue'))
@click.option('-v', '--verbose', is_flag=True, default=False, help='Enable verbose output (print debug statements).')
@click.option("-l", "--log-file", help="Save a verbose log to a file.")
@click.pass_context
def sfaira_cli(ctx, verbose, log_file):
    """
    Create and manage sfaira dataloaders.
    """
    # Set the base logger to output DEBUG
    log.setLevel(logging.DEBUG)

    # Set up logs to the console
    log.addHandler(
        rich.logging.RichHandler(
            level=logging.DEBUG if verbose else logging.INFO,
            console=rich.console.Console(file=sys.stderr),
            show_time=True,
            markup=True,
        )
    )

    # Set up logs to a file if we asked for one
    if log_file:
        log_fh = logging.FileHandler(log_file, encoding="utf-8")
        log_fh.setLevel(logging.DEBUG)
        log_fh.setFormatter(logging.Formatter("[%(asctime)s] %(name)-20s [%(levelname)-7s]  %(message)s"))
        log.addHandler(log_fh)


@sfaira_cli.command()
@click.option('--path-loader',
              default="sfaira/data/dataloaders/loaders/",
              type=click.Path(exists=True),
              help='Relative path from the current directory to the desired location of the dataloader.'
              )
@click.option('--path-data',
              default="sfaira/unit_tests/template_data/",
              type=click.Path(exists=False),
              help='Relative path from the current directory to the datafiles used by this dataloader.'
              )
@click.option('--doi', type=str, default=None, help="The doi of the paper you would like to create a dataloader for.")
def create_dataloader(path_loader, doi, path_data) -> None:
    """
    Interactively create a new sfaira dataloader.
    """
    if doi is None or re.match(r'\b10\.\d+/[\w.]+\b', doi):
        dataloader_creator = DataloaderCreator(path_loader, doi)
        dataloader_creator.create_dataloader()
        dataloader_creator.create_datadir(path_data)
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


@sfaira_cli.command()
@click.option('--path-loader',
              default="sfaira/data/dataloaders/loaders/",
              type=click.Path(exists=True),
              help='Relative path from the current directory to the desired location of the dataloader.'
              )
@click.option('--doi', type=str, default=None, help="The doi of the paper that the dataloader refers to.")
def validate_dataloader(path_loader, doi) -> None:
    """
    Verifies the dataloader against sfaira's requirements.

    PATH to the dataloader script.
    """
    if doi is None or re.match(r'\b10\.\d+/[\w.]+\b', doi):
        dataloader_validator = DataloaderValidator(path_loader, doi)
        dataloader_validator.validate()
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


@sfaira_cli.command()
@click.option('--path-loader',
              default="sfaira/data/dataloaders/loaders/",
              type=click.Path(exists=True),
              help='Relative path from the current directory to the location of the dataloader.'
              )
@click.option('--path-data',
              default="sfaira/unit_tests/template_data/",
              type=click.Path(exists=True),
              help='Relative path from the current directory to the datafiles used by this dataloader.'
              )
@click.option('--doi', type=str, default=None, help="The doi of the paper that the dataloader refers to.")
def annotate_dataloader(path_loader, path_data, doi) -> None:
    """
    Annotates a dataloader.

    PATH is the absolute path of the root of your sfaira clone.
    """
    if doi is None or re.match(r'\b10\.\d+/[\w.]+\b', doi):
        dataloader_validator = DataloaderValidator(path_loader, doi)
        dataloader_validator.validate()
        dataloader_annotater = DataloaderAnnotater()
        dataloader_annotater.annotate(path_loader, path_data, dataloader_validator.doi)
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


@sfaira_cli.command()
@click.option('--path-loader',
              default="sfaira/data/dataloaders/loaders/",
              type=click.Path(exists=True),
              help='Relative path from the current directory to the location of the dataloader.'
              )
@click.option('--path-data',
              default="sfaira/unit_tests/template_data/",
              type=click.Path(exists=True),
              help='Relative path from the current directory to the datafiles used by this dataloader.'
              )
@click.option('--doi', type=str, default=None, help="The doi of the paper that the dataloader refers to.")
def test_dataloader(path_loader, path_data, doi) -> None:
    """Runs a dataloader integration test.

    PATH is the absolute path of the root of your sfaira clone.
    """
    if doi is None or re.match(r'\b10\.\d+/[\w.]+\b', doi):
        dataloader_tester = DataloaderTester(path_loader, path_data, doi)
        dataloader_tester.test_dataloader()
    else:
        print('[bold red]The supplied DOI is malformed!')  # noqa: W605


if __name__ == "__main__":
    traceback.install()
    sys.exit(main())  # pragma: no cover
