import logging
import os
import sys

import click
import rich
import rich.logging
from rich import traceback
from rich import print

from sfaira.commands.annotate_dataloader import DataloaderAnnotater
from sfaira.commands.test_dataloader import DataloaderTester

from sfaira.commands.clean_dataloader import DataloaderCleaner
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
def create_dataloader() -> None:
    """
    Interactively create a new sfaira dataloader.
    """
    dataloader_creator = DataloaderCreator()
    dataloader_creator.create_dataloader()


@sfaira_cli.command()
@click.argument('path', type=click.Path(exists=True))
def clean_dataloader(path) -> None:
    """
    Clean a just written sfaira dataloader to adhere to sfaira's standards.

    PATH to the dataloader script.
    """
    dataloader_cleaner = DataloaderCleaner(path)
    dataloader_cleaner.clean_dataloader()


@sfaira_cli.command()
@click.argument('path', type=click.Path(exists=True))
def validate_dataloader(path) -> None:
    """
    Verifies the dataloader against sfaira's requirements.

    PATH to the dataloader script.
    """
    dataloader_validator = DataloaderValidator(path)
    dataloader_validator.validate()


@sfaira_cli.command()
@click.argument('path', type=click.Path(exists=True))
@click.option('--doi', type=str, default=None)
@click.option('--test-data', type=click.Path(exists=True))
def annotate_dataloader(path, doi, test_data) -> None:
    """
    Annotates a dataloader.

    PATH is the absolute path of the root of your sfaira clone.
    """
    dataloader_annotater = DataloaderAnnotater()
    dataloader_annotater.annotate(path, doi, test_data)


@sfaira_cli.command()
@click.argument('path', type=click.Path(exists=True))
@click.option('--test-data', type=click.Path(exists=True))
@click.option('--doi', type=str, default=None)
def test_dataloader(path, test_data, doi) -> None:
    """Runs a dataloader integration test.

    PATH is the absolute path of the root of your sfaira clone.
    """
    dataloader_tester = DataloaderTester(path, test_data, doi)
    dataloader_tester.test_dataloader()


if __name__ == "__main__":
    traceback.install()
    sys.exit(main())  # pragma: no cover
