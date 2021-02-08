import importlib
import logging
import os
import sys

import click
import rich
import rich.logging
from rich import traceback
from rich import print

import sfaira
from sfaira.commands.create_dataloader import DataloaderCreator
from sfaira.commands.templates.clean_dataloader import DataloaderCleaner
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
    Create state of the art projects from production ready templates.
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
def create_dataloader():
    dataloader_creator = DataloaderCreator()
    dataloader_creator.create_dataloader()


@sfaira_cli.command()
@click.argument('path', type=click.Path(exists=True))
def clean_dataloader(path):
    dataloader_cleaner = DataloaderCleaner(path)
    dataloader_cleaner.clean_dataloader()


@sfaira_cli.command()
def lint_dataloader():
    pass


@sfaira_cli.command()
def test_dataloader():
    pass


if __name__ == "__main__":
    traceback.install()
    sys.exit(main())  # pragma: no cover
