"""Command line interface (CLI) support for the aeolus entry-point."""
import click
from click_default_group import DefaultGroup

from ._version import version as __version__

__all__ = ["main"]
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(
    cls=DefaultGroup,
    default_if_no_args=True,
    invoke_without_command=True,
    context_settings=CONTEXT_SETTINGS,
)
@click.option(
    "-v",
    "--version",
    is_flag=True,
    help="Show the package version.",
)
def main(version: bool) -> None:
    """To get help for other commands, simply use "aeolus COMMAND --help"."""
    if version:
        click.echo("version ", nl=False)
        click.secho(f"{__version__}")
