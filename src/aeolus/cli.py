"""Command line interface (CLI) support for the aeolus entry-point."""

import click
from click_default_group import DefaultGroup

from ._version import version as __version__

__all__ = ["main"]

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(
    cls=DefaultGroup,
    default="pp",
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


@main.command(no_args_is_help=True, context_settings={"show_default": True})
@click.option(
    "-m",
    "--model",
    default="lfric",
    type=click.Choice(["lfric", "um"]),
    help="Model",
)
@click.option("-i", "--inpdir", type=str, help="Input directory")
@click.option("-o", "--outdir", type=str, help="Output directory")
@click.option(
    "-l", "--label", type=str, required=True, help="Simulation label"
)
@click.option(
    "-p", "--planet", type=str, required=True, help="Planet configuration"
)
@click.option(
    "-c", "--c_num", type=str, default="C48", help="Cubed Sphere Mesh Number"
)
@click.option(
    "--ref_cube",
    type=str,
    default="air_potential_temperature",
    help="Reference cube, to which coordinates all data will be regridded",
)
@click.option(
    "--level_height",
    type=str,
    default="uniform",
    help="Type of the vertical level height coordinate",
)
@click.option(
    "--top_height",
    type=float,
    default=40e3,
    help="If level_height=uniform, set the model top height.",
)
@click.option(
    "--time_prof",
    default="inst",
    help="Type of the time output",
    show_default=True,
    type=click.Choice(["inst", "mean"]),
)
@click.option(
    "--file_chunk_size",
    type=int,
    default=1000,
    help="Number of files to process in one go.",
)
@click.option(
    "--nlat",
    type=int,
    default=90,
    help="Number of latitudes in the target grid.",
)
@click.option(
    "--nlon",
    type=int,
    default=144,
    help="Number of longitudes in the target grid.",
)
def pp(
    model: str,
    inpdir: str,
    outdir: str,
    label: str,
    planet: str,
    c_num: str,
    ref_cube: str,
    level_height: str,
    top_height: float,
    time_prof: str,
    file_chunk_size: int,
    nlat: int,
    nlon: int,
) -> None:
    """Post-process model output."""
    from .postprocess import process_lfric  # noqa

    click.secho(f"Processing {model} output...")
    if model == "lfric":
        fname_out = process_lfric(
            inpdir,
            outdir,
            label,
            planet,
            c_num,
            ref_cube,
            level_height,
            top_height,
            time_prof,
            file_chunk_size,
            nlat,
            nlon,
        )
        click.secho(f"Saved to {fname_out}")
