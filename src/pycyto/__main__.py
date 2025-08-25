import typer
from typing_extensions import Annotated

app = typer.Typer()


@app.command()
def convert(
    path: Annotated[
        str,
        typer.Argument(help="Input file path (directory for MTX or direct MTX path)"),
    ],
    output: Annotated[str, typer.Argument(help="Output file path")],
    compress: Annotated[
        bool, typer.Option(help="Internally gzip compress the output h5ad (gzip)")
    ] = True,
    integer: Annotated[bool, typer.Option(help="Convert data to integer")] = False,
):
    """Convert an MTX directory into a sparse CSR h5ad file"""
    from .convert import convert_mtx_to_anndata, get_mtx_paths

    (mtx_path, feature_path, barcode_path) = get_mtx_paths(path)
    adata = convert_mtx_to_anndata(
        mtx_path, feature_path, barcode_path, dtype="int32" if integer else "float32"
    )
    adata.write(output, compression="gzip" if compress else None)


@app.command()
def pipeline(
    config_path: Annotated[
        str, typer.Argument(help="Path to pipeline configuration json")
    ],
    sequences_dir: Annotated[
        str,
        typer.Argument(
            help="Path to directory containing sequence files (compressed FASTQ or BQ/VBQ"
        ),
    ],
    force: Annotated[
        bool, typer.Option(help="Force overwrite of existing output files")
    ] = False,
    threads: Annotated[
        int, typer.Option(help="Number of threads to use for each cyto run")
    ] = 8,
):
    """Run a cyto pipeline over a collection of input files with sub-sample specification"""
    import subprocess
    import sys

    from .config import determine_cyto_runs, parse_config
    from .pipeline import initialize_pipeline

    sample_sheet = parse_config(config_path)
    cyto_runs = determine_cyto_runs(sample_sheet)
    commands = initialize_pipeline(
        cyto_runs, sequences_dir, force=force, threads=threads
    )
    for command in commands:
        print(f"{' '.join(command)}", file=sys.stderr)
        subprocess.run(command)


@app.command()
def aggregate(
    config_path: Annotated[
        str, typer.Argument(help="Path to aggregation configuration json")
    ],
    cyto_outdir: Annotated[
        str, typer.Argument(help="Path to directory containing cyto output files")
    ],
    outdir: Annotated[
        str, typer.Argument(help="Path to output directory to write aggregated files")
    ],
):
    """Aggregate cyto output files"""

    from .aggregate import aggregate_data
    from .config import parse_config

    config = parse_config(config_path)
    aggregate_data(config, cyto_outdir, outdir)


@app.command()
def version():
    from importlib.metadata import version

    version_number = version("pycyto")
    typer.echo(f"pycyto {version_number}")


def main():
    app()
