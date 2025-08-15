from typing_extensions import Annotated
import typer

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
    from .convert import get_mtx_paths, convert_mtx_to_anndata

    (mtx_path, feature_path, barcode_path) = get_mtx_paths(path)
    adata = convert_mtx_to_anndata(
        mtx_path, feature_path, barcode_path, dtype="int32" if integer else "float32"
    )
    adata.write(output, compression="gzip" if compress else None)


@app.command()
def version():
    from importlib.metadata import version

    version_number = version("pycyto")
    typer.echo(f"pycyto {version_number}")


def main():
    app()
