import os
import sys
import anndata as ad
import pandas as pd


def get_mtx_paths(input_path: str) -> tuple[str, str, str]:
    if input_path.endswith(".mtx") or input_path.endswith(".mtx.gz"):
        mtx_path = input_path
        feature_path = mtx_path.replace(".mtx", ".features.tsv")
        barcode_path = mtx_path.replace(".mtx", ".barcodes.tsv")

        if not os.path.exists(feature_path):
            print(f"Error: Expected Feature file {feature_path} does not exist.")
            sys.exit(1)

        if not os.path.exists(barcode_path):
            print(f"Error: Expected Barcode file {barcode_path} does not exist.")
            sys.exit(1)
    else:
        if os.path.exists(os.path.join(input_path, "matrix.mtx")):
            mtx_path = os.path.join(input_path, "matrix.mtx")
            feature_path = os.path.join(input_path, "features.tsv")
            barcode_path = os.path.join(input_path, "barcodes.tsv")

        elif os.path.exists(os.path.join(input_path, "matrix.mtx.gz")):
            mtx_path = os.path.join(input_path, "matrix.mtx.gz")
            feature_path = os.path.join(input_path, "features.tsv.gz")
            barcode_path = os.path.join(input_path, "barcodes.tsv.gz")

        else:
            mtx_path_uncompressed = os.path.join(input_path, "matrix.mtx")
            mtx_path_compressed = os.path.join(input_path, "matrix.mtx.gz")
            sys.exit(
                f"Error: Expected Matrix file {mtx_path_uncompressed} or {mtx_path_compressed} does not exist."
            )

        if not os.path.exists(feature_path):
            print(f"Error: Expected Feature file {feature_path} does not exist.")
            sys.exit(1)

        if not os.path.exists(barcode_path):
            print(f"Error: Expected Barcode file {barcode_path} does not exist.")
            sys.exit(1)

    return (mtx_path, feature_path, barcode_path)


def convert_mtx_to_anndata(
    mtx_path: str, feature_path: str, barcode_path: str, dtype: str = "int32"
) -> ad.AnnData:
    # note: mtx is gene x cell
    adata = ad.io.read_mtx(mtx_path, dtype=dtype)

    features = pd.read_csv(
        feature_path, sep="\t", header=None, index_col=0
    ).index.astype(str)
    features.name = "feature"
    adata.obs_names = features  # type: ignore

    barcodes = pd.read_csv(
        barcode_path, sep="\t", header=None, index_col=0
    ).index.astype(str)
    barcodes.name = "barcode"
    adata.var_names = barcodes  # type: ignore

    # convert to cell x gene
    adata = adata.T

    # ensure CSR format
    adata.X = adata.X.tocsr()

    return adata
