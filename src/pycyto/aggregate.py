import os
import re
import sys

import anndata as ad
import polars as pl


def _write_gex_h5ad(
    adata: ad.AnnData,
    sample_outdir: str,
    experiment: str,
    sample: str,
    compress: bool = False,
):
    adata.write_h5ad(
        os.path.join(sample_outdir, f"{experiment}_{sample}_gex.h5ad"),
        compression="gzip" if compress else None,
    )


def _write_assignments_tsv(
    assignments: pl.DataFrame,
    sample_outdir: str,
    experiment: str,
    sample: str,
):
    assignments.write_csv(
        os.path.join(sample_outdir, f"{experiment}_{sample}_assignments.tsv"),
        separator="\t",
    )


def _process_gex_crispr_set(
    gex_adata_list: list[ad.AnnData],
    assignments_list: list[pl.DataFrame],
    crispr_bcs: list[str],
    sample_outdir: str,
    experiment: str,
    sample: str,
    compress: bool = False,
):
    gex_adata = ad.concat(gex_adata_list)
    assignments = pl.concat(assignments_list, how="vertical_relaxed").unique()

    # rename CR -> BC for intersection
    if any(["CR" in cr for cr in crispr_bcs]):
        assignments = assignments.with_columns(
            match_barcode=pl.col("cell_id") + "-" + pl.col("lane_id").cast(pl.String)
        ).with_columns(pl.col("match_barcode").str.replace("CR", "BC"))
    else:
        assignments = assignments.with_columns(
            match_barcode=pl.col("cell_id") + "-" + pl.col("lane_id").cast(pl.String)
        )

    gex_adata.obs = gex_adata.obs.merge(  # type: ignore
        assignments.select(["match_barcode", "assignment", "counts", "moi"])
        .to_pandas()
        .set_index("match_barcode"),
        left_index=True,
        right_index=True,
        how="left",
    )

    # Write both modes
    _write_gex_h5ad(
        adata=gex_adata,
        sample_outdir=sample_outdir,
        experiment=experiment,
        sample=sample,
        compress=compress,
    )
    _write_assignments_tsv(
        assignments=assignments,
        sample_outdir=sample_outdir,
        experiment=experiment,
        sample=sample,
    )


def aggregate_data(
    config: pl.DataFrame, cyto_outdir: str, outdir: str, compress: bool = False
):
    unique_experiments = config["experiment"].unique().to_list()
    for e in unique_experiments:
        unique_samples = (
            config.filter(pl.col("experiment") == e)["sample"].unique().to_list()
        )
        for s in unique_samples:
            print(f"Processing experiment {e} sample {s}...", file=sys.stderr)

            subset = config.filter(pl.col("sample") == s, pl.col("experiment") == e)

            # identify necessary prefixes for output
            unique_prefixes = subset["expected_prefix"].unique().to_list()
            prefix_regex = re.compile(rf"^({'|'.join(unique_prefixes)}).*")

            # determine data regex
            crispr_regex = re.compile(r".+_CRISPR_Lane.+")
            gex_regex = re.compile(r".+_GEX_Lane.+")
            lane_regex = re.compile(r"_Lane(\d+)")

            gex_bcs = (
                subset.filter(pl.col("mode") == "gex")
                .select("bc_component")
                .to_series()
                .unique()
                .to_list()
            )
            crispr_bcs = (
                subset.filter(pl.col("mode") == "crispr")
                .select("bc_component")
                .to_series()
                .unique()
                .to_list()
            )

            gex_adata = []
            assignments = []
            n_matches = 0
            for root, _dirs, _files in os.walk(cyto_outdir):
                basename = os.path.basename(root)
                if prefix_regex.search(basename):
                    print(f"Processing {basename}...", file=sys.stderr)

                    lane_regex_match = lane_regex.search(basename)
                    if lane_regex_match:
                        lane_id = lane_regex_match.group(1)
                    else:
                        raise ValueError(f"Invalid basename: {basename}")

                    # process crispr data
                    if crispr_regex.match(basename):
                        expected_crispr_assignments_dir = os.path.join(
                            root, "assignments"
                        )

                        for crispr_bc in crispr_bcs:
                            expected_crispr_assignments_path = os.path.join(
                                expected_crispr_assignments_dir,
                                f"{crispr_bc}.assignments.tsv",
                            )
                            if os.path.exists(expected_crispr_assignments_path):
                                bc_assignments = pl.read_csv(
                                    expected_crispr_assignments_path,
                                    separator="\t",
                                ).with_columns(
                                    pl.lit(crispr_bc).alias("bc_idx"),
                                    pl.lit(lane_id).alias("lane_id"),
                                    pl.lit(e).alias("experiment"),
                                )
                                assignments.append(bc_assignments)
                            else:
                                print(
                                    f"Missing expected CRISPR assignments data for `{basename}` in {root} in path: {expected_crispr_assignments_path}",
                                    file=sys.stderr,
                                )

                    # process gex data
                    elif gex_regex.search(basename):
                        expected_gex_adata_dir = os.path.join(root, "counts")
                        for gex_bc in gex_bcs:
                            expected_gex_adata_path = os.path.join(
                                expected_gex_adata_dir, f"{gex_bc}.filt.h5ad"
                            )

                            if os.path.exists(expected_gex_adata_path):
                                bc_adata = ad.read_h5ad(expected_gex_adata_path)
                                bc_adata.obs["bc_idx"] = gex_bc
                                bc_adata.obs["lane_id"] = lane_id
                                bc_adata.obs["experiment"] = e
                                bc_adata.obs.index += "-" + bc_adata.obs[
                                    "lane_id"
                                ].astype(str)
                                gex_adata.append(bc_adata)
                            else:
                                print(
                                    f"Missing expected GEX data for `{gex_bc}` in {root} in path: {expected_gex_adata_path}",
                                    file=sys.stderr,
                                )

                    n_matches += 1

                    # finish on expected number of prefixes (don't recurse too deeply)
                    if n_matches == len(unique_prefixes):
                        break

            sample_outdir = os.path.join(outdir, s)
            os.makedirs(sample_outdir, exist_ok=True)

            # CRISPR + GEX case
            if len(gex_adata) > 0 and len(assignments) > 0:
                _process_gex_crispr_set(
                    gex_adata_list=gex_adata,
                    assignments_list=assignments,
                    crispr_bcs=crispr_bcs,
                    sample_outdir=sample_outdir,
                    experiment=e,
                    sample=s,
                    compress=compress,
                )

            elif len(gex_adata) > 0:
                print("Writing GEX data...", file=sys.stderr)
                gex_adata = ad.concat(gex_adata)
                _write_gex_h5ad(
                    adata=gex_adata,
                    sample_outdir=sample_outdir,
                    experiment=e,
                    sample=s,
                    compress=compress,
                )

            elif len(assignments) > 0:
                print("Writing assignments...", file=sys.stderr)
                assignments = pl.concat(assignments, how="vertical_relaxed").unique()
                _write_assignments_tsv(
                    assignments=assignments,
                    sample_outdir=sample_outdir,
                    experiment=e,
                    sample=s,
                )

            print(f"Found {n_matches} matches")
