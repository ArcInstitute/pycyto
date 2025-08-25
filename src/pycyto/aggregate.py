import os
import re
import sys

import anndata as ad
import polars as pl


def aggregate_data(config: pl.DataFrame, cyto_outdir: str, outdir: str):
    unique_samples = config["sample"].unique().to_list()

    for s in unique_samples:
        print(f"Processing sample {s}...", file=sys.stderr)

        subset = config.filter(pl.col("sample") == s)

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
            .to_list()
        )
        crispr_bcs = (
            subset.filter(pl.col("mode") == "crispr")
            .select("bc_component")
            .to_series()
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
                    expected_crispr_assignments_dir = os.path.join(root, "assignments")

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

        if len(gex_adata) > 0:
            print("Writing GEX data...", file=sys.stderr)
            gex_adata = ad.concat(gex_adata)
            gex_adata.write_h5ad(
                os.path.join(sample_outdir, f"{s}_gex.h5ad"), compression="gzip"
            )

        if len(assignments) > 0:
            print("Writing assignments...", file=sys.stderr)
            assignments = pl.concat(assignments, how="vertical_relaxed")
            assignments.write_csv(
                os.path.join(sample_outdir, f"{s}_assignments.tsv"), separator="\t"
            )

        print(f"Found {n_matches} matches")

    pass
