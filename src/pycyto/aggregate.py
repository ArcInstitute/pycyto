import logging
import os
import re

import anndata as ad
import anndata.experimental as ade
import polars as pl

# Set up logger for aggregation
logger = logging.getLogger("pycyto.aggregate")


def lazy_load_adata(path: str) -> ad.AnnData:
    """Lazy load the anndata object from the path - notably load in the obs and var as well."""
    adata = ade.read_lazy(path)
    adata.obs = adata.obs.to_memory()
    adata.var = adata.var.to_memory()
    return adata


def _write_h5ad(
    adata: ad.AnnData,
    sample_outdir: str,
    sample: str,
    compress: bool = False,
    mode: str = "gex",
):
    adata.obs_names_make_unique()  # always make unique
    output_path = os.path.join(sample_outdir, f"{sample}_{mode}.h5ad")
    adata.write_h5ad(
        output_path,
        compression="gzip" if compress else None,
    )
    logger.debug(
        f"Successfully wrote {mode} h5ad file: {output_path} (shape: {adata.shape})"
    )


def _write_assignments_parquet(
    assignments: pl.DataFrame,
    sample_outdir: str,
    sample: str,
):
    output_path = os.path.join(sample_outdir, f"{sample}_assignments.parquet")
    assignments.write_parquet(
        output_path,
        compression="zstd",
    )
    logger.debug(
        f"Successfully wrote assignments parquet file: {output_path} (shape: {assignments.shape})"
    )


def _write_reads_parquet(
    reads_df: pl.DataFrame,
    sample_outdir: str,
    sample: str,
):
    output_path = os.path.join(sample_outdir, f"{sample}_reads.parquet")
    reads_df.write_parquet(
        output_path,
        compression="zstd",
    )
    logger.debug(
        f"Successfully wrote reads parquet file: {output_path} (shape: {reads_df.shape})"
    )


def _filter_crispr_adata_to_gex_barcodes(
    gex_adata: ad.AnnData,
    crispr_adata: ad.AnnData,
) -> ad.AnnData:
    """Filters the CRISPR data to only include barcodes present in the GEX data.

    Creates a dummy column on each that captures all unique information.

    # already annotated
    index: (cell_barcode + flex_barcode + lane_id)

    # to create
    dummy = index + sample + experiment
    """
    gex_adata.obs["dummy"] = (
        gex_adata.obs.index
        + "-"
        + gex_adata.obs["sample"]
        + "-"
        + gex_adata.obs["experiment"]
    )
    crispr_adata.obs["dummy"] = (
        crispr_adata.obs.index
        + "-"
        + crispr_adata.obs["sample"]
        + "-"
        + crispr_adata.obs["experiment"]
    )
    mask = crispr_adata.obs["dummy"].isin(gex_adata.obs["dummy"])
    gex_adata.obs.drop(columns=["dummy"], inplace=True)  # type: ignore
    crispr_adata.obs.drop(columns=["dummy"], inplace=True)  # type: ignore
    return crispr_adata[mask]


def _process_gex_crispr_set(
    gex_adata_list: list[ad.AnnData],
    crispr_adata_list: list[ad.AnnData],
    assignments_list: list[pl.DataFrame],
    reads_list: list[pl.DataFrame],
    sample_outdir: str,
    sample: str,
    compress: bool = False,
):
    logger.debug(f"Concatenating {len(gex_adata_list)} GEX anndata objects")
    gex_adata = ad.concat(gex_adata_list, join="outer")
    logger.debug(f"Final GEX data shape: {gex_adata.shape}")

    logger.debug(f"Concatenating {len(crispr_adata_list)} CRISPR anndata objects")
    crispr_adata = ad.concat(crispr_adata_list, join="outer")
    logger.debug(f"Final CRISPR data shape: {crispr_adata.shape}")

    logger.debug(f"Concatenating {len(assignments_list)} assignment dataframes")
    assignments = pl.concat(assignments_list, how="vertical_relaxed").unique()
    logger.debug(f"Final assignments shape: {assignments.shape}")

    logger.debug(f"Concatenating {len(reads_list)} reads dataframes")
    reads_df = pl.concat(reads_list, how="vertical_relaxed").unique()
    logger.debug(f"Final reads shape: {reads_df.shape}")

    if assignments["cell"].str.contains("CR").any():
        logger.debug("Detected CR barcodes, converting to BC format for matching")
        assignments = assignments.with_columns(
            match_barcode=pl.col("cell") + "-" + pl.col("lane_id").cast(pl.String)
        ).with_columns(pl.col("match_barcode").str.replace("CR", "BC"))
        reads_df = reads_df.with_columns(
            match_barcode=pl.col("cell_id") + "-" + pl.col("lane_id").cast(pl.String)
        ).with_columns(pl.col("match_barcode").str.replace("CR", "BC"))
        crispr_adata.obs.index = crispr_adata.obs.index.str.replace("CR", "BC")
    else:
        logger.debug("Using standard barcode format for matching")
        assignments = assignments.with_columns(
            match_barcode=pl.col("cell") + "-" + pl.col("lane_id").cast(pl.String)
        )
        reads_df = reads_df.with_columns(
            match_barcode=pl.col("cell_id") + "-" + pl.col("lane_id").cast(pl.String)
        )

    logger.debug("Merging assignment data with GEX observations")
    assignment_data = (
        assignments.select(["match_barcode", "assignment", "umis", "moi"])
        .to_pandas()
        .set_index("match_barcode")
    )
    logger.debug(f"Assignment data shape for merge: {assignment_data.shape}")

    logger.debug("Merging reads statistics with GEX observations")
    reads_pivot = (
        reads_df.select(["match_barcode", "mode", "n_reads", "n_umis"])
        .pivot(index="match_barcode", on="mode", values=["n_reads", "n_umis"])
        .fill_null(0)
        .to_pandas()
        .set_index("match_barcode")
    )
    logger.debug(f"Reads pivot shape for merge: {reads_pivot.shape}")

    gex_adata.obs = gex_adata.obs.merge(  # type: ignore
        assignment_data,
        left_index=True,
        right_index=True,
        how="left",
    ).merge(
        reads_pivot,
        left_index=True,
        right_index=True,
        how="left",
    )

    # Filter crispr adata to filtered barcodes
    logger.debug(
        f"Filtering CRISPR data to match GEX barcodes (GEX: {gex_adata.shape[0]} cells, CRISPR: {crispr_adata.shape[0]} cells)"
    )
    filt_crispr_adata = _filter_crispr_adata_to_gex_barcodes(
        gex_adata=gex_adata,
        crispr_adata=crispr_adata,
    )
    logger.debug(f"Filtered CRISPR data shape: {filt_crispr_adata.shape}")

    # Write both modes
    logger.info("Writing GEX anndata...")
    _write_h5ad(
        adata=gex_adata,
        sample_outdir=sample_outdir,
        sample=sample,
        compress=compress,
        mode="gex",
    )
    logger.info("Writing CRISPR anndata...")
    _write_h5ad(
        adata=filt_crispr_adata,
        sample_outdir=sample_outdir,
        sample=sample,
        compress=compress,
        mode="crispr",
    )
    logger.info("Writing assignments data...")
    _write_assignments_parquet(
        assignments=assignments,
        sample_outdir=sample_outdir,
        sample=sample,
    )
    logger.info("Writing reads data...")
    _write_reads_parquet(
        reads_df=reads_df,
        sample_outdir=sample_outdir,
        sample=sample,
    )


def _load_assignments_for_experiment_sample(
    root: str,
    crispr_bcs: list[str],
    lane_id: str,
    experiment: str,
    sample: str,
) -> list[pl.DataFrame]:
    assignments_list = []
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
                pl.lit(sample).alias("sample"),
                pl.lit(experiment).alias("experiment"),
                pl.lit(lane_id).alias("lane_id"),
                pl.lit(crispr_bc).alias("bc_idx"),
            )
            assignments_list.append(bc_assignments)
        else:
            logger.warning(
                f"Missing expected CRISPR assignments data for `{crispr_bc}` in {root} in path: {expected_crispr_assignments_path}"
            )
    return assignments_list


def _load_gex_anndata_for_experiment_sample(
    root: str,
    gex_bcs: list[str],
    lane_id: str,
    experiment: str,
    sample: str,
) -> list[ad.AnnData]:
    gex_adata_list = []
    expected_gex_adata_dir = os.path.join(root, "counts")
    for gex_bc in gex_bcs:
        expected_gex_adata_path = os.path.join(
            expected_gex_adata_dir, f"{gex_bc}.filt.h5ad"
        )
        if os.path.exists(expected_gex_adata_path):
            bc_adata = lazy_load_adata(expected_gex_adata_path)
            bc_adata.obs["sample"] = sample
            bc_adata.obs["experiment"] = experiment
            bc_adata.obs["lane_id"] = lane_id
            bc_adata.obs["bc_idx"] = gex_bc
            bc_adata.obs.index += "-" + bc_adata.obs["lane_id"].astype(str)
            gex_adata_list.append(bc_adata)
        else:
            logger.warning(
                f"Missing expected GEX data for `{gex_bc}` in {root} in path: {expected_gex_adata_path}"
            )
    return gex_adata_list


def _load_crispr_anndata_for_experiment_sample(
    root: str,
    crispr_bcs: list[str],
    lane_id: str,
    experiment: str,
    sample: str,
) -> list[ad.AnnData]:
    crispr_adata_list = []
    expected_crispr_adata_dir = os.path.join(root, "counts")
    for crispr_bc in crispr_bcs:
        expected_crispr_adata_path = os.path.join(
            expected_crispr_adata_dir, f"{crispr_bc}.h5ad"
        )
        if os.path.exists(expected_crispr_adata_path):
            bc_adata = lazy_load_adata(expected_crispr_adata_path)
            bc_adata.obs["sample"] = sample
            bc_adata.obs["experiment"] = experiment
            bc_adata.obs["lane_id"] = lane_id
            bc_adata.obs["bc_idx"] = crispr_bc
            bc_adata.obs.index += "-" + bc_adata.obs["lane_id"].astype(str)
            crispr_adata_list.append(bc_adata)
        else:
            logger.warning(
                f"Missing expected CRISPR data for `{crispr_bc}` in {root} in path: {expected_crispr_adata_path}"
            )
    return crispr_adata_list


def _load_reads_for_experiment_sample(
    root: str, bcs: list[str], lane_id: str, experiment: str, sample: str, mode: str
) -> list[pl.DataFrame]:
    reads_list = []
    for bc in bcs:
        expected_reads_path = os.path.join(
            root, "stats", "reads", f"{bc}.reads.tsv.zst"
        )
        if os.path.exists(expected_reads_path):
            reads_df = (
                pl.read_csv(expected_reads_path, separator="\t", has_header=True)
                .with_columns(
                    pl.lit(bc).alias("bc_idx"),
                    pl.lit(lane_id).alias("lane_id"),
                    pl.lit(experiment).alias("experiment"),
                    pl.lit(sample).alias("sample"),
                    pl.lit(mode).alias("mode"),
                )
                .with_columns(cell_id=pl.col("barcode") + "-" + pl.col("bc_idx"))
            )
            reads_list.append(reads_df)
        else:
            logger.warning(
                f"Missing expected reads data for `{bc}` in {root} in path: {expected_reads_path}"
            )
    return reads_list


def aggregate_data(
    config: pl.DataFrame, cyto_outdir: str, outdir: str, compress: bool = False
):
    logger.info(f"Starting aggregation workflow with output directory: {outdir}")
    logger.debug(f"Compression enabled: {compress}")

    unique_samples = config["sample"].unique().to_list()
    logger.info(
        f"Found {len(unique_samples)} unique samples to process: {unique_samples}"
    )

    for s in unique_samples:
        unique_experiments = (
            config.filter(pl.col("sample") == s)["experiment"].unique().to_list()
        )

        gex_adata_list = []
        crispr_adata_list = []
        assignments_list = []
        reads_list = []

        for e in unique_experiments:
            logger.info(f"Processing sample '{s}' experiment '{e}'...")

            subset = config.filter(pl.col("sample") == s, pl.col("experiment") == e)

            # identify base prefixes (experiment + mode without specific lanes)
            base_prefixes = subset["expected_prefix"].unique().to_list()
            # Create regex to match any lane number for these base prefixes
            base_pattern = "|".join([re.escape(prefix) for prefix in base_prefixes])
            prefix_regex = re.compile(rf"^({base_pattern})\d+.*")

            # determine data regex
            crispr_regex = re.compile(r".+_CRISPR_Lane.+")
            gex_regex = re.compile(r".+_GEX_Lane.+")
            lane_regex = re.compile(r"_Lane(\d+)")

            gex_bcs = (
                subset.filter(pl.col("mode") == "gex")
                .select("bc_component")
                .to_series()
                .unique()
                .sort()
                .to_list()
            )
            crispr_bcs = (
                subset.filter(pl.col("mode") == "crispr")
                .select("bc_component")
                .to_series()
                .unique()
                .sort()
                .to_list()
            )
            if len(gex_bcs) > 0:
                logger.info(f"Expecting GEX Barcodes: {gex_bcs}")
            if len(crispr_bcs) > 0:
                logger.info(f"Expecting CRISPR Barcodes: {crispr_bcs}")

            # Discover all directories that match our experiment/mode patterns
            matched_directories = []
            for root, _dirs, _files in os.walk(cyto_outdir, followlinks=True):
                basename = os.path.basename(root)
                if prefix_regex.search(basename):
                    matched_directories.append((root, basename))

            logger.debug(
                f"Found {len(matched_directories)} matching directories for experiment '{e}'"
            )

            # Process all discovered directories
            for root, basename in matched_directories:
                logger.info(f"Processing directory: {basename}")

                lane_regex_match = lane_regex.search(basename)
                if lane_regex_match:
                    lane_id = lane_regex_match.group(1)
                else:
                    raise ValueError(f"Invalid basename: {basename}")

                # process crispr data
                if crispr_regex.match(basename):
                    # Load in assignments
                    logger.debug(f"Loading CRISPR assignments from {basename}")
                    local_assignments_list = _load_assignments_for_experiment_sample(
                        root=root,
                        crispr_bcs=crispr_bcs,
                        lane_id=lane_id,
                        experiment=e,
                        sample=s,
                    )
                    assignments_list.extend(local_assignments_list)
                    logger.debug(
                        f"Loaded {len(local_assignments_list)} assignment files from {basename}"
                    )

                    # Load in crispr anndata
                    logger.debug(f"Loading CRISPR anndata from {basename}")
                    local_crispr_adata_list = (
                        _load_crispr_anndata_for_experiment_sample(
                            root=root,
                            crispr_bcs=crispr_bcs,
                            lane_id=lane_id,
                            experiment=e,
                            sample=s,
                        )
                    )
                    crispr_adata_list.extend(local_crispr_adata_list)
                    logger.debug(
                        f"Loaded {len(local_crispr_adata_list)} CRISPR anndata files from {basename}"
                    )

                    # process barcode-level read statistics
                    logger.debug(f"Loading CRISPR read statistics from {basename}")
                    local_reads_list = _load_reads_for_experiment_sample(
                        root=root,
                        bcs=crispr_bcs,
                        lane_id=lane_id,
                        experiment=e,
                        sample=s,
                        mode="crispr",
                    )
                    reads_list.extend(local_reads_list)

                # process gex data
                elif gex_regex.search(basename):
                    logger.debug(f"Loading GEX anndata from {basename}")
                    local_gex_list = _load_gex_anndata_for_experiment_sample(
                        root=root,
                        gex_bcs=gex_bcs,
                        lane_id=lane_id,
                        experiment=e,
                        sample=s,
                    )
                    gex_adata_list.extend(local_gex_list)
                    logger.debug(
                        f"Loaded {len(local_gex_list)} GEX anndata files from {basename}"
                    )

                    # process barcode-level read statistics
                    logger.debug(f"Loading GEX read statistics from {basename}")
                    local_reads_list = _load_reads_for_experiment_sample(
                        root=root,
                        bcs=gex_bcs,
                        lane_id=lane_id,
                        experiment=e,
                        sample=s,
                        mode="gex",
                    )
                    reads_list.extend(local_reads_list)

        sample_outdir = os.path.join(outdir, s)
        os.makedirs(sample_outdir, exist_ok=True)
        logger.debug(f"Created output directory: {sample_outdir}")

        # CRISPR + GEX case
        if len(gex_adata_list) > 0 and len(assignments_list) > 0:
            logger.info(f"Processing combined GEX + CRISPR data for sample '{s}'")
            _process_gex_crispr_set(
                gex_adata_list=gex_adata_list,
                crispr_adata_list=crispr_adata_list,
                assignments_list=assignments_list,
                reads_list=reads_list,
                sample_outdir=sample_outdir,
                sample=s,
                compress=compress,
            )

        elif len(gex_adata_list) > 0:
            logger.info(f"Processing GEX-only data for sample '{s}'")
            logger.debug(f"Concatenating {len(gex_adata_list)} GEX anndata objects")
            gex_adata = ad.concat(gex_adata_list, join="outer")
            logger.debug(f"Final GEX data shape: {gex_adata.shape}")

            logger.info("Writing GEX data...")
            _write_h5ad(
                adata=gex_adata,
                sample_outdir=sample_outdir,
                sample=s,
                compress=compress,
                mode="gex",
            )

            logger.info("Writing reads data...")
            reads_df = pl.concat(reads_list, how="vertical_relaxed").unique()
            logger.debug(f"Reads data shape: {reads_df.shape}")
            _write_reads_parquet(
                reads_df=reads_df,
                sample_outdir=sample_outdir,
                sample=s,
            )

        elif len(assignments_list) > 0:
            logger.info(f"Processing CRISPR-only data for sample '{s}'")

            logger.info("Writing assignments data...")
            assignments = pl.concat(assignments_list, how="vertical_relaxed").unique()
            logger.debug(f"Assignments data shape: {assignments.shape}")
            _write_assignments_parquet(
                assignments=assignments,
                sample_outdir=sample_outdir,
                sample=s,
            )

            logger.info("Writing reads data...")
            reads_df = pl.concat(reads_list, how="vertical_relaxed").unique()
            logger.debug(f"Reads data shape: {reads_df.shape}")
            _write_reads_parquet(
                reads_df=reads_df,
                sample_outdir=sample_outdir,
                sample=s,
            )

            logger.info("Writing CRISPR anndata...")
            crispr_adata = ad.concat(crispr_adata_list, join="outer")
            logger.debug(f"CRISPR data shape: {crispr_adata.shape}")
            _write_h5ad(
                adata=crispr_adata,
                sample_outdir=sample_outdir,
                sample=s,
                compress=compress,
                mode="crispr",
            )

        else:
            logger.warning(f"No data found to process for sample '{s}'")

    logger.info(
        f"Aggregation workflow completed successfully. Processed {len(unique_samples)} samples."
    )
