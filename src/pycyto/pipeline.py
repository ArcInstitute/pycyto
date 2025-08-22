import os
import re
import shutil
from importlib import resources

import polars as pl

CYTO_HOME = os.path.join(os.environ["HOME"], ".cyto")
WHITELIST_PATH = os.path.join(CYTO_HOME, "737K-fixed-rna-profiling.txt.gz")
PROBE_SETS = [
    "probe-barcodes-fixed-rna-profiling-rna.txt",
    "probe-barcodes-fixed-rna-profiling-crispr.txt",
    "probe-barcodes-fixed-rna-profiling-ab.txt",
]
PROBE_SET_PATHS = [os.path.join(CYTO_HOME, probe_set) for probe_set in PROBE_SETS]
PROBE_SET_MAP = {"BC": 0, "CR": 1, "AB": 2}


def _init_10x_whitelist():
    if not os.path.exists(WHITELIST_PATH):
        path = resources.files("pycyto.vendor.barcodes").joinpath(
            "737K-fixed-rna-profiling.txt.gz"
        )
        assert path.is_file(), f"File not found: {path}"
        shutil.copyfile(path, WHITELIST_PATH)  # type: ignore


def _init_probesets():
    for probe_set, expected_path in zip(PROBE_SETS, PROBE_SET_PATHS):
        if not os.path.exists(expected_path):
            path = resources.files("pycyto.vendor.probe-sets").joinpath(probe_set)
            assert path.is_file(), f"File not found: {path}"
            shutil.copyfile(path, expected_path)  # type: ignore


def _build_homedir():
    if not os.path.exists(CYTO_HOME):
        os.makedirs(CYTO_HOME)
    _init_10x_whitelist()
    _init_probesets()


def _assert_interleaving(sequence_subset: list[str]):
    def _is_fastx(filepath):
        return re.match(r"_R[12]*.f*q(.gz|.zst)?$", filepath)

    fastq_subset = [f for f in sequence_subset if _is_fastx(f)]
    for a, b in zip(fastq_subset[::2], fastq_subset[1::2]):
        if "_R1" in a and "_R2" in b:
            continue
        raise ValueError(
            f"Fastq files are not interleaved R1,R2: {a}, {b}. Problem in file name encoding."
        )


def _identify_files(dir: str) -> list[str]:
    # Regex for _R1/_R2 files: .fastq/.fq with optional .gz/.zst compression
    fastq_pattern = re.compile(r".*_R[12]\.(fastq|fq)(\.gz|\.zst)?$")

    # Regex for .bq/.vbq files
    bq_pattern = re.compile(r".*\.(bq|vbq)$")

    matched_files = []

    for root, dirs, files in os.walk(dir):
        for file in files:
            if fastq_pattern.match(file) or bq_pattern.match(file):
                matched_files.append(os.path.join(root, file))

    return matched_files


def initialize_pipeline(
    cyto_runs: pl.DataFrame, sequences_dir: str, force: bool = False, threads: int = 8
) -> list[list[str]]:
    """Initialize the pipeline runner and return all commands.

    Args:
        cyto_runs (pl.DataFrame): DataFrame containing the cyto runs.
        sequences_dir (str): Directory containing the sequences.

    Returns:
        list[list[str]]: List of commands to run.
    """
    # Build the home directory if it doesn't exist already
    _build_homedir()

    files = _identify_files(sequences_dir)

    outdir = os.path.join(os.getcwd(), "pycyto_out")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        if force:
            shutil.rmtree(outdir)
            os.makedirs(outdir)
        else:
            raise ValueError(
                f"Output directory {outdir} already exists. Use --force to overwrite."
            )

    all_commands = []
    for entry in cyto_runs.to_dicts():
        probe_set_path = PROBE_SET_PATHS[PROBE_SET_MAP[entry["probe_set"]]]

        sequence_subset = [
            f for f in files if os.path.basename(f).startswith(entry["expected_prefix"])
        ]
        _assert_interleaving(sequence_subset)

        if len(sequence_subset) == 0:
            raise ValueError(f"No files found for entry {entry}")

        subdir = os.path.join(outdir, entry["expected_prefix"])
        command = [
            "cyto",
            "workflow",
            entry["mode"],
            "-T",
            str(threads),
            "-w",
            WHITELIST_PATH,
            "-p",
            probe_set_path,
            "-c",
            entry["feature_path"],
            "-o",
            subdir,
        ]
        if force:
            command.append("--force")
        if entry["mode"] == "crispr" and entry["probe_set"] == "CR":
            command.extend(
                [
                    "--offset",
                    "8",
                ]
            )
        command.extend(sequence_subset)

        all_commands.append(command)

    return all_commands
