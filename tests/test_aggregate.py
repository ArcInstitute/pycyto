"""Tests for the aggregate module.

Covers the cross-experiment same-sample fix for issue #46 plus regression
coverage of the existing single-experiment paths and CR/BC edge cases.

The fixtures build minimal cyto-output trees from scratch (no shared
infrastructure exists). Each test goes through `parse_config` -> `process_sample`
so that the directory-discovery and merge paths exercised mirror production.
"""

import json
import os
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp

# Add src to path for imports (mirror the pattern used by other tests).
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from pycyto.aggregate import process_sample
from pycyto.config import parse_config


def _write_mock_h5ad(path: Path, cell_sequences: list[str], flex_bc: str) -> None:
    """Write a minimal h5ad whose obs.index matches cyto's `{seq}-{flex_bc}` format."""
    obs_index = pd.Index([f"{seq}-{flex_bc}" for seq in cell_sequences])
    adata = ad.AnnData(
        X=sp.csr_matrix(np.zeros((len(cell_sequences), 1), dtype=np.float32)),
        obs=pd.DataFrame(index=obs_index),
        var=pd.DataFrame(index=pd.Index(["gene_1"])),
    )
    adata.write_h5ad(str(path))  # ty: ignore[invalid-argument-type]


def _write_mock_reads(
    path: Path,
    cell_sequences: list[str],
    n_reads: int,
    n_umis: int,
) -> None:
    """Write a mock reads TSV (zstd compressed).

    The `barcode` column is the raw cell sequence only; the loader appends the
    `bc_idx` (flex barcode) downstream to construct `cell_id`.
    """
    pl.DataFrame(
        {
            "barcode": cell_sequences,
            "n_reads": [n_reads] * len(cell_sequences),
            "n_umis": [n_umis] * len(cell_sequences),
        }
    ).write_csv(path, separator="\t", compression="zstd", include_header=True)


def _write_mock_assignments(
    path: Path,
    cell_sequences: list[str],
    flex_bc: str,
    assignment: str = "guide_1",
) -> None:
    """Write a mock assignments TSV.

    The `cell` column carries the flex-barcode suffix (matches what cyto emits
    and what `_process_gex_crispr_set` reads via `pl.col("cell").str.contains("CR")`).
    """
    cells = [f"{seq}-{flex_bc}" for seq in cell_sequences]
    pl.DataFrame(
        {
            "cell": cells,
            "assignment": [assignment] * len(cells),
            "umis": [5] * len(cells),
            "moi": [1] * len(cells),
        }
    ).write_csv(path, separator="\t", include_header=True)


def _write_probe_stub(path: Path) -> None:
    path.write_text("gene\tsequence\ngene_1\tATCG\n")


def _build_cyto_dir(
    root: Path,
    experiment: str,
    mode: str,
    lane: int,
    bc: str,
    cell_sequences: list[str],
    n_reads: int,
    n_umis: int,
    assignment: str = "guide_1",
) -> Path:
    """Build one cyto output directory for either GEX or CRISPR mode.

    Always creates `counts/` and `stats/reads/`. Only creates `assignments/`
    for the CRISPR mode.
    """
    mode_upper = mode.upper()
    dir_path = root / f"{experiment}_{mode_upper}_Lane{lane}"
    (dir_path / "counts").mkdir(parents=True, exist_ok=True)
    (dir_path / "stats" / "reads").mkdir(parents=True, exist_ok=True)

    suffix = ".filt.h5ad" if mode == "gex" else ".h5ad"
    _write_mock_h5ad(dir_path / "counts" / f"{bc}{suffix}", cell_sequences, bc)
    _write_mock_reads(
        dir_path / "stats" / "reads" / f"{bc}.reads.tsv.zst",
        cell_sequences,
        n_reads,
        n_umis,
    )

    if mode == "crispr":
        (dir_path / "assignments").mkdir(parents=True, exist_ok=True)
        _write_mock_assignments(
            dir_path / "assignments" / f"{bc}.assignments.tsv",
            cell_sequences,
            bc,
            assignment,
        )

    return dir_path


def _write_and_parse_config(
    tmp_path: Path,
    libraries: dict[str, str],
    samples: list[dict],
) -> pl.DataFrame:
    """Write a config JSON to tmp_path, materialize probe stubs, parse.

    `libraries` values must be ABSOLUTE paths so that `_pull_feature_path`'s
    `os.path.exists` check resolves regardless of test CWD.
    """
    for lib_path in libraries.values():
        _write_probe_stub(Path(lib_path))

    config_path = tmp_path / "config.json"
    with open(config_path, "w") as f:
        json.dump({"libraries": libraries, "samples": samples}, f)

    return parse_config(str(config_path))


# 16-character cell sequences. The first three overlap across experiments to
# create the (cell_barcode, flex_barcode, lane_id) collision exercised by #46.
_OVERLAP_SEQS = [
    "AAAACCCCGGGGTTTT",
    "ACGTACGTACGTACGT",
    "TTTTGGGGCCCCAAAA",
]
_EXP_A_ONLY = ["AAAATTTTGGGGCCCC", "CCCCAAAATTTTGGGG"]
_EXP_B_ONLY = ["GGGGCCCCAAAATTTT", "TTTTAAAACCCCGGGG"]


class TestMultiExperimentSameSample:
    """Reproduces issue #46: same `sample` under two experiments with overlapping cell barcodes."""

    def _make_two_experiment_fixture(
        self,
        tmp_path: Path,
        n_umis_a: int = 26,
        n_umis_b: int = 3,
    ) -> tuple[Path, Path, pl.DataFrame]:
        """Build two cyto trees (exp_A and exp_B) sharing sample 'sample01' with BC001+CR001 pairing."""
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        seqs_a = _OVERLAP_SEQS + _EXP_A_ONLY
        seqs_b = _OVERLAP_SEQS + _EXP_B_ONLY

        _build_cyto_dir(cyto_outdir, "exp_A", "gex", 1, "BC001", seqs_a, 100, n_umis_a)
        _build_cyto_dir(cyto_outdir, "exp_A", "crispr", 1, "CR001", seqs_a, 50, 7)
        _build_cyto_dir(cyto_outdir, "exp_B", "gex", 1, "BC001", seqs_b, 100, n_umis_b)
        _build_cyto_dir(cyto_outdir, "exp_B", "crispr", 1, "CR001", seqs_b, 50, 7)

        libraries = {
            "GEX": str(tmp_path / "gex_probes.tsv"),
            "GUIDES": str(tmp_path / "guides.tsv"),
        }
        samples = [
            {
                "experiment": "exp_A",
                "sample": "sample01",
                "mode": "gex+crispr",
                "features": "GEX+GUIDES",
                "barcodes": "BC001+CR001",
            },
            {
                "experiment": "exp_B",
                "sample": "sample01",
                "mode": "gex+crispr",
                "features": "GEX+GUIDES",
                "barcodes": "BC001+CR001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        return cyto_outdir, agg_outdir, config

    def test_gex_crispr_two_experiments_does_not_crash(self, tmp_path):
        """Reproduces issue #46. Fails on pre-fix code with ValueError from aggregate.py:236."""
        cyto_outdir, agg_outdir, config = self._make_two_experiment_fixture(tmp_path)

        # Must not raise ValueError: shape is inconsistent with obs
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        sample_dir = agg_outdir / "sample01"
        assert (sample_dir / "sample01_gex.h5ad").exists()
        assert (sample_dir / "sample01_crispr.h5ad").exists()
        assert (sample_dir / "sample01_assignments.parquet").exists()
        assert (sample_dir / "sample01_reads.parquet").exists()

    def test_gex_crispr_two_experiments_correct_cell_count(self, tmp_path):
        """5 cells per experiment with 3 overlapping sequences -> 10 distinct (seq, experiment) rows post-fix."""
        cyto_outdir, agg_outdir, config = self._make_two_experiment_fixture(tmp_path)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        assert gex.shape[0] == 10
        assert gex.obs["experiment"].value_counts().to_dict() == {
            "exp_A": 5,
            "exp_B": 5,
        }
        assert gex.obs.index.is_unique

        # GEX output is always BC-prefixed (CR appears nowhere in GEX output).
        pattern = re.compile(r"^[ACGTN]+-BC\d{3}-\d+-(exp_A|exp_B)$")
        for idx in gex.obs.index:
            assert pattern.match(idx) is not None, (
                f"obs.index entry {idx!r} does not match expected format"
            )

    def test_gex_crispr_two_experiments_correct_attribution(self, tmp_path):
        """Distinct n_umis per experiment must attribute to the correct rows after the merge."""
        cyto_outdir, agg_outdir, config = self._make_two_experiment_fixture(
            tmp_path, n_umis_a=26, n_umis_b=3
        )
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        exp_a_rows = gex.obs[gex.obs["experiment"] == "exp_A"]
        exp_b_rows = gex.obs[gex.obs["experiment"] == "exp_B"]

        assert (exp_a_rows["n_umis_gex"] == 26).all()
        assert (exp_b_rows["n_umis_gex"] == 3).all()
        assert (exp_a_rows["n_reads_gex"] == 100).all()
        assert (exp_b_rows["n_reads_gex"] == 100).all()
        assert (exp_a_rows["assignment"] == "guide_1").all()
        assert (exp_b_rows["assignment"] == "guide_1").all()

        # match_barcode column in both parquets carries experiment and has been CR->BC converted.
        assignments = pl.read_parquet(
            agg_outdir / "sample01" / "sample01_assignments.parquet"
        )
        reads = pl.read_parquet(agg_outdir / "sample01" / "sample01_reads.parquet")
        mb_pattern = re.compile(r"^[ACGTN]+-BC\d{3}-\d+-(exp_A|exp_B)$")
        for mb in assignments["match_barcode"].to_list():
            assert mb_pattern.match(mb) is not None, (
                f"assignments match_barcode {mb!r} unexpected"
            )
        for mb in reads["match_barcode"].to_list():
            assert mb_pattern.match(mb) is not None, (
                f"reads match_barcode {mb!r} unexpected"
            )

        # Reads parquet rows are correctly attributed per (mode, experiment) on the reads side.
        gex_a_reads = reads.filter(
            (pl.col("mode") == "gex") & (pl.col("experiment") == "exp_A")
        )
        gex_b_reads = reads.filter(
            (pl.col("mode") == "gex") & (pl.col("experiment") == "exp_B")
        )
        assert gex_a_reads["n_umis"].to_list() == [26] * 5
        assert gex_b_reads["n_umis"].to_list() == [3] * 5
        assert gex_a_reads["n_reads"].to_list() == [100] * 5
        assert gex_b_reads["n_reads"].to_list() == [100] * 5


class TestSingleExperimentRegression:
    """Regression coverage: the existing single-experiment path still works end-to-end."""

    def test_single_experiment_single_barcode_gex_crispr(self, tmp_path):
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        seqs = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"]
        _build_cyto_dir(cyto_outdir, "exp_only", "gex", 1, "BC001", seqs, 100, 26)
        _build_cyto_dir(cyto_outdir, "exp_only", "crispr", 1, "CR001", seqs, 50, 7)

        libraries = {
            "GEX": str(tmp_path / "gex_probes.tsv"),
            "GUIDES": str(tmp_path / "guides.tsv"),
        }
        samples = [
            {
                "experiment": "exp_only",
                "sample": "sample01",
                "mode": "gex+crispr",
                "features": "GEX+GUIDES",
                "barcodes": "BC001+CR001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        assert gex.shape[0] == 3
        assert (gex.obs["n_reads_gex"] == 100).all()
        assert (gex.obs["n_umis_gex"] == 26).all()
        assert (gex.obs["assignment"] == "guide_1").all()

        # End-to-end format lock-in.
        pattern = re.compile(r"^[ACGTN]+-BC001-1-exp_only$")
        for idx in gex.obs.index:
            assert pattern.match(idx) is not None, (
                f"obs.index entry {idx!r} does not match expected format"
            )


class TestGEXOnlyMultiExperiment:
    """The experiment suffix lands in obs.index for gex-only outputs too."""

    def test_gex_only_two_experiments_unique_obs_index(self, tmp_path):
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        # 3 cells per experiment, 2 overlapping sequences across experiments.
        seqs_a = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "AAAATTTTGGGGCCCC"]
        seqs_b = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "GGGGCCCCAAAATTTT"]

        _build_cyto_dir(cyto_outdir, "exp_A", "gex", 1, "BC001", seqs_a, 100, 26)
        _build_cyto_dir(cyto_outdir, "exp_B", "gex", 1, "BC001", seqs_b, 100, 3)

        libraries = {"GEX": str(tmp_path / "gex_probes.tsv")}
        samples = [
            {
                "experiment": "exp_A",
                "sample": "sample01",
                "mode": "gex",
                "features": "GEX",
                "barcodes": "BC001",
            },
            {
                "experiment": "exp_B",
                "sample": "sample01",
                "mode": "gex",
                "features": "GEX",
                "barcodes": "BC001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        # Pre-fix: obs_names_make_unique() silently renames the 2 collisions ->
        # shape would still be 6, but indices would have synthetic "-1" suffixes.
        # Post-fix: experiment suffix disambiguates -> all 6 indices unique and
        # the experiment column is preserved.
        assert gex.shape[0] == 6
        assert gex.obs.index.is_unique
        for idx in gex.obs.index:
            assert idx.endswith("-exp_A") or idx.endswith("-exp_B"), (
                f"obs.index entry {idx!r} missing experiment suffix"
            )
        assert gex.obs["experiment"].value_counts().to_dict() == {
            "exp_A": 3,
            "exp_B": 3,
        }


class TestCRISPROnlyMultiExperiment:
    """Symmetric to TestGEXOnlyMultiExperiment for the CRISPR-only single-modal path.

    Note: CRISPR-only never enters `_process_gex_crispr_set`, so the CR->BC
    conversion at line 165 does NOT run for this path. obs.index retains the
    raw `CR` prefix from the input file.
    """

    def test_crispr_only_two_experiments_unique_obs_index(self, tmp_path):
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        seqs_a = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "AAAATTTTGGGGCCCC"]
        seqs_b = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "GGGGCCCCAAAATTTT"]

        _build_cyto_dir(cyto_outdir, "exp_A", "crispr", 1, "CR001", seqs_a, 50, 7)
        _build_cyto_dir(cyto_outdir, "exp_B", "crispr", 1, "CR001", seqs_b, 50, 7)

        libraries = {"GUIDES": str(tmp_path / "guides.tsv")}
        samples = [
            {
                "experiment": "exp_A",
                "sample": "sample01",
                "mode": "crispr",
                "features": "GUIDES",
                "barcodes": "CR001",
            },
            {
                "experiment": "exp_B",
                "sample": "sample01",
                "mode": "crispr",
                "features": "GUIDES",
                "barcodes": "CR001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        crispr = ad.read_h5ad(agg_outdir / "sample01" / "sample01_crispr.h5ad")
        assert crispr.shape[0] == 6
        assert crispr.obs.index.is_unique
        # CRISPR-only path: raw CR prefix is retained.
        pattern = re.compile(r"^[ACGTN]+-CR\d{3}-\d+-(exp_A|exp_B)$")
        for idx in crispr.obs.index:
            assert pattern.match(idx) is not None, (
                f"obs.index entry {idx!r} does not match expected format"
            )
        assert crispr.obs["experiment"].value_counts().to_dict() == {
            "exp_A": 3,
            "exp_B": 3,
        }


class TestCREdgeCases:
    def test_cr_prefix_crispr_positive_conversion_with_experiment(self, tmp_path):
        """Positive path: CR->BC conversion succeeds and the experiment suffix survives."""
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        seqs = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT"]
        _build_cyto_dir(cyto_outdir, "exp_simple", "gex", 1, "BC001", seqs, 100, 26)
        _build_cyto_dir(cyto_outdir, "exp_simple", "crispr", 1, "CR001", seqs, 50, 7)

        libraries = {
            "GEX": str(tmp_path / "gex_probes.tsv"),
            "GUIDES": str(tmp_path / "guides.tsv"),
        }
        samples = [
            {
                "experiment": "exp_simple",
                "sample": "sample01",
                "mode": "gex+crispr",
                "features": "GEX+GUIDES",
                "barcodes": "BC001+CR001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        crispr = ad.read_h5ad(agg_outdir / "sample01" / "sample01_crispr.h5ad")

        gex_pattern = re.compile(r"^[ACGTN]+-BC001-1-exp_simple$")
        crispr_pattern = re.compile(r"^[ACGTN]+-BC001-1-exp_simple$")
        for idx in gex.obs.index:
            assert gex_pattern.match(idx) is not None
        for idx in crispr.obs.index:
            assert crispr_pattern.match(idx) is not None
        # Attribution proves the polars match_barcode join worked end-to-end.
        assert (gex.obs["n_umis_gex"] == 26).all()
        assert (gex.obs["n_umis_gex"] != 0).all()

    def test_cr_prefix_crispr_does_not_corrupt_experiment_named_CRISPR(self, tmp_path):
        """experiment="CRISPR_screen" must not be rewritten to "BCISPR_screen" by the CR->BC pass."""
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        seqs = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT"]
        _build_cyto_dir(cyto_outdir, "CRISPR_screen", "gex", 1, "BC001", seqs, 100, 26)
        _build_cyto_dir(cyto_outdir, "CRISPR_screen", "crispr", 1, "CR001", seqs, 50, 7)

        libraries = {
            "GEX": str(tmp_path / "gex_probes.tsv"),
            "GUIDES": str(tmp_path / "guides.tsv"),
        }
        samples = [
            {
                "experiment": "CRISPR_screen",
                "sample": "sample01",
                "mode": "gex+crispr",
                "features": "GEX+GUIDES",
                "barcodes": "BC001+CR001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        pattern = re.compile(r"^[ACGTN]+-BC001-1-CRISPR_screen$")
        for idx in gex.obs.index:
            assert pattern.match(idx) is not None, (
                f"obs.index entry {idx!r} -- experiment name corrupted by CR->BC pass"
            )
        assert (gex.obs["experiment"] == "CRISPR_screen").all()

        # Polars-side correctness: match_barcode in the parquets retains the literal experiment name.
        assignments = pl.read_parquet(
            agg_outdir / "sample01" / "sample01_assignments.parquet"
        )
        for mb in assignments["match_barcode"].to_list():
            assert mb.endswith("-CRISPR_screen"), (
                f"assignments match_barcode {mb!r} corrupted"
            )
            assert "BCISPR" not in mb

    def test_bc_prefix_crispr_skips_conversion(self, tmp_path):
        """CRISPR using a BC-prefix barcode takes the standard (non-CR-conversion) branch.

        Realistic configuration: GEX and CRISPR share the same flex barcode (`BC001`)
        so the merge keys align across modes. Assignments `cell` column carries `BC001`
        with no `CR` substring, so `_process_gex_crispr_set` takes the standard branch
        and obs.index retains the BC prefix end-to-end.
        """
        cyto_outdir = tmp_path / "cyto_outdir"
        cyto_outdir.mkdir()
        agg_outdir = tmp_path / "agg_outdir"
        agg_outdir.mkdir()

        seqs = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT"]
        _build_cyto_dir(cyto_outdir, "exp_bc", "gex", 1, "BC001", seqs, 100, 26)
        _build_cyto_dir(cyto_outdir, "exp_bc", "crispr", 1, "BC001", seqs, 50, 7)

        libraries = {
            "GEX": str(tmp_path / "gex_probes.tsv"),
            "GUIDES": str(tmp_path / "guides.tsv"),
        }
        samples = [
            {
                "experiment": "exp_bc",
                "sample": "sample01",
                "mode": "gex+crispr",
                "features": "GEX+GUIDES",
                "barcodes": "BC001+BC001",
            },
        ]
        config = _write_and_parse_config(tmp_path, libraries, samples)
        process_sample("sample01", config, str(cyto_outdir), str(agg_outdir))

        gex = ad.read_h5ad(agg_outdir / "sample01" / "sample01_gex.h5ad")
        crispr = ad.read_h5ad(agg_outdir / "sample01" / "sample01_crispr.h5ad")

        pattern = re.compile(r"^[ACGTN]+-BC001-1-exp_bc$")
        for idx in gex.obs.index:
            assert pattern.match(idx) is not None, (
                f"unexpected gex obs.index entry {idx!r}"
            )
        for idx in crispr.obs.index:
            assert pattern.match(idx) is not None, (
                f"unexpected crispr obs.index entry {idx!r}"
            )
        # Standard branch should have produced a real (non-NaN) join.
        assert (gex.obs["assignment"] == "guide_1").all()
        assert (gex.obs["n_umis_gex"] == 26).all()
