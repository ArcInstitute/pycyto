"""Regression test for multi-guide (MOI>1) assignment-column dtype inference.

Bug: ``_load_assignments_for_experiment_sample`` reads cyto's per-barcode
``<bc>.assignments.tsv`` with ``pl.read_csv``. For cells assigned more than one
guide (MOI>1), the per-guide columns ``assignment``, ``guide_ids_original``,
``umis``, ``fdr`` and ``log_odds`` hold pipe-delimited values (e.g.
``guide_ids_original = "129|262"``). Polars infers each column's dtype from a
bounded row sample (``infer_schema_length``); when every row in that sample is a
single-guide cell, the numeric-looking columns are inferred as ``i64`` / ``f64``
and the first later multi-guide row blows up with::

    ComputeError: could not parse `129|262` as dtype `i64`
                  at column 'guide_ids_original' (column number 7)

The fix reads those per-guide columns as ``String`` regardless of the sample, so
the dtype is deterministic across every per-barcode file (and uniform for the
later ``pl.concat``).
"""

import polars as pl

# Import the production constant so this test tracks the source of truth: if a
# new per-guide column is added there, these dtype assertions cover it too.
from pycyto.aggregate import (
    _ASSIGNMENT_PER_GUIDE_COLS,
    _load_assignments_for_experiment_sample,
)

# cyto's assignments.tsv schema (one row per CRISPR-detected cell).
HEADER = [
    "cell_id",
    "submatrix_id",
    "cell",
    "moi",
    "n_umi",
    "assignment",
    "guide_ids_original",
    "umis",
    "fdr",
    "log_odds",
    "tested",
]

PER_GUIDE_COLS = list(_ASSIGNMENT_PER_GUIDE_COLS)

# Push the first multi-guide row well past any plausible ``infer_schema_length``
# default (polars' is 100), so a naive ``read_csv`` infers the per-guide columns
# as numeric and then crashes on the pipe-delimited value -- faithfully
# reproducing the GCP failure. The String-dtype assertions below are independent
# of this window, so they keep guarding the fix even if polars' default changes.
N_SINGLE_GUIDE_ROWS = 1000


def _write_assignments_tsv(path, *, with_multiguide: bool) -> None:
    lines = ["\t".join(HEADER)]
    for i in range(N_SINGLE_GUIDE_ROWS):
        # single-guide cell: every per-guide column is a bare scalar.
        lines.append(
            "\t".join(
                [
                    str(1000 + i),  # cell_id
                    str(i),  # submatrix_id
                    f"CELL{i:05d}-D-A01",  # cell
                    "1",  # moi
                    "6",  # n_umi
                    "GENE_Protosp001_A",  # assignment
                    "7",  # guide_ids_original  -> inferred i64
                    "2",  # umis                -> inferred i64
                    "1.5e-06",  # fdr           -> inferred f64
                    "13.7",  # log_odds         -> inferred f64
                    "true",  # tested
                ]
            )
        )
    if with_multiguide:
        # MOI=3 cell appearing AFTER the inference window: pipe-delimited values.
        lines.append(
            "\t".join(
                [
                    "9999",
                    "210",
                    "CELLMULTI-D-A01",
                    "3",
                    "24",
                    "BLNK_Protosp104_B|APOE_Protosp168_B|PICALM_Protosp406_B",
                    "65|129|262",
                    "2|6|3",
                    "1.9e-06|3.4e-08|0.0002",
                    "6.1|17.1|8.2",
                    "true",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n")


def _write_bc(tmp_path, bc: str, *, with_multiguide: bool) -> None:
    assignments_dir = tmp_path / "assignments"
    assignments_dir.mkdir(exist_ok=True)
    _write_assignments_tsv(
        assignments_dir / f"{bc}.assignments.tsv", with_multiguide=with_multiguide
    )


def _load(tmp_path, bcs: list[str]) -> list[pl.DataFrame]:
    return _load_assignments_for_experiment_sample(
        root=str(tmp_path),
        crispr_bcs=bcs,
        lane_id="1",
        experiment="E1",
        sample="S1",
    )


def test_multiguide_assignments_load_without_crash(tmp_path):
    """A late multi-guide row must not crash the read; values survive intact."""
    _write_bc(tmp_path, "BC001", with_multiguide=True)
    (df,) = _load(tmp_path, ["BC001"])
    assert df.height == N_SINGLE_GUIDE_ROWS + 1
    # The pipe-delimited multi-guide values land on the right (multi-guide) row.
    multi_row = df.filter(pl.col("cell") == "CELLMULTI-D-A01")
    assert multi_row.height == 1
    assert multi_row["guide_ids_original"].item() == "65|129|262"
    assert multi_row["umis"].item() == "2|6|3"


def test_per_guide_columns_are_strings(tmp_path):
    """Per-guide columns are read as String, deterministically (no numeric infer)."""
    _write_bc(tmp_path, "BC001", with_multiguide=True)
    (df,) = _load(tmp_path, ["BC001"])
    for col in PER_GUIDE_COLS:
        assert df.schema[col] == pl.String, (
            f"{col} should be String, got {df.schema[col]}"
        )


def test_per_guide_columns_uniform_dtype_for_single_guide_only_file(tmp_path):
    """A barcode with NO multi-guide cell must still yield String per-guide columns.

    Otherwise it infers numeric dtypes and the later cross-barcode ``pl.concat``
    sees a mixed schema (some files i64, some String).
    """
    _write_bc(tmp_path, "BC001", with_multiguide=False)
    (df,) = _load(tmp_path, ["BC001"])
    for col in PER_GUIDE_COLS:
        assert df.schema[col] == pl.String, (
            f"{col} should be String, got {df.schema[col]}"
        )


def test_cross_file_concat_uniform_dtype(tmp_path):
    """The aggregation's cross-barcode concat must not hit a mixed schema.

    Mirrors ``_process_gex_crispr_set``'s ``pl.concat(..., how="vertical_relaxed")``
    over per-barcode frames: one barcode with no multi-guide cell (which would
    naively infer numeric) and one with multi-guide cells (which infers String).
    This guard is independent of polars' inference window.
    """
    _write_bc(tmp_path, "BC001", with_multiguide=False)
    _write_bc(tmp_path, "BC002", with_multiguide=True)
    dfs = _load(tmp_path, ["BC001", "BC002"])
    assert len(dfs) == 2
    combined = pl.concat(dfs, how="vertical_relaxed")
    for col in PER_GUIDE_COLS:
        assert combined.schema[col] == pl.String, (
            f"{col} should be String after concat, got {combined.schema[col]}"
        )
    assert "65|129|262" in combined["guide_ids_original"].to_list()
