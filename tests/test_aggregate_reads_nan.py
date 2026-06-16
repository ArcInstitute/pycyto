"""Regression tests for read-count NaNs on GEX-only cells in aggregation.

Bug: ``_process_gex_crispr_set`` built the per-cell metadata table anchored on
the CRISPR ``assignments`` table, which only contains cells detected in the
CRISPR modality. The final left-merge onto ``gex_adata.obs`` therefore left
GEX-only cells (no CRISPR mate) NaN in every merged column, including the
GEX-intrinsic ``n_reads_gex`` / ``n_umis_gex`` columns that actually exist in
the GEX reads table.
"""

import anndata as ad
import numpy as np
import polars as pl

from pycyto.aggregate import _process_gex_crispr_set

# ---------------------------------------------------------------------------
# Synthetic-fixture construction.
#
# We build a self-consistent scenario with two cells:
#   * ``cellBoth``  - present in GEX (adata + gex reads), CRISPR (adata +
#                     crispr reads), and the CRISPR assignments table.
#   * ``cellGexOnly`` - present ONLY in GEX (adata + gex reads). It has NO
#                     CRISPR mate: absent from assignments and crispr reads.
#
# ``match_barcode`` construction inside ``_process_gex_crispr_set`` (Flex-V1,
# non-CR path):
#   assignments.match_barcode = cell    + "-" + lane_id
#   reads_df.match_barcode    = cell_id + "-" + lane_id   (cell_id already
#                                                          built by the loader
#                                                          as barcode + "-" + bc_idx)
# The GEX adata obs index is built by the loader as <barcode> + "-" + lane_id.
# To make all three align we use a barcode that already embeds bc_idx so the
# adata index equals reads.match_barcode.
# ---------------------------------------------------------------------------

LANE_ID = "1"
BC_IDX = "TestGEX_BC1"

# adata index == <obs barcode> + "-" + lane_id, must equal reads match_barcode.
BARCODE_BOTH = "AAAACCCC-TestGEX_BC1"
BARCODE_GEX_ONLY = "GGGGTTTT-TestGEX_BC1"

MATCH_BOTH = f"{BARCODE_BOTH}-{LANE_ID}"
MATCH_GEX_ONLY = f"{BARCODE_GEX_ONLY}-{LANE_ID}"

# Known "real" GEX read values we expect to survive onto every GEX cell.
GEX_READS = {
    MATCH_BOTH: {"n_reads": 1000, "n_umis": 500},
    MATCH_GEX_ONLY: {"n_reads": 700, "n_umis": 350},
}
CRISPR_READS = {
    MATCH_BOTH: {"n_reads": 80, "n_umis": 40},
}
# CRISPR assignment for the "both" cell only.
ASSIGN_BOTH = {"assignment": "geneA|sgRNA1", "umis": 40, "moi": 1}


def _make_gex_adata() -> ad.AnnData:
    """GEX AnnData with both cells; index already in the post-load form."""
    obs_index = [BARCODE_BOTH, BARCODE_GEX_ONLY]
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
    )
    adata.obs_names = obs_index
    # index += "-" + lane_id, mimicking _load_gex_anndata_for_experiment_sample
    adata.obs["sample"] = "S1"
    adata.obs["experiment"] = "E1"
    adata.obs["lane_id"] = LANE_ID
    adata.obs["bc_idx"] = BC_IDX
    adata.obs.index = adata.obs.index + "-" + adata.obs["lane_id"].astype(str)
    return adata


def _make_crispr_adata() -> ad.AnnData:
    """CRISPR AnnData with only the 'both' cell."""
    adata = ad.AnnData(X=np.array([[5.0]], dtype=np.float32))
    adata.obs_names = [BARCODE_BOTH]
    adata.obs["sample"] = "S1"
    adata.obs["experiment"] = "E1"
    adata.obs["lane_id"] = LANE_ID
    adata.obs["bc_idx"] = BC_IDX
    adata.obs.index = adata.obs.index + "-" + adata.obs["lane_id"].astype(str)
    return adata


def _make_assignments() -> pl.DataFrame:
    """Assignments table (CRISPR modality) -- only the 'both' cell."""
    # `cell` + "-" + lane_id == match_barcode == BARCODE_BOTH + "-" + lane_id
    return pl.DataFrame(
        {
            "cell": [BARCODE_BOTH],
            "assignment": [ASSIGN_BOTH["assignment"]],
            "umis": [ASSIGN_BOTH["umis"]],
            "moi": [ASSIGN_BOTH["moi"]],
            "sample": ["S1"],
            "experiment": ["E1"],
            "lane_id": [LANE_ID],
            "bc_idx": [BC_IDX],
        }
    )


def _make_reads() -> pl.DataFrame:
    """Reads table holding both gex (per GEX cell) and crispr rows.

    ``cell_id`` is the loader-built ``barcode + "-" + bc_idx`` -- but our
    barcodes already embed the bc_idx, so cell_id is just the barcode here;
    ``match_barcode`` then becomes cell_id + "-" + lane_id, matching the adata
    index and the assignments match_barcode.
    """
    rows = []
    # gex reads: one per GEX cell (the bug source for the gex-only cell).
    for match, vals in GEX_READS.items():
        cell_id = match.rsplit("-", 1)[0]  # strip lane suffix back to cell_id
        rows.append(
            {
                "barcode": cell_id,
                "cell_id": cell_id,
                "n_reads": vals["n_reads"],
                "n_umis": vals["n_umis"],
                "bc_idx": BC_IDX,
                "lane_id": LANE_ID,
                "experiment": "E1",
                "sample": "S1",
                "mode": "gex",
            }
        )
    for match, vals in CRISPR_READS.items():
        cell_id = match.rsplit("-", 1)[0]
        rows.append(
            {
                "barcode": cell_id,
                "cell_id": cell_id,
                "n_reads": vals["n_reads"],
                "n_umis": vals["n_umis"],
                "bc_idx": BC_IDX,
                "lane_id": LANE_ID,
                "experiment": "E1",
                "sample": "S1",
                "mode": "crispr",
            }
        )
    return pl.DataFrame(rows)


def _run(tmp_path) -> ad.AnnData:
    """Run the function under test and read back the written GEX h5ad."""
    sample = "S1"
    sample_outdir = str(tmp_path)
    _process_gex_crispr_set(
        gex_adata_list=[_make_gex_adata()],
        crispr_adata_list=[_make_crispr_adata()],
        assignments_list=[_make_assignments()],
        reads_list=[_make_reads()],
        sample_outdir=sample_outdir,
        sample=sample,
        compress=False,
    )
    return ad.read_h5ad(f"{sample_outdir}/{sample}_gex.h5ad")


def test_gex_only_cell_keeps_gex_reads(tmp_path):
    """GEX-only cell must retain its real n_reads_gex / n_umis_gex (not NaN)."""
    gex = _run(tmp_path)
    obs = gex.obs

    assert MATCH_GEX_ONLY in obs.index, "GEX-only cell missing from obs"
    row = obs.loc[MATCH_GEX_ONLY]

    # The core bug assertion: these must be the real values, never NaN.
    assert not np.isnan(float(row["n_reads_gex"])), (
        "n_reads_gex is NaN for GEX-only cell (the bug)"
    )
    assert int(row["n_reads_gex"]) == GEX_READS[MATCH_GEX_ONLY]["n_reads"]
    assert int(row["n_umis_gex"]) == GEX_READS[MATCH_GEX_ONLY]["n_umis"]
    # No CRISPR mate -> crispr read counts fill to 0.
    assert int(row["n_reads_crispr"]) == 0
    assert int(row["n_umis_crispr"]) == 0


def test_both_modality_cell_unchanged(tmp_path):
    """Regression: a cell present in both modalities keeps all its values."""
    gex = _run(tmp_path)
    row = gex.obs.loc[MATCH_BOTH]

    assert int(row["n_reads_gex"]) == GEX_READS[MATCH_BOTH]["n_reads"]
    assert int(row["n_umis_gex"]) == GEX_READS[MATCH_BOTH]["n_umis"]
    assert int(row["n_reads_crispr"]) == CRISPR_READS[MATCH_BOTH]["n_reads"]
    assert int(row["n_umis_crispr"]) == CRISPR_READS[MATCH_BOTH]["n_umis"]
    assert row["assignment"] == ASSIGN_BOTH["assignment"]
    assert int(row["moi"]) == ASSIGN_BOTH["moi"]
    assert int(row["umis"]) == ASSIGN_BOTH["umis"]
