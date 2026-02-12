"""Unit tests for Flex-V2 384-barcode support."""

import pytest

from pycyto.config import (
    FLEX_V2_BARCODES,
    _assign_probeset,
    _detect_barcode_format,
    _expand_flex_v2_range,
    _get_flex_v2_prefix,
    _is_flex_v1_barcode,
    _is_flex_v2_barcode,
    _parse_flex_v2_component,
)


class TestFlexV2Detection:
    """Test Flex-V2 barcode detection functions."""

    def test_is_flex_v2_barcode_valid(self):
        """Test detection of valid Flex-V2 barcodes."""
        assert _is_flex_v2_barcode("A-A01")
        assert _is_flex_v2_barcode("A-H12")
        assert _is_flex_v2_barcode("B-C05")
        assert _is_flex_v2_barcode("D-H12")

    def test_is_flex_v2_barcode_valid_underscore(self):
        """Test detection of valid Flex-V2 barcodes with underscores."""
        assert _is_flex_v2_barcode("A_A01")
        assert _is_flex_v2_barcode("A_H12")
        assert _is_flex_v2_barcode("B_C05")
        assert _is_flex_v2_barcode("D_H12")

    def test_is_flex_v2_barcode_invalid_set(self):
        """Test rejection of invalid set prefixes."""
        assert not _is_flex_v2_barcode("E-A01")  # Invalid set
        assert not _is_flex_v2_barcode("Z-A01")

    def test_is_flex_v2_barcode_invalid_row(self):
        """Test rejection of invalid row letters."""
        assert not _is_flex_v2_barcode("A-I01")  # Invalid row (I is out of range)
        assert not _is_flex_v2_barcode("A-Z01")

    def test_is_flex_v2_barcode_invalid_column(self):
        """Test rejection of invalid column numbers."""
        assert not _is_flex_v2_barcode("A-A00")  # Column too low
        assert not _is_flex_v2_barcode("A-A13")  # Column too high
        assert not _is_flex_v2_barcode("A-A99")

    def test_is_flex_v2_barcode_not_flex_v1(self):
        """Test that Flex-V1 barcodes are not detected as Flex-V2."""
        assert not _is_flex_v2_barcode("BC001")
        assert not _is_flex_v2_barcode("CR016")
        assert not _is_flex_v2_barcode("AB008")

    def test_is_flex_v1_barcode_valid(self):
        """Test detection of valid Flex-V1 barcodes."""
        assert _is_flex_v1_barcode("BC001")
        assert _is_flex_v1_barcode("CR016")
        assert _is_flex_v1_barcode("AB008")

    def test_detect_format_v1(self):
        """Test format detection for Flex-V1."""
        assert _detect_barcode_format(["BC001", "BC002"]) == "flex-v1"
        assert _detect_barcode_format(["CR001"]) == "flex-v1"
        assert _detect_barcode_format(["AB016"]) == "flex-v1"

    def test_detect_format_v2(self):
        """Test format detection for Flex-V2."""
        assert _detect_barcode_format(["A-A01", "A-A02"]) == "flex-v2"
        assert _detect_barcode_format(["D-H12"]) == "flex-v2"

    def test_detect_format_empty(self):
        """Test that empty list raises error."""
        with pytest.raises(ValueError, match="Cannot detect format from empty"):
            _detect_barcode_format([])

    def test_detect_format_unknown(self):
        """Test that unknown format raises error."""
        with pytest.raises(ValueError, match="Unknown barcode format"):
            _detect_barcode_format(["INVALID"])

    def test_get_flex_v2_prefix(self):
        """Test extraction of Flex-V2 set prefix."""
        assert _get_flex_v2_prefix("A-A01") == "A"
        assert _get_flex_v2_prefix("B-C05") == "B"
        assert _get_flex_v2_prefix("C-H12") == "C"
        assert _get_flex_v2_prefix("D-A01") == "D"

    def test_get_flex_v2_prefix_invalid(self):
        """Test that invalid barcodes raise error."""
        with pytest.raises(ValueError, match="Invalid Flex-V2 barcode"):
            _get_flex_v2_prefix("BC001")


class TestFlexV2Expansion:
    """Test Flex-V2 barcode expansion functions."""

    def test_expand_flex_v2_range_same_row(self):
        """Test expansion within a single row."""
        result = _expand_flex_v2_range("A", "A", 1, "A", 12)
        assert len(result) == 12
        assert result[0] == "A-A01"
        assert result[-1] == "A-A12"

    def test_expand_flex_v2_range_partial_row(self):
        """Test expansion of partial row."""
        result = _expand_flex_v2_range("A", "A", 1, "A", 5)
        assert len(result) == 5
        assert result == ["A-A01", "A-A02", "A-A03", "A-A04", "A-A05"]

    def test_expand_flex_v2_range_multiple_rows(self):
        """Test expansion across multiple rows."""
        result = _expand_flex_v2_range("A", "A", 1, "B", 12)
        assert len(result) == 24  # 12 + 12
        assert result[0] == "A-A01"
        assert result[11] == "A-A12"
        assert result[12] == "A-B01"
        assert result[-1] == "A-B12"

    def test_expand_flex_v2_range_full_set(self):
        """Test expansion of full 96-well plate."""
        result = _expand_flex_v2_range("A", "A", 1, "H", 12)
        assert len(result) == 96  # 8 rows × 12 columns

    def test_expand_flex_v2_range_partial_start_end(self):
        """Test expansion with partial start and end rows."""
        result = _expand_flex_v2_range("B", "B", 5, "D", 3)
        # B05-B12 (8) + C01-C12 (12) + D01-D03 (3) = 23
        assert len(result) == 23
        assert result[0] == "B-B05"
        assert result[7] == "B-B12"
        assert result[8] == "B-C01"
        assert result[-1] == "B-D03"


class TestFlexV2Parsing:
    """Test Flex-V2 DSL parsing."""

    def test_parse_single_barcode(self):
        """Test parsing of single barcode."""
        result = _parse_flex_v2_component("A-A01")
        assert result == ["A-A01"]

    def test_parse_single_barcode_underscore(self):
        """Test parsing of single barcode with underscore (normalized to hyphen)."""
        result = _parse_flex_v2_component("A_A01")
        assert result == ["A-A01"]  # Should normalize to hyphen

    def test_parse_range_same_row(self):
        """Test parsing range within same row."""
        result = _parse_flex_v2_component("A-A01..A-A12")
        assert len(result) == 12
        assert result[0] == "A-A01"
        assert result[-1] == "A-A12"

    def test_parse_range_multiple_rows(self):
        """Test parsing range across multiple rows."""
        result = _parse_flex_v2_component("A-A01..A-H12")
        assert len(result) == 96

    def test_parse_selection(self):
        """Test parsing selection syntax."""
        result = _parse_flex_v2_component("A-A01|A-A05|A-C03")
        assert result == ["A-A01", "A-A05", "A-C03"]

    def test_parse_invalid_range_cross_set(self):
        """Test that cross-set ranges raise error."""
        with pytest.raises(
            ValueError, match="Cannot expand range across different sets"
        ):
            _parse_flex_v2_component("A-A01..B-A01")

    def test_parse_invalid_range_format(self):
        """Test that invalid range format raises error."""
        with pytest.raises(ValueError, match="Invalid Flex-V2 range notation"):
            _parse_flex_v2_component("A-A01..INVALID")

    def test_parse_invalid_single(self):
        """Test that invalid single barcode raises error."""
        with pytest.raises(ValueError, match="Invalid Flex-V2 barcode"):
            _parse_flex_v2_component("A-I01")


class TestProbeSetAssignment:
    """Test probe set assignment for Flex-V2."""

    def test_assign_probeset_flex_v1(self):
        """Test assignment for Flex-V1 barcodes."""
        assert _assign_probeset("BC001") == "BC"
        assert _assign_probeset("CR016") == "CR"
        assert _assign_probeset("AB008") == "AB"

    def test_assign_probeset_flex_v2(self):
        """Test assignment for Flex-V2 barcodes."""
        assert _assign_probeset("A-A01") == "A"
        assert _assign_probeset("B-C05") == "B"
        assert _assign_probeset("C-H12") == "C"
        assert _assign_probeset("D-A01") == "D"

    def test_assign_probeset_invalid(self):
        """Test that invalid barcode raises error."""
        with pytest.raises(ValueError, match="Unknown barcode format"):
            _assign_probeset("INVALID")


class TestFlexV2Constants:
    """Test Flex-V2 constant generation."""

    def test_flex_v2_barcodes_count(self):
        """Test that 384 Flex-V2 barcodes are generated."""
        assert len(FLEX_V2_BARCODES) == 384

    def test_flex_v2_barcodes_format(self):
        """Test that all generated barcodes follow correct format."""
        for barcode in FLEX_V2_BARCODES:
            assert _is_flex_v2_barcode(barcode), f"Invalid barcode: {barcode}"

    def test_flex_v2_barcodes_unique(self):
        """Test that all barcodes are unique."""
        assert len(set(FLEX_V2_BARCODES)) == 384

    def test_flex_v2_barcodes_coverage(self):
        """Test that all expected barcodes are present."""
        # Check first and last of each set
        assert "A-A01" in FLEX_V2_BARCODES
        assert "A-H12" in FLEX_V2_BARCODES
        assert "B-A01" in FLEX_V2_BARCODES
        assert "B-H12" in FLEX_V2_BARCODES
        assert "C-A01" in FLEX_V2_BARCODES
        assert "C-H12" in FLEX_V2_BARCODES
        assert "D-A01" in FLEX_V2_BARCODES
        assert "D-H12" in FLEX_V2_BARCODES
