"""Tests for barcode range parsing functionality."""

import os
import sys

import pytest

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from pycyto.config import (
    _expand_barcode_component,
    _expand_range,
    _expand_selection,
    _parse_barcodes,
    parse_config,
)


class TestExpandRange:
    """Test the _expand_range helper function."""

    def test_basic_ranges(self):
        """Test basic range expansion."""
        assert _expand_range("1..3") == [1, 2, 3]
        assert _expand_range("5..7") == [5, 6, 7]
        assert _expand_range("10..12") == [10, 11, 12]

    def test_single_element_range(self):
        """Test ranges with start == end."""
        assert _expand_range("1..1") == [1]
        assert _expand_range("5..5") == [5]

    def test_invalid_format_errors(self):
        """Test error handling for invalid formats."""
        with pytest.raises(ValueError, match="Invalid range format"):
            _expand_range("1.3")  # Wrong separator

        with pytest.raises(ValueError, match="Invalid range format"):
            _expand_range("1...3")  # Too many dots

        with pytest.raises(ValueError, match="Invalid range format"):
            _expand_range("abc..def")  # Non-numeric

    def test_invalid_order_error(self):
        """Test error when start > end."""
        with pytest.raises(ValueError, match="Invalid range \\(start > end\\)"):
            _expand_range("5..3")


class TestExpandSelection:
    """Test the _expand_selection helper function."""

    def test_single_number(self):
        """Test single number selection."""
        assert _expand_selection("1") == [1]
        assert _expand_selection("42") == [42]

    def test_pipe_separated_numbers(self):
        """Test pipe-separated individual numbers."""
        assert _expand_selection("1|3|5") == [1, 3, 5]
        assert _expand_selection("10|2|7") == [2, 7, 10]  # Should be sorted

    def test_ranges_in_selection(self):
        """Test ranges within selections."""
        assert _expand_selection("1..3|7|9..10") == [1, 2, 3, 7, 9, 10]
        assert _expand_selection("5..7") == [5, 6, 7]

    def test_mixed_complex_selection(self):
        """Test complex mixed selections."""
        assert _expand_selection("1|3..5|8|10..12") == [1, 3, 4, 5, 8, 10, 11, 12]

    def test_duplicate_removal(self):
        """Test that duplicates are removed."""
        assert _expand_selection("1|2|1|3|2") == [1, 2, 3]
        assert _expand_selection("1..3|2..4") == [1, 2, 3, 4]

    def test_invalid_selection_error(self):
        """Test error handling for invalid selections."""
        with pytest.raises(ValueError, match="Invalid selection format"):
            _expand_selection("abc")


class TestExpandBarcodeComponent:
    """Test the _expand_barcode_component helper function."""

    def test_explicit_barcode_backward_compatibility(self):
        """Test that explicit barcodes are returned as-is."""
        assert _expand_barcode_component("BC001") == ["BC001"]
        assert _expand_barcode_component("CR005") == ["CR005"]

    def test_range_expansion(self):
        """Test range expansion in barcode components."""
        assert _expand_barcode_component("BC1..3") == ["BC001", "BC002", "BC003"]
        assert _expand_barcode_component("CR5..7") == ["CR005", "CR006", "CR007"]

    def test_selection_expansion(self):
        """Test selection expansion in barcode components."""
        assert _expand_barcode_component("BC1|3|5") == ["BC001", "BC003", "BC005"]
        assert _expand_barcode_component("CR2|4|6") == ["CR002", "CR004", "CR006"]

    def test_mixed_range_selection(self):
        """Test mixed ranges and selections."""
        expected = ["BC001", "BC002", "BC003", "BC007", "BC009", "BC010"]
        assert _expand_barcode_component("BC1..3|7|9..10") == expected

    def test_zero_padding(self):
        """Test that numbers are zero-padded correctly."""
        assert _expand_barcode_component("BC1") == ["BC001"]
        assert _expand_barcode_component("CR12") == ["CR012"]
        assert _expand_barcode_component("AB1..2") == ["AB001", "AB002"]

    def test_invalid_component_error(self):
        """Test error handling for invalid components."""
        with pytest.raises(ValueError, match="Invalid barcode component format"):
            _expand_barcode_component("INVALIDFORMAT")


class TestParseBarcodes:
    """Test the main _parse_barcodes function."""

    def test_explicit_single_pair_backward_compatibility(self):
        """Test backward compatibility with explicit single pairs."""
        entry = {"barcodes": "BC001+CR001"}
        result = _parse_barcodes(entry, 2)
        expected = [["BC001", "CR001"]]
        assert result == expected

    def test_explicit_multiple_pairs_backward_compatibility(self):
        """Test backward compatibility with explicit multiple pairs."""
        entry = {"barcodes": "BC001+CR001|BC002+CR002"}
        result = _parse_barcodes(entry, 2)
        expected = [["BC001", "CR001"], ["BC002", "CR002"]]
        assert result == expected

    def test_range_syntax(self):
        """Test new range syntax."""
        entry = {"barcodes": "BC1..4+CR1..4"}
        result = _parse_barcodes(entry, 2)
        expected = [
            ["BC001", "CR001"],
            ["BC002", "CR002"],
            ["BC003", "CR003"],
            ["BC004", "CR004"],
        ]
        assert result == expected

    def test_selection_syntax(self):
        """Test new selection syntax."""
        entry = {"barcodes": "BC1|3|5+CR1|3|5"}
        result = _parse_barcodes(entry, 2)
        expected = [["BC001", "CR001"], ["BC003", "CR003"], ["BC005", "CR005"]]
        assert result == expected

    def test_mixed_range_selection_syntax(self):
        """Test mixed ranges and selections."""
        entry = {"barcodes": "BC1..3|7|9..10+CR1..3|7|9..10"}
        result = _parse_barcodes(entry, 2)
        expected = [
            ["BC001", "CR001"],
            ["BC002", "CR002"],
            ["BC003", "CR003"],
            ["BC007", "CR007"],
            ["BC009", "CR009"],
            ["BC010", "CR010"],
        ]
        assert result == expected

    def test_single_mode_library(self):
        """Test single-mode library (e.g., GEX only)."""
        entry = {"barcodes": "BC1..3"}
        result = _parse_barcodes(entry, 1)
        expected = [["BC001"], ["BC002"], ["BC003"]]
        assert result == expected

    def test_mixed_combination_types(self):
        """Test mixing explicit and range syntax in multiple combinations."""
        entry = {"barcodes": "BC1..2+CR1..2|BC005+CR005"}
        result = _parse_barcodes(entry, 2)
        expected = [["BC001", "CR001"], ["BC002", "CR002"], ["BC005", "CR005"]]
        assert result == expected

    def test_error_mismatched_lengths(self):
        """Test error when component lengths don't match."""
        entry = {"barcodes": "BC1..3+CR1..5"}
        with pytest.raises(ValueError, match="Mismatched component lengths"):
            _parse_barcodes(entry, 2)

    def test_error_wrong_number_components(self):
        """Test error when wrong number of components provided."""
        entry = {"barcodes": "BC1..3"}
        with pytest.raises(ValueError, match="Invalid number of barcodes found"):
            _parse_barcodes(entry, 2)  # Expecting 2 components, got 1

    def test_error_invalid_barcode_names(self):
        """Test error when invalid barcode names are generated."""
        entry = {"barcodes": "XY1..2+ZZ1..2"}  # Invalid prefixes
        with pytest.raises(ValueError, match="Invalid barcode found"):
            _parse_barcodes(entry, 2)


class TestParseConfigIntegration:
    """Test integration with parse_config function."""

    def create_test_config(self, samples):
        """Helper to create test config dictionary."""
        return {
            "libraries": {
                "GEX_PROBE_LIST": "./examples/gex_probes.tsv",
                "CRISPR_PROBE_LIST": "./examples/crispr_guides.tsv",
            },
            "samples": samples,
        }

    def test_range_syntax_integration(self, tmp_path):
        """Test that range syntax works in full config parsing."""
        config_data = self.create_test_config(
            [
                {
                    "experiment": "test_exp",
                    "mode": "gex+crispr",
                    "features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
                    "sample": "range_sample",
                    "barcodes": "BC1..4+CR1..4",
                }
            ]
        )

        # Create temporary config file
        config_file = tmp_path / "test_config.json"
        import json

        with open(config_file, "w") as f:
            json.dump(config_data, f)

        # Create temporary probe files
        gex_probes = tmp_path / "gex_probes.tsv"
        gex_probes.write_text("gene\tsequence\nGENE1\tATCG")

        crispr_probes = tmp_path / "crispr_guides.tsv"
        crispr_probes.write_text("guide\tsequence\nGUIDE1\tATCG")

        # Update config to use temporary files
        config_data["libraries"]["GEX_PROBE_LIST"] = str(gex_probes)
        config_data["libraries"]["CRISPR_PROBE_LIST"] = str(crispr_probes)
        with open(config_file, "w") as f:
            json.dump(config_data, f)

        # Parse config
        result = parse_config(str(config_file))

        # Verify results
        assert result.shape[0] == 8  # 4 pairs × 2 modes
        bc_components = result.filter(result["mode"] == "gex")["bc_component"].to_list()
        expected_bcs = ["BC001", "BC002", "BC003", "BC004"]
        assert sorted(bc_components) == expected_bcs

    def test_backward_compatibility_integration(self, tmp_path):
        """Test that explicit syntax still works in full config parsing."""
        config_data = self.create_test_config(
            [
                {
                    "experiment": "test_exp",
                    "mode": "gex+crispr",
                    "features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
                    "sample": "explicit_sample",
                    "barcodes": "BC001+CR001|BC003+CR003",
                }
            ]
        )

        # Create temporary files as before
        config_file = tmp_path / "test_config.json"
        gex_probes = tmp_path / "gex_probes.tsv"
        crispr_probes = tmp_path / "crispr_guides.tsv"

        gex_probes.write_text("gene\tsequence\nGENE1\tATCG")
        crispr_probes.write_text("guide\tsequence\nGUIDE1\tATCG")

        config_data["libraries"]["GEX_PROBE_LIST"] = str(gex_probes)
        config_data["libraries"]["CRISPR_PROBE_LIST"] = str(crispr_probes)

        import json

        with open(config_file, "w") as f:
            json.dump(config_data, f)

        # Parse config
        result = parse_config(str(config_file))

        # Verify results
        assert result.shape[0] == 4  # 2 pairs × 2 modes
        bc_components = result.filter(result["mode"] == "gex")["bc_component"].to_list()
        expected_bcs = ["BC001", "BC003"]
        assert sorted(bc_components) == expected_bcs
