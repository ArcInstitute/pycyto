"""Tests for config parsing functionality, especially lane removal refactoring."""

import json
import os
import sys

import pytest

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from pycyto.config import (
    EXPECTED_SAMPLE_KEYS,
    _parse_features,
    _parse_mode,
    _validate_json,
    _validate_keys,
    parse_config,
)


class TestConfigValidation:
    """Test config validation functions."""

    def test_expected_sample_keys_no_lane(self):
        """Test that 'lane' is not in expected sample keys."""
        assert "lane" not in EXPECTED_SAMPLE_KEYS
        expected_keys = ["experiment", "sample", "mode", "barcodes", "features"]
        assert set(EXPECTED_SAMPLE_KEYS) == set(expected_keys)

    def test_validate_json_success(self):
        """Test successful JSON validation."""
        valid_data = {
            "libraries": {"GEX_PROBE_LIST": "path/to/gex.tsv"},
            "samples": [],
        }
        # Should not raise any exception
        _validate_json(valid_data)

    def test_validate_json_missing_keys(self):
        """Test JSON validation with missing required keys."""
        with pytest.raises(ValueError, match="Missing key in data: libraries"):
            _validate_json({"samples": []})

        with pytest.raises(ValueError, match="Missing key in data: samples"):
            _validate_json({"libraries": {}})

    def test_validate_keys_success(self):
        """Test successful sample entry validation."""
        valid_entry = {
            "experiment": "test_exp",
            "sample": "test_sample",
            "mode": "gex",
            "barcodes": "BC001",
            "features": "GEX_PROBE_LIST",
        }
        # Should not raise any exception
        _validate_keys(valid_entry)

    def test_validate_keys_missing_lane_ok(self):
        """Test that missing 'lane' key doesn't cause validation error."""
        entry_without_lane = {
            "experiment": "test_exp",
            "sample": "test_sample",
            "mode": "gex",
            "barcodes": "BC001",
            "features": "GEX_PROBE_LIST",
            # Note: no 'lane' key, which should be fine now
        }
        # Should not raise any exception
        _validate_keys(entry_without_lane)

    def test_validate_keys_missing_required(self):
        """Test validation failure with missing required keys."""
        incomplete_entry = {
            "experiment": "test_exp",
            "sample": "test_sample",
            # Missing mode, barcodes, features
        }
        with pytest.raises(ValueError, match="Missing keys in entry"):
            _validate_keys(incomplete_entry)


class TestModeAndFeatureParsing:
    """Test mode and feature parsing functions."""

    def test_parse_single_mode(self):
        """Test parsing single mode."""
        entry = {"mode": "gex"}
        assert _parse_mode(entry) == ["gex"]

        entry = {"mode": "crispr"}
        assert _parse_mode(entry) == ["crispr"]

    def test_parse_multiple_modes(self):
        """Test parsing multiple modes."""
        entry = {"mode": "gex+crispr"}
        assert _parse_mode(entry) == ["gex", "crispr"]

        entry = {"mode": "gex+crispr+ab"}
        assert _parse_mode(entry) == ["gex", "crispr", "ab"]

    def test_parse_invalid_mode(self):
        """Test error handling for invalid modes."""
        entry = {"mode": "invalid_mode"}
        with pytest.raises(ValueError, match="Invalid mode invalid_mode found"):
            _parse_mode(entry)

        entry = {"mode": "gex+invalid_mode"}
        with pytest.raises(ValueError, match="Invalid mode found"):
            _parse_mode(entry)

    def test_parse_features_single(self):
        """Test parsing single feature."""
        entry = {"features": "GEX_PROBE_LIST"}
        known_features = ["GEX_PROBE_LIST", "CRISPR_PROBE_LIST"]
        result = _parse_features(entry, 1, known_features)
        assert result == ["GEX_PROBE_LIST"]

    def test_parse_features_multiple(self):
        """Test parsing multiple features."""
        entry = {"features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST"}
        known_features = ["GEX_PROBE_LIST", "CRISPR_PROBE_LIST"]
        result = _parse_features(entry, 2, known_features)
        assert result == ["GEX_PROBE_LIST", "CRISPR_PROBE_LIST"]

    def test_parse_features_wrong_count(self):
        """Test error when feature count doesn't match library count."""
        entry = {"features": "GEX_PROBE_LIST"}
        known_features = ["GEX_PROBE_LIST", "CRISPR_PROBE_LIST"]
        with pytest.raises(ValueError, match="Invalid number of features found"):
            _parse_features(entry, 2, known_features)  # Expecting 2, got 1

    def test_parse_features_unknown_feature(self):
        """Test error when using unknown feature."""
        entry = {"features": "UNKNOWN_FEATURE"}
        known_features = ["GEX_PROBE_LIST", "CRISPR_PROBE_LIST"]
        with pytest.raises(ValueError, match="Invalid feature found"):
            _parse_features(entry, 1, known_features)


class TestConfigParsingIntegration:
    """Test full config parsing integration."""

    def create_temp_config(self, config_data, tmp_path):
        """Helper to create temporary config and probe files."""
        # Create probe files
        gex_probes = tmp_path / "gex_probes.tsv"
        gex_probes.write_text("gene\tsequence\nGENE1\tATCG\nGENE2\tGCTA")

        crispr_probes = tmp_path / "crispr_guides.tsv"
        crispr_probes.write_text("guide\tsequence\nGUIDE1\tATCG\nGUIDE2\tGCTA")

        # Update config to use temp files
        config_data["libraries"]["GEX_PROBE_LIST"] = str(gex_probes)
        config_data["libraries"]["CRISPR_PROBE_LIST"] = str(crispr_probes)

        # Write config file
        config_file = tmp_path / "test_config.json"
        with open(config_file, "w") as f:
            json.dump(config_data, f, indent=2)

        return config_file

    def test_basic_config_without_lanes(self, tmp_path):
        """Test parsing a basic config without lane specifications."""
        config_data = {
            "libraries": {
                "GEX_PROBE_LIST": "./gex_probes.tsv",
                "CRISPR_PROBE_LIST": "./crispr_guides.tsv",
            },
            "samples": [
                {
                    "experiment": "test_exp",
                    "mode": "gex",
                    "features": "GEX_PROBE_LIST",
                    "sample": "single_mode_sample",
                    "barcodes": "BC001",
                }
            ],
        }

        config_file = self.create_temp_config(config_data, tmp_path)
        result = parse_config(str(config_file))

        # Verify basic structure
        assert result.shape[0] == 1  # 1 sample, 1 mode, 1 barcode
        assert "lane" not in result.columns  # Lane column should not exist
        assert list(result["sample"])[0] == "single_mode_sample"
        assert list(result["mode"])[0] == "gex"
        assert list(result["bc_component"])[0] == "BC001"

    def test_expected_prefix_without_lanes(self, tmp_path):
        """Test that expected_prefix is generated correctly without lanes."""
        config_data = {
            "libraries": {
                "GEX_PROBE_LIST": "./gex_probes.tsv",
                "CRISPR_PROBE_LIST": "./crispr_guides.tsv",
            },
            "samples": [
                {
                    "experiment": "test_experiment",
                    "mode": "gex+crispr",
                    "features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
                    "sample": "test_sample",
                    "barcodes": "BC001+CR001",
                }
            ],
        }

        config_file = self.create_temp_config(config_data, tmp_path)
        result = parse_config(str(config_file))

        # Check expected prefixes
        expected_prefixes = result["expected_prefix"].unique().sort().to_list()
        assert "test_experiment_GEX_Lane" in expected_prefixes
        assert "test_experiment_CRISPR_Lane" in expected_prefixes

        # Ensure no specific lane numbers in the prefix
        for prefix in expected_prefixes:
            assert not any(
                char.isdigit()
                for char in prefix.split("_Lane")[1]
                if len(prefix.split("_Lane")) > 1
            )

    def test_multiple_samples_and_experiments(self, tmp_path):
        """Test parsing config with multiple samples and experiments."""
        config_data = {
            "libraries": {
                "GEX_PROBE_LIST": "./gex_probes.tsv",
                "CRISPR_PROBE_LIST": "./crispr_guides.tsv",
            },
            "samples": [
                {
                    "experiment": "exp1",
                    "mode": "gex",
                    "features": "GEX_PROBE_LIST",
                    "sample": "sample1",
                    "barcodes": "BC1..3",
                },
                {
                    "experiment": "exp1",
                    "mode": "crispr",
                    "features": "CRISPR_PROBE_LIST",
                    "sample": "sample2",
                    "barcodes": "CR1..2",
                },
                {
                    "experiment": "exp2",
                    "mode": "gex+crispr",
                    "features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
                    "sample": "sample3",
                    "barcodes": "BC001+CR001",
                },
            ],
        }

        config_file = self.create_temp_config(config_data, tmp_path)
        result = parse_config(str(config_file))

        # Verify we have the right number of rows
        # exp1/sample1: 3 GEX barcodes = 3 rows
        # exp1/sample2: 2 CRISPR barcodes = 2 rows
        # exp2/sample3: 1 GEX + 1 CRISPR = 2 rows
        # Total: 7 rows
        assert result.shape[0] == 7

        # Verify unique experiments
        experiments = result["experiment"].unique().sort().to_list()
        assert experiments == ["exp1", "exp2"]

        # Verify unique samples
        samples = result["sample"].unique().sort().to_list()
        assert samples == ["sample1", "sample2", "sample3"]

    def test_complex_barcode_ranges_integration(self, tmp_path):
        """Test that complex barcode ranges work in full config parsing."""
        config_data = {
            "libraries": {
                "GEX_PROBE_LIST": "./gex_probes.tsv",
                "CRISPR_PROBE_LIST": "./crispr_guides.tsv",
            },
            "samples": [
                {
                    "experiment": "range_test",
                    "mode": "gex+crispr",
                    "features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
                    "sample": "complex_range",
                    "barcodes": "BC1..3|7|10..12+CR1..3|7|10..12",
                }
            ],
        }

        config_file = self.create_temp_config(config_data, tmp_path)
        result = parse_config(str(config_file))

        # Should have 7 pairs × 2 modes = 14 rows
        assert result.shape[0] == 14

        # Check GEX barcodes
        gex_barcodes = (
            result.filter(result["mode"] == "gex")["bc_component"]
            .unique()
            .sort()
            .to_list()
        )
        expected_gex = ["BC001", "BC002", "BC003", "BC007", "BC010", "BC011", "BC012"]
        assert gex_barcodes == expected_gex

        # Check CRISPR barcodes
        crispr_barcodes = (
            result.filter(result["mode"] == "crispr")["bc_component"]
            .unique()
            .sort()
            .to_list()
        )
        expected_crispr = [
            "CR001",
            "CR002",
            "CR003",
            "CR007",
            "CR010",
            "CR011",
            "CR012",
        ]
        assert crispr_barcodes == expected_crispr

    def test_config_with_missing_probe_files(self, tmp_path):
        """Test error handling when probe files don't exist."""
        config_data = {
            "libraries": {
                "GEX_PROBE_LIST": "./nonexistent_gex.tsv",
            },
            "samples": [
                {
                    "experiment": "test_exp",
                    "mode": "gex",
                    "features": "GEX_PROBE_LIST",
                    "sample": "test_sample",
                    "barcodes": "BC001",
                }
            ],
        }

        config_file = tmp_path / "bad_config.json"
        with open(config_file, "w") as f:
            json.dump(config_data, f)

        parse_config(str(config_file))

    def test_dataframe_structure_and_columns(self, tmp_path):
        """Test that the resulting dataframe has the correct structure."""
        config_data = {
            "libraries": {
                "GEX_PROBE_LIST": "./gex_probes.tsv",
            },
            "samples": [
                {
                    "experiment": "structure_test",
                    "mode": "gex",
                    "features": "GEX_PROBE_LIST",
                    "sample": "structure_sample",
                    "barcodes": "BC1..2",
                }
            ],
        }

        config_file = self.create_temp_config(config_data, tmp_path)
        result = parse_config(str(config_file))

        # Check all expected columns are present
        expected_columns = [
            "experiment",
            "sample",
            "mode",
            "bc_component",
            "bc_idx",
            "features",
            "probe_set",
            "feature_path",
            "expected_prefix",
        ]
        for col in expected_columns:
            assert col in result.columns

        # Ensure 'lane' is NOT in columns
        assert "lane" not in result.columns

        # Check data types and values for a specific row
        first_row = result.row(0)
        row_dict = dict(zip(result.columns, first_row))

        assert row_dict["experiment"] == "structure_test"
        assert row_dict["sample"] == "structure_sample"
        assert row_dict["mode"] == "gex"
        assert row_dict["bc_component"] in ["BC001", "BC002"]
        assert row_dict["bc_idx"] in [0, 1]
        assert row_dict["probe_set"] == "BC"
        assert row_dict["expected_prefix"] == "structure_test_GEX_Lane"
