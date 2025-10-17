# pycyto Test Suite

This directory contains comprehensive tests for the pycyto package, focusing on the configuration parsing and barcode range functionality.

## Test Structure

### `test_barcode_parsing.py`
Tests the new barcode range and selection parsing functionality introduced to streamline barcode configuration.

**Key Features Tested:**
- **Range Syntax**: `BC1..8+CR1..8` → 8 barcode pairs
- **Selection Syntax**: `BC1|3|5+CR1|3|5` → 3 specific pairs
- **Mixed Syntax**: `BC1..4|7|12..14+CR1..4|7|12..14` → ranges + selections
- **Backward Compatibility**: Existing explicit syntax like `BC001+CR001|BC002+CR002`
- **Error Handling**: Mismatched lengths, invalid formats, etc.

**Test Classes:**
- `TestExpandRange`: Basic range expansion (`1..3` → `[1,2,3]`)
- `TestExpandSelection`: Selection parsing (`1|3|5..7` → `[1,3,5,6,7]`)
- `TestExpandBarcodeComponent`: Full barcode expansion (`BC1..3` → `['BC001','BC002','BC003']`)
- `TestParseBarcodes`: Complete barcode parsing logic
- `TestParseConfigIntegration`: End-to-end config file parsing

### `test_config_parsing.py`
Tests the configuration parsing system, especially the lane removal refactoring that eliminated the need to specify lane numbers in config files.

**Key Features Tested:**
- **Lane Removal**: Configs no longer require `"lane": "1|2|3|4"` specifications
- **Config Validation**: Proper validation of required vs optional fields
- **Mode Parsing**: Single (`"gex"`) and multi-mode (`"gex+crispr"`) handling
- **Feature Parsing**: Library feature validation and counting
- **Integration**: Full config parsing with temporary files

**Test Classes:**
- `TestConfigValidation`: Basic config structure and key validation
- `TestModeAndFeatureParsing`: Mode and feature string parsing
- `TestConfigParsingIntegration`: End-to-end config parsing scenarios

## Running Tests

### Run All Tests
```bash
uv run pytest tests/ -v
```

### Run Specific Test Module
```bash
uv run pytest tests/test_barcode_parsing.py -v
uv run pytest tests/test_config_parsing.py -v
```

### Run Specific Test Class
```bash
uv run pytest tests/test_barcode_parsing.py::TestExpandRange -v
```

### Run Specific Test
```bash
uv run pytest tests/test_barcode_parsing.py::TestExpandRange::test_basic_ranges -v
```

## Test Coverage

The test suite covers:

### ✅ Barcode Functionality
- Range expansion with inclusive bounds (`1..8` includes both 1 and 8)
- Selection parsing with pipe delimiters (`1|3|5`)
- Mixed range and selection syntax (`1..4|7|12..14`)
- Zero-padding of barcode numbers (`1` → `001`)
- Backward compatibility with explicit syntax
- Error handling for invalid formats and mismatched lengths

### ✅ Config Parsing
- Removal of lane requirements from config files
- Dynamic lane discovery during aggregation
- Mode and feature validation
- Expected prefix generation without specific lane numbers
- Integration with temporary config files
- Proper error handling for missing files and invalid configurations

### ✅ Integration Testing
- Full config parsing with temporary probe files
- Complex barcode range scenarios in real configs
- Verification of dataframe structure and column presence
- Backward compatibility with existing config formats

## Example Usage Scenarios Tested

### New Barcode Range Syntax
```json
{
  "barcodes": "BC1..8+CR1..8"           // 8 contiguous pairs
  "barcodes": "BC1|3|5+CR1|3|5"         // 3 specific pairs
  "barcodes": "BC1..4|7+CR1..4|7"       // Mixed ranges and selections
  "barcodes": "BC1..16"                 // Single mode, 16 barcodes
}
```

### Simplified Config Structure (No Lanes)
```json
{
  "samples": [
    {
      "experiment": "my_experiment",
      "mode": "gex+crispr",
      "features": "GEX_PROBE_LIST+CRISPR_PROBE_LIST",
      "sample": "my_sample",
      "barcodes": "BC1..8+CR1..8"
      // Note: No "lane" field required!
    }
  ]
}
```

## Test Philosophy

These tests follow several key principles:

1. **Comprehensive Coverage**: Test both happy paths and error conditions
2. **Backward Compatibility**: Ensure existing configs continue to work
3. **Integration Focus**: Test real-world usage scenarios, not just isolated functions
4. **Clear Error Messages**: Verify that error handling provides helpful feedback
5. **Temporary Files**: Use pytest's `tmp_path` fixture for safe file-based testing

## Future Scaling

The test suite is designed to handle future scaling to 96+ barcodes:
- Range syntax naturally scales: `BC1..96+CR1..96`
- Tests verify large ranges work correctly
- No hardcoded limits in test expectations
