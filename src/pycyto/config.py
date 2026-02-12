import json
import os
import re

import polars as pl

KNOWN_LIBMODES = ["gex", "crispr", "ab"]

# Flex-V1 barcodes (16-plex)
FLEX_V1_PROBE_SETS = ["BC", "CR", "AB"]
FLEX_V1_BARCODES = [
    f"{name}0{i:02d}" for i in range(1, 17) for name in FLEX_V1_PROBE_SETS
]

# Flex-V2 barcodes (384-plex)
# Pattern: [ABCD]-[ABCDEFGH][01-12] or [ABCD]_[ABCDEFGH][01-12]
FLEX_V2_SETS = ["A", "B", "C", "D"]
FLEX_V2_ROWS = ["A", "B", "C", "D", "E", "F", "G", "H"]
FLEX_V2_COLS = [f"{i:02d}" for i in range(1, 13)]
FLEX_V2_BARCODES = [
    f"{set_}-{row}{col}"
    for set_ in FLEX_V2_SETS
    for row in FLEX_V2_ROWS
    for col in FLEX_V2_COLS
]
# Also include underscore variants for validation
FLEX_V2_BARCODES_UNDERSCORE = [bc.replace("-", "_") for bc in FLEX_V2_BARCODES]

# Combined list for validation
KNOWN_PROBE_SET = FLEX_V1_PROBE_SETS + FLEX_V2_SETS
KNOWN_BARCODES = FLEX_V1_BARCODES + FLEX_V2_BARCODES + FLEX_V2_BARCODES_UNDERSCORE
EXPECTED_SAMPLE_KEYS = [
    "experiment",
    "sample",
    "mode",
    "barcodes",
    "features",
]
EXPECTED_KEYS = [
    "libraries",
    "samples",
]


def _is_flex_v1_barcode(barcode: str) -> bool:
    """Check if barcode follows Flex-V1 format: BC001, CR001, AB001, etc."""
    return re.match(r"^(BC|CR|AB)0\d{2}$", barcode) is not None


def _is_flex_v2_barcode(barcode: str) -> bool:
    """Check if barcode follows Flex-V2 format: A-A01, B-C05, D-H12, etc.

    Also accepts underscore separator: A_A01, B_C05, D_H12
    """
    return re.match(r"^[ABCD][-_][ABCDEFGH](0[1-9]|1[0-2])$", barcode) is not None


def _detect_barcode_format(barcodes: list[str]) -> str:
    """Detect the barcode format from a list of barcodes."""
    if not barcodes:
        raise ValueError("Cannot detect format from empty barcode list")

    sample = barcodes[0]
    if _is_flex_v1_barcode(sample):
        return "flex-v1"
    elif _is_flex_v2_barcode(sample):
        return "flex-v2"
    else:
        raise ValueError(f"Unknown barcode format: {sample}")


def _normalize_flex_v2_barcode(barcode: str) -> str:
    """Normalize Flex-V2 barcode by replacing underscore with hyphen.

    A_A01 → A-A01
    A-A01 → A-A01 (unchanged)
    """
    if "_" in barcode and _is_flex_v2_barcode(barcode):
        return barcode.replace("_", "-")
    return barcode


def _get_flex_v2_prefix(barcode: str) -> str:
    """Extract the set prefix from a Flex-V2 barcode: A-A01 → A, A_A01 → A"""
    match = re.match(r"^([ABCD])[-_]", barcode)
    if match:
        return match.group(1)
    raise ValueError(f"Invalid Flex-V2 barcode: {barcode}")


def _expand_range(range_str: str) -> list[int]:
    """Expand a range string like '5..7' into [5, 6, 7]."""
    if ".." not in range_str:
        raise ValueError(f"Invalid range format: {range_str}")

    parts = range_str.split("..")
    if len(parts) != 2:
        raise ValueError(f"Invalid range format: {range_str}")

    try:
        start = int(parts[0])
        end = int(parts[1])
    except ValueError:
        raise ValueError(f"Invalid range format: {range_str}")

    if start > end:
        raise ValueError(f"Invalid range (start > end): {range_str}")

    return list(range(start, end + 1))


def _expand_selection(selection: str) -> list[int]:
    """Expand a selection string like '1|3|5..7|12' into [1, 3, 5, 6, 7, 12]."""
    if "|" not in selection and ".." not in selection:
        # Single number
        try:
            return [int(selection)]
        except ValueError:
            raise ValueError(f"Invalid selection format: {selection}")

    result = []
    parts = selection.split("|")

    for part in parts:
        if ".." in part:
            result.extend(_expand_range(part))
        else:
            try:
                result.append(int(part))
            except ValueError:
                raise ValueError(f"Invalid selection format: {selection}")

    # Remove duplicates and sort
    return sorted(list(set(result)))


def _expand_flex_v2_range(
    prefix: str, start_row: str, start_col: int, end_row: str, end_col: int
) -> list[str]:
    """
    Expand Flex-V2 plate range notation.

    Examples:
        _expand_flex_v2_range("A", "A", 1, "A", 12) → ["A-A01", "A-A02", ..., "A-A12"]
        _expand_flex_v2_range("A", "A", 1, "H", 12) → ["A-A01", ..., "A-A12", "A-B01", ..., "A-H12"]
    """
    rows = FLEX_V2_ROWS
    start_row_idx = rows.index(start_row)
    end_row_idx = rows.index(end_row)

    barcodes = []
    for row_idx in range(start_row_idx, end_row_idx + 1):
        row = rows[row_idx]

        # Determine column range
        if row == start_row and row == end_row:
            # Same row: use specified column range
            col_range = range(start_col, end_col + 1)
        elif row == start_row:
            # First row: start_col to 12
            col_range = range(start_col, 13)
        elif row == end_row:
            # Last row: 1 to end_col
            col_range = range(1, end_col + 1)
        else:
            # Middle rows: all columns (1-12)
            col_range = range(1, 13)

        for col in col_range:
            barcodes.append(f"{prefix}-{row}{col:02d}")

    return barcodes


def _parse_flex_v2_component(component: str) -> list[str]:
    """
    Parse Flex-V2 barcode component.

    Supports:
        - Single: "A-A01" or "A_A01" → ["A-A01"]
        - Range: "A-A01..A-A12" → ["A-A01", "A-A02", ..., "A-A12"]
        - Range: "A-A01..A-H12" → all 96 barcodes in set A
        - Selection: "A-A01|A-A05|A-A09" → ["A-A01", "A-A05", "A-A09"]

    Accepts both hyphen (-) and underscore (_) as separators, normalizes to hyphen.
    """
    # Normalize underscores to hyphens
    component = component.replace("_", "-")

    # Check if it contains range notation (..)
    if ".." in component:
        match = re.match(
            r"^([ABCD])-([ABCDEFGH])(\d{2})\.\.([ABCD])-([ABCDEFGH])(\d{2})$", component
        )
        if not match:
            raise ValueError(f"Invalid Flex-V2 range notation: {component}")

        start_prefix = match.group(1)
        start_row = match.group(2)
        start_col = int(match.group(3))
        end_prefix = match.group(4)
        end_row = match.group(5)
        end_col = int(match.group(6))

        if start_prefix != end_prefix:
            raise ValueError(f"Cannot expand range across different sets: {component}")

        return _expand_flex_v2_range(
            start_prefix, start_row, start_col, end_row, end_col
        )

    # Check if it contains selection notation (|)
    elif "|" in component:
        parts = component.split("|")
        barcodes = []
        for part in parts:
            barcodes.extend(_parse_flex_v2_component(part.strip()))
        return barcodes

    # Single barcode
    else:
        if not _is_flex_v2_barcode(component):
            raise ValueError(f"Invalid Flex-V2 barcode: {component}")
        return [component]


def _expand_barcode_component(component: str) -> list[str]:
    """Expand a barcode component like 'BC1..8' into ['BC001', 'BC002', ..., 'BC008'].

    Also handles Flex-V2 format with hyphen (A-A01) or underscore (A_A01).
    """
    # Check if this is already an explicit barcode (backward compatibility)
    if component in KNOWN_BARCODES:
        return [component]

    # Check if it's Flex-V2 format (contains hyphen or underscore like A-A01 or A_A01)
    if ("-" in component or "_" in component) and component[0] in FLEX_V2_SETS:
        return _parse_flex_v2_component(component)

    # Otherwise, use existing Flex-V1 logic
    # Extract prefix and numeric part
    # Find where the letters end and numbers begin
    i = 0
    while i < len(component) and component[i].isalpha():
        i += 1

    if i == len(component):
        raise ValueError(f"Invalid barcode component format: {component}")

    prefix = component[:i]
    numeric_part = component[i:]

    # Expand the numeric part
    numbers = _expand_selection(numeric_part)

    # Generate barcode names with zero-padding
    result = []
    for num in numbers:
        barcode = f"{prefix}{num:03d}"
        result.append(barcode)

    return result


def _validate_json(data: dict):
    for key in EXPECTED_KEYS:
        if key not in data:
            raise ValueError(f"Missing key in data: {key}")


def _validate_keys(entry: dict):
    if not all(key in entry for key in EXPECTED_SAMPLE_KEYS):
        raise ValueError(f"Missing keys in entry: {entry}")


def _parse_mode(entry: dict) -> list[str]:
    libmode = entry["mode"]
    if "+" in libmode:
        modes = libmode.split("+")
        if not all(mode in KNOWN_LIBMODES for mode in modes):
            raise ValueError(f"Invalid mode found: {libmode}")
        return modes
    else:
        if libmode not in KNOWN_LIBMODES:
            raise ValueError(f"Invalid mode {libmode} found: {libmode}")
        return [libmode]


def _validate_component_barcode(barcode: str):
    if barcode not in KNOWN_BARCODES:
        raise ValueError(f"Invalid barcode found in barcodes: {barcode}")


def _parse_features(entry: dict, nlibs: int, known_features: list[str]) -> list[str]:
    if "+" in entry["features"]:
        features = entry["features"].split("+")
    else:
        features = [entry["features"]]

    if len(features) != nlibs:
        raise ValueError(
            f"Invalid number of features found in features: {entry['features']}. Expected {nlibs} features."
        )

    for f in features:
        if f not in known_features:
            raise ValueError(
                f"Invalid feature found in features: {f}. Missing from provided features: {known_features}"
            )

    return features


def _parse_barcodes(entry: dict, nlib: int) -> list[list[str]]:
    """Parse and validate barcodes in a configuration entry.

    The number of paired barcodes must match the number of libraries.
    Supports both explicit format (BC001+CR001) and range format (BC1..8+CR1..8).
    """
    barcodes = entry["barcodes"]

    # Determine if this is multiple combinations or a single combination with selections
    # Multiple combinations: each part should have the right number of '+' separators
    # Single combination: may have '|' within components for selections
    combinations = []

    if "|" in barcodes:
        # Try splitting on '|' and check if each part looks like a complete combination
        potential_combinations = barcodes.split("|")

        # Check if all parts have the expected number of '+' separators (nlib-1)
        all_complete = all(
            part.count("+") == nlib - 1 for part in potential_combinations
        )

        if all_complete:
            # This is multiple combinations
            combinations = potential_combinations
        else:
            # This is a single combination with '|' used for selections within components
            combinations = [barcodes]
    else:
        combinations = [barcodes]

    pairings = []
    for combination in combinations:
        # Split into components by mode (BC+CR, etc.)
        components = combination.split("+")

        # Special case for Flex-V2: if we have nlib > 1 but only 1 component,
        # and that component is a Flex-V2 barcode, duplicate it for all libraries
        # (Flex-V2 uses single naming scheme for multi-modal)
        if len(components) == 1 and nlib > 1:
            expanded = _expand_barcode_component(components[0])
            if expanded and _is_flex_v2_barcode(expanded[0]):
                # This is a Flex-V2 barcode - use same barcode for all modes
                components = [components[0]] * nlib

        if len(components) != nlib:
            raise ValueError(
                f"Invalid number of barcodes found in barcode combination: {combination}. Expected {nlib} barcodes."
            )

        # Expand each component into explicit barcode lists
        expanded_components = []
        for component in components:
            expanded = _expand_barcode_component(component)
            expanded_components.append(expanded)

        # Validate all components have same length
        lengths = [len(comp) for comp in expanded_components]
        if len(set(lengths)) > 1:
            raise ValueError(
                f"Mismatched component lengths in combination '{combination}': {lengths}"
            )

        # Generate pairings by index
        for i in range(lengths[0]):
            pairing = [comp[i] for comp in expanded_components]
            # Validate each barcode in the pairing
            for barcode in pairing:
                _validate_component_barcode(barcode)
            pairings.append(pairing)

    return pairings


def _pull_feature_path(feature: str, libraries: dict) -> str:
    path = libraries[feature]
    if not os.path.exists(path):
        print(f"Warning: Feature path not found: {path}")
    return path


def _assign_probeset(barcode: str) -> str:
    """
    Assign probe set based on barcode prefix.

    Flex-V1: BC001 → "BC", CR001 → "CR", AB001 → "AB"
    Flex-V2: A-A01 → "A", B-C05 → "B", D-H12 → "D"
    """
    if barcode.startswith("BC"):
        return "BC"
    elif barcode.startswith("CR"):
        return "CR"
    elif barcode.startswith("AB"):
        return "AB"
    elif _is_flex_v2_barcode(barcode):
        return _get_flex_v2_prefix(barcode)
    else:
        raise ValueError(f"Unknown barcode format: {barcode}")


def parse_config(config_path: str):
    """Parse and validate a configuration json file."""
    with open(config_path, "r") as f:
        config = json.load(f)
    _validate_json(config)

    dataframe = []
    for entry in config["samples"]:
        _validate_keys(entry)
        libmode = _parse_mode(entry)
        nlib = len(libmode)
        barcodes = _parse_barcodes(entry, nlib)
        features = _parse_features(entry, nlib, config["libraries"].keys())
        for bc_idx, bc in enumerate(barcodes):
            for mode, bc_component, mode_feature in zip(libmode, bc, features):
                dataframe.append(
                    {
                        "experiment": entry["experiment"],
                        "sample": entry["sample"],
                        "mode": mode,
                        "bc_component": bc_component,
                        "bc_idx": bc_idx,
                        "features": mode_feature,
                        "probe_set": _assign_probeset(bc_component),
                        "feature_path": _pull_feature_path(
                            mode_feature, config["libraries"]
                        ),
                    }
                )

    return pl.DataFrame(dataframe).with_columns(
        expected_prefix=(
            pl.col("experiment") + "_" + pl.col("mode").str.to_uppercase() + "_Lane"
        )
    )


def determine_cyto_runs(sample_sheet: pl.DataFrame) -> pl.DataFrame:
    """Determine the expected cyto run names based on the sample sheet.

    Args:
        sample_sheet: A dataframe containing the sample sheet information.

    Returns:
        A dataframe containing the expected cyto run names.
    """
    return (
        sample_sheet.select(
            ["experiment", "mode", "features", "probe_set", "feature_path"]
        )
        .unique()
        .with_columns(
            (pl.col("experiment") + "_" + pl.col("mode").str.to_uppercase()).alias(
                "expected_prefix"
            )
        )
    )
