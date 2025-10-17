import json
import os

import polars as pl

KNOWN_LIBMODES = ["gex", "crispr", "ab"]
KNOWN_PROBE_SET = ["BC", "CR", "AB"]
KNOWN_BARCODES = [f"{name}0{i:02d}" for i in range(1, 17) for name in KNOWN_PROBE_SET]
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


def _expand_barcode_component(component: str) -> list[str]:
    """Expand a barcode component like 'BC1..8' into ['BC001', 'BC002', ..., 'BC008']."""
    # Check if this is already an explicit barcode (backward compatibility)
    if component in KNOWN_BARCODES:
        return [component]

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
    if os.path.exists(path):
        return path
    else:
        raise FileNotFoundError(f"Feature path not found: {path}")


def _assign_probeset(barcode: str) -> str:
    if barcode.startswith("BC"):
        return "BC"
    elif barcode.startswith("CR"):
        return "CR"
    elif barcode.startswith("AB"):
        return "AB"
    else:
        raise ValueError(f"Invalid barcode format: {barcode}")


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
