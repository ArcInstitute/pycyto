import json
import os

import polars as pl

KNOWN_LIBMODES = ["gex", "crispr", "ab"]
KNOWN_PROBE_SET = ["BC", "CR", "AB"]
KNOWN_BARCODES = [f"{name}0{i:02d}" for i in range(1, 17) for name in KNOWN_PROBE_SET]
EXPECTED_SAMPLE_KEYS = [
    "name",
    "subname",
    "libmode",
    "lanes",
    "barcodes",
]
EXPECTED_KEYS = [
    "libraries",
    "samples",
]


def _validate_json(data: dict):
    for key in EXPECTED_KEYS:
        if key not in data:
            raise ValueError(f"Missing key in data: {key}")


def _validate_keys(entry: dict):
    if not all(key in entry for key in EXPECTED_SAMPLE_KEYS):
        raise ValueError(f"Missing keys in entry: {entry}")


def _parse_libmode(entry: dict) -> list[str]:
    libmode = entry["libmode"]
    if "+" in libmode:
        modes = libmode.split("+")
        if not all(mode in KNOWN_LIBMODES for mode in modes):
            raise ValueError(f"Invalid mode found in libmode: {libmode}")
        return modes
    else:
        if libmode not in KNOWN_LIBMODES:
            raise ValueError(f"Invalid mode {libmode} found in libmode: {libmode}")
        return [libmode]


def _parse_gem_lanes(entry: dict) -> list[int]:
    gem_lanes = entry["lanes"]
    if "|" in gem_lanes:
        lanes = gem_lanes.split("|")
        if not all(lane.isdigit() for lane in lanes):
            raise ValueError(f"Invalid lane found in gem_lanes: {gem_lanes}")
        return [int(lane) for lane in lanes]
    else:
        if not gem_lanes.isdigit():
            raise ValueError(f"Invalid lane found in gem_lanes: {gem_lanes}")
        return [int(gem_lanes)]


def _validate_component_barcode(barcode: str):
    if barcode not in KNOWN_BARCODES:
        raise ValueError(f"Invalid barcode found in barcodes: {barcode}")


def _parse_features(entry: dict, nlibs: int, known_features: list[str]) -> list[str]:
    if "+" in entry["features"]:
        features = entry["features"].split("+")
        if len(features) != nlibs:
            raise ValueError(
                f"Invalid number of features found in features: {entry['features']}. Expected {nlibs} features."
            )
    else:
        features = [entry["features"]]

    for f in features:
        if f not in known_features:
            raise ValueError(
                f"Invalid feature found in features: {f}. Missing from provided features: {known_features}"
            )

    return features


def _parse_barcodes(entry: dict, nlib: int) -> list[list[str]]:
    """Parse and validate barcodes in a configuration entry.

    The number of paired barcodes must match the number of libraries.
    """
    barcodes = entry["barcodes"]
    if "|" in barcodes:
        combinations = barcodes.split("|")
        pairings = [c.split("+") for c in combinations]
    else:
        pairings = barcodes.split("+")

    for p in pairings:
        if len(p) != nlib:
            raise ValueError(
                f"Invalid number of barcodes found in barcode pair: {p}. Expected {nlib} barcodes."
            )
        for component in p:
            _validate_component_barcode(component)
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
        libmode = _parse_libmode(entry)
        nlib = len(libmode)
        gem_lanes = _parse_gem_lanes(entry)
        barcodes = _parse_barcodes(entry, nlib)
        features = _parse_features(entry, nlib, config["libraries"].keys())
        for lane in gem_lanes:
            for bc_idx, bc in enumerate(barcodes):
                for mode, bc_component, mode_feature in zip(libmode, bc, features):
                    dataframe.append(
                        {
                            "name": entry["name"],
                            "subname": entry["subname"],
                            "mode": mode,
                            "lane": lane,
                            "bc_component": bc_component,
                            "bc_idx": bc_idx,
                            "features": mode_feature,
                            "probe_set": _assign_probeset(bc_component),
                            "feature_path": _pull_feature_path(
                                mode_feature, config["libraries"]
                            ),
                        }
                    )

    return pl.DataFrame(dataframe)


def determine_cyto_runs(sample_sheet: pl.DataFrame) -> pl.DataFrame:
    """Determine the expected cyto run names based on the sample sheet.

    Args:
        sample_sheet: A dataframe containing the sample sheet information.

    Returns:
        A dataframe containing the expected cyto run names.
    """
    return (
        sample_sheet.select(
            ["name", "mode", "lane", "features", "probe_set", "feature_path"]
        )
        .unique()
        .with_columns(
            (
                pl.col("name")
                + "_Lane"
                + pl.col("lane").cast(pl.String)
                + "_"
                + pl.col("mode").str.to_uppercase()
            ).alias("expected_prefix")
        )
    )
