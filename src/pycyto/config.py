import json
import polars as pl

KNOWN_LIBMODES = ["gex", "crispr", "ab"]
KNOWN_BARCODES = [f"{name}0{i:02d}" for i in range(1, 17) for name in ["BC", "CR", "AB"]]
EXPECTED_KEYS = [
    "name",
    "subname",
    "libmode",
    "lanes",
    "barcodes",
]

def _validate_keys(entry: dict):
    if not all(key in entry for key in EXPECTED_KEYS):
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
    if not barcode in KNOWN_BARCODES:
        raise ValueError(f"Invalid barcode found in barcodes: {barcode}")

def _parse_features(entry: dict, nlibs: int) -> list[str]:
    if "+" in entry["features"]:
        features = entry["features"].split("+")
        if len(features) != nlibs:
            raise ValueError(f"Invalid number of features found in features: {entry['features']}. Expected {nlibs} features.")
        return features
    else:
        return [entry["features"]]

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
            raise ValueError(f"Invalid number of barcodes found in barcode pair: {p}. Expected {nlib} barcodes.")
        for component in p:
            _validate_component_barcode(component)
    return pairings

def parse_config(config_path: str):
    """Parse and validate a configuration json file."""
    with open(config_path, "r") as f:
        config = json.load(f)

    dataframe = []
    for entry in config:
        _validate_keys(entry)
        libmode = _parse_libmode(entry)
        nlib = len(libmode)
        gem_lanes = _parse_gem_lanes(entry)
        barcodes = _parse_barcodes(entry, nlib)
        features = _parse_features(entry, nlib)
        for lane in gem_lanes:
            for bc_idx, bc in enumerate(barcodes):
                for mode, bc_component, mode_feature in zip(libmode, bc, features):
                    dataframe.append({
                        "name": entry["name"],
                        "subname": entry["subname"],
                        "mode": mode,
                        "lane": lane,
                        "bc_component": bc_component,
                        "bc_idx": bc_idx,
                        "features": mode_feature,
                    })

    return pl.DataFrame(dataframe)

def determine_cyto_runs(sample_sheet: pl.DataFrame) -> pl.DataFrame:
    """Determine the expected cyto run names based on the sample sheet.

    Args:
        sample_sheet: A dataframe containing the sample sheet information.

    Returns:
        A dataframe containing the expected cyto run names.
    """
    return sample_sheet.select(["name", "mode", "lane", "features"]).unique().with_columns(
        (
            pl.col("name")
            + "__" + pl.col("mode")
            + "__W" + pl.col("lane").cast(pl.String)
        ).alias("expected_prefix")
    )
