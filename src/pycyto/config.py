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
        name = entry["name"]
        subname = entry["subname"]
        libmode = _parse_libmode(entry)
        gem_lanes = _parse_gem_lanes(entry)
        barcodes = _parse_barcodes(entry, len(libmode))
        for lane in gem_lanes:
            for bc_idx, bc in enumerate(barcodes):
                for mode, bc_component in zip(libmode, bc):
                    dataframe.append({
                        "name": name,
                        "subname": subname,
                        "mode": mode,
                        "lane": lane,
                        "bc_component": bc_component,
                        "bc_idx": bc_idx
                    })

    return pl.DataFrame(dataframe)
