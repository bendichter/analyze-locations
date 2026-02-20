"""Scan NWB files on the DANDI Archive to collect unique values used for:
  - ImagingPlane.location              (optical physiology)
  - electrodes.location                (extracellular electrophysiology)
  - IntracellularElectrode.location    (intracellular electrophysiology)

Usage:
    python scan_locations.py [--max-dandisets N] [--max-assets N] [--output FILE]
"""

import argparse
import json
import os
from collections import Counter, defaultdict

import h5py
import remfile
import requests
from tqdm import tqdm

DANDI_API = "https://api.dandiarchive.org/api"


def iter_dandisets():
    """Yield all dandiset metadata dicts from the DANDI API."""
    url = f"{DANDI_API}/dandisets/?page_size=200&ordering=-modified"
    while url:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        yield from data["results"]
        url = data.get("next")


def get_nwb_assets(dandiset_id, version="draft", max_assets=5):
    """Return up to *max_assets* NWB asset dicts for a dandiset."""
    assets = []
    url = (
        f"{DANDI_API}/dandisets/{dandiset_id}/versions/{version}"
        f"/assets/?page_size=50&glob=*.nwb"
    )
    while url and len(assets) < max_assets:
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            break
        data = resp.json()
        for asset in data["results"]:
            if asset["path"].endswith(".nwb"):
                assets.append(asset)
                if len(assets) >= max_assets:
                    break
        url = data.get("next")
    return assets


def get_download_url(dandiset_id, asset_id, version="draft"):
    """Build the DANDI download URL for an asset (remfile follows redirects)."""
    return (
        f"{DANDI_API}/dandisets/{dandiset_id}/versions/{version}"
        f"/assets/{asset_id}/download/"
    )


def _read_scalar_or_array(dataset):
    """Read an HDF5 dataset and return a list of decoded strings."""
    raw = dataset[()]
    if hasattr(raw, "__iter__") and not isinstance(raw, (str, bytes)):
        return [v.decode("utf-8", errors="replace") if isinstance(v, bytes) else v
                for v in raw]
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")
    return [raw]


def extract_locations(url):
    """Open an NWB file via HTTP streaming and return location values.

    Returns (imaging_locations, electrode_locations, icephys_locations).
    """
    imaging_locations = []
    electrode_locations = []
    icephys_locations = []

    rf = remfile.File(url)
    with h5py.File(rf, "r") as f:
        # --- ImagingPlane objects ---
        if "general/optophysiology" in f:
            opto = f["general/optophysiology"]
            for name in opto:
                plane = opto[name]
                if isinstance(plane, h5py.Group) and "location" in plane:
                    val = plane["location"][()]
                    if isinstance(val, bytes):
                        val = val.decode("utf-8", errors="replace")
                    imaging_locations.append(val)

        # --- Extracellular electrodes table ---
        if "general/extracellular_ephys/electrodes" in f:
            electrodes = f["general/extracellular_ephys/electrodes"]
            if "location" in electrodes:
                electrode_locations = _read_scalar_or_array(electrodes["location"])

        # --- IntracellularElectrode objects ---
        if "general/intracellular_ephys" in f:
            icephys = f["general/intracellular_ephys"]
            for name in icephys:
                item = icephys[name]
                if isinstance(item, h5py.Group) and "location" in item:
                    val = item["location"][()]
                    if isinstance(val, bytes):
                        val = val.decode("utf-8", errors="replace")
                    icephys_locations.append(val)

    return imaging_locations, electrode_locations, icephys_locations


CACHE_FILE = "scan_cache.jsonl"


def load_cache():
    """Load per-dandiset results from the JSONL cache file."""
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                entry = json.loads(line)
                cache[entry["dandiset_id"]] = entry
    return cache


def append_cache(entry):
    """Append a single dandiset result to the JSONL cache file."""
    with open(CACHE_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-dandisets", type=int, default=None,
        help="Limit the number of dandisets to scan (default: all)",
    )
    parser.add_argument(
        "--max-assets", type=int, default=3,
        help="Max NWB files to check per dandiset (default: 3)",
    )
    parser.add_argument(
        "--output", type=str, default="location_results.json",
        help="Output JSON file (default: location_results.json)",
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Ignore existing cache and rescan everything",
    )
    args = parser.parse_args()

    # Load cache
    cache = {} if args.no_cache else load_cache()
    if cache:
        print(f"Loaded cache with {len(cache)} dandisets already scanned")

    # Counters: value -> count of files containing it
    imaging_counter = Counter()
    electrode_counter = Counter()
    icephys_counter = Counter()

    # Per-dandiset tracking for provenance
    imaging_by_dandiset = defaultdict(set)
    electrode_by_dandiset = defaultdict(set)
    icephys_by_dandiset = defaultdict(set)

    n_files_scanned = 0
    n_files_with_imaging = 0
    n_files_with_electrodes = 0
    n_files_with_icephys = 0
    n_errors = 0

    # Collect dandiset list so tqdm knows the total
    print("Fetching dandiset list …")
    all_dandisets = list(iter_dandisets())
    if args.max_dandisets:
        all_dandisets = all_dandisets[: args.max_dandisets]

    for ds in tqdm(all_dandisets, desc="Dandisets", unit="ds"):
        ds_id = ds["identifier"]

        # Use cached result if available (re-scan if missing icephys field)
        if ds_id in cache and "icephys_locations" in cache[ds_id]:
            entry = cache[ds_id]
            n_files_scanned += entry["files_scanned"]
            n_files_with_imaging += entry["files_with_imaging"]
            n_files_with_electrodes += entry["files_with_electrodes"]
            n_files_with_icephys += entry.get("files_with_icephys", 0)
            n_errors += entry["errors"]
            for loc in entry["imaging_locations"]:
                imaging_counter[loc] += entry["imaging_locations"][loc]
                imaging_by_dandiset[loc].add(ds_id)
            for loc in entry["electrode_locations"]:
                electrode_counter[loc] += entry["electrode_locations"][loc]
                electrode_by_dandiset[loc].add(ds_id)
            for loc in entry.get("icephys_locations", {}):
                icephys_counter[loc] += entry["icephys_locations"][loc]
                icephys_by_dandiset[loc].add(ds_id)
            continue

        assets = get_nwb_assets(ds_id, max_assets=args.max_assets)

        ds_imaging = Counter()
        ds_electrode = Counter()
        ds_icephys = Counter()
        ds_files = 0
        ds_img_files = 0
        ds_elec_files = 0
        ds_ice_files = 0
        ds_errors = 0

        for asset in assets:
            asset_id = asset["asset_id"]
            path = asset["path"]

            url = get_download_url(ds_id, asset_id)

            try:
                img_locs, elec_locs, ice_locs = extract_locations(url)
                ds_files += 1
            except Exception as exc:
                ds_errors += 1
                tqdm.write(f"ERROR {ds_id}/{path}: {exc}")
                continue

            if img_locs:
                ds_img_files += 1
                for loc in img_locs:
                    ds_imaging[loc] += 1

            if elec_locs:
                ds_elec_files += 1
                for loc in set(elec_locs):
                    ds_electrode[loc] += 1

            if ice_locs:
                ds_ice_files += 1
                for loc in ice_locs:
                    ds_icephys[loc] += 1

        # Save to cache
        entry = {
            "dandiset_id": ds_id,
            "files_scanned": ds_files,
            "files_with_imaging": ds_img_files,
            "files_with_electrodes": ds_elec_files,
            "files_with_icephys": ds_ice_files,
            "errors": ds_errors,
            "imaging_locations": dict(ds_imaging),
            "electrode_locations": dict(ds_electrode),
            "icephys_locations": dict(ds_icephys),
        }
        append_cache(entry)
        cache[ds_id] = entry

        # Accumulate into global counters
        n_files_scanned += ds_files
        n_files_with_imaging += ds_img_files
        n_files_with_electrodes += ds_elec_files
        n_files_with_icephys += ds_ice_files
        n_errors += ds_errors
        for loc, count in ds_imaging.items():
            imaging_counter[loc] += count
            imaging_by_dandiset[loc].add(ds_id)
        for loc, count in ds_electrode.items():
            electrode_counter[loc] += count
            electrode_by_dandiset[loc].add(ds_id)
        for loc, count in ds_icephys.items():
            icephys_counter[loc] += count
            icephys_by_dandiset[loc].add(ds_id)

    # Build JSON-serializable results
    n_dandisets = len(all_dandisets)
    results = {
        "summary": {
            "dandisets_scanned": n_dandisets,
            "files_scanned": n_files_scanned,
            "files_with_imaging_plane": n_files_with_imaging,
            "files_with_electrodes": n_files_with_electrodes,
            "files_with_icephys": n_files_with_icephys,
            "errors": n_errors,
        },
        "imaging_plane_locations": {
            loc: {
                "file_count": count,
                "dandiset_count": len(imaging_by_dandiset[loc]),
                "dandisets": sorted(imaging_by_dandiset[loc]),
            }
            for loc, count in imaging_counter.most_common()
        },
        "electrode_locations": {
            loc: {
                "file_count": count,
                "dandiset_count": len(electrode_by_dandiset[loc]),
                "dandisets": sorted(electrode_by_dandiset[loc]),
            }
            for loc, count in electrode_counter.most_common()
        },
        "icephys_locations": {
            loc: {
                "file_count": count,
                "dandiset_count": len(icephys_by_dandiset[loc]),
                "dandisets": sorted(icephys_by_dandiset[loc]),
            }
            for loc, count in icephys_counter.most_common()
        },
    }

    with open(args.output, "w") as fout:
        json.dump(results, fout, indent=2)

    # Print summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Dandisets scanned:          {n_dandisets}")
    print(f"NWB files scanned:          {n_files_scanned}")
    print(f"Files with ImagingPlane:     {n_files_with_imaging}")
    print(f"Files with electrodes:       {n_files_with_electrodes}")
    print(f"Files with icephys:          {n_files_with_icephys}")
    print(f"Errors:                      {n_errors}")
    print()

    for label, counter, by_ds in [
        ("ImagingPlane.location", imaging_counter, imaging_by_dandiset),
        ("electrodes.location", electrode_counter, electrode_by_dandiset),
        ("IntracellularElectrode.location", icephys_counter, icephys_by_dandiset),
    ]:
        if counter:
            print(f"{label} — {len(counter)} unique values:")
            for loc, count in counter.most_common(50):
                ds_count = len(by_ds[loc])
                print(f"  {count:5d} files  {ds_count:4d} dandisets  │ {loc}")
            if len(counter) > 50:
                print(f"  ... and {len(counter) - 50} more (see {args.output})")
            print()

    print(f"Full results saved to {args.output}")


if __name__ == "__main__":
    main()
