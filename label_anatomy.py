"""Tag DANDI assets with brain region ontology terms from the Allen CCF.

Scans NWB files on DANDI, matches brain region location values against the
Allen Mouse Brain Atlas Common Coordinate Framework, and updates each asset's
``about`` metadata field with ``Anatomy`` entries using MBAO identifiers.

Usage:
    python label_anatomy.py [options]

    --dry-run              Show what would change (default behavior)
    --apply                Actually write metadata changes
    --dandiset ID [ID ...] Target specific dandiset(s)
    --max-assets N         Max NWB files per dandiset
    --max-dandisets N      Max dandisets to process
    --refresh-allen        Force re-fetch Allen CCF data
    --no-cache             Ignore label progress cache
    --output FILE          Summary report file (default: label_results.json)
"""

import argparse
import ast
import json
import os
import re
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import h5py
import remfile
import requests
from tqdm import tqdm

from scan_locations import extract_locations, get_download_url, _read_scalar_or_array

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DANDI_API = "https://api.dandiarchive.org/api"
ALLEN_STRUCTURE_GRAPH_URL = (
    "http://api.brain-map.org/api/v2/structure_graph_download/1.json"
)
ALLEN_MAPPING_FILE = "allen_ccf_mapping.json"
LABEL_CACHE_FILE = "label_cache.jsonl"
SCAN_CACHE_FILE = "scan_cache.jsonl"
MBAO_PREFIX = "https://purl.brain-bican.org/ontology/mbao/MBA_"
MOUSE_TAXON = "http://purl.obolibrary.org/obo/NCBITaxon_10090"

TRIVIAL_LOCATIONS = {
    "unknown", "none", "", " ", "n/a", "void", "unspecific",
    "na", "not applicable", "other", "nan",
}

MAX_RETRIES = 3
RETRY_BACKOFF = 2  # seconds, multiplied by attempt number

# ---------------------------------------------------------------------------
# Step 1: Allen CCF Mapping
# ---------------------------------------------------------------------------


def _flatten_allen_tree(node, out):
    """Recursively flatten the Allen structure graph tree into *out*."""
    out.append({
        "id": node["id"],
        "acronym": node["acronym"],
        "name": node["name"],
    })
    for child in node.get("children", []):
        _flatten_allen_tree(child, out)


def fetch_allen_mapping():
    """Download the Allen structure graph and return a flat list of structures."""
    print("Fetching Allen CCF structure graph …")
    resp = requests.get(ALLEN_STRUCTURE_GRAPH_URL, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    # The response has {"success": true, "msg": [...]} where msg[0] is the root
    root = data["msg"][0]
    structures = []
    _flatten_allen_tree(root, structures)
    print(f"  Fetched {len(structures)} structures")
    return structures


def load_or_fetch_allen_mapping(refresh=False):
    """Load the Allen CCF mapping from cache or fetch from API.

    Returns the list of structure dicts.
    """
    if not refresh and os.path.exists(ALLEN_MAPPING_FILE):
        with open(ALLEN_MAPPING_FILE) as f:
            structures = json.load(f)
        print(f"Loaded {len(structures)} Allen CCF structures from {ALLEN_MAPPING_FILE}")
        return structures

    structures = fetch_allen_mapping()
    with open(ALLEN_MAPPING_FILE, "w") as f:
        json.dump(structures, f, indent=2)
    print(f"Cached Allen CCF mapping to {ALLEN_MAPPING_FILE}")
    return structures


def build_lookup_dicts(structures):
    """Build lookup dictionaries from the Allen CCF structures list.

    Returns (by_acronym, by_name, by_acronym_lower, by_name_lower).
    """
    by_acronym = {}
    by_name = {}
    by_acronym_lower = {}
    by_name_lower = {}

    for s in structures:
        by_acronym[s["acronym"]] = s
        by_name[s["name"]] = s
        acr_lower = s["acronym"].lower()
        name_lower = s["name"].lower()
        # First match wins for case-insensitive (shouldn't have collisions)
        if acr_lower not in by_acronym_lower:
            by_acronym_lower[acr_lower] = s
        if name_lower not in by_name_lower:
            by_name_lower[name_lower] = s

    return by_acronym, by_name, by_acronym_lower, by_name_lower


# ---------------------------------------------------------------------------
# Step 2: Location Matching
# ---------------------------------------------------------------------------


def _extract_area(location):
    """Try to extract a brain-area value from structured location strings.

    Handles two common formats found in NWB files:
      - "{'area': 'VISp', 'depth': '20'}"   (Python dict repr)
      - "area: VISp,depth: 175"              (key: value pairs)

    Returns the area string if found, or None.
    """
    # Python dict repr:  {'area': 'VISp', 'depth': '20'}
    if location.startswith("{"):
        try:
            d = ast.literal_eval(location)
            if isinstance(d, dict) and "area" in d:
                return str(d["area"]).strip()
        except (ValueError, SyntaxError):
            pass

    # Comma-separated key: value pairs:  area: VISp,depth: 175
    m = re.match(r"area:\s*([^,]+)", location)
    if m:
        return m.group(1).strip()

    return None


def _match_single(loc, by_acronym, by_name, by_acronym_lower, by_name_lower):
    """Match a single token against Allen CCF structures. Returns a structure dict or None."""
    if loc in by_acronym:
        return by_acronym[loc]
    if loc in by_name:
        return by_name[loc]
    loc_lower = loc.lower()
    if loc_lower in by_acronym_lower:
        return by_acronym_lower[loc_lower]
    if loc_lower in by_name_lower:
        return by_name_lower[loc_lower]
    return None


def match_location(location, by_acronym, by_name, by_acronym_lower, by_name_lower):
    """Match a location string against Allen CCF structures.

    Returns a list of matched structure dicts (may be empty).
    Handles plain values, structured strings, and comma-separated lists.
    """
    loc = location.strip()
    if loc.lower() in TRIVIAL_LOCATIONS:
        return []

    # Try as a single value first
    result = _match_single(loc, by_acronym, by_name, by_acronym_lower, by_name_lower)
    if result:
        return [result]

    # Try extracting an area value from structured strings
    area = _extract_area(loc)
    if area:
        return match_location(area, by_acronym, by_name, by_acronym_lower, by_name_lower)

    # Try comma-separated list (e.g. "VISp,VISrl,VISlm,VISal")
    if "," in loc:
        parts = [p.strip() for p in loc.split(",")]
        # Only treat as a list if at least one part matches — avoids splitting
        # structured strings like "area: VISp,depth: 175" (already handled above)
        matches = []
        for part in parts:
            if part and part.lower() not in TRIVIAL_LOCATIONS:
                s = _match_single(part, by_acronym, by_name, by_acronym_lower, by_name_lower)
                if s:
                    matches.append(s)
        if matches:
            return matches

    return []


def structure_to_anatomy(structure):
    """Convert an Allen CCF structure dict to a DANDI Anatomy entry."""
    return {
        "schemaKey": "Anatomy",
        "name": structure["name"],
        "identifier": f"{MBAO_PREFIX}{structure['id']}",
    }


# ---------------------------------------------------------------------------
# Step 3: Scan Cache Pre-filtering
# ---------------------------------------------------------------------------


def load_scan_cache():
    """Load scan_cache.jsonl and return dict of dandiset_id -> entry."""
    cache = {}
    if not os.path.exists(SCAN_CACHE_FILE):
        return cache
    with open(SCAN_CACHE_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            cache[entry["dandiset_id"]] = entry
    return cache


def get_candidate_dandisets(scan_cache, lookups):
    """Filter scan_cache to dandisets that have at least one CCF-matching location.

    Returns dict of dandiset_id -> set of matching location strings.
    """
    by_acronym, by_name, by_acronym_lower, by_name_lower = lookups
    candidates = {}

    for ds_id, entry in scan_cache.items():
        all_locations = set()
        for loc_dict in [
            entry.get("imaging_locations", {}),
            entry.get("electrode_locations", {}),
            entry.get("icephys_locations", {}),
        ]:
            all_locations.update(loc_dict.keys())

        matching = set()
        for loc in all_locations:
            if match_location(loc, by_acronym, by_name, by_acronym_lower, by_name_lower):  # non-empty list = truthy
                matching.add(loc)

        if matching:
            candidates[ds_id] = matching

    return candidates


# ---------------------------------------------------------------------------
# Label Cache (per-asset progress tracking)
# ---------------------------------------------------------------------------


def load_label_cache():
    """Load label_cache.jsonl and return dict of (dandiset_id, asset_id) -> entry."""
    cache = {}
    if not os.path.exists(LABEL_CACHE_FILE):
        return cache
    with open(LABEL_CACHE_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            key = (entry["dandiset_id"], entry["asset_id"])
            cache[key] = entry
    return cache


def append_label_cache(entry):
    """Append a single asset result to label_cache.jsonl."""
    with open(LABEL_CACHE_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")


# ---------------------------------------------------------------------------
# DANDI API helpers
# ---------------------------------------------------------------------------


def _request_with_retry(method, url, **kwargs):
    """Make an HTTP request with retry on 429 and 5xx errors."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = method(url, **kwargs)
            if resp.status_code == 429 or resp.status_code >= 500:
                if attempt < MAX_RETRIES:
                    wait = RETRY_BACKOFF * attempt
                    tqdm.write(f"  Retry {attempt}/{MAX_RETRIES} after {resp.status_code}, waiting {wait}s …")
                    time.sleep(wait)
                    continue
            return resp
        except requests.RequestException as exc:
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF * attempt
                tqdm.write(f"  Retry {attempt}/{MAX_RETRIES} after error: {exc}, waiting {wait}s …")
                time.sleep(wait)
            else:
                raise
    return resp  # Should not reach here, but just in case


def check_species_mouse(dandiset_id):
    """Check if a dandiset's assetsSummary.species includes Mus musculus."""
    url = f"{DANDI_API}/dandisets/{dandiset_id}/versions/draft/"
    resp = _request_with_retry(requests.get, url, timeout=30)
    if resp.status_code != 200:
        return False
    data = resp.json()
    species_list = data.get("asset_count") and data  # just check response is valid
    # Navigate to assetsSummary.species
    assets_summary = data.get("metadata", data).get("assetsSummary", {})
    species = assets_summary.get("species", [])
    for sp in species:
        identifier = sp.get("identifier", "")
        if "NCBITaxon_10090" in identifier or "10090" in identifier:
            return True
    return False


def get_nwb_assets_paged(dandiset_id, version="draft", max_assets=None):
    """Yield NWB asset dicts for a dandiset, with pagination."""
    url = (
        f"{DANDI_API}/dandisets/{dandiset_id}/versions/{version}"
        f"/assets/?page_size=100&glob=*.nwb"
    )
    count = 0
    while url:
        resp = _request_with_retry(requests.get, url, timeout=30)
        if resp.status_code != 200:
            break
        data = resp.json()
        for asset in data["results"]:
            if asset["path"].endswith(".nwb"):
                yield asset
                count += 1
                if max_assets and count >= max_assets:
                    return
        url = data.get("next")


def get_asset_metadata(dandiset_id, asset_id, version="draft"):
    """Get raw metadata for a specific asset."""
    url = f"{DANDI_API}/dandisets/{dandiset_id}/versions/{version}/assets/{asset_id}/"
    resp = _request_with_retry(requests.get, url, timeout=30)
    resp.raise_for_status()
    return resp.json()


def set_asset_metadata(dandiset_id, asset_id, metadata, api_key, version="draft"):
    """Update raw metadata for a specific asset.

    Uses the DANDI API PUT endpoint to replace asset metadata.
    """
    url = f"{DANDI_API}/dandisets/{dandiset_id}/versions/{version}/assets/{asset_id}/"
    headers = {
        "Authorization": f"token {api_key}",
        "Content-Type": "application/json",
    }
    resp = _request_with_retry(
        requests.put, url, headers=headers, json={"metadata": metadata}, timeout=30
    )
    resp.raise_for_status()
    return resp.json()


# ---------------------------------------------------------------------------
# Step 4: Asset Processing
# ---------------------------------------------------------------------------


def process_asset(dandiset_id, asset, lookups, apply=False, api_key=None):
    """Process a single NWB asset: extract locations, match, optionally update.

    Returns a result dict for the label cache.
    """
    asset_id = asset["asset_id"]
    path = asset["path"]
    by_acronym, by_name, by_acronym_lower, by_name_lower = lookups

    result = {
        "dandiset_id": dandiset_id,
        "asset_id": asset_id,
        "path": path,
        "status": "unknown",
        "matched_locations": {},
        "unmatched_locations": [],
        "anatomy_entries": [],
        "error": None,
    }

    # Extract locations from NWB file
    try:
        download_url = get_download_url(dandiset_id, asset_id)
        img_locs, elec_locs, ice_locs = extract_locations(download_url)
    except Exception as exc:
        result["status"] = "error"
        result["error"] = str(exc)
        return result

    # Collect unique locations
    all_locations = set()
    for loc in img_locs:
        all_locations.add(loc)
    for loc in elec_locs:
        all_locations.add(loc)
    for loc in ice_locs:
        all_locations.add(loc)

    if not all_locations:
        result["status"] = "skipped_no_locations"
        return result

    # Match against Allen CCF
    anatomy_entries = []
    seen_ids = set()
    for loc in sorted(all_locations):
        structures = match_location(loc, by_acronym, by_name, by_acronym_lower, by_name_lower)
        if structures:
            result["matched_locations"][loc] = [
                {"id": s["id"], "acronym": s["acronym"], "name": s["name"]}
                for s in structures
            ]
            for structure in structures:
                if structure["id"] not in seen_ids:
                    seen_ids.add(structure["id"])
                    anatomy_entries.append(structure_to_anatomy(structure))
        else:
            if loc.strip().lower() not in TRIVIAL_LOCATIONS:
                result["unmatched_locations"].append(loc)

    if not anatomy_entries:
        result["status"] = "skipped_no_match"
        return result

    result["anatomy_entries"] = anatomy_entries

    # Read-modify-write asset metadata
    try:
        metadata = get_asset_metadata(dandiset_id, asset_id)
    except Exception as exc:
        result["status"] = "error"
        result["error"] = f"Failed to get metadata: {exc}"
        return result

    existing_about = metadata.get("about") or []
    existing_ids = {e.get("identifier") for e in existing_about}

    new_entries = [
        e for e in anatomy_entries if e["identifier"] not in existing_ids
    ]

    if not new_entries:
        result["status"] = "skipped_no_change"
        return result

    result["new_entries"] = new_entries

    if apply:
        merged_about = existing_about + new_entries
        metadata["about"] = merged_about
        try:
            set_asset_metadata(dandiset_id, asset_id, metadata, api_key)
            result["status"] = "updated"
        except Exception as exc:
            result["status"] = "error"
            result["error"] = f"Failed to set metadata: {exc}"
    else:
        result["status"] = "would_update"

    return result


# ---------------------------------------------------------------------------
# Step 5: CLI & Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Tag DANDI assets with Allen CCF brain region ontology terms."
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=True,
        help="Show what would change (default behavior)",
    )
    parser.add_argument(
        "--apply", action="store_true",
        help="Actually write metadata changes",
    )
    parser.add_argument(
        "--dandiset", nargs="+", metavar="ID",
        help="Target specific dandiset(s)",
    )
    parser.add_argument(
        "--max-assets", type=int, default=None,
        help="Max NWB files per dandiset",
    )
    parser.add_argument(
        "--max-dandisets", type=int, default=None,
        help="Max dandisets to process",
    )
    parser.add_argument(
        "--refresh-allen", action="store_true",
        help="Force re-fetch Allen CCF data",
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Ignore label progress cache",
    )
    parser.add_argument(
        "--output", type=str, default="label_results.json",
        help="Summary report file (default: label_results.json)",
    )
    parser.add_argument(
        "--workers", type=int, default=8,
        help="Number of parallel workers for NWB streaming (default: 8)",
    )
    args = parser.parse_args()

    apply = args.apply
    api_key = os.environ.get("DANDI_API_KEY")

    # Auth check
    if apply and not api_key:
        print("ERROR: --apply requires DANDI_API_KEY environment variable", file=sys.stderr)
        sys.exit(1)

    # Step 1: Load Allen CCF mapping
    structures = load_or_fetch_allen_mapping(refresh=args.refresh_allen)
    lookups = build_lookup_dicts(structures)
    by_acronym, by_name, by_acronym_lower, by_name_lower = lookups
    print(f"  {len(by_acronym)} acronyms, {len(by_name)} names")

    # Step 2: Identify candidate dandisets from scan_cache
    if args.dandiset:
        # User specified specific dandisets — use those directly
        target_dandisets = args.dandiset
        print(f"\nTargeting {len(target_dandisets)} specified dandiset(s): {', '.join(target_dandisets)}")
    else:
        scan_cache = load_scan_cache()
        if not scan_cache:
            print(f"\nNo {SCAN_CACHE_FILE} found. Run scan_locations.py first or specify --dandiset.")
            sys.exit(1)
        print(f"\nLoaded {len(scan_cache)} dandisets from {SCAN_CACHE_FILE}")

        candidates = get_candidate_dandisets(scan_cache, lookups)
        print(f"Found {len(candidates)} dandisets with CCF-matching locations")

        target_dandisets = sorted(candidates.keys())

    if args.max_dandisets:
        target_dandisets = target_dandisets[:args.max_dandisets]

    # Load label cache
    label_cache = {} if args.no_cache else load_label_cache()
    if label_cache:
        print(f"Loaded label cache with {len(label_cache)} asset results")

    # Step 3: Process each dandiset
    summary = {
        "dandisets_processed": 0,
        "dandisets_skipped_species": 0,
        "assets_processed": 0,
        "assets_cached": 0,
        "assets_updated": 0,
        "assets_would_update": 0,
        "assets_no_change": 0,
        "assets_no_match": 0,
        "assets_no_locations": 0,
        "assets_error": 0,
    }
    all_results = []
    targeted_mode = bool(args.dandiset)
    cache_lock = threading.Lock()

    def _log_result(result, prefix=""):
        """Print details for a processed asset and update summary counters."""
        status = result["status"]
        tqdm.write(f"  {prefix}-> {status}")
        if result.get("matched_locations"):
            for loc, infos in result["matched_locations"].items():
                names = ", ".join(f"{i['name']} (MBA_{i['id']})" for i in infos)
                tqdm.write(f"     {prefix}matched: {loc!r} -> {names}")
        if result.get("unmatched_locations"):
            for loc in result["unmatched_locations"]:
                tqdm.write(f"     {prefix}unmatched: {loc!r}")
        if result.get("new_entries"):
            for entry in result["new_entries"]:
                tqdm.write(f"     {prefix}+ {entry['name']} ({entry['identifier']})")
        if result.get("error"):
            tqdm.write(f"     {prefix}ERROR: {result['error']}")

    def _record_result(result):
        """Update summary counters (call under cache_lock)."""
        status = result["status"]
        summary["assets_processed"] += 1
        if status == "updated":
            summary["assets_updated"] += 1
        elif status == "would_update":
            summary["assets_would_update"] += 1
        elif status == "skipped_no_change":
            summary["assets_no_change"] += 1
        elif status == "skipped_no_match":
            summary["assets_no_match"] += 1
        elif status == "skipped_no_locations":
            summary["assets_no_locations"] += 1
        elif status == "error":
            summary["assets_error"] += 1

    def _process_and_record(ds_id, asset, pbar):
        """Process one asset, log, and persist. Thread-safe."""
        asset_id = asset["asset_id"]
        path = asset["path"]
        cache_key = (ds_id, asset_id)

        with cache_lock:
            if cache_key in label_cache:
                cached = label_cache[cache_key]
                tqdm.write(f"  [cached] {ds_id}/{path}: {cached['status']}")
                all_results.append(cached)
                summary["assets_cached"] += 1
                pbar.update(1)
                return

        result = process_asset(ds_id, asset, lookups, apply=apply, api_key=api_key)

        with cache_lock:
            tqdm.write(f"  {ds_id}/{path}")
            _log_result(result, prefix="  ")
            _record_result(result)
            append_label_cache(result)
            label_cache[cache_key] = result
            all_results.append(result)
            pbar.update(1)

    # Collect all (ds_id, asset) work items, filtering by species
    work_items = []  # list of (ds_id, asset)

    if targeted_mode:
        print("Checking species and collecting assets …")
        for ds_id in target_dandisets:
            if not check_species_mouse(ds_id):
                print(f"  {ds_id}: skipping — not a mouse dataset")
                summary["dandisets_skipped_species"] += 1
                continue
            summary["dandisets_processed"] += 1
            for asset in get_nwb_assets_paged(ds_id, max_assets=args.max_assets):
                work_items.append((ds_id, asset))
    else:
        for ds_id in tqdm(target_dandisets, desc="Filtering", unit="ds"):
            if not check_species_mouse(ds_id):
                tqdm.write(f"  {ds_id}: skipping — not a mouse dataset")
                summary["dandisets_skipped_species"] += 1
                continue
            summary["dandisets_processed"] += 1
            for asset in get_nwb_assets_paged(ds_id, max_assets=args.max_assets):
                work_items.append((ds_id, asset))

    print(f"\n{len(work_items)} assets to process across {summary['dandisets_processed']} dandisets "
          f"({args.workers} workers)")

    with tqdm(total=len(work_items), desc="Assets", unit="asset") as pbar:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(_process_and_record, ds_id, asset, pbar)
                for ds_id, asset in work_items
            ]
            for future in as_completed(futures):
                exc = future.exception()
                if exc:
                    tqdm.write(f"  Worker error: {exc}")

    # Print summary
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    mode = "APPLY" if apply else "DRY RUN"
    print(f"Mode:                    {mode}")
    print(f"Dandisets processed:     {summary['dandisets_processed']}")
    print(f"Dandisets skipped (spp): {summary['dandisets_skipped_species']}")
    print(f"Assets processed:        {summary['assets_processed']}")
    print(f"Assets from cache:       {summary['assets_cached']}")
    if apply:
        print(f"Assets updated:          {summary['assets_updated']}")
    else:
        print(f"Assets would update:     {summary['assets_would_update']}")
    print(f"Assets no change:        {summary['assets_no_change']}")
    print(f"Assets no CCF match:     {summary['assets_no_match']}")
    print(f"Assets no locations:     {summary['assets_no_locations']}")
    print(f"Assets with errors:      {summary['assets_error']}")

    # Save report
    report = {
        "summary": summary,
        "results": all_results,
    }
    with open(args.output, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nReport saved to {args.output}")


if __name__ == "__main__":
    main()
