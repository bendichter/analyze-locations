# analyze-locations

Tools to scan NWB files on the [DANDI Archive](https://dandiarchive.org) for brain region location values and tag assets with ontology terms from the [Allen Mouse Brain Atlas](http://atlas.brain-map.org/).

## Scripts

### `scan_locations.py`

Scans NWB files across DANDI to catalogue values used in:

- `ImagingPlane.location` (optical physiology)
- `electrodes.location` (extracellular electrophysiology)
- `IntracellularElectrode.location` (intracellular electrophysiology)

```bash
python scan_locations.py [--max-dandisets N] [--max-assets N] [--output FILE] [--no-cache]
```

Results are cached per-dandiset in `scan_cache.jsonl` for resumability. See [FINDINGS.md](FINDINGS.md) for a detailed analysis of the results.

### `label_anatomy.py`

Matches location values against the Allen Common Coordinate Framework (1,327 structures) and updates each asset's `about` metadata field with [`Anatomy`](https://www.dandiarchive.org/handbook/135_metadata/#anatomy) entries using [MBAO identifiers](https://bioregistry.io/registry/mba) (e.g. `MBA_385` for Primary visual area).

```bash
# Dry run (default) — show what would change
python label_anatomy.py

# Target specific dandiset(s)
python label_anatomy.py --dandiset 000017

# Actually write metadata (requires DANDI_API_KEY)
DANDI_API_KEY=your_key python label_anatomy.py --apply --dandiset 000017

# Limit scope
python label_anatomy.py --max-dandisets 10 --max-assets 5
```

**Options:**

| Flag | Description |
|---|---|
| `--dry-run` | Show what would change (default) |
| `--apply` | Write metadata changes (requires `DANDI_API_KEY`) |
| `--dandiset ID [ID ...]` | Target specific dandiset(s) |
| `--max-assets N` | Max NWB files per dandiset |
| `--max-dandisets N` | Max dandisets to process |
| `--refresh-allen` | Force re-fetch Allen CCF data |
| `--no-cache` | Ignore label progress cache |
| `--output FILE` | Summary report file (default: `label_results.json`) |

**How it works:**

1. Downloads and caches the Allen CCF structure graph (1,327 structures with IDs, acronyms, and names)
2. Pre-filters dandisets using `scan_cache.jsonl` to find those with CCF-matching locations
3. Checks species — only processes mouse datasets (`NCBITaxon_10090`)
4. For each NWB asset: extracts locations, matches against Allen CCF (exact then case-insensitive, with extraction from structured strings like `{'area': 'VISp', 'depth': '20'}`)
5. Merges matched `Anatomy` entries into the asset's `about` field (deduplicating by identifier)
6. Tracks progress in `label_cache.jsonl` for resumability

## Installation

```bash
pip install h5py remfile requests tqdm
```

## Findings

See [FINDINGS.md](FINDINGS.md) for a detailed analysis of location values across 786 dandisets and 4,218 NWB files, including Allen CCF match rates, common unmatched patterns, and recommendations.
