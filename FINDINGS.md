# Brain Region Location Values in NWB Files on DANDI

## Overview

We scanned **786 dandisets** on the [DANDI Archive](https://dandiarchive.org),
reading up to 10 NWB files per dandiset (**4,218 files total, 0 errors**), to
catalogue the values used for:

- **`ImagingPlane.location`** (optical physiology)
- **`electrodes.location`** column (extracellular electrophysiology)
- **`IntracellularElectrode.location`** (intracellular electrophysiology)

The goal was to understand how the community currently populates these fields
and how many values already conform to a controlled vocabulary such as the
[Allen Common Coordinate Framework (CCF)](http://atlas.brain-map.org/).

## Summary Statistics

| Metric | ImagingPlane | Electrodes | Intracellular | Combined |
|---|---|---|---|---|
| Files with field present | 880 | 1,510 | 45 | 2,435 |
| Unique location values | 131 | 664 | 15 | 794 |
| Values matching Allen CCF | 14 | 363 | 0 | 368 |
| % of unique values matched | 10.7% | 54.7% | 0.0% | **46.3%** |

## Allen CCF Matching

We compared every unique location value against the 1,327 structures in the
Allen Mouse Brain Atlas ontology (ontology ID 1, fetched from the Allen
Institute API). Matching was done on exact acronym, exact full name, and
case-insensitive variants of both.

**368 out of 794 unique values (46.3%) map directly to the Allen CCF.**

The electrode table has much better coverage (54.7%) than the imaging plane
field (10.7%), likely because Neuropixels and other high-density probe
experiments tend to use Allen CCF acronyms for per-channel registration.
Intracellular electrode locations have 0% match rate because they overwhelmingly
use free-text descriptions with embedded metadata.

### Top Allen CCF-Matched Values

| Value | Files | Dandisets | Allen CCF Structure |
|---|---|---|---|
| `CA1` | 144 | 25 | Field CA1 |
| `root` | 71 | 10 | root (top of hierarchy) |
| `MB` | 45 | 12 | Midbrain |
| `CA3` | 36 | 10 | Field CA3 |
| `DG-mo` | 36 | 8 | Dentate gyrus, molecular layer |
| `TH` | 32 | 10 | Thalamus |
| `VISp` | 32 | 8 | Primary visual area |
| `DG-sg` | 29 | 8 | Dentate gyrus, granule cell layer |
| `LP` | 26 | 10 | Lateral posterior nucleus of thalamus |
| `MOs5` | 26 | 5 | Secondary motor area, layer 5 |

## Values That Do NOT Match the Allen CCF

The 426 unmatched values (53.7%) fall into several categories:

### 1. Placeholders and Missing Data

These are clearly placeholder or null-like values and carry no anatomical
information.

| Value | Files | Dandisets |
|---|---|---|
| `unknown` | 482 | 67 |
| `none` / `None` | 134 | 17 |
| `""` (empty string) | 11 | 4 |
| `" "` (whitespace) | 32 | 5 |
| `N/A` / `n.a.` | 25 | 4 |
| `void` | 3 | 1 |
| `unspecific` | 20 | 2 |

### 2. Common Abbreviations Not in the Allen CCF

These are widely-used neuroscience abbreviations that do not appear in the
Allen ontology, either because they are informal, from a different atlas, or
refer to composite/functional areas.

| Value | Files | Dandisets | Likely Meaning |
|---|---|---|---|
| `M1` | 40 | 8 | Primary motor cortex |
| `PFC` | 37 | 5 | Prefrontal cortex |
| `ALM` | 35 | 6 | Anterior lateral motor cortex |
| `V1` | 74 | 4 | Primary visual area (= `VISp` in Allen CCF) |
| `MEC` | 3 | 1 | Medial entorhinal cortex |
| `NAc` | 5 | 3 | Nucleus accumbens (= `ACB` in Allen CCF) |
| `DLS` / `DMS` | 6 | 2 | Dorsolateral / dorsomedial striatum |

### 3. Human Brain Regions

Human clinical electrode recordings (e.g. epilepsy monitoring) use region
names from human atlases that have no Allen CCF equivalent.

| Value | Files | Dandisets |
|---|---|---|
| `hippocampus_left` / `hippocampus_right` | 70 | 5 |
| `amygdala_left` / `amygdala_right` | 73 | 5 |
| `dorsal_anterior_cingulate_cortex_right` | 19 | 3 |
| `pre_supplementary_motor_area_left` | 17 | 3 |
| `Hipp, Right Hippocampus cHipp, caudal hippocampus` | 20 | 2 |
| `MTG, Right Middle Temporal Gyrus A21c, caudal area 21` | 18 | 2 |

### 4. Non-Mouse Species

Many NWB files come from non-mouse organisms where the Allen CCF does not
apply.

| Value | Files | Dandisets | Species |
|---|---|---|---|
| `Hindbrain RSSystem and ncMLF` | 350 | 1 | Zebrafish |
| `Ventral Midbrain` / `Dorsal Midbrain` | 540 | 1 | Zebrafish |
| `Whole Worm` / `head` / `Head` | 142 | 9 | C. elegans / Drosophila |
| `right antennal lobe` | 20 | 2 | Drosophila |
| `lumbar spinal cord` | 10 | 1 | Various |

### 5. Free-Text Descriptions and Composite Values

Some labs embed extra metadata (depth, hemisphere, coordinates) into the
location string rather than using separate fields. This pattern is especially
prevalent in intracellular electrophysiology data, where **all 15 unique
values** are either free-text or structured strings.

| Value | Files | Dandisets |
|---|---|---|
| `area: VISp,depth: 175` | 7 | 2 |
| `{'area': 'VISp', 'depth': '175'}` | 4 | 2 |
| `VISp,VISrl,VISlm,VISal` | 80 | 1 |
| `Left ALM, Right ALM` | 3 | 1 |
| `AP:-4.2mm,ML:+3.0mm,DV:-2.2mm` | 20 | 2 |
| `Primary motor cortex (M1), area 4, arm area` | 3 | 1 |
| `{"subject_id": 462458, ...}` (JSON blob) | 24 | 1 |

### 6. Errors

| Value | Files | Dandisets | Issue |
|---|---|---|---|
| `bpy.context.scene.location` | 2 | 2 | Blender Python expression, likely a conversion bug |

## IntracellularElectrode.location

Intracellular electrophysiology data was found in **45 files** across **4
dandisets** (000005, 000013, 000025, 001750/001751/001752). The 15 unique
location values show two distinct patterns:

1. **`Soma`** (100 files, 3 dandisets) — the most common value, referring to
   somatic whole-cell recordings. This is a cell compartment, not a brain
   region.

2. **Structured free-text** (13 values, 1 dandiset each) — semicolon-delimited
   key-value strings encoding brain region, subregion, cortical layer,
   hemisphere, coordinates, and depth all in a single field:
   - `brain_region: barrel cortex; brain_subregion: N/A; cortical_layer: 4; hemisphere: left; coordinate_ref: bregma; coordinate_ap: 1.80; coordinate_ml: 3.50; coordinate_dv: 0.48`
   - `brain_region: barrel; brain_subregion: C2; cortical_layer: 5; hemisphere: left; brain_location_full_name: S1, barrel field, Paxinos; depth: 637.00`

3. **JSON blob** (1 file, 1 dandiset) — a full JSON object embedded in the
   location string with keys like `Anatomical_Location`, `Brain_Hemisphere`,
   `Brain_Area`, and cell position coordinates.

None of these values match the Allen CCF directly, highlighting the need for
structured metadata fields for intracellular recordings.

## Automated Labeling Dry Run

Using `label_anatomy.py`, we performed a dry run across all candidate
dandisets — those with at least one Allen CCF-matching location value in
`scan_cache.jsonl`. The tool streamed up to 10 NWB files per dandiset,
matched location values against the Allen CCF, and identified which assets
could be tagged with MBAO `Anatomy` entries.

### Candidate Filtering

Of the 786 scanned dandisets, **68** had at least one location value matching
the Allen CCF. Of those, **48** are mouse datasets (`NCBITaxon_10090`) and
**20** were skipped as non-mouse species.

### Asset-Level Results

| Status | Assets |
|---|---|
| Would update (has new anatomy terms) | 411 |
| No CCF match (unmatched locations only) | 24 |
| No location data in file | 15 |
| Errors | 0 |
| **Total assets processed** | **450** |

**91% of assets** (411 / 450) across all 48 mouse dandisets would receive at
least one new `Anatomy` entry.

### Dandiset-Level Resolution Rates

Among the 48 mouse dandisets processed, we checked what proportion have at
least one asset with at least one resolvable location, broken down by
location type:

| Location type | Dandisets with field | Dandisets with resolvable location | % |
|---|---|---|---|
| `ImagingPlane.location` | 24 | 23 | **96%** |
| `electrodes.location` | 25 | 25 | **100%** |
| `IntracellularElectrode.location` | 0 | 0 | N/A |
| **Any location type** | **48** | **48** | **100%** |

Every mouse dandiset with CCF-candidate locations had at least one
resolvable asset. Electrode locations resolved at 100%, and imaging plane
locations at 96% (the one unresolved dandiset uses non-standard location
strings). No mouse intracellular electrophysiology dandisets had
CCF-matching locations in the scan cache.

### Top Matched Structures

| Structure | Assets | Dandisets |
|---|---|---|
| Field CA1 | 153 | 23 |
| Primary visual area (VISp) | 97 | 16 |
| root | 81 | 10 |
| Midbrain (MB) | 55 | 12 |
| Field CA3 | 46 | 11 |
| Dentate gyrus, molecular layer (DG-mo) | 40 | 8 |
| Thalamus (TH) | 37 | 11 |
| Striatum (STR) | 37 | 6 |
| Lateral posterior nucleus (LP) | 34 | 11 |
| Dentate gyrus, granule cell layer (DG-sg) | 33 | 8 |

### Notable Unmatched Values

These are non-trivial location strings in mouse dandisets that did not match
the Allen CCF, even after structured-string extraction and comma splitting:

| Value | Assets | Dandisets | Notes |
|---|---|---|---|
| `RSC` | 11 | 2 | Informal abbreviation for retrosplenial cortex (= `RSP` in Allen CCF) |
| `ca1-pyr`, `ca1-so`, `ca1-sr`, `ca1-slm` | 10 each | 1 | CA1 sublayers — not in Allen CCF |
| `dg-mol`, `dg-val gc`, `dg-val hil` | 10 each | 1 | DG sublayers — not in Allen CCF |
| `PFC` | 10 | 1 | Prefrontal cortex — composite, not a single Allen structure |
| `No Area` | 10 | 1 | Placeholder |
| `Entorhinal area medial part dorsal zone` | 10 | 1 | Close to Allen name but not exact |
| `Above brain` | 10 | 1 | Non-anatomical (electrode outside brain) |
| `MZMG` | 4 | 1 | Non-standard abbreviation |
| `PPC, right hemisphere` | 2 | 1 | Posterior parietal cortex with hemisphere suffix |
| `PM, right hemisphere` | 8 | 1 | Ambiguous — matched as principal mammillary tract (likely incorrect; probably means posterior medial cortex) |

## Recommendations

1. **Adopt a controlled vocabulary.** For mouse data, the Allen CCF ontology
   provides excellent coverage. Common aliases (`V1` -> `VISp`, `M1` -> `MOp`,
   `NAc` -> `ACB`) could be mapped automatically during validation.

2. **Support multiple atlases.** The Allen CCF only covers mouse. Human data
   needs a different ontology (e.g. the
   [Brainnetome Atlas](https://atlas.brainnetome.org/) or
   [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) parcellations). NWB could
   allow specifying which atlas a location value refers to.

3. **Separate concerns.** Depth, hemisphere, coordinates, and region should be
   stored in separate fields rather than packed into the location string. This
   is especially critical for intracellular electrophysiology, where complex
   structured strings are used as a workaround for missing metadata fields.

4. **Validate on upload.** DANDI could warn when a location value does not
   match any known atlas, reducing placeholder values like `unknown` (present
   in 67 dandisets) and catching errors like the Blender expression.

5. **Distinguish cell compartment from brain region.** For intracellular
   recordings, `Soma` (the most common value) refers to the recording site on
   the cell, not the brain region. These are separate concepts and could
   benefit from distinct fields.

## Reproducibility

The location scan was performed on 2025-02-19 using `scan_locations.py` in
this repository. Results are cached in `scan_cache.jsonl` (per-dandiset) and
aggregated in `location_results_all.json`. The Allen CCF ontology was fetched
from `http://api.brain-map.org/api/v2/structure_graph_download/1.json`
(1,327 structures).

The labeling dry run was performed on 2026-02-20 using `label_anatomy.py`
with `--max-assets 10 --no-cache`. Full results are in
`label_results_full.json`.
