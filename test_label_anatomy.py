"""Tests for label_anatomy.py."""

import json
import os
import tempfile

import pytest

from label_anatomy import (
    _extract_area,
    _flatten_allen_tree,
    _match_single,
    build_lookup_dicts,
    get_candidate_dandisets,
    load_label_cache,
    append_label_cache,
    match_location,
    structure_to_anatomy,
    LABEL_CACHE_FILE,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_STRUCTURES = [
    {"id": 385, "acronym": "VISp", "name": "Primary visual area"},
    {"id": 394, "acronym": "VISam", "name": "Anteromedial visual area"},
    {"id": 409, "acronym": "VISrl", "name": "Rostrolateral visual area"},
    {"id": 402, "acronym": "VISal", "name": "Anterolateral visual area"},
    {"id": 549, "acronym": "TH", "name": "Thalamus"},
    {"id": 382, "acronym": "CA1", "name": "Field CA1"},
    {"id": 463, "acronym": "CA3", "name": "Field CA3"},
    {"id": 726, "acronym": "DG", "name": "Dentate gyrus"},
    {"id": 997, "acronym": "root", "name": "root"},
    {"id": 313, "acronym": "MB", "name": "Midbrain"},
    {"id": 242, "acronym": "LS", "name": "Lateral septal nucleus"},
    {"id": 500, "acronym": "RSP", "name": "Retrosplenial area"},
    {"id": 31, "acronym": "ACA", "name": "Anterior cingulate area"},
    {"id": 44, "acronym": "VISlm", "name": "Lateral-medial visual area"},  # fictitious id for testing
]


@pytest.fixture
def lookups():
    return build_lookup_dicts(SAMPLE_STRUCTURES)


@pytest.fixture
def lookup_args(lookups):
    """Unpack lookups tuple for passing to match functions."""
    return lookups


# ---------------------------------------------------------------------------
# _flatten_allen_tree
# ---------------------------------------------------------------------------


class TestFlattenAllenTree:
    def test_single_node(self):
        node = {"id": 1, "acronym": "root", "name": "root"}
        out = []
        _flatten_allen_tree(node, out)
        assert out == [{"id": 1, "acronym": "root", "name": "root"}]

    def test_nested(self):
        tree = {
            "id": 1, "acronym": "root", "name": "root",
            "children": [
                {"id": 2, "acronym": "A", "name": "Area A", "children": [
                    {"id": 3, "acronym": "A1", "name": "Area A1"},
                ]},
                {"id": 4, "acronym": "B", "name": "Area B"},
            ],
        }
        out = []
        _flatten_allen_tree(tree, out)
        assert len(out) == 4
        assert [s["id"] for s in out] == [1, 2, 3, 4]

    def test_extra_fields_stripped(self):
        node = {"id": 1, "acronym": "X", "name": "X area", "color": "red", "depth": 0}
        out = []
        _flatten_allen_tree(node, out)
        assert set(out[0].keys()) == {"id", "acronym", "name"}


# ---------------------------------------------------------------------------
# build_lookup_dicts
# ---------------------------------------------------------------------------


class TestBuildLookupDicts:
    def test_all_dicts_populated(self, lookups):
        by_acronym, by_name, by_acronym_lower, by_name_lower = lookups
        assert "VISp" in by_acronym
        assert "Primary visual area" in by_name
        assert "visp" in by_acronym_lower
        assert "primary visual area" in by_name_lower

    def test_sizes(self, lookups):
        by_acronym, by_name, by_acronym_lower, by_name_lower = lookups
        assert len(by_acronym) == len(SAMPLE_STRUCTURES)
        assert len(by_name) == len(SAMPLE_STRUCTURES)


# ---------------------------------------------------------------------------
# _extract_area
# ---------------------------------------------------------------------------


class TestExtractArea:
    def test_dict_repr(self):
        assert _extract_area("{'area': 'VISp', 'depth': '20'}") == "VISp"

    def test_dict_repr_other_area(self):
        assert _extract_area("{'area': 'RSP', 'depth': '350'}") == "RSP"

    def test_key_value_pairs(self):
        assert _extract_area("area: VISp,depth: 175") == "VISp"

    def test_key_value_with_spaces(self):
        assert _extract_area("area:  VISam ,depth: 100") == "VISam"

    def test_no_area_key(self):
        assert _extract_area("{'depth': '20'}") is None

    def test_plain_value(self):
        assert _extract_area("VISp") is None

    def test_invalid_dict(self):
        assert _extract_area("{not valid python}") is None

    def test_empty(self):
        assert _extract_area("") is None


# ---------------------------------------------------------------------------
# _match_single
# ---------------------------------------------------------------------------


class TestMatchSingle:
    def test_exact_acronym(self, lookup_args):
        s = _match_single("CA1", *lookup_args)
        assert s is not None
        assert s["id"] == 382

    def test_exact_name(self, lookup_args):
        s = _match_single("Thalamus", *lookup_args)
        assert s is not None
        assert s["id"] == 549

    def test_case_insensitive_acronym(self, lookup_args):
        s = _match_single("visp", *lookup_args)
        assert s is not None
        assert s["id"] == 385

    def test_case_insensitive_name(self, lookup_args):
        s = _match_single("thalamus", *lookup_args)
        assert s is not None
        assert s["id"] == 549

    def test_no_match(self, lookup_args):
        assert _match_single("V1", *lookup_args) is None

    def test_no_match_trivial(self, lookup_args):
        # _match_single does NOT filter trivials — that's match_location's job
        # "unknown" isn't an acronym or name, so it returns None anyway
        assert _match_single("unknown", *lookup_args) is None


# ---------------------------------------------------------------------------
# match_location
# ---------------------------------------------------------------------------


class TestMatchLocation:
    def test_exact_acronym(self, lookup_args):
        results = match_location("CA1", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 382

    def test_exact_name(self, lookup_args):
        results = match_location("Primary visual area", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 385

    def test_case_insensitive_acronym(self, lookup_args):
        results = match_location("visp", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 385

    def test_case_insensitive_name(self, lookup_args):
        results = match_location("thalamus", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 549

    def test_whitespace_stripped(self, lookup_args):
        results = match_location("  CA1  ", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 382

    # --- Trivial / no-match ---

    def test_unknown(self, lookup_args):
        assert match_location("unknown", *lookup_args) == []

    def test_none_string(self, lookup_args):
        assert match_location("none", *lookup_args) == []

    def test_empty(self, lookup_args):
        assert match_location("", *lookup_args) == []

    def test_whitespace_only(self, lookup_args):
        assert match_location(" ", *lookup_args) == []

    def test_na(self, lookup_args):
        assert match_location("N/A", *lookup_args) == []

    def test_void(self, lookup_args):
        assert match_location("void", *lookup_args) == []

    def test_unspecific(self, lookup_args):
        assert match_location("unspecific", *lookup_args) == []

    def test_v1_no_match(self, lookup_args):
        assert match_location("V1", *lookup_args) == []

    # --- Structured strings ---

    def test_dict_repr(self, lookup_args):
        results = match_location("{'area': 'VISp', 'depth': '20'}", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 385

    def test_key_value_pairs(self, lookup_args):
        results = match_location("area: VISp,depth: 175", *lookup_args)
        assert len(results) == 1
        assert results[0]["id"] == 385

    def test_dict_repr_rsp(self, lookup_args):
        results = match_location("{'area': 'RSP', 'depth': '350'}", *lookup_args)
        assert len(results) == 1
        assert results[0]["acronym"] == "RSP"

    # --- Comma-separated ---

    def test_comma_separated_all_match(self, lookup_args):
        results = match_location("VISp,VISrl,VISal", *lookup_args)
        ids = {s["id"] for s in results}
        assert ids == {385, 409, 402}

    def test_comma_separated_partial_match(self, lookup_args):
        results = match_location("VISp,VISrl,VISlm,VISal", *lookup_args)
        ids = {s["id"] for s in results}
        # VISlm is in our sample structures (id 44)
        assert 385 in ids
        assert 409 in ids
        assert 402 in ids

    def test_comma_separated_with_spaces(self, lookup_args):
        results = match_location("CA1, CA3, DG", *lookup_args)
        ids = {s["id"] for s in results}
        assert ids == {382, 463, 726}

    def test_comma_separated_none_match(self, lookup_args):
        assert match_location("FOO,BAR,BAZ", *lookup_args) == []

    def test_comma_separated_not_confused_with_key_value(self, lookup_args):
        # "area: VISp,depth: 175" should match via _extract_area, not comma split
        results = match_location("area: VISp,depth: 175", *lookup_args)
        assert len(results) == 1
        assert results[0]["acronym"] == "VISp"


# ---------------------------------------------------------------------------
# structure_to_anatomy
# ---------------------------------------------------------------------------


class TestStructureToAnatomy:
    def test_format(self):
        s = {"id": 549, "acronym": "TH", "name": "Thalamus"}
        entry = structure_to_anatomy(s)
        assert entry == {
            "schemaKey": "Anatomy",
            "name": "Thalamus",
            "identifier": "https://purl.brain-bican.org/ontology/mbao/MBA_549",
        }

    def test_large_id(self):
        s = {"id": 312782546, "acronym": "VISa", "name": "Anterior area"}
        entry = structure_to_anatomy(s)
        assert entry["identifier"] == "https://purl.brain-bican.org/ontology/mbao/MBA_312782546"


# ---------------------------------------------------------------------------
# get_candidate_dandisets
# ---------------------------------------------------------------------------


class TestGetCandidateDandisets:
    def test_matches(self, lookups):
        scan_cache = {
            "000001": {
                "dandiset_id": "000001",
                "imaging_locations": {"VISp": 5},
                "electrode_locations": {},
                "icephys_locations": {},
            },
            "000002": {
                "dandiset_id": "000002",
                "imaging_locations": {},
                "electrode_locations": {"unknown": 3},
                "icephys_locations": {},
            },
            "000003": {
                "dandiset_id": "000003",
                "imaging_locations": {},
                "electrode_locations": {"CA1": 2, "V1": 1},
                "icephys_locations": {},
            },
        }
        candidates = get_candidate_dandisets(scan_cache, lookups)
        assert "000001" in candidates
        assert "000002" not in candidates  # only "unknown"
        assert "000003" in candidates

    def test_empty_cache(self, lookups):
        assert get_candidate_dandisets({}, lookups) == {}

    def test_no_locations(self, lookups):
        scan_cache = {
            "000001": {
                "dandiset_id": "000001",
                "imaging_locations": {},
                "electrode_locations": {},
                "icephys_locations": {},
            },
        }
        assert get_candidate_dandisets(scan_cache, lookups) == {}

    def test_structured_location_matches(self, lookups):
        scan_cache = {
            "000001": {
                "dandiset_id": "000001",
                "imaging_locations": {"{'area': 'VISp', 'depth': '20'}": 3},
                "electrode_locations": {},
                "icephys_locations": {},
            },
        }
        candidates = get_candidate_dandisets(scan_cache, lookups)
        assert "000001" in candidates

    def test_comma_separated_matches(self, lookups):
        scan_cache = {
            "000001": {
                "dandiset_id": "000001",
                "imaging_locations": {"VISp,VISrl,VISal": 2},
                "electrode_locations": {},
                "icephys_locations": {},
            },
        }
        candidates = get_candidate_dandisets(scan_cache, lookups)
        assert "000001" in candidates


# ---------------------------------------------------------------------------
# Label cache I/O
# ---------------------------------------------------------------------------


class TestLabelCache:
    def test_load_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr("label_anatomy.LABEL_CACHE_FILE", str(tmp_path / "nonexistent.jsonl"))
        cache = load_label_cache()
        assert cache == {}

    def test_round_trip(self, tmp_path, monkeypatch):
        cache_file = str(tmp_path / "label_cache.jsonl")
        monkeypatch.setattr("label_anatomy.LABEL_CACHE_FILE", cache_file)

        entry = {
            "dandiset_id": "000001",
            "asset_id": "abc123",
            "status": "would_update",
        }
        append_label_cache(entry)
        append_label_cache({
            "dandiset_id": "000002",
            "asset_id": "def456",
            "status": "skipped_no_match",
        })

        cache = load_label_cache()
        assert len(cache) == 2
        assert ("000001", "abc123") in cache
        assert cache[("000001", "abc123")]["status"] == "would_update"
        assert ("000002", "def456") in cache
