"""Offline tests for the cBioPortal fetch module -- no network, no real data."""

from __future__ import annotations

import gzip
import json

import pytest

from brca.data.cbioportal import RAW_FILES, save_raw, to_wide, update_checksum_file

RECORDS = [
    {"patientId": "P2", "clinicalAttributeId": "OS_MONTHS", "value": "12.5"},
    {"patientId": "P1", "clinicalAttributeId": "OS_MONTHS", "value": "40.0"},
    {"patientId": "P1", "clinicalAttributeId": "OS_STATUS", "value": "0:LIVING"},
]
PAYLOADS = {
    "study": {"studyId": "toy", "name": "Toy", "publicStudy": True},
    "attributes": [{"clinicalAttributeId": "OS_MONTHS"}],
    "patient": RECORDS,
    "sample": [],
}


def test_to_wide_one_row_per_id_and_missing_stays_missing():
    wide = to_wide(RECORDS, "patientId")
    assert list(wide.index) == ["P1", "P2"]
    assert wide.loc["P1", "OS_STATUS"] == "0:LIVING"
    assert wide["OS_STATUS"].isna().sum() == 1  # P2 never reported it


def test_to_wide_rejects_duplicate_attributes():
    with pytest.raises(ValueError, match="duplicate"):
        to_wide([*RECORDS, RECORDS[0]], "patientId")


def test_save_raw_roundtrips_and_refuses_to_overwrite(tmp_path):
    checksums = save_raw(PAYLOADS, tmp_path)
    assert set(checksums) == set(RAW_FILES.values())
    with gzip.open(tmp_path / RAW_FILES["patient"]) as fh:
        assert json.load(fh) == RECORDS
    assert json.loads((tmp_path / "PROVENANCE.json").read_text())["public_study"] is True
    with pytest.raises(FileExistsError, match="immutable"):
        save_raw(PAYLOADS, tmp_path)


def test_save_raw_is_byte_deterministic(tmp_path):
    a = save_raw(PAYLOADS, tmp_path / "a")
    b = save_raw(PAYLOADS, tmp_path / "b")
    assert a == b


def test_checksum_file_merges_in_sha256sum_format(tmp_path):
    f = tmp_path / "CHECKSUMS.sha256"
    f.write_text("aaa  other/file.gz\n")
    update_checksum_file(f, {"x.gz": "bbb"}, "snap")
    assert f.read_text() == "aaa  other/file.gz\nbbb  snap/x.gz\n"
