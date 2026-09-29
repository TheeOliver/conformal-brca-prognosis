"""Fetch the public METABRIC clinical tables from the cBioPortal REST API.

The study ``brca_metabric`` (Pereira 2016 / Curtis 2012 / Rueda 2019) is flagged
``publicStudy: true`` on cBioPortal. Only clinical tables are fetched -- the thesis
needs clinical predictors, PAM50 and survival, not expression or raw genomics
(those remain controlled-access at EGA ``EGAS00000000083``).

Raw responses are stored exactly as returned (gzipped JSON) and never edited;
reshaping happens downstream in stage 01. ``data/`` is gitignored and must stay
so: public availability does not make this repo a redistribution point.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import urllib.parse
import urllib.request
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

API_BASE = "https://www.cbioportal.org/api"
USER_AGENT = "conformal-brca-prognosis/0.1 (academic thesis; clinical tables only)"
PAGE_SIZE = 10_000_000  # the API returns everything in one page at this size

RAW_FILES = {
    "study": "study.json.gz",
    "attributes": "clinical_attributes.json.gz",
    "patient": "clinical_patient.json.gz",
    "sample": "clinical_sample.json.gz",
}


def _get(path: str, params: dict[str, str | int] | None = None, timeout: int = 120) -> object:
    url = f"{API_BASE}/{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def fetch_study(study_id: str) -> dict[str, object]:
    """Return the raw API payloads for one study's metadata and clinical tables."""
    clinical = f"studies/{study_id}/clinical-data"
    common = {"projection": "SUMMARY", "pageSize": PAGE_SIZE}
    return {
        "study": _get(f"studies/{study_id}"),
        "attributes": _get(f"studies/{study_id}/clinical-attributes"),
        "patient": _get(clinical, {"clinicalDataType": "PATIENT", **common}),
        "sample": _get(clinical, {"clinicalDataType": "SAMPLE", **common}),
    }


def to_wide(records: Iterable[dict], id_field: str) -> pd.DataFrame:
    """Pivot cBioPortal long-format clinical records to one row per ``id_field``.

    Values stay as the strings the API returned; typing is stage 01's job.
    """
    columns = [id_field, "clinicalAttributeId", "value"]
    df = pd.DataFrame.from_records(list(records), columns=columns)
    if df.duplicated([id_field, "clinicalAttributeId"]).any():
        raise ValueError(f"duplicate ({id_field}, attribute) pairs in API response")
    wide = df.pivot(index=id_field, columns="clinicalAttributeId", values="value")
    wide.columns.name = None
    return wide.sort_index()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_raw(payloads: dict[str, object], out_dir: Path, *, force: bool = False) -> dict[str, str]:
    """Write payloads as gzipped JSON plus ``PROVENANCE.json``; return sha256 per file.

    ``data/raw`` is immutable: refuses to overwrite an existing snapshot unless
    ``force`` -- a silently replaced snapshot would make earlier results
    untraceable.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    existing = [out_dir / name for name in RAW_FILES.values() if (out_dir / name).exists()]
    if existing and not force:
        raise FileExistsError(
            f"raw snapshot already exists in {out_dir}; data/raw is immutable. "
            "Pass --force only if you intend to replace it, and record why in STATUS.md."
        )
    checksums = {}
    for key, name in RAW_FILES.items():
        path = out_dir / name
        # mtime=0 keeps the gzip bytes deterministic for identical content
        with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            gz.write(json.dumps(payloads[key], sort_keys=True, indent=1).encode())
        checksums[name] = _sha256(path)

    study = payloads["study"]
    provenance = {
        "source": "cBioPortal REST API",
        "api_base": API_BASE,
        "study_id": study.get("studyId"),
        "study_name": study.get("name"),
        "citation": study.get("citation"),
        "pmid": study.get("pmid"),
        "public_study": study.get("publicStudy"),
        "cbioportal_import_date": study.get("importDate"),
        "downloaded_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "record_counts": {k: len(payloads[k]) for k in ("attributes", "patient", "sample")},
        "files_sha256": checksums,
    }
    (out_dir / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return checksums


def update_checksum_file(checksum_file: Path, checksums: dict[str, str], rel_dir: str) -> None:
    """Merge entries into ``data/raw/CHECKSUMS.sha256`` (``sha256sum -c`` format)."""
    lines = {}
    if checksum_file.exists():
        for line in checksum_file.read_text().splitlines():
            if line.strip():
                digest, name = line.split(maxsplit=1)
                lines[name.strip()] = digest
    for name, digest in checksums.items():
        lines[f"{rel_dir}/{name}"] = digest
    checksum_file.write_text("".join(f"{d}  {n}\n" for n, d in sorted(lines.items())))
