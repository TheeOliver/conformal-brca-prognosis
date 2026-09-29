"""Access to the single source of truth for every analysis constant.

No seed, fraction, alpha or threshold may be written as a literal in ``src/`` or
``scripts/``. See ``.claude/rules/reproducibility.md``.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "default.yaml"


def _namespace(obj: Any) -> Any:
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: _namespace(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [_namespace(v) for v in obj]
    return obj


def load_config(path: Path | str | None = None) -> SimpleNamespace:
    """Load the project configuration as a dot-accessible namespace."""
    path = Path(path) if path is not None else DEFAULT_CONFIG_PATH
    with path.open() as fh:
        return _namespace(yaml.safe_load(fh))
