---
description: Seeds, configuration and run manifests - how any reported number stays traceable to the code that made it.
paths:
  - "src/brca/**"
  - "scripts/**"
---

# Reproducibility

- **One seed, one config.** `config/default.yaml` holds every seed, fraction, alpha and
  threshold. No analysis constant is written as a literal in `src/` or `scripts/` — if you
  find yourself typing `0.05` or `42`, it belongs in the config.
- **No global RNG state.** Never `np.random.seed(...)`, `np.random.rand(...)` or
  `random.random()`. Always `rng = np.random.default_rng(config.seed)` and pass `rng` down.
  Ruff's `NPY002` enforces this; do not silence it.
- Derive child seeds deterministically (`rng.spawn(...)` or an explicit offset per model), so
  adding a fourth model does not shift the draws of the first three.
- **Every stage writes a manifest** to `outputs/manifests/<stage>_<timestamp>.json`: git SHA,
  whether the tree was dirty, the config hash, resolved package versions, the SLURM job ID if
  any, and wall-clock. A result without a manifest cannot be defended.
- Refuse to write into `outputs/` from a dirty working tree unless explicitly overridden —
  a figure produced from uncommitted code is untraceable.
- Pin dependencies through `uv.lock`, and commit it. `uv sync` must rebuild the exact
  environment; never `pip install` into the venv by hand.
- Stages are idempotent: rerunning one with the same config and seed overwrites its outputs
  with byte-identical content. If it does not, something is reading unseeded randomness.
