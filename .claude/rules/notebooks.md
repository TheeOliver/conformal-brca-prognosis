---
description: What notebooks are for here, and why the reference notebooks must never be copied from.
paths:
  - "notebooks/**"
---

# Notebooks

Notebooks are for **exploration and nothing else**. Nothing in the thesis depends on one
having been run.

- No analysis logic is defined in a notebook. Import it: `from brca.models import ...`. If a
  function is worth keeping, move it into `src/brca/` and add a test — then import it back.
- A notebook never writes to `outputs/` and never appears in a `make` target. Results come
  from `scripts/01…05`, which are reproducible and manifested.
- Strip outputs before committing. `nbstripout` runs automatically as a PostToolUse hook, but
  check anyway — cell outputs can carry patient-level rows straight into git history.
- Notebooks may read `data/processed/`, but they respect the split like everything else:
  **never explore the test set.** Looking at test data informally still burns it.

## The reference notebooks are lossy

`references/notebooks_reference/*.md` are notebook JSON rendered to markdown and split on
commas. They are **damaged**: dicts are cut mid-literal (`SAMPLE_KWARGS = {"chains": 4`),
function arguments are dropped (`np.percentile(y["futime"]`), and base64 image blobs are
inlined. Read them for the approach and the API surface; **never copy code out of them** — it
will look correct and fail, or worse, run and be wrong.

Their filenames contain spaces and em-dashes — **always quote the path** in shell commands.
