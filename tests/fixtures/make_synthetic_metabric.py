"""Generate a synthetic METABRIC-shaped cohort.

Contains no real patient. It exists so the pipeline and the test suite run before
the access-controlled dataset is available, and so no test ever reads patient data.
See ``.claude/rules/data-and-privacy.md``.

The generating process is a crude Weibull AFT with administrative and random
censoring. It reproduces the *shape* of the data, not its biology -- numbers from
it are never presented as results.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT = Path(__file__).with_name("synthetic_metabric.csv")
N_PATIENTS = 400
SEED = 8927
ADMIN_CENSOR_MONTHS = 300.0

_SUBTYPES = ("LumA", "LumB", "Her2", "Basal", "Normal")
_SUBTYPE_PREVALENCE = (0.35, 0.25, 0.12, 0.18, 0.10)
# Time-scale offsets: positive = longer survival. Ordering mirrors the published
# prognostic ordering so a Kaplan-Meier sanity check is meaningful.
_SUBTYPE_OFFSET = {"LumA": 0.35, "LumB": 0.05, "Her2": -0.25, "Basal": -0.40, "Normal": 0.20}


def make_cohort(n: int = N_PATIENTS, seed: int = SEED) -> pd.DataFrame:
    """Build a synthetic cohort. Deterministic given ``seed``."""
    rng = np.random.default_rng(seed)

    pam50 = rng.choice(_SUBTYPES, size=n, p=_SUBTYPE_PREVALENCE)
    age = rng.normal(61.0, 13.0, n).clip(25.0, 95.0)
    size = rng.gamma(4.0, 6.0, n).clip(1.0, 200.0)
    grade = rng.choice([1.0, 2.0, 3.0], size=n, p=[0.15, 0.45, 0.40])
    nodes = rng.poisson(1.8, n).astype(float)

    er = rng.choice(["positive", "negative"], size=n, p=[0.75, 0.25])
    pr = np.where(rng.random(n) < 0.85, er, np.where(er == "positive", "negative", "positive"))
    her2 = rng.choice(["positive", "negative"], size=n, p=[0.15, 0.85])

    log_scale = (
        4.6
        - 0.012 * (age - 61.0)
        - 0.006 * (size - 24.0)
        - 0.25 * (grade - 2.0)
        - 0.05 * nodes
        + np.array([_SUBTYPE_OFFSET[s] for s in pam50])
    )
    true_time = np.exp(log_scale) * rng.weibull(1.3, n)
    random_censor = rng.exponential(200.0, n)

    observed = np.minimum(np.minimum(true_time, random_censor), ADMIN_CENSOR_MONTHS)
    event = true_time <= np.minimum(random_censor, ADMIN_CENSOR_MONTHS)

    return pd.DataFrame(
        {
            "patient_id": [f"SYN-{i:04d}" for i in range(n)],
            "survival_months": np.round(np.maximum(observed, 0.1), 2),
            "event": event,
            "age_at_diagnosis": np.round(age, 1),
            "tumor_size": np.round(size, 1),
            "tumor_grade": grade,
            "lymph_nodes_positive": nodes,
            "er_status": er,
            "pr_status": pr,
            "her2_status": her2,
            "pam50_subtype": pam50,
        }
    )


def main() -> None:
    cohort = make_cohort()
    cohort.to_csv(OUTPUT, index=False)
    print(f"wrote {OUTPUT} ({len(cohort)} rows, {cohort['event'].mean():.1%} events)")


if __name__ == "__main__":
    main()
