---
description: The survival outcome contract - dtypes, units and censoring - shared by every model and metric.
paths:
  - "src/brca/**"
  - "scripts/**"
  - "tests/**"
---

# Survival labels

One representation, everywhere:

```python
from sksurv.util import Surv

y = Surv.from_arrays(
    event=df["event"].to_numpy(bool),
    time=df["survival_months"].to_numpy(float),
)
```

- `event` is **bool**, not 0/1 int: scikit-survival rejects int, and PyMC silently builds a
  nonsense likelihood from a float.
- `time` is **months**, float, strictly positive. Years or days corrupt every downstream metric.
- **Never drop censored rows.** Right censoring is information — the patient was event-free up
  to that time. Dropping them, or treating censoring times as events, biases survival downward.
- Never impute a survival time; never censor administratively without recording the horizon.
- `time <= 0` is a data error: fail loudly, do not clip. Keep `y` and the feature frame
  index-aligned — reorder them together or not at all.
