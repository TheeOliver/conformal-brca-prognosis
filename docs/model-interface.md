# Model interface

> Pulled in by: `add-survival-model` skill.

Every model implements this contract, so the pipeline, the conformal layer and the evaluation
table can treat them interchangeably.

```python
class SurvivalModel(Protocol):
    def fit(self, X: pd.DataFrame, y: np.ndarray) -> "SurvivalModel": ...

    def predict_risk(self, X: pd.DataFrame) -> np.ndarray:
        """Relative risk, higher = worse. Shape (n,). Used for the C-index."""

    def predict_survival_function(self, X: pd.DataFrame, times: np.ndarray) -> np.ndarray:
        """S(t | x) on the SHARED grid. Shape (n, len(times)), values in [0, 1],
        monotone non-increasing along axis 1."""

    def save(self, path: Path) -> None: ...

    @classmethod
    def load(cls, path: Path) -> "SurvivalModel": ...
```

## Why `predict_survival_function` carries the weight

Both the Brier score and the conformal layer consume it. It must be correct **on the shared
time grid**, not on whatever internal grid the estimator happens to use — a `StepFunction` from
scikit-survival must be evaluated at `times`, not sampled at its own knots.

Two invariants deserve a test, because both fail silently:

- monotone non-increasing in `t`;
- values within `[0, 1]`.

## Per-model notes

| Model | Implementation | Watch for |
| --- | --- | --- |
| Cox PH | `sksurv.linear_model.CoxPHSurvivalAnalysis` | check proportional hazards; report ties handling |
| Bayesian Weibull AFT | PyMC, `pm.Censored` | must also expose diagnostics and persist InferenceData |
| Random Survival Forest | `sksurv.ensemble` | untuned by design; state that in every comparison |

A Bayesian model must expose its posterior, not only a point prediction — collapsing it to a
mean discards exactly what RQ4 asks about.

## Registration

1. Module at `src/brca/models/<name>.py`, one model per file.
2. Register in `src/brca/models/__init__.py`.
3. Config block named after the model in `config/default.yaml`.
4. It must run for **both** feature sets, or the RQ2 comparison has a hole.
