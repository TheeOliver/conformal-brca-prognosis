"""Bayesian uncertainty illustrations for explicitly hypothetical training profiles."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def hypothetical_profiles(
    training: pd.DataFrame, config: Any
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Combine training medians/modes and age quartiles; never select patient rows.

    Subtype profiles alter only the subtype of the median/mode reference. They
    illustrate the fitted prognostic predictions and are not causal contrasts.
    """
    if len(training) < config.eda.minimum_group_size:
        raise ValueError("Too few training observations for aggregate hypothetical profiles")
    features = [*config.feature_sets.clinical, *config.feature_sets.clinical_molecular]
    numeric = set(config.preprocessing.numeric_columns)
    reference = {}
    for column in features:
        if column in numeric:
            value = training[column].median()
            if not np.isfinite(value):
                raise ValueError("A numeric profile predictor has no finite training median")
            reference[column] = float(value)
        else:
            modes = training[column].mode(dropna=True)
            value = modes.iloc[0] if len(modes) else config.preprocessing.unknown_category
            reference[column] = value.item() if isinstance(value, np.generic) else value
    quartiles = list(config.eda.quantiles[1:-1])
    if len(quartiles) != 3 or not (0 < quartiles[0] < quartiles[1] < quartiles[2] < 1):
        raise ValueError("EDA quantiles must contain ordered lower/median/upper quartiles")
    profile_rows = [reference.copy()]
    definitions = [
        {
            "profile_id": "reference",
            "label": "Median age and modal predictors",
            "construction": "Training numeric medians and categorical modes",
            "hypothetical": True,
        }
    ]
    for name, quantile in zip(
        ("age_lower_quartile", "age_upper_quartile"), quartiles[::2], strict=True
    ):
        age = float(training["age_at_diagnosis"].quantile(quantile))
        profile_rows.append(reference | {"age_at_diagnosis": age})
        definitions.append(
            {
                "profile_id": name,
                "label": f"Age P{100 * quantile:g}: {age:.1f} years; other predictors at reference",
                "construction": "Training age quantile with numeric medians and categorical modes",
                "age_quantile": quantile,
                "hypothetical": True,
            }
        )
    counts = training["pam50_subtype"].dropna().value_counts()
    eligible = sorted(
        str(level) for level, n in counts.items() if n >= config.eda.minimum_group_size
    )
    for j, level in enumerate(eligible):
        if level == str(reference["pam50_subtype"]):
            definitions[0]["subtype_training_count"] = int(counts[level])
            continue
        profile_rows.append(reference | {"pam50_subtype": level})
        definitions.append(
            {
                "profile_id": f"subtype_{j}",
                "label": f"{level} subtype; age and other predictors at reference",
                "construction": (
                    "Training-observed subtype with numeric medians and categorical modes"
                ),
                "subtype_training_count": int(counts[level]),
                "hypothetical": True,
            }
        )
    frame = pd.DataFrame(profile_rows, columns=features)
    frame.index = pd.Index([row["profile_id"] for row in definitions], name="hypothetical_profile")
    return frame, definitions


def hypothetical_uncertainty_report(
    profiles: pd.DataFrame,
    definitions: list[dict[str, Any]],
    models: dict[str, Any],
    preprocessors: dict[str, Any],
    config: Any,
) -> dict[str, Any]:
    """Compute curve credible bands and individual predictive time intervals once."""
    if len(definitions) != len(profiles) or any(not row.get("hypothetical") for row in definitions):
        raise ValueError("Every illustration must carry an aligned hypothetical-profile definition")
    times = np.asarray(config.bayesian.ppc_time_grid_months, dtype=float)
    if (
        times.ndim != 1
        or not np.isfinite(times).all()
        or (times <= 0).any()
        or (np.diff(times) <= 0).any()
    ):
        raise ValueError(
            "The configured posterior prediction grid must be increasing positive months"
        )
    alphas = np.asarray(config.conformal.alphas, dtype=float)
    if not np.isfinite(alphas).all() or ((alphas <= 0) | (alphas >= 1)).any():
        raise ValueError("Predictive interval alpha values must be within (0, 1)")
    probabilities = np.concatenate(([0.5], np.ravel(np.column_stack((alphas / 2, 1 - alphas / 2)))))
    result = {
        "schema_version": 1,
        "endpoint": config.outcome.primary_endpoint,
        "model": "bayesian_weibull_aft",
        "profile_source": "hypothetical_combinations_of_training_aggregates",
        "patient_records_included": False,
        "calibration_accessed": False,
        "test_accessed": False,
        "time_grid_months": times.tolist(),
        "credible_probability": float(config.bayesian.predictive_interval_probability),
        "credible_band_type": "pointwise_parameter_uncertainty_in_survival_function",
        "predictive_interval_type": (
            "posterior_mixture_event_time_including_future_outcome_variation"
        ),
        "predictive_time_horizon": "unrestricted_model_extrapolation",
        "caption_template": (
            "Bayesian Weibull predictions for the explicitly hypothetical profile {label}, "
            "assembled from training aggregates. Left: posterior mean survival and pointwise "
            f"{100 * config.bayesian.predictive_interval_probability:g}% credible bands. "
            "Right: unrestricted posterior predictive event-time intervals at the labelled "
            "probabilities, with medians marked. Blue solid circles denote clinical predictors; "
            "orange dashed squares denote clinical plus PAM50. Curve uncertainty and individual "
            "event-time variation are different quantities; neither panel measures coverage."
        ),
        "interpretation": (
            "Hypothetical profiles are assembled from training aggregates and do not represent "
            "identified patients. Curve bands describe uncertainty in predicted survival "
            "probabilities; time intervals include future outcome variation. These are model-based "
            "uncertainty illustrations, not measured coverage or causal comparisons. Net DSS "
            "interpretation retains the competing-death censoring assumptions."
        ),
        "profiles": [],
    }
    for j, definition in enumerate(definitions):
        result["profiles"].append(
            {
                **definition,
                "features": profiles.iloc[j].to_dict(),
                "models": {},
            }
        )
    for feature in ("clinical", "clinical_molecular"):
        X = preprocessors[feature].transform(profiles)
        model = models[feature]
        bands = model.predict_survival_bands(X, times)
        quantiles = model.predict_quantiles(X, probabilities)
        if (
            quantiles.shape != (len(profiles), len(probabilities))
            or not np.isfinite(quantiles).all()
        ):
            raise ValueError("Hypothetical-profile time quantiles must be finite and aligned")
        if (quantiles <= 0).any():
            raise ValueError("Hypothetical-profile predictive survival times must be positive")
        for name in ("mean", "lower", "upper"):
            values = np.asarray(bands[name])
            if values.shape != (len(profiles), len(times)) or not np.isfinite(values).all():
                raise ValueError("Hypothetical-profile credible bands must be finite and aligned")
            if ((values < 0) | (values > 1)).any():
                raise ValueError("Hypothetical-profile survival probabilities must be in [0, 1]")
        if (bands["lower"] > bands["upper"]).any():
            raise ValueError("Hypothetical-profile credible bounds are reversed")
        for j in range(len(profiles)):
            predictive = []
            for a, alpha in enumerate(alphas):
                lower, upper = quantiles[j, 1 + 2 * a : 3 + 2 * a]
                if lower > upper:
                    raise ValueError("Hypothetical-profile predictive time bounds are reversed")
                predictive.append(
                    {
                        "alpha": float(alpha),
                        "nominal_probability": float(1 - alpha),
                        "lower_months": float(lower),
                        "upper_months": float(upper),
                        "width_months": float(upper - lower),
                    }
                )
            result["profiles"][j]["models"][feature] = {
                "survival_mean": bands["mean"][j].tolist(),
                "survival_credible_lower": bands["lower"][j].tolist(),
                "survival_credible_upper": bands["upper"][j].tolist(),
                "predictive_median_months": float(quantiles[j, 0]),
                "predictive_intervals": predictive,
            }
    return result
