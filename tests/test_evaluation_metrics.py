"""IPCW evaluation support, API equivalence and genuinely paired bootstrap checks."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from sksurv.metrics import (
    brier_score,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    integrated_brier_score,
)
from sksurv.util import Surv

from brca.evaluation.metrics import (
    EvaluationSupportError,
    _capped_test_labels,
    evaluate_models,
    training_time_grid,
)


def evaluation_config(**overrides):
    values = dict(
        time_grid_percentiles=[20, 40, 60, 80],
        maximum_horizon_months=90.0,
        confidence_level=0.95,
        bootstrap_replicates=20,
        minimum_bootstrap_success_fraction=0.9,
        primary_metric="uno_c",
    )
    return SimpleNamespace(**(values | overrides))


@pytest.fixture
def evaluation_data():
    y_train = Surv.from_arrays(np.arange(301) % 4 != 0, np.linspace(0.5, 150.5, 301))
    rng = np.random.default_rng(7557)
    y_test = Surv.from_arrays(rng.random(150) > 0.3, np.linspace(0.25, 200, 150))
    settings = evaluation_config()
    times = training_time_grid(y_train, settings)
    risk = -y_test["time"] / 80 + rng.normal(scale=0.3, size=len(y_test))
    survival = np.exp(-np.exp(risk[:, None]) * times[None, :] / 80)
    predictions = {
        "cox__clinical": {"risk": risk, "survival": survival},
        "cox__clinical_molecular": {"risk": risk, "survival": survival},
    }
    return y_train, y_test, predictions, times, settings


def test_training_grid_depends_only_on_training_and_config():
    y_train = Surv.from_arrays(
        [True, False, True, False, True, False], [1.0, 2.0, 3.0, 4.0, 5.0, 8.0]
    )
    config = evaluation_config(time_grid_percentiles=[20, 50, 80], maximum_horizon_months=4.0)
    assert_array_equal(training_time_grid(y_train, config), [2.0, 3.5, 4.0])
    assert_array_equal(training_time_grid(y_train, {"evaluation": vars(config)}), [2.0, 3.5, 4.0])
    assert_array_equal(
        training_time_grid(y_train, SimpleNamespace(evaluation=config)), [2.0, 3.5, 4.0]
    )
    config.maximum_horizon_months = 100
    assert training_time_grid(y_train, config)[-1] < 5.0


def test_full_sample_metrics_equal_public_sksurv_and_all_intervals_present(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    original_y = y_test.copy()
    result = evaluate_models(
        y_train, y_test, predictions, times, config, np.random.default_rng(332)
    )
    assert_array_equal(y_test, original_y)
    capped, count = _capped_test_labels(y_test, times[-1])
    clinical = predictions["cox__clinical"]
    expected_auc, expected_mean_auc = cumulative_dynamic_auc(
        y_train, capped, 1 - clinical["survival"], times
    )
    actual = result["models"]["cox__clinical"]
    assert_allclose(
        actual["uno_c"]["estimate"],
        concordance_index_ipcw(
            y_train, capped, clinical["risk"], tau=np.nextafter(times[-1], np.inf)
        )[0],
    )
    assert_allclose(
        actual["ibs"]["estimate"],
        integrated_brier_score(y_train, capped, clinical["survival"], times),
    )
    assert_allclose(actual["mean_auc"]["estimate"], expected_mean_auc)
    assert_allclose(
        [row["estimate"] for row in actual["brier"]],
        brier_score(y_train, capped, clinical["survival"], times)[1],
    )
    assert_allclose([row["estimate"] for row in actual["auc"]], expected_auc)
    for row in [
        actual["uno_c"],
        actual["ibs"],
        actual["mean_auc"],
        *actual["brier"],
        *actual["auc"],
    ]:
        assert row["status"] == "ok"
        assert 0 <= row["ci_lower"] <= row["ci_upper"] <= 1
        assert row["n_bootstrap_success"] + row["n_bootstrap_failed"] == config.bootstrap_replicates
    assert result["reportable"]
    assert result["n_test_followup_capped"] == count
    json.dumps(result, allow_nan=False)


def test_identical_models_have_exactly_zero_paired_intervals(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    result = evaluate_models(y_train, y_test, predictions, times, config, np.random.default_rng(97))
    comparison = result["paired_feature_comparisons"]["cox"]
    for row in [
        comparison["uno_c"],
        comparison["ibs"],
        comparison["mean_auc"],
        *comparison["brier"],
        *comparison["auc"],
    ]:
        assert row["estimate"] == row["ci_lower"] == row["ci_upper"] == 0.0


@pytest.mark.parametrize("stratified", [False, True])
def test_paired_brier_ci_matches_manual_same_patient_resamples(evaluation_data, stratified):
    y_train, y_test, predictions, times, config = evaluation_data
    predictions["cox__clinical_molecular"] = {
        "risk": predictions["cox__clinical"]["risk"] + 0.1,
        "survival": predictions["cox__clinical"]["survival"] ** 1.5,
    }
    strata = np.arange(len(y_test)) % 3 == 0 if stratified else None
    result = evaluate_models(
        y_train,
        y_test,
        predictions,
        times,
        config,
        np.random.default_rng(64),
        bootstrap_strata=strata,
    )
    rng = np.random.default_rng(64)
    capped, _ = _capped_test_labels(y_test, times[-1])
    differences = []
    for _ in range(config.bootstrap_replicates):
        if stratified:
            groups = [np.flatnonzero(strata == level) for level in [False, True]]
            indices = np.concatenate([rng.choice(group, size=len(group)) for group in groups])
        else:
            indices = rng.integers(0, len(y_test), size=len(y_test))
        base = brier_score(
            y_train, capped[indices], predictions["cox__clinical"]["survival"][indices], times
        )[1]
        molecular = brier_score(
            y_train,
            capped[indices],
            predictions["cox__clinical_molecular"]["survival"][indices],
            times,
        )[1]
        differences.append(molecular - base)
    expected = np.quantile(differences, [0.025, 0.975], axis=0)
    actual = result["paired_feature_comparisons"]["cox"]["brier"]
    assert_allclose([row["ci_lower"] for row in actual], expected[0])
    assert_allclose([row["ci_upper"] for row in actual], expected[1])


def test_bootstrap_is_reproducible_and_model_order_independent(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    one = evaluate_models(y_train, y_test, predictions, times, config, np.random.default_rng(82))
    two = evaluate_models(
        y_train,
        y_test,
        dict(reversed(list(predictions.items()))),
        times,
        config,
        np.random.default_rng(82),
    )
    assert one == two


def test_long_followup_can_be_evaluated_at_prespecified_horizon(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    y_test["time"][-1] = y_train["time"].max() * 10
    result = evaluate_models(y_train, y_test, predictions, times, config, np.random.default_rng(65))
    assert result["evaluation_horizon_months"] == times[-1]
    assert result["test_label_administrative_censoring_months"] > times[-1]
    assert result["minimum_censoring_survival_on_evaluation_support"] > 0


def test_grid_is_never_moved_using_heldout_outcomes(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    y_test["time"] = np.linspace(times[0] + 1, 200, len(y_test))
    with pytest.raises(EvaluationSupportError, match="grid was not changed") as caught:
        evaluate_models(y_train, y_test, predictions, times, config, np.random.default_rng(1))
    assert caught.value.diagnostics["n_events_by_first_grid_time"] == 0
    assert caught.value.diagnostics["time_grid_months"] == times.tolist()


def test_all_censored_heldout_sample_fails_with_aggregate_diagnostic(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    y_test["event"] = False
    with pytest.raises(EvaluationSupportError, match="no observed events") as caught:
        evaluate_models(y_train, y_test, predictions, times, config, np.random.default_rng(1))
    assert caught.value.diagnostics == {"n_test": len(y_test)}


def test_undefined_bootstrap_estimates_are_counted_and_ci_is_gated():
    y_train = Surv.from_arrays(
        [True, False, True, False, True, False], [1.0, 2.0, 3.0, 4.0, 5.0, 8.0]
    )
    config = evaluation_config(time_grid_percentiles=[20, 50], bootstrap_replicates=60)
    times = training_time_grid(y_train, config)
    y_test = Surv.from_arrays([True, False], [1.0, 7.0])
    predictions = {
        "cox__clinical": {
            "risk": np.array([1.0, 0.0]),
            "survival": np.array([[0.5, 0.2], [0.9, 0.8]]),
        }
    }
    result = evaluate_models(
        y_train, y_test, predictions, times, config, np.random.default_rng(493)
    )
    assert result["status"] == "bootstrap_gate_failed"
    assert not result["reportable"]
    uno = result["models"]["cox__clinical"]["uno_c"]
    assert uno["n_bootstrap_failed"] > 0
    assert uno["ci_lower"] is None and uno["ci_upper"] is None
    json.dumps(result, allow_nan=False)


def test_invalid_predictions_and_nonprotocol_grid_are_rejected(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    with pytest.raises(ValueError, match="training-only time grid"):
        evaluate_models(y_train, y_test, predictions, times + 0.1, config, np.random.default_rng(1))
    predictions["cox__clinical"]["survival"][0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        evaluate_models(y_train, y_test, predictions, times, config, np.random.default_rng(1))


def test_training_with_no_supported_time_range_is_rejected():
    y_train = Surv.from_arrays([True, False, False], [1.0, 1.0, 1.0])
    with pytest.raises(EvaluationSupportError, match="two distinct"):
        training_time_grid(y_train, evaluation_config())


def test_bootstrap_preserves_primary_strata_even_for_all_event_os_labels():
    y_train = Surv.from_arrays(
        [True, False, True, False, True, False], [1.0, 2.0, 3.0, 4.0, 5.0, 8.0]
    )
    config = evaluation_config(time_grid_percentiles=[20, 50], bootstrap_replicates=60)
    times = training_time_grid(y_train, config)
    # Both OS events are observed, but the primary DSS split strata differ.
    y_test = Surv.from_arrays([True, True], [1.0, 7.0])
    predictions = {
        "cox__clinical": {
            "risk": np.array([1.0, 0.0]),
            "survival": np.array([[0.5, 0.2], [0.9, 0.8]]),
        }
    }
    result = evaluate_models(
        y_train,
        y_test,
        predictions,
        times,
        config,
        np.random.default_rng(7),
        bootstrap_strata=np.array([True, False]),
    )
    assert result["bootstrap"]["sampling_design"] == "fixed_primary_event_strata"
    assert result["bootstrap"]["stratum_counts"] == {"False": 1, "True": 1}
    model = result["models"]["cox__clinical"]
    for row in [model["uno_c"], model["ibs"], model["mean_auc"], *model["brier"], *model["auc"]]:
        assert row["n_bootstrap_success"] == 60
        assert row["n_bootstrap_failed"] == 0
        assert row["ci_lower"] == row["ci_upper"]
    assert result["status"] == "ok"


def test_passing_gate_with_undefined_draws_is_conditionally_approximate():
    y_train = Surv.from_arrays(
        [True, False, True, False, True, False], [1.0, 2.0, 3.0, 4.0, 5.0, 8.0]
    )
    config = evaluation_config(time_grid_percentiles=[20, 50], bootstrap_replicates=10)
    times = training_time_grid(y_train, config)
    y_test = Surv.from_arrays([True, False], [1.0, 7.0])
    predictions = {
        "cox__clinical": {
            "risk": np.array([1.0, 0.0]),
            "survival": np.array([[0.5, 0.2], [0.9, 0.8]]),
        }
    }
    predictions["cox__clinical_molecular"] = predictions["cox__clinical"]

    class OneUndefinedDraw:
        calls = 0

        def integers(self, low, high, size):
            self.calls += 1
            return np.zeros(size, dtype=int) if self.calls == 1 else np.arange(size)

    result = evaluate_models(y_train, y_test, predictions, times, config, OneUndefinedDraw())
    assert result["status"] == "conditional_approximate"
    assert result["reportable"]
    assert result["bootstrap"]["sampling_design"] == "iid_approximation_for_stratified_split"
    for model in [*result["models"].values(), *result["paired_feature_comparisons"].values()]:
        assert model["status"] == "conditional_approximate"
        assert model["reportable"]
        assert model["uno_c"]["status"] == "conditional_approximate"
        assert model["uno_c"]["conditional_on_defined_replicates"]
        assert model["uno_c"]["n_bootstrap_failed"] == 1
        assert model["uno_c"]["ci_lower"] is not None


def test_bootstrap_strata_must_be_aligned_primary_event_booleans(evaluation_data):
    y_train, y_test, predictions, times, config = evaluation_data
    for strata in [np.array([True]), np.arange(len(y_test))]:
        with pytest.raises(ValueError, match="aligned boolean"):
            evaluate_models(
                y_train,
                y_test,
                predictions,
                times,
                config,
                np.random.default_rng(8),
                bootstrap_strata=strata,
            )
