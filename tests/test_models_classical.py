"""Classical model contracts on generated data; no real patient data are loaded."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from sklearn.exceptions import ConvergenceWarning, NotFittedError
from sksurv.preprocessing import OneHotEncoder
from sksurv.util import Surv

from brca.models.base import step_survival_on_grid, step_survival_quantiles
from brca.models.cox import CoxPHModel, proportional_hazards_test
from brca.models.rsf import RandomSurvivalForestModel


def cox_config(**overrides):
    values = dict(
        alpha=0.0,
        ties="breslow",
        n_iter=100,
        tol=1e-9,
        ph_diagnostic_transform="log",
        ph_diagnostic_alpha=0.05,
    )
    return SimpleNamespace(**(values | overrides))


def rsf_config(**overrides):
    values = dict(
        n_estimators=25, min_samples_split=10, min_samples_leaf=5, max_features="sqrt", n_jobs=1
    )
    return SimpleNamespace(**(values | overrides))


@pytest.fixture(scope="module")
def cohort():
    data = pd.read_csv(Path(__file__).parent / "fixtures" / "synthetic_metabric.csv")
    X = data.drop(columns=["patient_id", "survival_months", "event"])
    for column in X.select_dtypes(include=["object", "str"]).columns:
        X[column] = X[column].astype("category")
    X = OneHotEncoder().fit_transform(X)
    X = (X - X.mean()) / X.std()
    y = Surv.from_arrays(data.event.to_numpy(bool), data.survival_months.to_numpy(float))
    return X, y


@pytest.fixture(params=[CoxPHModel, RandomSurvivalForestModel])
def model(request):
    config = cox_config() if request.param is CoxPHModel else rsf_config()
    return request.param(config, rng=np.random.default_rng(724))


def test_fit_shared_grid_quantiles_and_persistence(model, cohort, tmp_path):
    X, y = cohort
    assert model.fit(X, y) is model
    times = np.concatenate(([0.0], np.linspace(1, y["time"].max() * 2, 45)))
    predicted = model.predict_survival_function(X.iloc[:5], times)
    assert predicted.shape == (5, len(times))
    assert np.isfinite(predicted).all()
    assert ((predicted >= 0) & (predicted <= 1)).all()
    assert (np.diff(predicted, axis=1) <= 0).all()
    assert_array_equal(predicted[:, 0], 1)
    assert_allclose(predicted[:, -1], predicted[:, -2])
    quantiles = model.predict_quantiles(X.iloc[:5], np.array([0.1, 0.5, 0.99]))
    assert quantiles.shape == (5, 3)
    assert (quantiles[:, 1:] >= quantiles[:, :-1]).all()
    artifact = tmp_path / "model.pkl"
    model.save(artifact)
    restored = type(model).load(artifact)
    assert_array_equal(restored.predict_survival_function(X.iloc[:5], times), predicted)
    assert_array_equal(restored.predict_risk(X.iloc[:5]), model.predict_risk(X.iloc[:5]))
    assert_array_equal(restored.predict_quantiles(X.iloc[:5], [0.1, 0.5, 0.99]), quantiles)


def test_same_seed_reproduces_predictions(model, cohort):
    X, y = cohort
    second = type(model)(model.config, rng=np.random.default_rng(724))
    model.fit(X, y)
    second.fit(X, y)
    assert_array_equal(model.predict_risk(X.iloc[:10]), second.predict_risk(X.iloc[:10]))


def test_risk_ordering_on_known_synthetic_signal(model):
    rng = np.random.default_rng(5281)
    X = pd.DataFrame({"risk_covariate": rng.normal(size=350)})
    event_time = rng.exponential(60 * np.exp(-X.risk_covariate))
    censor_time = rng.exponential(200, len(X))
    y = Surv.from_arrays(event_time <= censor_time, np.minimum(event_time, censor_time))
    model.fit(X, y)
    risk = model.predict_risk(pd.DataFrame({"risk_covariate": [-1.5, 1.5]}))
    assert risk[1] > risk[0]


def test_feature_contract_and_unfitted_errors(model, cohort, tmp_path):
    X, y = cohort
    with pytest.raises(NotFittedError):
        model.predict_risk(X)
    with pytest.raises(NotFittedError):
        model.save(tmp_path / "unfitted.pkl")
    model.fit(X, y)
    with pytest.raises(ValueError, match="columns and order"):
        model.predict_risk(X.iloc[:, ::-1])
    with pytest.raises(ValueError, match="columns and order"):
        model.predict_risk(X.rename(columns={X.columns[0]: "wrong_name"}))
    invalid = X.copy()
    invalid.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        model.predict_risk(invalid)


def test_integer_event_labels_and_invalid_time_rejected(model, cohort):
    X, y = cohort
    malformed = np.empty(len(y), dtype=[("event", "i4"), ("time", "f8")])
    malformed["event"], malformed["time"] = y["event"], y["time"]
    with pytest.raises(ValueError, match="boolean event"):
        model.fit(X, malformed)
    invalid = y.copy()
    invalid["time"][0] = 0
    with pytest.raises(ValueError, match="strictly positive"):
        model.fit(X, invalid)
    invalid["time"][0] = np.inf
    with pytest.raises(ValueError, match="finite"):
        model.fit(X, invalid)
    invalid = y.copy()
    invalid["event"] = False
    with pytest.raises(ValueError, match="at least one observed event"):
        model.fit(X, invalid)


def test_step_curve_boundaries_and_unidentified_quantiles():
    knots = np.array([2.0, 4.0, 8.0])
    survival = np.array([[0.9, 0.7, 0.6], [1.0, 0.5, 0.2]])
    assert_array_equal(
        step_survival_on_grid(knots, survival, [0, 1, 2, 3, 4, 9]),
        [[1, 1, 0.9, 0.9, 0.7, 0.6], [1, 1, 1, 1, 0.5, 0.2]],
    )
    assert_array_equal(
        step_survival_quantiles(knots, survival, [0.2, 0.5, 0.9]),
        [[4, np.inf, np.inf], [4, 4, np.inf]],
    )
    with pytest.raises(ValueError, match="nonnegative"):
        step_survival_on_grid(knots, survival, [-1])
    with pytest.raises(ValueError, match="finite"):
        step_survival_on_grid(knots, survival, [np.inf])
    with pytest.raises(ValueError, match="strictly within"):
        step_survival_quantiles(knots, survival, [0, 1])


def test_cox_constant_columns_removed_from_train_only(cohort):
    X, y = cohort
    X = X.assign(unseen_category=0.0, constant=1.0)
    model = CoxPHModel(cox_config(), rng=np.random.default_rng(72)).fit(X, y)
    assert model.constant_features_ == ("unseen_category", "constant")
    assert "unseen_category" not in model.ph_diagnostic_["columns"]
    changed = X.iloc[:2].copy()
    changed["unseen_category"] = 1.0
    assert_array_equal(model.predict_risk(changed), model.predict_risk(X.iloc[:2]))
    assert model.ph_diagnostic_["training_only"]
    assert model.ph_diagnostic_["n_training"] == len(X)
    json.dumps(model.ph_diagnostic_, allow_nan=False)


def test_penalized_ph_inference_unavailable(cohort):
    X, y = cohort
    model = CoxPHModel(cox_config(alpha=0.1), rng=np.random.default_rng(21)).fit(X, y)
    assert model.ph_diagnostic_["status"] == "unavailable"
    assert "unpenalized" in model.ph_diagnostic_["reason"]


def test_cox_convergence_warning_is_an_error(cohort):
    X, y = cohort
    with pytest.raises(ConvergenceWarning):
        CoxPHModel(cox_config(n_iter=1), rng=np.random.default_rng(21)).fit(X, y)


def test_rsf_rejects_unbounded_cpu_request(cohort):
    X, y = cohort
    with pytest.raises(ValueError, match="positive allocated CPU"):
        RandomSurvivalForestModel(rsf_config(n_jobs=-1), rng=np.random.default_rng(1)).fit(X, y)


def test_loading_wrong_model_type_rejected(cohort, tmp_path):
    X, y = cohort
    artifact = tmp_path / "cox.pkl"
    CoxPHModel(cox_config(), rng=np.random.default_rng(1)).fit(X, y).save(artifact)
    with pytest.raises(ValueError, match="model class"):
        RandomSurvivalForestModel.load(artifact)


@pytest.mark.parametrize("ties", ["breslow", "efron"])
def test_ph_score_matches_finite_difference_likelihood(ties):
    """Independent likelihood derivatives verify tie handling and nuisance projection."""
    X = pd.DataFrame({"x": [-1, 1, 0.5, -0.8, 0.3, 1.3, -0.4, 0.7]})
    y = Surv.from_arrays(
        np.array([True, True, False, True, True, False, True, True]),
        np.array([1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0]),
    )
    fitted = CoxPHModel(cox_config(ties=ties), rng=np.random.default_rng(21)).fit(X, y)
    event_times = np.unique(y["time"][y["event"]])
    center = np.log(y["time"][y["event"]]).mean()

    def likelihood(theta):
        value = 0.0
        for time in event_times:
            event = y["event"] & (y["time"] == time)
            risk = y["time"] >= time
            g = np.log(time) - center
            linear = X.x.to_numpy() * (theta[0] + g * theta[1])
            value += linear[event].sum()
            count = event.sum()
            for k in range(count):
                fraction = k / count if ties == "efron" else 0
                value -= np.log(np.exp(linear[risk]).sum() - fraction * np.exp(linear[event]).sum())
        return value

    theta = np.array([fitted.estimator_.coef_[0], 0.0])
    step = np.eye(2) * 1e-4
    score = np.array([(likelihood(theta + h) - likelihood(theta - h)) / 2e-4 for h in step])
    hessian = np.empty((2, 2))
    for j in range(2):
        for k in range(2):
            hessian[j, k] = (
                -(
                    likelihood(theta + step[j] + step[k])
                    - likelihood(theta + step[j] - step[k])
                    - likelihood(theta - step[j] + step[k])
                    + likelihood(theta - step[j] - step[k])
                )
                / 4e-8
            )
    expected = (score[1] - hessian[1, 0] / hessian[0, 0] * score[0]) ** 2 / (
        hessian[1, 1] - hessian[1, 0] ** 2 / hessian[0, 0]
    )
    assert_allclose(fitted.ph_diagnostic_["global"]["statistic"], expected, rtol=2e-5)


def test_ph_diagnostic_detects_large_time_varying_association():
    rng = np.random.default_rng(3659)
    X = pd.DataFrame({"x": rng.binomial(1, 0.5, size=900).astype(float)})
    exponential = rng.exponential(size=len(X))
    early_rate = 0.02 * np.exp(1.2 * X.x)
    late_rate = 0.02 * np.exp(-1.2 * X.x)
    event_time = np.where(
        exponential < 20 * early_rate,
        exponential / early_rate,
        20 + (exponential - 20 * early_rate) / late_rate,
    )
    y = Surv.from_arrays(event_time <= 180, np.minimum(event_time, 180))
    model = CoxPHModel(cox_config(), rng=np.random.default_rng(83)).fit(X, y)
    assert model.ph_diagnostic_["global"]["p_value"] < 0.001


def test_ph_diagnostic_unidentified_time_slope_is_unavailable():
    X = pd.DataFrame({"x": [-1.0, 0.0, 1.0, 2.0]})
    y = Surv.from_arrays([True, True, False, False], [1.0, 1.0, 2.0, 3.0])
    result = proportional_hazards_test(
        X, y, np.zeros(len(X)), ties="breslow", transform="log", significance_level=0.05
    )
    assert result["status"] == "unavailable"
    assert result["reason"] == "singular_interaction_information"
