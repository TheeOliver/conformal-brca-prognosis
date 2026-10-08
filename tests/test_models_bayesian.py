"""Synthetic-only checks of censored likelihood, posterior prediction and MCMC."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.stats import lognorm, weibull_min

from brca.models.bayesian import (
    BayesianAFTModel,
    BayesianDiagnosticsError,
    BayesianInitializationError,
    _validate_y,
)


def config():
    return SimpleNamespace(
        seed=8927,
        bayesian=SimpleNamespace(
            prior_predictive_draws=250,
            prior_time_histogram_bins=40,
            predictive_draws=300,
            predictive_interval_probability=0.90,
            ppc_time_grid_months=[12, 24, 60, 120, 180],
            prior_intercept_mu=np.log(120),
            prior_intercept_sigma=0.5,
            prior_coefficient_sigma=0.35,
            prior_log_shape_mu=np.log(1.5),
            prior_log_shape_sigma=0.35,
            prior_lognormal_sigma_scale=0.75,
            sensitivity_coefficient_sigma=0.7,
            sensitivity_intercept_sigma=0.75,
            prior_plausible_time_bounds_months=[0.01, 1200],
            prior_max_extreme_fraction=0.20,
            quantile_bisection_iterations=80,
        ),
        mcmc=SimpleNamespace(
            init="adapt_diag",
            initialization_jitter_fraction=0.1,
            chains=4,
            tune=700,
            draws=1000,
            target_accept=0.95,
            diagnostics=SimpleNamespace(
                max_rhat=1.01, min_ess_bulk=400, min_ess_tail=400, max_divergences=0
            ),
        ),
    )


def fixed_model(family="weibull"):
    model = BayesianAFTModel(config(), np.random.default_rng(8927), family=family)
    model.feature_names_ = ["x"]
    model.diagnostics_ = {"passed": True}
    model._posterior = {
        "intercept": np.array([np.log(120)]),
        "beta": np.array([[-0.5]]),
        "shape" if family == "weibull" else "sigma": np.array(
            [1.5 if family == "weibull" else 0.8]
        ),
    }
    return model


@pytest.mark.parametrize("family", ["weibull", "lognormal"])
def test_analytic_survival_quantiles_and_risk_order(family):
    model = fixed_model(family)
    X = pd.DataFrame({"x": [-1, 0, 1]})
    times = np.array([0, 12, 60, 120, 240, 1200.0])
    survival = model.predict_survival_function(X, times)
    assert survival.shape == (3, len(times))
    assert ((survival >= 0) & (survival <= 1)).all()
    assert (np.diff(survival, axis=1) <= 0).all()
    assert model.predict_risk(X)[-1] > model.predict_risk(X)[0]
    probabilities = np.array([0.05, 0.5, 0.95])
    quantiles = model.predict_quantiles(X, probabilities)
    scale = 120 * np.exp(-0.5 * X["x"].to_numpy())
    distribution = weibull_min(c=1.5) if family == "weibull" else lognorm(s=0.8)
    np.testing.assert_allclose(quantiles, scale[:, None] * distribution.ppf(probabilities))
    np.testing.assert_allclose(survival, distribution.sf(times[None, :] / scale[:, None]))


@pytest.mark.parametrize("family", ["weibull", "lognormal"])
def test_mixture_quantiles_invert_posterior_predictive_cdf(family):
    model = fixed_model(family)
    model._posterior = {
        "intercept": np.log([60, 240]),
        "beta": np.array([[-0.7], [0.2]]),
        "shape" if family == "weibull" else "sigma": np.array([0.8, 1.9]),
    }
    X = pd.DataFrame({"x": [0.5]})
    probabilities = np.array([0.05, 0.5, 0.95])
    quantiles = model.predict_quantiles(X, probabilities)
    cdf = 1 - model.predict_survival_function(X, quantiles[0])[0]
    np.testing.assert_allclose(cdf, probabilities, atol=1e-12)
    bands = model.predict_survival_bands(X, np.array([12, 60, 120]))
    assert (bands["lower"] <= bands["upper"]).all()


def test_predictions_reject_wrong_columns_and_unaccepted_fit():
    model = fixed_model()
    with pytest.raises(ValueError, match="columns"):
        model.predict_risk(pd.DataFrame({"wrong": [0]}))
    model.diagnostics_ = {"passed": False}
    with pytest.raises(RuntimeError, match="gate"):
        model.predict_quantiles(pd.DataFrame({"x": [0]}), [0.5])


def test_survival_labels_are_strict():
    y = np.array([(True, 12.0)], dtype=[("event", "?"), ("time", "f8")])
    _validate_y(y, 1)
    for bad in [
        np.array([(1, 12.0)], dtype=[("event", "i4"), ("time", "f8")]),
        np.array([(True, 0.0)], dtype=y.dtype),
        np.array([(True, np.nan)], dtype=y.dtype),
    ]:
        with pytest.raises(ValueError):
            _validate_y(bad, 1)


def test_initializer_rejects_finite_logp_with_nonfinite_gradient():
    model = fixed_model()
    model.model_ = SimpleNamespace(
        compile_logp=lambda: lambda point: -123.0,
        compile_dlogp=lambda: lambda point: np.array([np.nan]),
    )
    with pytest.raises(BayesianInitializationError, match="posterior sampling was not started"):
        model._initial_points(np.random.default_rng(8927))


def test_censored_likelihood_matches_density_and_survival():
    import pymc as pm

    values = np.array([50.0, 80.0])
    distribution = pm.Weibull.dist(alpha=1.5, beta=120.0)
    censored = pm.Censored.dist(distribution, lower=None, upper=np.array([np.inf, 80.0]))
    actual = pm.logp(censored, values).eval()
    expected = [weibull_min.logpdf(50, 1.5, scale=120), weibull_min.logsf(80, 1.5, scale=120)]
    np.testing.assert_allclose(actual, expected)


def test_diagnostics_gate_and_roundtrip(tmp_path):
    import xarray as xr

    model = fixed_model()
    rng = np.random.default_rng(13)
    dims = ("chain", "draw")
    model.idata = xr.DataTree.from_dict(
        {
            "posterior": xr.Dataset(
                {
                    "intercept": (dims, rng.normal(np.log(120), 0.1, (4, 2000))),
                    "beta": ((*dims, "feature"), rng.normal(-0.5, 0.1, (4, 2000, 1))),
                    "shape": (dims, rng.lognormal(np.log(1.5), 0.1, (4, 2000))),
                }
            ),
            "sample_stats": xr.Dataset({"diverging": (dims, np.zeros((4, 2000), dtype=bool))}),
        }
    )
    assert model.check_diagnostics()["passed"]
    model._posterior = model._draws(model.idata, "posterior", config().bayesian.predictive_draws)
    full_beta = model.idata["posterior"]["beta"].values.ravel()
    summary = model.coefficient_summary()["x"]
    np.testing.assert_allclose(summary["mean_log_time_ratio"], full_beta.mean())
    np.testing.assert_allclose(
        summary["time_ratio_interval"], np.exp(np.quantile(full_beta, [0.05, 0.95]))
    )
    model.save(tmp_path / "fit.nc")
    loaded = BayesianAFTModel.load(tmp_path / "fit.nc")
    X = pd.DataFrame({"x": [-1, 0, 1]})
    np.testing.assert_array_equal(
        model.predict_quantiles(X, [0.1, 0.9]), loaded.predict_quantiles(X, [0.1, 0.9])
    )
    model.idata["sample_stats"]["diverging"][0, 0] = True
    with pytest.raises(BayesianDiagnosticsError, match="divergences"):
        model.check_diagnostics()


@pytest.mark.slow
@pytest.mark.parametrize("family", ["weibull", "lognormal"])
def test_synthetic_sampling_workflow(family, tmp_path):
    """Real NUTS sampling must pass the production gate on an identifiable toy fit."""
    rng = np.random.default_rng(89)
    X = pd.DataFrame({"x": rng.normal(size=160)})
    scale = 120 * np.exp(-0.5 * X["x"].to_numpy())
    event_time = scale * (
        rng.weibull(1.5, len(X)) if family == "weibull" else np.exp(rng.normal(0, 0.8, len(X)))
    )
    censoring = rng.uniform(100, 300, len(X))
    y = np.empty(len(X), dtype=[("event", "?"), ("time", "f8")])
    y["event"] = event_time <= censoring
    y["time"] = np.minimum(event_time, censoring)
    model = BayesianAFTModel(config(), np.random.default_rng(8927), family).fit(X, y)
    assert model.diagnostics_["passed"]
    assert model.workflow_["prior_predictive"]["passed"]
    assert model.workflow_["training_ppc"]["split"] == "train"
    assert model.coefficient_summary()["x"]["mean_log_time_ratio"] < 0
    model.save(tmp_path / f"{family}.nc")
    print(f"SYNTHETIC {family}: {model.diagnostics_}")
    if family == "weibull":
        repeat = BayesianAFTModel(config(), np.random.default_rng(8927), family).fit(X, y)
        np.testing.assert_array_equal(model.predict_risk(X), repeat.predict_risk(X))
