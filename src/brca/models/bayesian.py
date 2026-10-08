"""Censored Bayesian AFT models with a training-only model-checking workflow.

All prediction integrates the stored posterior analytically; scoring new covariates
never constructs or refits a PyMC model. InferenceData uses xarray.DataTree in PyMC 6.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri

from brca.compute_guard import require_compute_node


class BayesianDiagnosticsError(RuntimeError):
    """A sampled fit failed the configured acceptance gate."""


class PriorPredictiveError(ValueError):
    """The prior predicts too much mass outside the configured time bounds."""


class BayesianInitializationError(ValueError):
    """A proposed chain start has a nonfinite likelihood or gradient."""


def _namespace(value: Any) -> Any:
    if isinstance(value, dict):
        return SimpleNamespace(**{key: _namespace(item) for key, item in value.items()})
    if isinstance(value, list):
        return [_namespace(item) for item in value]
    return value


def _plain(value: Any) -> Any:
    if isinstance(value, SimpleNamespace):
        return {key: _plain(item) for key, item in vars(value).items()}
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return value


def _dataset(idata: Any, group: str) -> Any:
    node = idata[group]
    return node.to_dataset() if hasattr(node, "to_dataset") else node


def _validate_y(y: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    if y.dtype.names != ("event", "time") or y.dtype["event"].kind != "b":
        raise ValueError("y must have boolean event and floating time fields, in that order")
    if y.dtype["time"].kind != "f" or len(y) != n:
        raise ValueError("y must have floating times and match the number of feature rows")
    if not np.isfinite(y["time"]).all() or np.any(y["time"] <= 0):
        raise ValueError("Survival times must be finite and strictly positive months")
    if not np.any(y["event"]):
        raise ValueError("At least one observed event is required for an AFT fit")
    return y["event"], y["time"]


class BayesianAFTModel:
    """Weibull primary model or prespecified lognormal family sensitivity.

    X must already be encoded and scaled by a preprocessing fit on training data.
    Positive coefficients imply longer characteristic survival, hence lower risk.
    """

    def __init__(
        self,
        config: SimpleNamespace,
        rng: np.random.Generator,
        family: str = "weibull",
        prior_variant: str = "default",
    ) -> None:
        if family not in {"weibull", "lognormal"}:
            raise ValueError("family must be weibull or lognormal")
        if prior_variant not in {"default", "sensitivity"}:
            raise ValueError("prior_variant must be default or sensitivity")
        self.config = config
        self.rng = rng
        self.family = family
        self.prior_variant = prior_variant
        self.idata = None
        self.diagnostics_: dict[str, Any] = {}
        self.workflow_: dict[str, Any] = {}
        self.feature_names_: list[str] = []
        self._posterior: dict[str, np.ndarray] | None = None

    def _x(self, X: pd.DataFrame, *, fitting: bool = False) -> np.ndarray:
        if not isinstance(X, pd.DataFrame) or X.shape[1] == 0 or len(X) == 0:
            raise ValueError("X must be a nonempty numeric DataFrame")
        if X.columns.has_duplicates:
            raise ValueError("Feature names must be unique")
        if fitting:
            self.feature_names_ = list(X.columns)
        elif list(X.columns) != self.feature_names_:
            raise ValueError("Prediction columns must match fitted columns in the same order")
        matrix = X.to_numpy(dtype=float)
        if not np.isfinite(matrix).all():
            raise ValueError("X must contain finite, already imputed numeric values")
        return matrix

    def _build_model(self, X: np.ndarray, event: np.ndarray, time: np.ndarray) -> Any:
        import pymc as pm
        import pytensor.tensor as pt

        cfg = self.config.bayesian
        sensitivity = self.prior_variant == "sensitivity"
        intercept_sd = cfg.sensitivity_intercept_sigma if sensitivity else cfg.prior_intercept_sigma
        beta_sd = cfg.sensitivity_coefficient_sigma if sensitivity else cfg.prior_coefficient_sigma
        with pm.Model(coords={"feature": self.feature_names_}) as model:
            covariates = pm.Data("covariates", X, dims=("observation", "feature"))
            intercept = pm.Normal("intercept", mu=cfg.prior_intercept_mu, sigma=intercept_sd)
            beta = pm.Normal("beta", mu=0, sigma=beta_sd, dims="feature")
            location = intercept + pt.dot(covariates, beta)
            if self.family == "weibull":
                shape = pm.LogNormal(
                    "shape", mu=cfg.prior_log_shape_mu, sigma=cfg.prior_log_shape_sigma
                )
                distribution = pm.Weibull.dist(alpha=shape, beta=pt.exp(location))
            else:
                sigma = pm.HalfNormal("sigma", sigma=cfg.prior_lognormal_sigma_scale)
                distribution = pm.LogNormal.dist(mu=location, sigma=sigma)
            # Uncensored events have an infinite upper bound and contribute their
            # density. Censored rows sit exactly at their upper bound and contribute
            # P(T >= observed follow-up), never an event-time density.
            upper = pm.Data("censoring_upper", np.where(event, np.inf, time))
            pm.Censored("observed_time", distribution, lower=None, upper=upper, observed=time)
        return model

    @property
    def _parameters(self) -> list[str]:
        return ["intercept", "beta", "shape" if self.family == "weibull" else "sigma"]

    def _draws(self, idata: Any, group: str, limit: int | None = None) -> dict[str, np.ndarray]:
        ds = _dataset(idata, group)
        result = {}
        for name in self._parameters:
            values = ds[name].transpose("chain", "draw", ...).values
            result[name] = values.reshape((-1, *values.shape[2:]))
        n = len(result["intercept"])
        if limit is not None and n > limit:
            # A fixed evenly spaced subset makes prediction reproducible and bounded
            # in memory, without consuming or changing any sampling RNG stream.
            indices = np.linspace(0, n - 1, limit, dtype=int)
            result = {name: values[indices] for name, values in result.items()}
        return result

    def _event_time_draws(
        self, X: np.ndarray, draws: dict[str, np.ndarray], rng: np.random.Generator
    ) -> np.ndarray:
        location = draws["intercept"][:, None] + draws["beta"] @ X.T
        if self.family == "weibull":
            return np.exp(location) * rng.weibull(draws["shape"][:, None], location.shape)
        return np.exp(location + draws["sigma"][:, None] * rng.normal(size=location.shape))

    def _initial_points(self, rng: np.random.Generator) -> list[dict[str, Any]]:
        """Bound chain starts near prior centers and check numerical derivatives.

        PyMC's default unit jitter can create extreme linear predictors on a
        standardized but skewed covariate. A finite censored log probability
        alone does not ensure its autodifferentiated gradient is finite.
        """
        cfg = self.config.bayesian
        sampler = self.config.mcmc
        fraction = float(sampler.initialization_jitter_fraction)
        if sampler.init != "adapt_diag" or not np.isfinite(fraction) or not 0 <= fraction <= 1:
            raise ValueError("Initialization requires adapt_diag and a jitter fraction in [0, 1]")
        sensitivity = self.prior_variant == "sensitivity"
        intercept_sd = cfg.sensitivity_intercept_sigma if sensitivity else cfg.prior_intercept_sigma
        beta_sd = cfg.sensitivity_coefficient_sigma if sensitivity else cfg.prior_coefficient_sigma
        logp = self.model_.compile_logp()
        gradient = self.model_.compile_dlogp()
        starts = []
        checks = []
        for chain in range(sampler.chains):
            point = {
                "intercept": cfg.prior_intercept_mu
                + rng.uniform(-fraction, fraction) * intercept_sd,
                "beta": rng.uniform(-fraction, fraction, len(self.feature_names_)) * beta_sd,
            }
            if self.family == "weibull":
                point["shape"] = np.exp(
                    cfg.prior_log_shape_mu
                    + rng.uniform(-fraction, fraction) * cfg.prior_log_shape_sigma
                )
            else:
                point["sigma"] = cfg.prior_lognormal_sigma_scale * np.exp(
                    rng.uniform(-fraction, fraction)
                )
            transformed = {
                (f"{name}_log__" if name in {"shape", "sigma"} else name): np.asarray(
                    np.log(value) if name in {"shape", "sigma"} else value
                )
                for name, value in point.items()
            }
            probability = float(logp(transformed))
            derivative = gradient(transformed)
            if not np.isfinite(probability) or not np.isfinite(derivative).all():
                raise BayesianInitializationError(
                    f"Chain {chain} start has nonfinite log probability or gradients; "
                    "posterior sampling was not started."
                )
            checks.append(
                {
                    "chain": chain,
                    "log_probability": probability,
                    "max_absolute_gradient": float(np.max(np.abs(derivative))),
                    "finite_gradients": True,
                }
            )
            starts.append(point)
        self.workflow_["initialization"] = {
            "method": sampler.init,
            "jitter_fraction": fraction,
            "chain_checks": checks,
        }
        return starts

    def _survival_draws(
        self, X: np.ndarray, times: np.ndarray, draws: dict[str, np.ndarray]
    ) -> np.ndarray:
        location = draws["intercept"][:, None, None] + (draws["beta"] @ X.T)[:, :, None]
        with np.errstate(divide="ignore", over="ignore"):
            log_time = np.log(times)[None, None, :]
            if self.family == "weibull":
                return np.exp(-np.exp(draws["shape"][:, None, None] * (log_time - location)))
            return ndtr((location - log_time) / draws["sigma"][:, None, None])

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> BayesianAFTModel:
        """Check priors, sample, enforce diagnostics, then run a training-only PPC."""
        require_compute_node("Bayesian AFT fitting")
        import pymc as pm

        self.idata = None
        self._posterior = None
        self.diagnostics_ = {}
        self.workflow_ = {}
        matrix = self._x(X, fitting=True)
        event, time = _validate_y(y, len(matrix))
        self.workflow_["constant_training_features"] = [
            name for j, name in enumerate(self.feature_names_) if np.ptp(matrix[:, j]) == 0
        ]
        cfg = self.config.bayesian
        sampler = self.config.mcmc
        streams = self.rng.spawn(5)
        self.model_ = self._build_model(matrix, event, time)
        with self.model_:
            prior = pm.sample_prior_predictive(
                draws=cfg.prior_predictive_draws,
                var_names=self._parameters,
                random_seed=streams[0],
            )
        prior_draws = self._draws(prior, "prior")
        prior_times = self._event_time_draws(matrix, prior_draws, streams[1])
        bounds = np.asarray(cfg.prior_plausible_time_bounds_months, dtype=float)
        extreme = (
            (~np.isfinite(prior_times)) | (prior_times < bounds[0]) | (prior_times > bounds[1])
        )
        tail = (1 - cfg.predictive_interval_probability) / 2
        histogram_edges = np.geomspace(*bounds, cfg.prior_time_histogram_bins + 1)
        histogram_counts, _ = np.histogram(prior_times, bins=histogram_edges)
        self.workflow_["prior_predictive"] = {
            "draws": cfg.prior_predictive_draws,
            "time_quantiles_months": np.quantile(prior_times, [tail, 0.5, 1 - tail]).tolist(),
            "time_histogram_edges_months": histogram_edges.tolist(),
            "time_histogram_probability": (histogram_counts / prior_times.size).tolist(),
            "plausible_bounds_months": bounds.tolist(),
            "extreme_fraction": float(np.mean(extreme)),
            "max_extreme_fraction": cfg.prior_max_extreme_fraction,
            "passed": bool(np.mean(extreme) <= cfg.prior_max_extreme_fraction),
        }
        if not self.workflow_["prior_predictive"]["passed"]:
            raise PriorPredictiveError(
                json.dumps(self.workflow_["prior_predictive"], sort_keys=True)
            )
        initial_points = self._initial_points(streams[4])
        with self.model_:
            self.idata = pm.sample(
                init=sampler.init,
                initvals=initial_points,
                chains=sampler.chains,
                cores=min(sampler.chains, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
                tune=sampler.tune,
                draws=sampler.draws,
                target_accept=sampler.target_accept,
                random_seed=streams[2],
                progressbar=False,
                return_inferencedata=True,
                blas_cores=min(sampler.chains, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
            )
        self.idata["prior"] = prior["prior"]
        self.check_diagnostics()
        self._posterior = self._draws(self.idata, "posterior", cfg.predictive_draws)
        self.workflow_["training_ppc"] = self._ppc(matrix, y, prior_draws, streams[3])
        self.workflow_["coefficients"] = self.coefficient_summary()
        self.workflow_["family"] = self.family
        self.workflow_["prior_variant"] = self.prior_variant
        return self

    def check_diagnostics(self) -> dict[str, Any]:
        """Raise for every failed gate, including non-finite diagnostic estimates."""
        import arviz as az

        if self.idata is None:
            raise RuntimeError("Model has not been sampled")
        ds = _dataset(self.idata, "posterior")[self._parameters]
        rhat_ds = az.rhat(ds)
        bulk_ds = az.ess(ds, method="bulk")
        tail_ds = az.ess(ds, method="tail")
        rhat = rhat_ds.to_array().values
        bulk = bulk_ds.to_array().values
        tail = tail_ds.to_array().values
        stats = _dataset(self.idata, "sample_stats")
        values = {
            "max_rhat": float(np.max(rhat)),
            "min_ess_bulk": float(np.min(bulk)),
            "min_ess_tail": float(np.min(tail)),
            "divergences": int(np.asarray(stats["diverging"]).sum()),
        }
        gate = self.config.mcmc.diagnostics
        checks = {
            "rhat": bool(np.isfinite(rhat).all() and values["max_rhat"] <= gate.max_rhat),
            "ess_bulk": bool(
                np.isfinite(bulk).all() and values["min_ess_bulk"] >= gate.min_ess_bulk
            ),
            "ess_tail": bool(
                np.isfinite(tail).all() and values["min_ess_tail"] >= gate.min_ess_tail
            ),
            "divergences": values["divergences"] <= gate.max_divergences,
        }
        self.diagnostics_ = {
            **values,
            "checks": checks,
            "passed": all(checks.values()),
            "thresholds": _plain(gate),
            "parameters": {
                name: {
                    "rhat": rhat_ds[name].values.tolist(),
                    "ess_bulk": bulk_ds[name].values.tolist(),
                    "ess_tail": tail_ds[name].values.tolist(),
                }
                for name in self._parameters
            },
        }
        if not self.diagnostics_["passed"]:
            raise BayesianDiagnosticsError(json.dumps(self.diagnostics_, sort_keys=True))
        return self.diagnostics_

    def _fitted(self) -> dict[str, np.ndarray]:
        if self._posterior is None or not self.diagnostics_.get("passed", False):
            raise RuntimeError("A fitted posterior passing the diagnostics gate is required")
        return self._posterior

    def predict_risk(self, X: pd.DataFrame) -> np.ndarray:
        draws = self._fitted()
        matrix = self._x(X)
        return -(draws["intercept"].mean() + matrix @ draws["beta"].mean(axis=0))

    def predict_survival_function(self, X: pd.DataFrame, times: np.ndarray) -> np.ndarray:
        draws = self._fitted()
        matrix = self._x(X)
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or np.any(times < 0) or not np.isfinite(times).all():
            raise ValueError("times must be a finite nonnegative one-dimensional array")
        # Loop over times: a full draws x patients x grid tensor can exceed RAM
        # when bootstrap evaluation uses dense grids.
        return np.column_stack(
            [
                self._survival_draws(matrix, np.array([time]), draws).mean(axis=0)[:, 0]
                for time in times
            ]
        )

    def predict_survival_bands(self, X: pd.DataFrame, times: np.ndarray) -> dict[str, np.ndarray]:
        """Pointwise parameter uncertainty in S(t|x), not individual time intervals."""
        matrix = self._x(X)
        tail = (1 - self.config.bayesian.predictive_interval_probability) / 2
        draws = self._survival_draws(matrix, np.asarray(times), self._fitted())
        return {
            "mean": draws.mean(axis=0),
            "lower": np.quantile(draws, tail, axis=0),
            "upper": np.quantile(draws, 1 - tail, axis=0),
        }

    def predict_quantiles(self, X: pd.DataFrame, probabilities: np.ndarray) -> np.ndarray:
        """Quantiles of the posterior mixture of event times, including future noise."""
        draws = self._fitted()
        matrix = self._x(X)
        probabilities = np.atleast_1d(np.asarray(probabilities, dtype=float))
        if probabilities.ndim != 1 or np.any((probabilities <= 0) | (probabilities >= 1)):
            raise ValueError("probabilities must lie strictly between zero and one")
        if not np.isfinite(probabilities).all():
            raise ValueError("probabilities must be finite")
        location = draws["intercept"][:, None] + draws["beta"] @ matrix.T
        result = np.empty((len(matrix), len(probabilities)))
        for j, probability in enumerate(probabilities):
            if self.family == "weibull":
                component_log_quantiles = (
                    location + np.log(-np.log1p(-probability)) / draws["shape"][:, None]
                )
            else:
                component_log_quantiles = location + draws["sigma"][:, None] * ndtri(probability)
            lower = component_log_quantiles.min(axis=0)
            upper = component_log_quantiles.max(axis=0)
            # A mixture p-quantile lies between the smallest and largest component
            # p-quantiles. Bisection on log time avoids an arbitrary horizon or cap.
            for _ in range(self.config.bayesian.quantile_bisection_iterations):
                middle = (lower + upper) / 2
                if self.family == "weibull":
                    with np.errstate(over="ignore"):
                        cdf = -np.expm1(-np.exp(draws["shape"][:, None] * (middle - location)))
                else:
                    cdf = ndtr((middle - location) / draws["sigma"][:, None])
                below = cdf.mean(axis=0) < probability
                lower = np.where(below, middle, lower)
                upper = np.where(below, upper, middle)
            result[:, j] = np.exp((lower + upper) / 2)
        return result

    def coefficient_summary(self) -> dict[str, Any]:
        self._fitted()
        # Coefficient intervals use every retained draw. The bounded subset is
        # exclusively a prediction-cost setting, not a reduction of inference.
        draws = self._draws(self.idata, "posterior")
        tail = (1 - self.config.bayesian.predictive_interval_probability) / 2
        return {
            name: {
                "mean_log_time_ratio": float(draws["beta"][:, j].mean()),
                "time_ratio_interval": np.exp(
                    np.quantile(draws["beta"][:, j], [tail, 1 - tail])
                ).tolist(),
            }
            for j, name in enumerate(self.feature_names_)
        }

    def _ppc(
        self,
        X: np.ndarray,
        y: np.ndarray,
        prior: dict[str, np.ndarray],
        rng: np.random.Generator,
    ) -> dict[str, Any]:
        from sksurv.nonparametric import kaplan_meier_estimator

        cfg = self.config.bayesian
        times = np.asarray(cfg.ppc_time_grid_months, dtype=float)
        times = times[times < y["time"].max()]
        km_times, km, km_ci = kaplan_meier_estimator(
            y["event"],
            y["time"],
            conf_level=cfg.predictive_interval_probability,
            conf_type="log-log",
        )
        indices = np.searchsorted(km_times, times, side="right")
        observed = np.r_[1.0, km][indices]
        observed_ci = np.column_stack([np.ones(2), km_ci])[:, indices]
        # Curves marginalize over the empirical TRAIN covariates, retaining one
        # curve per posterior draw for epistemic uncertainty.
        curves = self._survival_draws(X, times, self._fitted()).mean(axis=1)
        prior_curves = self._survival_draws(X, times, prior).mean(axis=1)
        tail = (1 - cfg.predictive_interval_probability) / 2
        event_times = self._event_time_draws(X, self._fitted(), rng)
        c_times, c_survival = kaplan_meier_estimator(y["event"], y["time"], reverse=True)
        censor_cdf = 1 - c_survival
        c_indices = np.searchsorted(censor_cdf, rng.uniform(size=event_times.shape), side="left")
        # Unidentified censoring-tail mass is left at infinity, not made into
        # artificial deaths at the longest follow-up.
        censor_times = np.r_[c_times, np.inf][c_indices]
        simulated_event_fraction = (event_times <= censor_times).mean(axis=1)
        return {
            "split": "train",
            "n": len(y),
            "times_months": times.tolist(),
            "km_survival": observed.tolist(),
            "km_lower": observed_ci[0].tolist(),
            "km_upper": observed_ci[1].tolist(),
            "posterior_survival_mean": curves.mean(axis=0).tolist(),
            "posterior_survival_lower": np.quantile(curves, tail, axis=0).tolist(),
            "posterior_survival_upper": np.quantile(curves, 1 - tail, axis=0).tolist(),
            "prior_survival_lower": np.quantile(prior_curves, tail, axis=0).tolist(),
            "prior_survival_upper": np.quantile(prior_curves, 1 - tail, axis=0).tolist(),
            "absolute_survival_discrepancy": np.abs(curves.mean(axis=0) - observed).tolist(),
            "observed_event_fraction": float(y["event"].mean()),
            "replicated_event_fraction_interval": np.quantile(
                simulated_event_fraction, [tail, 1 - tail]
            ).tolist(),
            "censoring_assumption": (
                "Replicated censoring uses training reverse-KM and marginal independent "
                "censoring; this descriptive check does not establish conditional independence."
            ),
            "predictive_time_quantiles_months": np.quantile(
                event_times, [tail, 0.5, 1 - tail]
            ).tolist(),
            "prior_times_summary": self.workflow_["prior_predictive"],
        }

    def save(self, path: str | Path) -> None:
        """Persist posterior/prior draws plus aggregate diagnostics and metadata."""
        if self.idata is None:
            raise RuntimeError("No sampled posterior to save")
        path = Path(path)
        path = path if path.suffix == ".nc" else path.with_suffix(".nc")
        path.parent.mkdir(parents=True, exist_ok=True)
        self.idata.to_netcdf(path, engine="h5netcdf")
        metadata = {
            "config": _plain(self.config),
            "family": self.family,
            "prior_variant": self.prior_variant,
            "feature_names": self.feature_names_,
            "diagnostics": self.diagnostics_,
            "workflow": self.workflow_,
        }
        path.with_suffix(".json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> BayesianAFTModel:
        import xarray as xr

        path = Path(path)
        path = path if path.suffix == ".nc" else path.with_suffix(".nc")
        metadata = json.loads(path.with_suffix(".json").read_text())
        config = _namespace(metadata["config"])
        model = cls(
            config,
            np.random.default_rng(config.seed),
            metadata["family"],
            metadata["prior_variant"],
        )
        with xr.open_datatree(path, engine="h5netcdf") as stored:
            model.idata = stored.load()
        model.feature_names_ = metadata["feature_names"]
        model.diagnostics_ = metadata["diagnostics"]
        model.workflow_ = metadata["workflow"]
        model.check_diagnostics()
        model._posterior = model._draws(model.idata, "posterior", config.bayesian.predictive_draws)
        return model
