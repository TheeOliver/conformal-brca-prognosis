.DEFAULT_GOAL := help
.PHONY: help setup hooks check download lint format test test-all data split fit conformal eval eda uncertainty export all clean-outputs

CONFIG ?= config/default.yaml
RUN    ?= uv run
SPLIT_ARGS ?=
STAGE_ARGS ?=

help: ## Show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

setup: hooks ## Create the venv, resolve deps, install git hooks
	uv sync

hooks: ## Install the agent-independent git hooks (.githooks/)
	git config core.hooksPath .githooks

check: lint test ## Definition of done for any agent: lint + fast tests

lint: ## Ruff check (no writes)
	$(RUN) ruff check .
	$(RUN) ruff format --check .

format: ## Ruff autoformat + autofix
	$(RUN) ruff format .
	$(RUN) ruff check --fix .

test: ## Fast tests only -- this is what the Stop hook gate runs
	$(RUN) pytest -m "not slow"

test-all: ## Every test including the slow MCMC ones (submit with sbatch)
	$(RUN) pytest

# --- Pipeline stages ---------------------------------------------------------
# Each stage is resumable and writes a run manifest to outputs/manifests/.
# Stages 03-05 are COMPUTE: submit them with sbatch, do not run on the login
# node. See docs/slurm-guide.md. A PreToolUse hook enforces this for Claude.

download: ## 00 - fetch public METABRIC clinical tables (network only)
	$(RUN) python scripts/00_download_data.py --config $(CONFIG)

data: ## 01 - load raw METABRIC, validate schema, write processed table
	$(RUN) python scripts/01_prepare_data.py --config $(CONFIG) $(STAGE_ARGS)

split: ## 02 - draw the ONE train/calibration/test split and freeze it
	$(RUN) python scripts/02_make_splits.py --config $(CONFIG) $(SPLIT_ARGS)

fit: ## 03 - fit Cox, Bayesian AFT and RSF on train, both feature sets [COMPUTE]
	$(RUN) python scripts/03_fit_models.py --config $(CONFIG) $(STAGE_ARGS)

conformal: ## 04 - calibrate conformal intervals on the calibration set [COMPUTE]
	$(RUN) python scripts/04_conformal.py --config $(CONFIG) $(STAGE_ARGS)

eval: ## 05 - score on test ONCE; write outputs/metrics/*.json [COMPUTE]
	$(RUN) python scripts/05_evaluate.py --config $(CONFIG) $(STAGE_ARGS)

eda: ## Training-only EDA [SLURM]
	$(RUN) python scripts/eda.py --config $(CONFIG) $(STAGE_ARGS)

uncertainty: ## Hypothetical-profile Bayesian uncertainty [SLURM]
	$(RUN) python scripts/training_uncertainty.py --config $(CONFIG) $(STAGE_ARGS)

export: ## Export saved aggregate metrics to PDF/PNG/LaTeX [SLURM]
	$(RUN) python scripts/06_export_results.py --config $(CONFIG) $(STAGE_ARGS)

.NOTPARALLEL: all
all: data split eda fit uncertainty conformal eval export ## Full pipeline [COMPUTE -- use sbatch]

clean-outputs: ## Delete generated artefacts (never touches data/)
	rm -rf outputs/figures/* outputs/tables/* outputs/metrics/* outputs/models/* outputs/manifests/*
