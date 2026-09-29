You are setting up the complete .claude/ directory for a new project. Work in plan mode first — show me the full plan before creating any files.

## The Project / Deliverable
Trustworthy prognostic modelling for breast cancer using Bayesian analysis
and conformal prediction

The goal is to investigate how an individual prognosis can be made for
breast cancer patients that is not only accurate, but also trustworthy.
The project connects the courses you've taken with the direction you want
to pursue for your degree, with an emphasis on uncertainty quantification
and trustworthy AI. The main idea is to understand that a good model
shouldn't just say what it expects to happen, but also how much we can
trust that prediction.

Data
You'll be working with the METABRIC dataset, a large dataset of breast
cancer patients containing: time to clinical event and censoring
information; age; tumour size and grade; lymph nodes; ER/PR/HER2 status;
other clinical characteristics; molecular information, including PAM50
molecular subtype.
You'll compare two types of information:
1. Clinical model – clinical predictors only.
2. Clinical + molecular model – the same clinical predictors + PAM50.
This way you'll investigate whether the molecular information adds extra
value to the prognosis.
Simona did a thesis on this dataset with some of these models (from the
Statistical Modelling course), so you can consult her about some details.

Main research questions
1. How well can clinical characteristics predict survival or disease
   recurrence?
2. Does adding the PAM50 molecular subtype improve the prognosis?
3. How do a classical statistical model, a Bayesian model, and a machine
   learning model differ?
4. What uncertainty about the individual prognosis does the Bayesian
   model give?
5. Can conformal prediction create prediction intervals that actually
   achieve, for example, 80%, 90%, or 95% empirical coverage?
6. Is the model with the best predictive performance also the model with
   the most useful and best-calibrated uncertainty?

Methods
You'll use three main models:

1. Cox proportional hazards model
   A classical biostatistics model, for which you'll need the following
   concepts: survival function; hazard function; right censoring; hazard
   ratio; proportional hazards assumption.

2. Bayesian Weibull / AFT model
   A Bayesian model for survival analysis — prior distribution;
   likelihood with censoring; posterior distribution; posterior
   predictive distribution; credible and predictive intervals; MCMC
   diagnostics; posterior predictive checks.

3. Random Survival Forest
   A machine-learning model that allows for non-linear relationships and
   will serve as the ML comparison. Detailed hyperparameter optimization
   isn't needed — that's not the goal.

4. Conformal prediction – the main new topic, which I hope you haven't
   worked with before.
   First you need to understand how the split for conformal prediction
   works:
   training set -> calibration set -> test set.
   Then you'll use an existing conformal method adapted for survival data
   with censoring. You need to clearly understand the difference between:
   Bayesian uncertainty — tells us what the uncertainty is according to
   the chosen probabilistic model and the prior distribution.
   Conformal prediction — checks whether the prediction intervals
   actually achieve their declared coverage on new data.
   This is one of the central ideas of the thesis.

Evaluation
Predictive performance — you'll use: C-index; Brier score; Integrated
Brier Score; time-dependent AUC as needed.

Bayesian uncertainty — you'll analyse: posterior distributions;
individual posterior predictions; predictive intervals; prior
sensitivity; posterior predictive checks.

Conformal uncertainty — for nominal coverage levels of 80%, 90%, and 95%,
you'll compute: empirical coverage; prediction interval width.

It's important to understand that a wide interval can have good coverage
but still be practically uninformative.
That's why you also need to analyse: validity — whether the interval
covers as much as it promises; efficiency — how narrow and informative
it is.

What I expect will be new and/or useful for you for your master's studies
Biostatistics: survival analysis; right censoring; Kaplan-Meier; Cox
model; hazard ratio; AFT models; prognostic modelling.
Bioinformatics: gene expression; molecular breast-cancer subtypes; PAM50;
batch effects; feature-selection leakage.
Trustworthy ML: calibration; conformal prediction; exchangeability;
marginal coverage; interval efficiency; discrimination vs. calibration;
distribution shift.

You also need to understand the difference between prediction and
causality.
If a variable is a good predictor, that doesn't mean it has a causal
effect on the outcome. That's why this project is about prognostic
modelling, not causality.

Experimental design: You need to split METABRIC into training,
calibration, and test sets for fitting the models / conformal
prediction / final evaluation, respectively.
On the same splits, you'll compare the models: clinical vs. clinical +
PAM50, to see whether the molecular information improves: prediction,
calibration, and uncertainty.

In the end, you should be able to clearly answer three different
questions:
1. What is most likely to happen to the patient? (Statistical and ML
   prognostic modelling)
2. How uncertain is the prognosis according to the probabilistic model
   itself? (Bayesian predictive uncertainty)
3. Does the prediction interval actually achieve the confidence it
   declares? (Conformal prediction)

Proposed thesis structure
Introduction and motivation
Survival analysis and censoring
Bayesian survival model
Conformal prediction and trustworthy ML
METABRIC and PAM50
Methodology and experimental design
Results
Discussion
Conclusion and future work

This is enough for now — if you'd like to continue further, there are
more directions for future work: external validation on another dataset;
TCGA-BRCA; distribution shift; recalibration; weighted conformal
prediction; transfer learning; using PAM50 gene-expression features
directly.

Angelopoulos & Bates – A Gentle Introduction to Conformal Prediction
Candès, Lei & Ren – Conformalized Survival Analysis
Qin et al. – Conformal Predictive Intervals in Survival Analysis
PyMC tutorials for Bayesian Weibull/AFT survival analysis
scikit-survival tutorials for Cox, Random Survival Forest, C-index, and
Brier score
the original METABRIC paper for the clinical and molecular context of
the data

## Instructions

Analyze the deliverable above, then design and create the full .claude/ structure. 
For each file you create, explain WHY it belongs in that layer 
(CLAUDE.md vs rule vs skill vs hook vs subagent) before writing it.

Work through each layer in order:

### 1. CLAUDE.md (root)
Create a lean CLAUDE.md under 100 lines covering only:
- Project summary (1-2 lines)
- Build, test, lint, and run commands
- Directory structure and what lives where
- Hard conventions Claude cannot infer from code
- Pointers to rules/skills for detail (not the detail itself)
Do NOT include anything Claude can figure out on its own.

### 2. .claude/rules/
For each distinct concern in this project (e.g. API design, database, 
frontend components, testing, security), create a separate rule file.
- Use path frontmatter to scope rules to the files they apply to
- Each file should cover exactly one topic
- Keep each file under 30 lines

### 3. .claude/skills/
Identify every repeatable multi-step workflow this project needs 
(e.g. "add a new API endpoint", "cut a release", "run a code review").
Create a SKILL.md for each one with:
- Frontmatter: name, description (written so Claude knows when to auto-invoke it)
- Step-by-step instructions
- References to any docs/ files it should read during execution

### 4. .claude/agents/
Identify any review or specialist tasks that need an isolated fresh context
(e.g. security reviewer, code quality reviewer, test coverage checker).
Create agent files for each one.

### 5. settings.json hooks
Identify what MUST happen deterministically — not "Claude should do X" but 
"X must happen regardless of what Claude decides":
- PostToolUse: formatting, linting after file edits
- Stop: test gates, quality checks before declaring done
- PreToolUse: any dangerous commands that need validation
Write the hooks config for settings.json.

### 6. docs/ reference files
Identify what detailed reference material should live in docs/ 
(too long for CLAUDE.md, only needed for specific tasks).
Create stubs for each file and note which skill or rule should reference them.

### Output format
For each file, show:
CREATING: <path>
REASON: <why this layer, not another>
<file content>

Start with the plan. Wait for my approval before creating files.