## Multi-Metric Computational Prioritization of Candidate Ligands for Interleukin-8 (IL-8)

## Project overview

This project explores the computational discovery and prioritization of candidate ligands for interleukin-8 (IL-8). It brings together selected docking, predicted affinity, protein–ligand interaction, and binding-site metrics to compare candidates using a consensus-ranking workflow.

The Python analysis normalizes selected metrics, combines them into group and weighted consensus scores, compares rankings, and estimates how often compounds appear in the top five when the weights are varied. The purpose is to help prioritize candidates for further investigation.

**These results are computational predictions, not evidence that a compound binds IL-8, inhibits its activity, or is a validated drug. Experimental validation is required to assess those claims.**

This public repository contains code and documentation only. Input data, images, rankings, plots, and other thesis materials are intentionally kept private and are not included.

## What this repository contains

- `scripts/model.py` — primary analysis script for validating the input table, normalizing selected metrics, calculating group and consensus scores, assessing ranking stability, and generating local outputs.
- `scripts/make_panels.py` — optional script for making compound panels from local ranking, cluster data, and images.
- `requirements.txt` — Python package dependencies.
- `.gitignore` — excludes private input data, images, generated results, and common local files.

The numbered scripts `scripts/step1.py` through `scripts/step10.py` are earlier exploratory scripts. They refer to `data/test_data.csv`, which is not included. The documented workflow uses `scripts/model.py`.

## Workflow

The current documented workflow focuses on computational comparison and prioritization:

1. Provide the project input data locally.
2. Normalize selected metrics so that higher normalized values represent more favorable values.
3. Calculate scores for metric groups.
4. Calculate a weighted consensus score and rank candidates.
5. Compare weighted ranking with mean group ranking and leave-one-group-out rankings.
6. Estimate top-five ranking stability under randomized weights.
7. Generate tables and plots locally for review.

This repository does not include or claim to reproduce any earlier candidate-generation or experimental drug-discovery stages that are not documented in its code.

## Ranking approach

The primary model uses the following groups and weights:

| Group | Metrics | Weight |
|---|---|---:|
| Docking | `vina_kcal`, `cbdock_score` | 0.25 |
| HADDOCK | `haddock_score` | 0.25 |
| Affinity | `prodigy_dg_kcal` | 0.15 |
| Interactions | `plip_saltbridge_no`, `plip_hydrogen_no`, `plip_hydrophobic_no` | 0.15 |
| Site | `vina_contacts_4A`, `vina_min_dist_pocket` | 0.20 |

Each metric is scaled between 0 and 1, with 1 representing the most favorable value among the compounds in the input table. Lower values are treated as favorable for docking scores, HADDOCK score, predicted binding free energy, and minimum pocket distance. Higher values are treated as favorable for interaction counts and Vina contacts within 4 Å.

The model averages metrics within each group and calculates a weighted consensus score. If a compound is missing a group score, the available group weights are rescaled for that compound. Ranking stability is estimated using 1,000 randomized weight sets with a fixed random seed.

The ranking depends on the selected input compounds, metrics, and weights. Adding or removing compounds can change normalized scores and rankings. These scores are screening heuristics, not probabilities of binding or efficacy.

## Requirements

- Python 3
- Packages listed in `requirements.txt`

Install dependencies from the repository root:

```bash
python -m pip install -r requirements.txt
```

## Running the analysis

The input files are private and are not included in this repository. To run the model, place authorized copies of the required files at these paths on your local machine:

- `data/main_data.csv`
- `data/haddock_clusters.csv`

Then run:

```bash
python scripts/model.py
```

The script generates files in the local `results/` directory, including:

- `ranking.csv`
- `bar_chart.png`
- `heatmap.png`
- `haddock_vs_irmsd.png` when the cluster input contains the required columns

These outputs are ignored by Git and should not be committed to this public repository.

To generate compound panels, first run the model and make the required compound images available locally in `images/`. Then run:

```bash
python scripts/make_panels.py
```

The panel script expects image names based on each compound ID with `_haddock`, `_pymol`, and `_2d` suffixes. It saves generated panels under `results/panels/`. The private images and generated panels are not included in this repository.

Run the commands from the repository root because the scripts use paths relative to that directory.

## Private files and publication

The repository is configured as code-only. The following directories are excluded from Git:

- `data/` — input datasets
- `images/` — molecular and interaction images
- `results/` — rankings, plots, and panels

Do not remove these exclusions unless you have checked that the files may be shared and that publishing them is consistent with your thesis requirements. Keep sensitive or unpublished research results private.

## Data provenance and citations

The data and images used by the local workflow are not included here. Before sharing any research materials publicly, document their sources, how they were generated or processed, relevant software and database versions, and applicable permissions.

Cite the relevant target biology, computational tools, methods, and datasets in your thesis or in repository documentation when their details and sharing permissions are confirmed. No data sources or reuse permissions are asserted by this code-only repository.

## Limitations

- Docking and predicted affinity scores are estimates and may not correspond to experimental binding.
- Interaction counts and calculated molecular properties depend on the input structures and analysis methods.
- Consensus rankings depend on the selected metrics, normalization range, and weights.
- A high computational rank does not establish safety, efficacy, selectivity, or drug candidacy.
- Experimental validation is needed to evaluate biological activity.
