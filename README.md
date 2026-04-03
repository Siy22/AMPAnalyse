# AMPAnalyse 🧬

**An in silico pipeline for antimicrobial peptide (AMP) design via Structure-Activity Relationship (SAR) analysis and intelligent variant generation.**

AMPAnalyse is a two-stage Python toolkit built for the rational design of optimized antimicrobial peptides. Starting from mutagenesis data, it extracts data-driven design principles and uses them to generate ranked peptide variants with improved safety and/or activity profiles.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline](#pipeline)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Stage 1 — Pattern Analyzer](#stage-1--pattern-analyzer)
  - [Stage 2 — Variant Generator](#stage-2--variant-generator)
- [Input Format](#input-format)
- [Output Files](#output-files)
- [Configuration](#configuration)
- [Example Workflow](#example-workflow)
- [Project Context](#project-context)
- [License](#license)

---

## Overview

Antimicrobial peptides (AMPs) are promising candidates against drug-resistant pathogens, but optimizing them — maximizing antimicrobial potency while minimizing hemolytic toxicity — is a challenging multidimensional problem. AMPAnalyse addresses this by:

1. Mining single-point mutagenesis data (e.g., from PYAMPA) to uncover which amino acid substitutions improve fitness.
2. Using the extracted principles to intelligently generate and rank peptide variants at user-specified positions.

---

## Pipeline

```
PYAMPA Mutagenesis CSV
        │
        ▼
┌─────────────────────────────┐
│  pyampa_pattern_analyzer.py │  ← Stage 1: SAR Extraction
└─────────────────────────────┘
        │
        ▼
  design_principles.json
        │
        ▼
┌──────────────────────────────────┐
│  peptide_variant_generator.py    │  ← Stage 2: Variant Generation
└──────────────────────────────────┘
        │
        ▼
  <sequence>_dual_variants.csv
```

---

## Features

### Stage 1 — PYAMPA Pattern Analyzer (`pyampa_pattern_analyzer.py`)
- Flexible CSV ingestion with automatic column name matching
- Data-driven thresholds (median-based for AMP/hemolysis, zero-baseline for fitness)
- Amino acid property enrichment (type, size, charge, hydrophobicity, aromaticity)
- Amino acid preference analysis ranked by activity, safety, and fitness
- Property-change analysis (type substitutions, hydrophobicity shifts, aromatic effects)
- Position-specific mutation profiling
- Extraction of 6 actionable SAR design principles
- JSON export of principles and amino acid rankings

### Stage 2 — Smart Peptide Variant Generator (`peptide_variant_generator.py`)
- Loads SAR principles from Stage 1 output
- User-defined target positions with range support (e.g., `4,5,7` or `1-3,10`)
- Three optimization strategies: **Safety**, **Balanced**, **Activity**
- Generates single, double, triple, and quadruple mutation combinations
- Property prediction for HC50 (hemolytic concentration) and MIC (minimum inhibitory concentration)
- Composite scoring with mutation penalty to prefer simpler variants
- Dual output: combined multi-mutation variants and unmasked single-mutation variants
- CSV export of ranked candidates

---

## Requirements

- Python 3.8+
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)

Install dependencies:

```bash
pip install pandas numpy
```

---

## Installation

```bash
git clone https://github.com/<your-username>/AMPAnalyse.git
cd AMPAnalyse
pip install pandas numpy
```

---

## Usage

### Stage 1 — Pattern Analyzer

Edit the configuration section at the top of `pyampa_pattern_analyzer.py`:

```python
PYAMPA_FILE    = "your_mutagenesis_data.csv"
OUTPUT_REPORT  = "pyampa_analysis_report.txt"
OUTPUT_PRINCIPLES = "design_principles_YOURSEQ.json"
```

Then run:

```bash
python pyampa_pattern_analyzer.py
```

The script will print a full SAR report to the terminal and save a JSON file of design principles.

---

### Stage 2 — Variant Generator

Edit the configuration at the top of `peptide_variant_generator.py` to point to your Stage 1 output:

```python
PRINCIPLES_FILE = "principles/design_principles_YOURSEQ.json"
OUTPUT_FILE     = "YOURSEQ_variants.csv"
MAX_MUTATIONS   = 2       # Maximum simultaneous mutations (1–4)
TOP_N_CANDIDATES = 10     # Number of top variants to report
```

Then run:

```bash
python peptide_variant_generator.py
```

The script is interactive and will prompt you for:

| Prompt | Description |
|---|---|
| Peptide sequence | Your parent AMP sequence (e.g., `KRRVSWWWMEERR`) |
| Baseline HC50 (μM) | Hemolytic concentration of the parent peptide |
| Baseline MIC (μg/mL) | Minimum inhibitory concentration of the parent peptide |
| Target positions | Positions to mutate (e.g., `4,5,7` or `1-3`) |
| Strategy | `1` Safety · `2` Balanced · `3` Activity |

---

## Input Format

Stage 1 expects a **CSV file** with the following columns (alternative names are accepted automatically):

| Standard Column | Accepted Aliases |
|---|---|
| `Position` | `Pos` |
| `Original_AA` | `Original`, `OriginalAA` |
| `Mutated_AA` | `Mutated`, `MutatedAA` |
| `AMP_Prob` | `AMP_Probab`, `AMP_Probability`, `amp_prob` |
| `Hemo_Prob` | `Hemolytic_P`, `Hemolytic_Prob`, `Hemolytic_Probability`, `hemo_prob` |
| `Fitness` | `Fitness_Score`, `fitness_score` |

Example CSV row:

```
Position,Original_AA,Mutated_AA,AMP_Prob,Hemo_Prob,Fitness
1,K,A,0.812,0.341,0.124
```

---

## Output Files

| File | Description |
|---|---|
| `design_principles_<SEQ>.json` | SAR principles, amino acid rankings, and analysis thresholds |
| `<SEQUENCE>_dual_variants.csv` | Ranked variant candidates (multi-mutation + single-mutation groups) |

The variants CSV includes: `rank`, `sequence`, `mutations`, `num_mutations`, `predicted_hc50`, `predicted_mic`, `hc50_improvement_%`, `mic_improvement_%`, `score`, `hemo_modifier`, `amp_modifier`, `group`.

---

## Configuration

### Analysis Thresholds (Stage 1)

| Parameter | Default | Description |
|---|---|---|
| `GOOD_AMP_THRESHOLD` | Median (data-driven) | Minimum AMP probability to be considered active |
| `GOOD_HEMO_THRESHOLD` | Median (data-driven) | Maximum hemolysis probability to be considered safe |
| `GOOD_FITNESS_THRESHOLD` | `0.0` | Zero-baseline: positive means better than wildtype |

### Generation Parameters (Stage 2)

| Parameter | Default | Description |
|---|---|---|
| `MAX_MUTATIONS` | `2` | Maximum simultaneous mutations per variant |
| `TOP_N_CANDIDATES` | `10` | Number of top-ranked variants to display and export |

---

## Example Workflow

```bash
# Stage 1: Analyze mutagenesis data for KRRVSWWWMEERR
python pyampa_pattern_analyzer.py
# → Outputs: design_principles_KRRVSWWWMEERR.json

# Stage 2: Generate optimized variants
python peptide_variant_generator.py

# Interactive prompts:
# Sequence:       KRRVSWWWMEERR
# HC50 (μM):      37.146
# MIC (μg/mL):    0.714
# Positions:      5,6,7       ← hemolytic hotspots identified by HemoPI2
# Strategy:       1 (Safety)

# → Outputs: KRRVSWWWMEERR_dual_variants.csv
```

After generating variants, the recommended next steps are:
1. Review the CSV and shortlist 2–3 candidates
2. Validate predictions with external tools (e.g., HemoPI2, ToxinPred, APD3)
3. Progress shortlisted candidates to experimental testing

---

## Project Context

AMPAnalyse was developed as part of a **Final Year Project** in an in silico AMP design pipeline. The tool is designed to bridge computational mutagenesis screening (via tools like PYAMPA) with rational variant design, providing interpretable SAR principles alongside quantitative variant rankings.

---

## License

This project is currently unlicensed. All rights reserved by the author. If you would like to use or build upon this work, please open an issue or contact the repository owner.
