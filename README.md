# De novo Assembly for Microbial AA Profiling

A computational workflow to infer the amino acid composition of active microbial communities directly from metatranscriptomic sequencing data using de novo transcriptome assembly.

This pipeline assembles microbial transcripts, predicts protein coding sequences, estimates transcript abundance, and calculates amino acid composition across the inferred proteome. It was developed to characterize the de novo amino acid profile of rumen microbial protein, but it can be adapted to other microbial ecosystems.

---

## Overview

See the wiki [here](https://github.com/GRMC-Cornell/de-novo-transcriptome-assembly-for-microbial-AA-profile/wiki#welcome-to-the-de-novo-transcriptome-assembly-for-microbial-aa-profile-wiki)!

Traditional estimates of microbial protein composition rely on cultured organisms or experimentally determined reference values. This workflow instead reconstructs the expressed proteome directly from RNA sequencing data, allowing amino acid composition to be estimated from the transcriptionally active microbial community.

The pipeline performs the following steps:

1. Assemble transcripts from metatranscriptomic reads.

2. Predict protein coding sequences.

3. Quantify transcript abundance.

4. Weight protein sequences according to transcript abundance.

5. Calculate amino acid composition for the inferred microbial proteome.

6. Generate amino acid profiles for downstream statistical analyses.

---

## Repository Contents

| File | Description |

| `Script.py` | Main pipeline execution script |

| `Funcs.py` | Helper functions used throughout the workflow |

| `Counts.py` | Amino acid counting and abundance calculations |

| `README.md` | Documentation |

| `LICENSE` | MIT License |

---

## Requirements

Python 3.9 or newer

Python packages:

- pandas

- numpy

- Biopython

Install dependencies with:

```bash

pip install pandas numpy biopython

```

---

## Input

The workflow expects predicted protein sequences and transcript abundance estimates generated from a de novo transcriptome assembly.

Upstream tools include:

- Trinity

- TransDecoder

- Salmon

Input files should include:

- Predicted protein FASTA file (.faa)

- Transcript abundance table (.csv)

---

## Output

The pipeline generates:

- Amino acid counts

- Relative amino acid abundances

- Normalized amino acid composition

- Summary tables for downstream statistical analyses

---

## Workflow

```text

Raw RNA-seq reads

        │

        ▼

De novo transcriptome assembly

        │

        ▼

Protein prediction

        │

        ▼

Transcript abundance estimation

        │

        ▼

Weighted protein abundance

        │

        ▼

Amino acid composition

        │

        ▼

Statistical analyses

```

---

## Citation

If you use this repository, please cite this repo until the manuscript is published.

> *Citation will be added upon publication.*
