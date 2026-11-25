# De novo Assembly for Microbial AA Profiling

A computational workflow to infer the amino acid composition of active microbial communities directly from metatranscriptomic sequencing data using de novo transcriptome assembly.

This pipeline assembles microbial transcripts, predicts protein-coding sequences, estimates transcript abundance as transcripts per million (TPM), and calculates amino acid composition by weighting predicted protein AA residue counts by their parent transcript TPM. It was developed to characterize the de novo amino acid profile of rumen microbial protein, but it can be adapted to other microbial ecosystems.

---

## Overview

See the wiki [here](https://github.com/GRMC-Cornell/de-novo-transcriptome-assembly-for-microbial-AA-profile/wiki#welcome-to-the-de-novo-transcriptome-assembly-for-microbial-aa-profile-wiki)!

Traditional estimates of microbial protein composition rely on cultured organisms or experimentally determined reference values. This workflow predicts protein sequences from assembled microbial transcripts and estimates their amino acid composition using transcript TPM as an abundance weight.

The pipeline performs the following steps:

1. Assemble transcripts from metatranscriptomic reads.

2. Predict protein coding sequences.

3. Quantify transcript abundance separately for each sample using Salmon TPM values.

4. Multiply each predicted protein’s AA residue counts by its parent transcript TPM, incorporating protein length through the residue counts.

5. Sum the TPM-weighted AA residue counts to calculate residue percentages, and apply AA-specific residue molecular weights to calculate anhydrous AA mass percentages.

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

- Predicted protein sequences in FASTA format (.pep or .faa).

- A transcript-level TPM matrix, with transcript IDs in rows and samples in columns.

- transcript-to-taxonomic-group mapping table for domain-specific AA profiles.

---

## Output

The pipeline generates:

- TPM-weighted AA residue totals.

- Amino acid residue percentages.

- Anhydrous AA mass percentages calculated using AA-specific residue molecular weights.

- Overall and, when taxonomy is supplied, domain-specific AA profiles.

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

Transcript abundance estimation (Salmon TPM)

        │

        ▼

Weighted protein abundance (AA residue counts × parent transcript TPM)

        │

        ▼

Amino acid composition (AA residue percentages and anhydrous AA mass percentages)

        │

        ▼

Statistical analyses

```

---

## Citation

If you use this repository, please cite this repo until the manuscript is published.

> *Arshad, M. A., M. A. Majersik, K. Nichols, B. C. Weimer, R. Profeta, M. I. Marcondes, R. S. Bisinotto, J. E. P. Santos, E. S. Ribeiro, F. S. Lima, and H. F. Monteiro. 2026. A molecular framework for rumen microbial protein: de novo characterization of the amino acid profile. Manuscript under revision (R1), Journal of Dairy Science.*
