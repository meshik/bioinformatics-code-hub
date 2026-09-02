# Bioinformatics Code Hub
A community collection of bioinformatic analyses for anyone to learn and adapt.

Need help with an analysis? You can open an issue:
- Request a specific analysis
- Request an analysis expanded for your needs
- Request help with downloading/processing a specific dataset

Good at bioinformatics? Consider contributing your code!

## Existing Analyses:

### Transcriptomics
- Bulk RNA-seq
  - Differential expression:
    - [DESeq2 - R](transcriptomics/bulk/DESeq2/r.ipynb)
    - [PyDESeq2](transcriptomics/bulk/DESeq2/py.ipynb)
  - Counts generation:
    - [Kallisto](transcriptomics/kallisto/bash.ipynb)
- Single-cell RNA-seq
  - Cell-type annotation:
    - [CellTypist](transcriptomics/single-cell/cell%20annotations/cell%20typist/py.ipynb)
  - Differential expression
  - [TF analysis](transcriptomics/single-cell/transcription%20factor%20analysis/r.ipynb)
  - [Pseudotime / trajectory analysis](transcriptomics/single-cell/pseudotime%20analysis/r.ipynb)
  - [Ligand-receptor interactions](transcriptomics/single-cell/ligand-receptor%20analysis/r.ipynb)
- Spatial transcriptomics
  - [Visualize clusters and genes on spatial coordinates](transcriptomics/spatial-transcriptomics/basic%20analysis/r.ipynb)

### Epigenomics
- DNA methylation
  - Illumina array ([450k](epigenomics/dna_methylation/r.ipynb), 850k, EPIC)

### Proteomics
- Differential Abundance:
  - [proDA](proteomics/differential_abundance/r.ipynb)

## Datasets used

<!-- DATASETS:START -->
<!-- Generated from datasets.yml by scripts/render_datasets.py. Do not edit manually. -->

| Data type | Dataset | Used in |
| --- | --- | --- |
| Bulk RNA-seq | [GSE164073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164073) — SARS-CoV-2 infection of human cornea, limbus, and sclera | [DESeq2](transcriptomics/bulk/DESeq2/r.ipynb), [PyDESeq2](transcriptomics/bulk/DESeq2/py.ipynb) |
| Bulk RNA-seq | [GSE229677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229677) — Mouse endocardium | [Kallisto data source](transcriptomics/kallisto/data/README.md) |
| Single-cell RNA-seq | [GSE109816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109816) — Human adult normal heart | [Single-cell analyses](transcriptomics/single-cell/data/README.md) |
| Single-cell RNA-seq | [E-MTAB-13632](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13632) — Human organogenesis model | [FASTQ to counts with STARsolo](transcriptomics/single-cell/fastq-to-counts/10x-gene-expression-starsolo/bash.ipynb) |
| Spatial transcriptomics | [mbvhhf8m62 v2](https://data.mendeley.com/datasets/mbvhhf8m62/2) — Developmental heart - filtered and unfiltered count matrices and meta tables (Asp, M. (2021), Mendeley Data, V2, doi: 10.17632/mbvhhf8m62.2) — License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [Spatial transcriptomics analysis](transcriptomics/spatial-transcriptomics/basic%20analysis/r.ipynb) |
| DNA methylation | [GSE47915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47915) — Benign and tumor prostate tissue on Illumina 450K arrays | [Illumina 450K analysis](epigenomics/dna_methylation/r.ipynb) |
| Proteomics | [UbiLength](https://doi.org/10.1016/j.molcel.2017.01.004) — Interactors of linear ubiquitin baits of different lengths (Zhang et al., 2017; distributed with DEP 1.32.0) | [Differential abundance with proDA](proteomics/differential_abundance/r.ipynb) |
<!-- DATASETS:END -->

### Third-party data

The MIT License applies only to original code and documentation in this repository. Third-party datasets are downloaded from their original sources, are not redistributed here, and remain subject to their original licenses and repository terms. The linked analyses process or transform source data as educational reanalyses.

 
Contributions are welcome ❤️



This project is licensed under the [MIT License](./LICENSE).
