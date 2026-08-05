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

Datasets used:
   - Bulk RNA-Seq:
     - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164073 
   - Single-cell RNA-Seq:
     - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109816
   - Spatial transcriptomics:
     - https://data.mendeley.com/datasets/mbvhhf8m62/2
   - Transcriptomics:
     - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229677

 
Contributions are welcome ❤️



This project is licensed under the [MIT License](./LICENSE).
