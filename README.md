# Bioinformatics Code Hub
A community collection of bioinformatic analyses for anyone to learn and adapt.

Need help with an analysis? You can open an issue:
- Request a specific analysis
- Request an analysis expanded for your needs
- Request help with downloading/processing a specific dataset

Good at bioinformatics? Consider contributing your code!

## Planned Analyses:
### Genomics
- Copy-number variation
- Variant calling
  - Small variants (GATK)
  - Structural variants (Manta)
  - Haplotype phasing (WhatsHap)
- Variant annotation (VEP, ANNOVAR)
- GWAS
  - QC & imputation (PLINK, Eagle)
  - Association testing
- Population genetics
  - PCA / ADMIXTURE
  - Selection scans (iHS, Fst)

### Transcriptomics
- Bulk RNA-seq
  - QC
  - Alignment
  - Differential expression:
    - [DESeq2 - R](transcriptomics/bulk/DESeq2/r.ipynb)
    - [PyDESeq2](transcriptomics/bulk/DESeq2/py.ipynb)
    - edgeR
  - Counts generation:
    - [Kallisto](transcriptomics/kallisto/bash.ipynb)
  - Isoform analysis (IsoformSwitchAnalyzeR)
  - Co-expression networks (WGCNA)
  - GSEA / pathway analysis
- Single-cell RNA-seq
  - Preprocessing (Cell Ranger / Alevin-fry)
  - Normalisation & batch correction (SCTransform / Harmony)
  - Doublet detection (Scrublet)
  - Clustering & visualisation
  - Cell-type annotation:
    - [CellTypist](transcriptomics/single-cell/cell%20annotations/cell%20typist/py.ipynb)
    - SingleR
    - scmap
  - Differential expression
  - [TF analysis](transcriptomics/single-cell/transcription%20factor%20analysis/r.ipynb)
  - [Pseudotime / trajectory analysis](transcriptomics/single-cell/pseudotime%20analysis/r.ipynb)
  - [Ligand-receptor interactions](transcriptomics/single-cell/ligand-receptor%20analysis/r.ipynb)
- Spatial transcriptomics
  - [Visualize clusters and genes on spatial coordinates](transcriptomics/spatial-transcriptomics/basic%20analysis/r.ipynb)
  - Spot deconvolution (RCTD)
  - Image-based QC
- Long-read transcriptomics (Iso-seq, Nanopore)

### Epigenomics
- DNA methylation (Bismark)
- ATAC-seq / chromatin accessibility
- ChIP-seq
- Hi-C / 3D genome
- Single-cell epigenomics (snmC-seq, scATAC)

### Proteomics
- LC-MS preprocessing (MaxQuant)
- Differential protein abundance
- PTM site localisation
- Spectral library generation
- Protein–protein interaction networks (STRING / Cytoscape)

### Metabolomics
- LC-MS untargeted workflows (MS-Dial / XCMS)
- NMR spectroscopy
- Pathway mapping (Mummichog)

### Metagenomics
- 16S/18S rRNA amplicon pipelines (QIIME 2)
- Shotgun metagenome assembly & binning (MetaBAT)
- Functional profiling (HUMAnN)

### Structural Bioinformatics & Molecular Simulation
- Homology modelling (MODELLER)
- Docking (AutoDock Vina)
- Molecular dynamics
  - MDVerse
  - GROMACS
  - OpenMM

### Imaging & Computer Vision
- Bright-field / fluorescence segmentation (Cellpose, Stardist)
- Cell tracking (DeepCell)
- Medical imaging (MRI, CT)

### Machine Learning & Statistics
- ML for omics
- Deep learning (keras-tensor) for sequence data
- AutoML notebooks (auto-sklearn)

### Multi-Omics Integration
- MOFA / DIABLO integrative analysis
- Network-based integration

### Utilities & QC
- File-format conversions (BAM ⇆ CRAM)
- Reference genome downloads (NCBI, Ensembl)
- Workflow management (Snakemake, Nextflow)

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
