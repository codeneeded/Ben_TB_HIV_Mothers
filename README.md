# scRNA-seq + TCR of HIV-Positive Mothers Stratified by Latent TB Infection (IGRA)

Single-cell transcriptomic and TCR repertoire analysis of peripheral immune cells from **HIV-positive mothers**, stratified by **latent TB infection (LTBI)** status as determined by the **Interferon-γ Release Assay (IGRA⁺ vs IGRA⁻)**. The analysis profiles how *Mycobacterium tuberculosis* co-infection reshapes the immune landscape of women living with HIV — including cell-type composition, per-cell-type transcriptional programs, module-level activation states, HIV reservoir–favoring signatures, and clonal architecture of the αβ T-cell repertoire.

The data are generated on the **Parse Biosciences** combinatorial-barcoding scRNA-seq platform (10 samples) with paired αβ TCR information processed through scRepertoire.

---

## Cohort and platform

- **10 mothers** (sample IDs `1273088`, `1273505`, `6056780`, `6060100`, `6060150`, `6060180`, `6060280`, `6065400`, `6072090`, `6072210`), all HIV-positive, partitioned into IGRA⁺ and IGRA⁻ groups (mapping in [`Annotation/IGRA_annotation.csv`](Annotation/IGRA_annotation.csv)).
- **Platform:** Parse Biosciences GEX (combined sublibrary outputs read with `ReadParseBio()`).
- **Modalities:** Gene expression + αβ TCR.

---

## Repository structure

```
Ben_TB_HIV_Mothers/
├── Scripts/                                  # End-to-end numbered R pipeline
│   ├── 1A. Create_Seurat_Parse_raw.r         # Read raw Parse outputs → merged Seurat
│   ├── 1B. Create_Seurat_Parse_Parse_Filtered.r
│   ├── 2.Basic_QC.R                          # nFeature / nCount / mito / ribo / heme / platelet / complexity
│   ├── 3. Cell_Cycle_Effect-Doublet_Removal.R
│   ├── 4.Integration.R                       # CCA + Harmony + MNN, side-by-side
│   ├── 5.Annotation.R                        # Manual + Azimuth-assisted cluster annotation
│   ├── 6. Labels+DGE.R                       # Per-cluster IGRA⁺ vs IGRA⁻ DGE + pseudobulk
│   ├── 7.Module_Scoring.R                    # Themed module scores + composite indices
│   ├── 8.Pathway_Analysis_EnrichR.R          # EnrichR over per-cell-type DEGs
│   └── 9.TCR.R                               # scRepertoire workflow
│
├── QC_Plots/
│   ├── RAW_QC/                               # Pre-filter exploration: cells per sample/IGRA, scatter, histograms
│   ├── PARSEF_QC/                            # Post-Parse-filter equivalents
│   ├── Cell_Cycle/                           # Cycle scoring, ridge plots, PCA by phase, Azimuth UMAP
│   ├── Doublets/                             # Per-sample doublet diagnostic plots (10 mothers)
│   └── FINAL_QC/                             # Post-QC cells per sample / IGRA, feature distributions
│
├── Integration_Plots/
│   ├── CCA_Clustering_Grid.png
│   ├── Harmony_Clustering_Grid.png
│   └── MNN_Clustering_Grid.png               # Three integration methods compared
│
├── Annotation/
│   ├── Annotated_UMAP.png
│   ├── ClusterDistribution_Absolute.png
│   ├── ClusterDistribution_Proportion.png
│   ├── ClusterDistribution_IGRA_PosVsNeg.png # Composition shift by IGRA status
│   ├── ClusterFrequency_Sample_IGRA_PosVsNeg.png
│   ├── IGRA_annotation.csv                   # Sample → IGRA status mapping
│   ├── FeaturePlots/                         # ~140 marker FeaturePlots
│   ├── VlnPlots/                             # Pre-annotation marker violins
│   └── Post-Annotation/
│       ├── VlnPlots/                         # Marker violins on annotated clusters
│       └── VlnPlots_Extra/                   # Additional markers (AP-1, KLFs, exhaustion, cycling)
│
├── DGE/                                      # IGRA⁺ vs IGRA⁻ within each annotated cell type
│   ├── DGE_<celltype>_IGRApos_vs_neg.csv
│   ├── Pseudobulk_DGE_IGRApos_vs_neg.csv
│   └── Volcano_Plots/                        # Per-cell-type volcanoes + custom-label variants
│
├── Pathway_Analysis_EnrichR/IGRAPosvsNeg/
│   ├── CSVs/                                 # Per-cell-type EnrichR XLSX exports
│   └── Plots/                                # EnrichR pathway / TF plots
│
├── Module_Scoring/                           # Themed AddModuleScore() + composite indices
│   ├── Autophagy/                            # AMPK, mTORC, Initiation, Conjugation, Cargo, Fusion,
│   │                                         # Lysosome, Mitophagy, Xenophagy, cGAS-STING, CLEAR-TFs
│   ├── ISG/                                  # Type I & Type II interferon signatures
│   ├── IL10_STAT3/                           # IL10R-JAK1-TYK2-STAT3 axis, M2-like, Treg overlap
│   ├── TGFb/                                 # SMAD2/3-SMAD4, ligand/receptor, noncanonical MAPK/PI3K, transcriptional targets
│   ├── Myeloid/                              # Inflammatory, Regulatory, DC maturation
│   ├── NK/                                   # Activation, Effector, Inhibition
│   ├── Bcell/                                # GC-like, Atypical (ABC), Plasmablast
│   ├── CD8/                                  # Effector, Memory, Exhaustion
│   ├── CD4_Reservoir/, Reservoir_CD4/        # Reservoir Memory, Tfh Core, Tph Discriminator
│   ├── Reservoir_Favoring/                   # Apoptosis, Cell-cycle E2F/G2M, Chromatin repression,
│   │                                         # Glycolysis, OxPhos, MYC targets, NF-κB, STAT reactivation
│   ├── Composites/                           # Per-cell composite indices (see below)
│   │   ├── MS_RPI.png                        # Reactivation Propensity Index
│   │   ├── MS_MAI.png                        # mTOR–MAPK Autophagy Index
│   │   ├── MS_IL10_composite.png, MS_TGFb_composite.png
│   │   ├── IL10_vs_TGFb_composites_sample_scatter*.png
│   │   └── Sample_Level/                     # RPI/MAI per-sample (overall + per-cluster) CSVs + plots
│   ├── Module_Summary_*.png                  # One summary per module family
│   ├── Module_Gene_Coverage.tsv              # How many panel genes were detected in the data
│   └── Module_Stats_byCluster.csv            # Per-cluster module statistics
│
└── TCR/                                      # scRepertoire outputs
    ├── CD3_Composition/                      # TRA/TRB V- and J-gene heatmaps, kmer composition,
    │                                         # %AA usage, positional entropy
    ├── Clonal_Diversity/                     # Diversity (TARA), Morisita / raw clonal overlap
    ├── Clonal_Visualizations/                # Abundance, homeostasis, length, strict-clone counts,
    │                                         # unique-clone IGRA⁺ vs IGRA⁻
    └── Seurat_Plots/
        ├── Clonal_Overlay/                   # Clonal sizes overlaid on UMAP, by sample / by condition
        ├── Clones_per_Cluster/               # Clonal occupancy per annotated cluster
        └── Hyperexpansion/                   # Hyperexpanded clones overall and split by IGRA status
```

---

## Pipeline overview

The full analysis lives in `Scripts/`, run in numerical order.

### 1 – Object construction (`1A`, `1B`)
Per-sample Parse Bio outputs (`DGE_unfiltered/` and `DGE_filtered/`) are loaded with `ReadParseBio()` and merged into a single Seurat object. Sample identity is propagated from the Parse folder name into `orig.ident` / `sample`.

### 2 – Basic QC (`2.Basic_QC.R`)
Cell-level QC over `nFeature_RNA`, `nCount_RNA`, percent mitochondrial / ribosomal / hemoglobin / platelet, and a complexity (`log10GenesPerUMI`) score. Both pre- and post-Parse-filter views are kept (`QC_Plots/RAW_QC/`, `QC_Plots/PARSEF_QC/`) so the effect of Parse's own cell calling is visible.

### 3 – Cell cycle and doublets (`3. Cell_Cycle_Effect-Doublet_Removal.R`)
Seurat `CellCycleScoring` + per-sample doublet detection (one diagnostic plot per donor in `QC_Plots/Doublets/`), with a filtered, doublet-cleaned object saved for downstream steps.

### 4 – Integration (`4.Integration.R`)
The cleaned object is split by sample and integrated three different ways with `IntegrateLayers()`:
- **CCA** (`CCAIntegration`)
- **Harmony** (`HarmonyIntegration`)
- **MNN** (fast MNN via SeuratWrappers)

Side-by-side clustering grids in `Integration_Plots/` allow direct visual comparison so the choice of integration method is auditable.

### 5 – Annotation (`5.Annotation.R`)
Clusters are annotated using a combination of Azimuth reference mapping and manual marker review (`FeaturePlots/`, `VlnPlots/`). The resulting cell types span the major peripheral compartments:

> CD4⁺ naïve / TCM / Tfh / Treg / Activated · CD8⁺ TCM / EM / TEMRA · dnT · NK CD56-bright · NK Cytotoxic · CD14⁺ Monocytes · CD16⁺ NCM · cDC2 · DC-like APCs · moDC · moDC precursor · Naïve B · Memory B · Intermediate B · GC-like B

Composition shifts between IGRA⁺ and IGRA⁻ are summarized in `Annotation/ClusterDistribution_IGRA_PosVsNeg.png` and `ClusterFrequency_Sample_IGRA_PosVsNeg.png`.

### 6 – Differential expression (`6. Labels+DGE.R`)
For each annotated cell type, **IGRA⁺ vs IGRA⁻** is tested at single-cell resolution and exported as `DGE/DGE_<celltype>_IGRApos_vs_neg.csv`, with matched volcanoes in `DGE/Volcano_Plots/`. A **pseudobulk** DESeq2 analysis (`Pseudobulk_DGE_IGRApos_vs_neg.csv`) is run alongside as a more conservative donor-level test that controls for sample identity.

### 7 – Module scoring (`7.Module_Scoring.R`)
Curated gene panels are scored cell-by-cell with `AddModuleScore` and grouped into thematic families (Autophagy, ISG, IL10/STAT3, TGFβ, Myeloid, NK, B-cell, CD8, CD4 Reservoir, Reservoir-Favoring metabolism/transcription). Two **composite indices** are derived from these submodules and z-scored within cluster:

- **RPI — Reactivation Propensity Index.** A signed combination of HIV-reservoir-permissive features (Reservoir_Favoring submodules + glycolysis / OxPhos contribution) minus reservoir-suppressive features. Higher RPI = a transcriptional state more permissive to viral reactivation.
- **MAI — mTOR–MAPK Autophagy Index.** `z(pro-autophagy + noncanonical MAPK)` minus `z(mTORC)`. Higher MAI = more autophagy-favoring metabolic posture.

Both per-cell and per-sample (overall and per-cluster) summaries are written for RPI and MAI under `Module_Scoring/Composites/` and `Composites/Sample_Level/`, alongside an IL10 vs TGFβ composite scatter for examining the regulatory axis at the donor level.

### 8 – Pathway enrichment (`8.Pathway_Analysis_EnrichR.R`)
Per-cell-type DEG lists are passed to **enrichR**; XLSX exports per cell type live in `Pathway_Analysis_EnrichR/IGRAPosvsNeg/CSVs/`, with summary plots in the matching `Plots/` directory.

### 9 – TCR repertoire (`9.TCR.R`)
The αβ TCR analysis uses **scRepertoire**:
- **CD3 composition** — V- and J-gene usage heatmaps for TRA and TRB, k-mer composition, percent amino-acid composition, and positional entropy.
- **Clonal diversity / overlap** — TARA-style diversity, Morisita and raw overlap matrices comparing donors.
- **Clonal visualization** — abundance, homeostasis, length distributions, strict-definition clone counts, IGRA⁺ vs IGRA⁻ unique-clone comparisons.
- **Seurat overlays** — clonal sizes overlaid on the UMAP by sample and by condition; per-cluster occupancy; hyperexpansion analyses overall and split by IGRA status.

---

## Dependencies

**Core single-cell**
- [`Seurat`](https://satijalab.org/seurat/) v5 — data structure, integration, DGE, module scoring
- [`SeuratWrappers`](https://github.com/satijalab/seurat-wrappers), [`SeuratExtend`](https://github.com/huayc09/SeuratExtend), [`scCustomize`](https://samuel-marsh.github.io/scCustomize/)
- [`harmony`](https://github.com/immunogenomics/harmony) — batch correction
- [`Azimuth`](https://azimuth.hubmapconsortium.org/) — reference-based annotation
- [`SingleCellExperiment`](https://bioconductor.org/packages/SingleCellExperiment/), [`scater`](https://bioconductor.org/packages/scater/)

**Repertoire**
- [`scRepertoire`](https://www.borch.dev/uploads/screpertoire/) — TCR clonotype analysis

**Enrichment**
- [`enrichR`](https://CRAN.R-project.org/package=enrichR)

**Visualization & utilities**
- `ggplot2`, `patchwork`, `cowplot`, `ggrepel`, `gridExtra`
- [`ComplexHeatmap`](https://bioconductor.org/packages/ComplexHeatmap/), `circlize`
- `RColorBrewer`, `Polychrome`, `scales`
- `tidyverse`, `data.table`, `Matrix`, `hdf5r`, `readxl`, `Cairo`, `igraph`

---

## Key questions the analysis addresses

| Question | Where to look |
|---|---|
| Does LTBI (IGRA⁺) shift the peripheral immune composition in HIV⁺ mothers? | `Annotation/ClusterDistribution_IGRA_PosVsNeg.png`, `ClusterFrequency_Sample_IGRA_PosVsNeg.png` |
| Which cell types respond most transcriptionally to LTBI? | `DGE/Volcano_Plots/`, `DGE/Pseudobulk_DGE_IGRApos_vs_neg.csv` |
| What pathways and TFs drive those responses? | `Pathway_Analysis_EnrichR/IGRAPosvsNeg/` |
| Is the HIV reservoir–permissive program elevated in IGRA⁺ donors? | `Module_Scoring/Reservoir_Favoring/`, `Composites/MS_RPI.png`, `Composites/Sample_Level/RPI_*` |
| Does autophagic / metabolic posture differ by LTBI status? | `Module_Scoring/Autophagy/`, `Composites/MS_MAI.png`, `Composites/Sample_Level/MAI_*` |
| Are clonal expansion / hyperexpansion patterns different in IGRA⁺ vs IGRA⁻? | `TCR/Clonal_Visualizations/Unique_Clone_TRAB_IgraposvsNeg.png`, `TCR/Seurat_Plots/Hyperexpansion/Clonal_Hyperexpansion_by_IGRA_Status.png` |

---

## Reproducing the analysis

1. Clone the repo and open `Ben_TB_HIV_Mothers.Rproj` in RStudio.
2. Update the absolute paths near the top of each script (the scripts currently point at `/home/akshay-iyer/Documents/Ben_TB_HIV_Mothers/...`).
3. Place the per-sample Parse Bio output folders (`sample_<id>/DGE_unfiltered/` and `DGE_filtered/`) under the directory referenced in `1A`/`1B`.
4. Run scripts 1–9 in order. Intermediate Seurat objects are saved to `saved_R_data/` (gitignored) so later scripts pick up where earlier ones left off.
