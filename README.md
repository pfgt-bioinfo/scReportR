# scReportR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
<!-- badges: end -->

`scReportR` is a personal R package for standardised single-cell RNA-seq
analysis. It provides a modular workflow — QC, preprocessing, normalisation,
integration, clustering, and annotation — driven by a single `params.yml`
configuration file and rendered as a self-contained Quarto HTML report.

## Key features

- **One config file** — all parameters (QC thresholds, normalisation method,
  integration covariates, marker lists, colours) live in `params.yml`
- **qs2 cache** — each expensive step is cached; re-running picks up where
  you left off
- **Mixed library types** — designed for projects combining OCM, 3'v4, and
  sorted populations from the same patient
- **Clean separation** — computation functions in the package, display logic
  in `analysis.qmd`; plots never auto-insert, the notebook decides layout
- **Flexible metadata** — optional columns (tissue, population) are handled
  gracefully when absent

## Installation

```r
# Install from GitHub
devtools::install_github("pfgt-bioinfo/scReportR")
```

**Dependencies** — the following Bioconductor packages must be installed
separately:

```r
BiocManager::install(c(
  "SingleCellExperiment", "scDblFinder",
  "celda", "celldex", "SingleR", "BiocParallel"
))
```

Optional:

```r
# Pan-Human Azimuth (cloud annotation)
devtools::install_github("satijalab/AzimuthAPI")

# Cluster stability visualisation
install.packages("clustree")
```

## Quick start

```r
library(scReportR)

# 1. Initialise a new project
init_sc_report(
  project_name = "Patient01",
  path         = "~/analyses",
  author       = "Your Name",
  samples      = c("CD45neg_T1", "CD45pos_T1", "immune_T1")
)
# → Opens Patient01.Rproj in RStudio

# 2. Edit params.yml  (sample paths, QC thresholds, colours, ...)

# 3. Render analysis.qmd section by section, or all at once:
quarto::quarto_render("analysis.qmd")
```

## Workflow

```
params.yml
    │
    ▼
sc_config()          # load parameters, set up directories & palettes
    │
    ├── sc_load()            # 10x data loading + QC metrics
    ├── sc_decontX()         # ambient RNA decontamination
    ├── sc_doublets()        # doublet detection (scDblFinder)
    ├── sc_filter()          # per-sample QC filtering
    ├── sc_merge()           # merge into Seurat v5 object
    │
    ├── sc_normalize()       # SCTransform v2 or LogNormalize + PCA
    ├── sc_integrate()       # Harmony / RPCA / CCA + UMAP
    │
    ├── sc_cluster()         # Leiden clustering (multi-resolution)
    ├── sc_find_markers()    # FindAllMarkers (SCT-aware)
    ├── sc_module_scores()   # per-cluster module scores
    │
    ├── sc_singler()         # SingleR annotation (celldex references)
    ├── sc_panhuman_azimuth()# Pan-Human Azimuth (cloud API)
    └── sc_annotate_manual() # manual cluster-to-celltype mapping
```

Each step saves its result to `cache/` as a `.qs2` file. Re-running a
notebook section loads from cache — no recomputation unless `force = TRUE`
or the cache file is deleted.

## Project structure

```
MyProject/
├── MyProject.Rproj
├── params.yml          ← edit this for each project
├── analysis.qmd        ← render this
├── data/               ← 10x output folders
├── cache/              ← qs2 cache files (auto-generated)
└── output/             ← figures + metadata CSV
```

## Reusing across projects

Copy `params.yml` and `analysis.qmd` from a previous project, update the
`project` block and `samples` list, and you are ready to run. The package
code never needs to be modified for standard analyses.

To free disk space after finalising:

```r
sc_cache_finalize(cfg, dry_run = TRUE)   # preview
sc_cache_finalize(cfg, dry_run = FALSE)  # confirm
```

## Documentation

Full function reference: <https://pfgt-bioinfo.github.io/scReportR/>
