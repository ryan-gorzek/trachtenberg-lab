# Trachtenberg Lab

Code from my PhD in Josh Trachtenberg’s lab (UCLA) spanning systems, anatomy, and transcriptomics of visual cortex across development and evolution. This repo contains analysis notebooks, figure-generation scripts, and lightweight utilities used to process raw data through to publication-quality figures.

![plot](https://github.com/ryan-gorzek/trachtenberg-lab/blob/main/plots/TLab_Summary.png)

Also be sure to check out my other self-contained projects conducted in the lab:
- [comparatome](https://github.com/ryan-gorzek/comparatome) — R package for cross-species single-cell and spatial transcriptomic analysis.
- [intergenome](https://github.com/ryan-gorzek/intergenome) — Nextflow pipeline for STARsolo alignment and intergenic read quantification.
- [electroviz](https://github.com/ryan-gorzek/electroviz) — Python toolkit for extracellular electrophysiology data visualization and analysis.
- [celltype-tuning-V1](https://github.com/ryan-gorzek/celltype-tuning-V1) — MATLAB workflow for visually evoked neural activity-based ML classification of cortical cell types.

## Repository Layout
* `retinotopy/` — functional imaging and map estimation (stimulus logs → ΔF/F processing → phase maps and visual field fits).
* `tuning/` — neuronal response characterization (OSI/DSI, SF/TF tuning, contrast response; per-cell and population summaries).
* `tracing/` — anatomical projection quantification (registration, ROI masks, axon density metrics, laminar profiles).
* `transcriptomics/` — single-nucleus/cell and spatial transcriptomics (QC, normalization, DE, clustering, integration, gene-set analyses).
* `modeling/` — statistical and computational modeling of neural circuits (theoretical connections between simulated cell types).
* `README.md`, `LICENSE`, `.gitignore`.

Structure and top-level folders match the public repo listing. ([GitHub][1])

## Related Preprints

Use and citation of this code should refer to the relevant manuscript(s):

* Gorzek R, Trachtenberg JT. 2025. **Comparative transcriptomics reveals shifts in cortical architecture at the metatherian/eutherian transition.** *bioRxiv* doi:10.1101/2025.10.03.680397. ([BioRxiv][2])
  Figure-generation notebooks and cross-species transcriptomic tools live under `transcriptomics/V1` and `transcriptomics/Stereo-seq`. These are currently being curated into a manuscript-specific repository and R package, respectively.

* Yoo J, Xie F, Butrus S, Xu R, Tan Z, **Gorzek R**, Mirshahidi P, Tring E, Suresh S, Kim J, Fleishman G, Tan L, Ringach D, Trachtenberg JT, Xu X, Zipursky SL, Shekhar K, Jain S. 2025. **Establishing a continuum of cell types in the visual cortex.** *bioRxiv* doi:10.1101/2025.09.22.677893. ([BioRxiv][3])
  Supporting scripts for V1 feedforward connections to LM, RL, and AM/PM in NR vs. DR are located in `tracing/`.

## Environments

The codebase mixes Python, R, and MATLAB. Minimal working versions:

* Python ≥ 3.9 with: `numpy`, `pandas`, `scanpy`, `anndata`, `scikit-learn`, `matplotlib`, `seaborn`, `statsmodels`.
* R ≥ 4.3 with: `Seurat`, `SeuratObject`, `Matrix`, `dplyr`, `ggplot2`, `patchwork`, `matrixStats`, `glmGamPoi`.
* MATLAB ≥ R2023a with Image Processing Toolbox (for legacy figure scripts).

## License

MIT. See `LICENSE`.
