# GAN_example

Reproducible example for training the original conditional single-cell GAN and
scoring blastoid cells with its critic after natural E4.5 calibration.

The calibrated score is **not** the raw WGAN critic output. For natural E4.5
condition `c`, the workflow computes

```text
W1(query critic outputs, natural-calibration critic outputs) /
robust SD of the natural-calibration critic outputs
```

Lower values indicate a critic-output distribution closer to the corresponding
natural E4.5 condition.

## Repository layout

```text
GAN_example/
├── GAN_scoring_example.ipynb          # main, recommended notebook
├── original_training_notebook.ipynb   # original GAN training/generation code
├── data/                              # natural E4.5 and blastoid inputs
├── checkpoints/                       # original critic and generator weights
├── prepared/                          # gene order; matrices are generated here
├── scripts/                           # preprocessing, scoring and plotting
├── results/example/                   # small precomputed example results/figures
├── requirements.txt
├── environment.yml
└── run_pipeline.ps1                   # optional command-line reproduction
```

## Quick start

1. Clone the repository with Git LFS enabled because `.rds` and `.pth` files
   are binary assets.
2. Create the Python environment:

   ```bash
   conda env create -f environment.yml
   conda activate sccgan-example
   ```

3. Install R packages `Seurat`, `Matrix`, and `ggplot2` in R 4.4 or a compatible
   release.
4. Open `GAN_scoring_example.ipynb` from the repository root.
5. Set `RUN_PREPROCESS = True`, run through the preprocessing cell, then set
   `RUN_SCORING = True` and run the scoring cell. Set `RUN_PLOTS = True` to
   regenerate figures.

On Windows, the entire workflow can also be run with:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_pipeline.ps1
```

## Cell-type-specific scoring

- EPI-like cells are compared with natural `Pre_EPI`.
- PE-like cells are compared separately with natural `Pre_VE` and `PaE`.
- TE-like cells are compared separately with natural `mTE` and `pTE`.
- The mismatch challenge evaluates each annotated lineage under all five critic
  labels.

Natural E4.5 is split by biological cell ID into a 70% calibration subset and a
30% internal holdout subset. This split tests the scoring calibration. If
`E450_expression.csv` was derived from the same natural dataset, it must not be
described as an external GAN-training holdout.

## Data and public release

The `data/` folder contains the files used in the local analysis. Before making
the GitHub repository public, verify redistribution permissions for every RDS
file and add the corresponding publication/accession citations to
`data/README.md`. If redistribution is not permitted, remove those files and
replace them with accession/download instructions while retaining the expected
filenames.

## Reproducibility notes

- The original network architecture and supplied `critic.pth` are unchanged.
- All datasets use the same 2,000 genes, fixed-depth downsampling, and
  log-normalization.
- Raw critic values are unbounded and must not be interpreted as probabilities.
- Bootstrap distributions and full calibrated outputs are generated under
  `results/calibrated/`; compact example summaries are stored in
  `results/example/`.

# GAN-scoring
