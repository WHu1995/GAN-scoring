# GAN_scoring

Example for training the conditional single-cell GAN and
scoring blastoid cells with its critic after natural E4.5 calibration.

Lower values indicate a critic-output distribution closer to the corresponding
natural E4.5 condition.

## Repository layout

```text
GAN_example/
├── GAN_scoring_example.ipynb          # main, recommended notebook
├── GAN_training.ipynb   # original GAN training/generation code
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


