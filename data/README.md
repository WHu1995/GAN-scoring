# Input data

Expected files:

| File | Role |
|---|---|
| `E450_expression.csv` | Original 2,000-gene GAN training matrix |
| `E4.5.rds` | Natural E4.5 reference Seurat object |
| `Jorgensen_iGata4_blastoids.rds` | iGata4 blastoid query |
| `Xu_TPS_blastoids.rds` | Xu TPS blastoid query |
| `Li_EPSC_blastoids.rds` | Li EPSC blastoid query |
| `Bing_TBLC_blastoids.rds` | Bing TBLC blastoid query |
| `Sozen_EPSC_TSC.rds` | Sozen EPSC + TSC query |
| `Rivron_ESC_TSC_blastoids.rds` | Rivron ESC + TSC query |

Each RDS file is expected to contain a Seurat object with RNA counts and a
metadata column named `cell_type`. Query annotations should use `EPI-like`,
`PE-like`, and `TE-like`. Natural reference annotations should use `Pre_EPI`,
`Pre_VE`/`PE`, `PaE`, `mTE`, and `pTE` as applicable.

Before public release, add the source publication, accession, preprocessing
history, and redistribution license for each dataset. If a dataset cannot be
redistributed, replace the file with documented download instructions.

