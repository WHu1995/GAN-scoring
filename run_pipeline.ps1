param(
    [string]$Python = "python",
    [string]$Rscript = "C:\Program Files\R\R-4.4.2\bin\Rscript.exe",
    [int]$Bootstrap = 2000
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path

& $Rscript `
    (Join-Path $ProjectRoot "scripts\prepare_inputs.R") `
    (Join-Path $ProjectRoot "data") `
    (Join-Path $ProjectRoot "prepared") `
    "1000" "4" "2" "20260930"

& $Python `
    (Join-Path $ProjectRoot "scripts\score_critic.py") `
    --data-dir (Join-Path $ProjectRoot "prepared") `
    --checkpoint (Join-Path $ProjectRoot "checkpoints\critic.pth") `
    --training-csv (Join-Path $ProjectRoot "data\E450_expression.csv") `
    --output-dir (Join-Path $ProjectRoot "results\calibrated") `
    --seed 20261006 `
    --bootstrap $Bootstrap

foreach ($PlotScript in @(
    "plot_summary.R",
    "plot_four_conditions.R",
    "plot_mismatch_heatmap.R"
)) {
    & $Rscript `
        (Join-Path $ProjectRoot "scripts\$PlotScript") `
        (Join-Path $ProjectRoot "results\calibrated")
}

Write-Host "Completed. Results: $(Join-Path $ProjectRoot 'results\calibrated')"

