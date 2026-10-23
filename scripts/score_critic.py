"""Calibrate the original sccGAN critic against natural E4.5 distributions.

This script intentionally keeps the network architecture and ``critic.pth``
weights from ``sccGAN code.ipynb`` unchanged.  The unconstrained WGAN output is
not treated as a probability or a distance.  Instead, for each natural E4.5
condition, it reports a robustly standardised Wasserstein-1 distance between
the query critic-output distribution and a natural calibration distribution.

The input matrices in ``prepared`` were generated using one common gene
order and fixed-depth preprocessing.  Replicate downsamplings are averaged per
biological cell before calculating distribution distances.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn


LABEL_ORDER = ["PaE", "Pre_EPI", "Pre_VE", "mTE", "pTE"]
LINEAGE_CONDITIONS = {
    "EPI": ["Pre_EPI"],
    "PE": ["Pre_VE", "PaE"],
    "TE": ["mTE", "pTE"],
}


def wasserstein_distance(first: np.ndarray, second: np.ndarray) -> float:
    """Exact one-dimensional empirical Wasserstein-1 distance."""
    first = np.sort(np.asarray(first, dtype=float))
    second = np.sort(np.asarray(second, dtype=float))
    support = np.sort(np.concatenate([first, second]))
    if len(support) < 2:
        return 0.0
    first_cdf = np.searchsorted(first, support[:-1], side="right") / len(first)
    second_cdf = np.searchsorted(second, support[:-1], side="right") / len(second)
    return float(np.sum(np.abs(first_cdf - second_cdf) * np.diff(support)))


class OriginalCritic(nn.Module):
    """Architecture copied verbatim from the original notebook."""

    def __init__(self, input_dim: int = 2000, label_dim: int = 5):
        super().__init__()
        layers: list[nn.Module] = []
        current_dim = input_dim + label_dim
        for width in [1024, 512, 256, 128]:
            layers.extend([nn.Linear(current_dim, width), nn.ReLU()])
            current_dim = width
        layers.append(nn.Linear(current_dim, 1))
        self.model = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
        return self.model(torch.cat([x, labels], dim=1))


def load_matrix(data_dir: Path, stem: str, n_genes: int) -> tuple[np.ndarray, pd.DataFrame]:
    meta = pd.read_csv(data_dir / f"{stem}.meta.csv")
    flat = np.memmap(data_dir / f"{stem}.f32", dtype="float32", mode="r")
    expected = len(meta) * n_genes
    if flat.size != expected:
        raise ValueError(f"{stem}: expected {expected} values, found {flat.size}")
    return np.asarray(flat.reshape(len(meta), n_genes)), meta


def all_scores(
    critic: nn.Module,
    matrix: np.ndarray,
    batch_size: int = 512,
) -> np.ndarray:
    result: list[np.ndarray] = []
    critic.eval()
    with torch.no_grad():
        for start in range(0, len(matrix), batch_size):
            x = torch.from_numpy(np.asarray(matrix[start : start + batch_size]).copy())
            rows = []
            for label_index in range(len(LABEL_ORDER)):
                labels = torch.zeros((len(x), len(LABEL_ORDER)), dtype=torch.float32)
                labels[:, label_index] = 1.0
                rows.append(critic(x, labels).squeeze(1).cpu().numpy())
            result.append(np.column_stack(rows))
    return np.vstack(result)


def aggregate_cells(meta: pd.DataFrame, scores: np.ndarray) -> pd.DataFrame:
    scored = meta.copy()
    for index, condition in enumerate(LABEL_ORDER):
        scored[f"critic_{condition}"] = scores[:, index]
    group_columns = [
        column
        for column in ["dataset", "cell_id", "cell_type_raw", "condition", "lineage"]
        if column in scored.columns
    ]
    score_columns = [f"critic_{condition}" for condition in LABEL_ORDER]
    return scored.groupby(group_columns, dropna=False, observed=True, sort=False)[score_columns].mean().reset_index()


def stratified_natural_split(cells: pd.DataFrame, seed: int) -> tuple[set[str], set[str]]:
    rng = np.random.default_rng(seed)
    calibration: set[str] = set()
    test: set[str] = set()
    for _, group in cells.groupby("condition", observed=True, sort=False):
        ids = group["cell_id"].astype(str).unique().to_numpy(copy=True)
        rng.shuffle(ids)
        cut = max(2, int(np.floor(0.70 * len(ids))))
        cut = min(cut, len(ids) - 2)
        calibration.update(ids[:cut])
        test.update(ids[cut:])
    return calibration, test


def robust_scale(values: np.ndarray) -> float:
    q25, q75 = np.quantile(values, [0.25, 0.75])
    scale = float((q75 - q25) / 1.349)
    if not np.isfinite(scale) or scale < 1e-6:
        scale = float(np.std(values, ddof=1))
    return max(scale, 1e-6)


def bootstrap_w1(
    query: np.ndarray,
    reference: np.ndarray,
    scale: float,
    rng: np.random.Generator,
    iterations: int,
) -> np.ndarray:
    values = np.empty(iterations, dtype=float)
    for iteration in range(iterations):
        q = rng.choice(query, size=len(query), replace=True)
        r = rng.choice(reference, size=len(reference), replace=True)
        values[iteration] = wasserstein_distance(q, r) / scale
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("prepared"))
    parser.add_argument(
        "--checkpoint", type=Path, default=Path("checkpoints/critic.pth")
    )
    parser.add_argument(
        "--training-csv", type=Path, default=Path("data/E450_expression.csv")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/calibrated")
    )
    parser.add_argument("--seed", type=int, default=20261006)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    genes = [line.strip() for line in (args.data_dir / "genes.txt").read_text().splitlines() if line.strip()]
    training = pd.read_csv(args.training_csv, nrows=1)
    training_genes = [column for column in training.columns if column not in {training.columns[0], "cell_type"}]
    gene_map = pd.read_csv(args.data_dir / "gene_name_map.csv")
    if (
        genes != gene_map["canonical_gene"].astype(str).tolist()
        or training_genes != gene_map["model_gene_name"].astype(str).tolist()
    ):
        raise ValueError("Prepared genes are not in the exact mapped order used by the original critic")

    critic = OriginalCritic(len(genes), len(LABEL_ORDER))
    state = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    critic.load_state_dict(state)

    natural_matrix, natural_meta = load_matrix(args.data_dir, "E4_5_reference", len(genes))
    natural_cells = aggregate_cells(natural_meta, all_scores(critic, natural_matrix))
    calibration_ids, test_ids = stratified_natural_split(natural_cells, args.seed)
    natural_calibration = natural_cells[natural_cells["cell_id"].astype(str).isin(calibration_ids)].copy()
    natural_test = natural_cells[natural_cells["cell_id"].astype(str).isin(test_ids)].copy()

    references: dict[str, np.ndarray] = {}
    scales: dict[str, float] = {}
    for condition in LABEL_ORDER:
        values = natural_calibration.loc[
            natural_calibration["condition"].eq(condition), f"critic_{condition}"
        ].to_numpy(float)
        references[condition] = values
        scales[condition] = robust_scale(values)

    # Natural E4.5 label-swap challenge.  The same held-out cells are kept
    # fixed while their condition input is deliberately changed.  A paired
    # bootstrap reports the excess distance relative to the correct label.
    natural_mismatch_rows: list[dict[str, object]] = []
    natural_mismatch_boot_rows: list[dict[str, object]] = []
    natural_rng = np.random.default_rng(args.seed + 11)
    for true_condition in LABEL_ORDER:
        true_cells = natural_test[natural_test["condition"].eq(true_condition)].reset_index(drop=True)
        if len(true_cells) < 2:
            continue
        point_distances = {
            supplied_condition: wasserstein_distance(
                true_cells[f"critic_{supplied_condition}"].to_numpy(float),
                references[supplied_condition],
            )
            / scales[supplied_condition]
            for supplied_condition in LABEL_ORDER
        }
        paired_penalties = {condition: [] for condition in LABEL_ORDER}
        paired_distances = {condition: [] for condition in LABEL_ORDER}
        for iteration in range(args.bootstrap):
            query_indices = natural_rng.integers(0, len(true_cells), size=len(true_cells))
            matched_query = true_cells[f"critic_{true_condition}"].to_numpy(float)[query_indices]
            matched_reference = natural_rng.choice(
                references[true_condition], size=len(references[true_condition]), replace=True
            )
            matched_distance = (
                wasserstein_distance(matched_query, matched_reference) / scales[true_condition]
            )
            for supplied_condition in LABEL_ORDER:
                if supplied_condition == true_condition:
                    supplied_distance = matched_distance
                    penalty = 0.0
                else:
                    supplied_query = true_cells[
                        f"critic_{supplied_condition}"
                    ].to_numpy(float)[query_indices]
                    supplied_reference = natural_rng.choice(
                        references[supplied_condition],
                        size=len(references[supplied_condition]),
                        replace=True,
                    )
                    supplied_distance = (
                        wasserstein_distance(supplied_query, supplied_reference)
                        / scales[supplied_condition]
                    )
                    penalty = supplied_distance - matched_distance
                paired_distances[supplied_condition].append(float(supplied_distance))
                paired_penalties[supplied_condition].append(float(penalty))
                natural_mismatch_boot_rows.append(
                    {
                        "true_condition": true_condition,
                        "supplied_condition": supplied_condition,
                        "iteration": iteration,
                        "supplied_distance": float(supplied_distance),
                        "mismatch_penalty": float(penalty),
                    }
                )

        matched_point = point_distances[true_condition]
        for supplied_condition in LABEL_ORDER:
            penalty = point_distances[supplied_condition] - matched_point
            boot_penalty = np.asarray(paired_penalties[supplied_condition], dtype=float)
            natural_mismatch_rows.append(
                {
                    "true_condition": true_condition,
                    "true_lineage": next(
                        lineage
                        for lineage, conditions in LINEAGE_CONDITIONS.items()
                        if true_condition in conditions
                    ),
                    "supplied_condition": supplied_condition,
                    "supplied_lineage": next(
                        lineage
                        for lineage, conditions in LINEAGE_CONDITIONS.items()
                        if supplied_condition in conditions
                    ),
                    "is_correct_condition": supplied_condition == true_condition,
                    "is_correct_lineage": any(
                        true_condition in conditions and supplied_condition in conditions
                        for conditions in LINEAGE_CONDITIONS.values()
                    ),
                    "n_cells": len(true_cells),
                    "matched_distance": matched_point,
                    "supplied_distance": point_distances[supplied_condition],
                    "mismatch_penalty": penalty,
                    "ci_low": float(np.quantile(boot_penalty, 0.025)),
                    "ci_high": float(np.quantile(boot_penalty, 0.975)),
                    "probability_penalty_positive": float(np.mean(boot_penalty > 0)),
                }
            )

    pd.DataFrame(natural_mismatch_rows).to_csv(
        args.output_dir / "natural_e45_celltype_mismatch_matrix.csv", index=False
    )
    pd.DataFrame(natural_mismatch_boot_rows).to_csv(
        args.output_dir / "natural_e45_celltype_mismatch_bootstrap.csv", index=False
    )

    stems = [
        path.name.removesuffix(".meta.csv")
        for path in sorted(args.data_dir.glob("*.meta.csv"))
        if path.name != "E4_5_reference.meta.csv"
    ]
    datasets: dict[str, pd.DataFrame] = {"Natural E4.5 holdout": natural_test}
    for stem in stems:
        matrix, meta = load_matrix(args.data_dir, stem, len(genes))
        datasets[stem] = aggregate_cells(meta, all_scores(critic, matrix))

    rng = np.random.default_rng(args.seed + 1)

    # Complete cell-type mismatch challenge: every annotated blastoid lineage
    # is scored under every natural condition.  Distances stay on the same
    # natural-calibrated scale as the primary score above.
    mismatch_rows: list[dict[str, object]] = []
    mismatch_boot_rows: list[dict[str, object]] = []
    for dataset, cells in datasets.items():
        if dataset == "Natural E4.5 holdout":
            continue
        for annotated_lineage in LINEAGE_CONDITIONS:
            lineage_cells = cells[cells["lineage"].eq(annotated_lineage)]
            if lineage_cells.empty:
                continue
            for supplied_condition in LABEL_ORDER:
                query = lineage_cells[f"critic_{supplied_condition}"].to_numpy(float)
                reference = references[supplied_condition]
                estimate = wasserstein_distance(query, reference) / scales[supplied_condition]
                boot = bootstrap_w1(
                    query, reference, scales[supplied_condition], rng, args.bootstrap
                )
                mismatch_rows.append(
                    {
                        "dataset": dataset,
                        "annotated_lineage": annotated_lineage,
                        "supplied_condition": supplied_condition,
                        "supplied_lineage": next(
                            lineage
                            for lineage, conditions in LINEAGE_CONDITIONS.items()
                            if supplied_condition in conditions
                        ),
                        "within_annotated_lineage": supplied_condition
                        in LINEAGE_CONDITIONS[annotated_lineage],
                        "n_cells": len(query),
                        "distance": estimate,
                        "ci_low": float(np.quantile(boot, 0.025)),
                        "ci_high": float(np.quantile(boot, 0.975)),
                    }
                )
                for iteration, value in enumerate(boot):
                    mismatch_boot_rows.append(
                        {
                            "dataset": dataset,
                            "annotated_lineage": annotated_lineage,
                            "supplied_condition": supplied_condition,
                            "iteration": iteration,
                            "distance": float(value),
                        }
                    )

    mismatch = pd.DataFrame(mismatch_rows)
    matched_min = (
        mismatch[mismatch["within_annotated_lineage"]]
        .groupby(["dataset", "annotated_lineage"], observed=True)["distance"]
        .min()
        .rename("best_matched_distance")
    )
    mismatch = mismatch.join(matched_min, on=["dataset", "annotated_lineage"])
    mismatch["distance_minus_best_matched"] = (
        mismatch["distance"] - mismatch["best_matched_distance"]
    )
    mismatch.to_csv(args.output_dir / "blastoid_celltype_mismatch_matrix.csv", index=False)
    pd.DataFrame(mismatch_boot_rows).to_csv(
        args.output_dir / "blastoid_celltype_mismatch_bootstrap.csv", index=False
    )

    rows: list[dict[str, object]] = []
    boot_rows: list[dict[str, object]] = []
    condition_boot_rows: list[dict[str, object]] = []
    raw_rows: list[pd.DataFrame] = []

    for dataset, cells in datasets.items():
        raw_piece = cells.copy()
        raw_piece.insert(0, "scored_dataset", dataset)
        raw_rows.append(raw_piece)
        for lineage, conditions in LINEAGE_CONDITIONS.items():
            lineage_cells = cells[cells["lineage"].eq(lineage)]
            if lineage_cells.empty:
                continue
            condition_results = []
            for condition in conditions:
                if dataset == "Natural E4.5 holdout":
                    query = lineage_cells.loc[
                        lineage_cells["condition"].eq(condition), f"critic_{condition}"
                    ].to_numpy(float)
                else:
                    query = lineage_cells[f"critic_{condition}"].to_numpy(float)
                if len(query) < 2:
                    continue
                reference = references[condition]
                estimate = wasserstein_distance(query, reference) / scales[condition]
                boot = bootstrap_w1(query, reference, scales[condition], rng, args.bootstrap)
                condition_results.append((condition, estimate, boot, len(query)))

                rows.append(
                    {
                        "dataset": dataset,
                        "level": "condition",
                        "lineage": lineage,
                        "condition": condition,
                        "n_cells": len(query),
                        "distance": estimate,
                        "ci_low": float(np.quantile(boot, 0.025)),
                        "ci_high": float(np.quantile(boot, 0.975)),
                    }
                )
                for iteration, value in enumerate(boot):
                    condition_boot_rows.append(
                        {
                            "dataset": dataset,
                            "lineage": lineage,
                            "condition": condition,
                            "iteration": iteration,
                            "distance": float(value),
                        }
                    )

            if not condition_results:
                continue
            # A blastoid lineage is allowed to match its closest natural stage.
            # For natural holdout, average the stage-specific holdout checks.
            if dataset == "Natural E4.5 holdout":
                chosen = condition_results
                lineage_estimate = float(np.mean([item[1] for item in chosen]))
                lineage_boot = np.mean(np.vstack([item[2] for item in chosen]), axis=0)
                chosen_condition = "+".join(item[0] for item in chosen)
            else:
                best_index = int(np.argmin([item[1] for item in condition_results]))
                chosen_condition, lineage_estimate, lineage_boot, _ = condition_results[best_index]
            rows.append(
                {
                    "dataset": dataset,
                    "level": "lineage",
                    "lineage": lineage,
                    "condition": chosen_condition,
                    "n_cells": len(lineage_cells),
                    "distance": lineage_estimate,
                    "ci_low": float(np.quantile(lineage_boot, 0.025)),
                    "ci_high": float(np.quantile(lineage_boot, 0.975)),
                }
            )
            for iteration, value in enumerate(lineage_boot):
                boot_rows.append(
                    {
                        "dataset": dataset,
                        "lineage": lineage,
                        "iteration": iteration,
                        "distance": float(value),
                    }
                )

    scores = pd.DataFrame(rows)
    lineage_scores = scores[scores["level"].eq("lineage")].copy()
    overall_rows = []
    overall_boot_rows = []
    bootstrap = pd.DataFrame(boot_rows)
    for dataset, group in lineage_scores.groupby("dataset", sort=False):
        # Only complete EPI+PE+TE models receive an overall score.
        if set(group["lineage"]) != set(LINEAGE_CONDITIONS):
            continue
        estimate = float(group.set_index("lineage").loc[list(LINEAGE_CONDITIONS), "distance"].mean())
        wide = bootstrap[bootstrap["dataset"].eq(dataset)].pivot(
            index="iteration", columns="lineage", values="distance"
        )
        overall_boot = wide[list(LINEAGE_CONDITIONS)].mean(axis=1).to_numpy(float)
        overall_rows.append(
            {
                "dataset": dataset,
                "distance": estimate,
                "ci_low": float(np.quantile(overall_boot, 0.025)),
                "ci_high": float(np.quantile(overall_boot, 0.975)),
            }
        )
        for iteration, value in enumerate(overall_boot):
            overall_boot_rows.append(
                {"dataset": dataset, "iteration": iteration, "distance": float(value)}
            )

    pd.concat(raw_rows, ignore_index=True).to_csv(args.output_dir / "cell_critic_outputs.csv", index=False)
    scores.to_csv(args.output_dir / "condition_and_lineage_distances.csv", index=False)
    pd.DataFrame(condition_boot_rows).to_csv(
        args.output_dir / "condition_bootstrap.csv", index=False
    )
    bootstrap.to_csv(args.output_dir / "lineage_bootstrap.csv", index=False)
    overall = pd.DataFrame(overall_rows).sort_values("distance")
    overall.to_csv(args.output_dir / "overall_distances.csv", index=False)
    pd.DataFrame(overall_boot_rows).to_csv(args.output_dir / "overall_bootstrap.csv", index=False)
    pd.DataFrame(
        [
            {
                "condition": condition,
                "natural_calibration_n": len(references[condition]),
                "natural_critic_median": float(np.median(references[condition])),
                "natural_critic_scale": scales[condition],
            }
            for condition in LABEL_ORDER
        ]
    ).to_csv(args.output_dir / "natural_calibration.csv", index=False)

    print(overall.to_string(index=False))
    print(f"Wrote calibrated original-critic results to {args.output_dir}")


if __name__ == "__main__":
    main()
