'''
Permutation tests for the original two-dimensional Poincare embeddings.

Run the epoch and subsample-count combinations in the paper:
    python scripts/poincare_embedding.py --subsamples 10 15 20 --subsample-size 200 --permutations 10000

To use embeddings produced by train_poincare.py:
    python scripts/poincare_embedding.py --embeddings outputs/poincare_training --output outputs/poincare_retrained
'''

import argparse
from pathlib import Path
import sys

import gudhi as gd
import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ApproxPH as ph

SEED = 20211121
EPOCHS = (50, 150, 200)
DIMENSION = 2


def poincare_distances(points):
    margin = 1 - np.sum(points ** 2, axis=1)
    if not np.isfinite(points).all() or np.any(margin <= 0):
        raise ValueError("Poincare vectors must be finite and strictly inside the unit ball.")
    argument = 1 + 2 * cdist(points, points, metric="sqeuclidean") / np.outer(margin, margin)
    distances = np.arccosh(np.maximum(argument, 1))
    np.fill_diagonal(distances, 0)
    return distances


def diagram_bank(points, n, count, seed, prefix, output):
    rng = np.random.default_rng(seed)
    diagrams = []
    for index in range(count):
        indices = rng.choice(len(points), n, replace=False)
        destination = output / f"{prefix}_n{n}_diagram_{index + 1}.npy"
        if destination.exists():
            diagram = np.load(destination)
        else:
            distances = poincare_distances(points[indices])
            rips = gd.RipsComplex(distance_matrix=distances, sparse=0.3)
            tree = rips.create_simplex_tree(max_dimension=2)
            tree.persistence(min_persistence=0.01)
            diagram = tree.persistence_intervals_in_dimension(1)
            diagram = diagram[np.isfinite(diagram).all(axis=1) & (diagram[:, 1] > diagram[:, 0])]
            np.save(destination, diagram)
        diagrams.append(diagram)
        print(f"{prefix}: {index + 1}/{count} diagrams ready", flush=True)
    return diagrams


def permutation_test(first, second, costs, scale, permutations, destination):
    pooled = np.concatenate([first, second])
    split = len(first)

    def statistic(left, right):
        return scale * np.sqrt(max(0, ph.wass_dist(left.mean(axis=0), right.mean(axis=0), costs, reg=None)))

    observed = statistic(first, second)
    statistics = list(np.load(destination)) if destination.exists() else []
    rng = np.random.default_rng(SEED)
    cached = len(statistics)
    for index in range(permutations):
        order = rng.permutation(len(pooled))
        if index < cached:
            continue
        statistics.append(statistic(pooled[order[:split]], pooled[order[split:]]))
        if (index + 1) % 100 == 0 or index + 1 == permutations:
            np.save(destination, statistics)
    statistics = np.asarray(statistics[:permutations])
    p_value = (1 + np.count_nonzero(statistics >= observed)) / (permutations + 1)
    return observed, p_value


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--subsamples", type=int, nargs="+", default=[10, 15, 20])
    parser.add_argument("--subsample-size", type=int, default=200)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--embeddings", type=Path, default=ROOT / "data/poincare")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/poincare_embedding")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    band_path = args.output / "band_points.npy"
    if not band_path.exists():
        rng = np.random.default_rng(SEED)
        theta = rng.uniform(0, 2 * np.pi, 10000)
        radius = rng.uniform(0.98, 0.999, 10000)
        np.save(band_path, radius[:, None] * np.column_stack([np.cos(theta), np.sin(theta)]))
    band = np.load(band_path)
    count = max(args.subsamples)
    band_diagrams = diagram_bank(band, args.subsample_size, count, SEED + 1, "band", args.output)
    costs = ph.dist_mat(ph.mesh_gen(), 2)
    for epoch in EPOCHS:
        points = np.load(args.embeddings / f"embedding_dim{DIMENSION}_epoch{epoch}.npy")
        prefix = f"embedding_dim{DIMENSION}_epoch{epoch}"
        snapshot = args.output / f"{prefix}_points.npy"
        np.save(snapshot, points)
        diagrams = diagram_bank(points, args.subsample_size, count, SEED + epoch, prefix, args.output)
        pooled = diagrams + band_diagrams
        scale = max(1.0, max((diagram[:, 1].max() for diagram in pooled if len(diagram)), default=1.0))
        measures = np.array([np.maximum(ph.diag_to_mesr(diagram / scale, 1)[0] - ph.float_error, 0)
                             for diagram in pooled])
        for first_count in args.subsamples:
            for second_count in args.subsamples:
                first, second = measures[:first_count], measures[count:count + second_count]
                name = f"epoch{epoch}_n{args.subsample_size}_B{first_count}_{second_count}_bank{count}"
                np.save(args.output / f"{name}_mpm.npy", [first.mean(axis=0), second.mean(axis=0)])
                observed, p_value = permutation_test(first, second, costs, scale, args.permutations,
                                                    args.output / f"{name}_permutations.npy")
                np.savez_compressed(args.output / f"{name}_result.npz", statistic=observed, p_value=p_value,
                                    scale=scale, permutations=args.permutations, subsamples=[first_count, second_count],
                                    subsample_size=args.subsample_size, epoch=epoch, seed=SEED,
                                    grid_units=ph.nb_units, min_persistence=0.01, sparse=0.3)
                print(f"epoch={epoch}, B0={first_count}, B1={second_count}: p={p_value:.6g}", flush=True)


if __name__ == "__main__":
    main()
