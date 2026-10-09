'''
Experiment for parameter tuning.

Run from the repository directory:
    python scripts/parameter_tuning.py --datasets torus sphere --max-subsamples 300

To select particular subsample sizes and transport powers:
    python scripts/parameter_tuning.py --datasets torus --subsample-sizes 100 200 500 --powers 2 5 7 10
'''

import argparse
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ApproxPH as ph

SEED = 20211121
TOLERANCE = 1e-3
CUTOFFS = {"torus": 0.9, "sphere": 0.45}


def compute_mean(points, n, count, cutoff, rng, destination):
    cached = destination.exists()
    mean = np.zeros(ph.mat_size)
    for _ in range(count):
        indices = rng.choice(len(points), n, replace=False)
        if cached:
            continue
        diagram = ph.get_PD(points[indices], max_edge_length=cutoff, min_persistence=0.01)
        diagram[np.isposinf(diagram[:, 1]), 1] = cutoff
        measure, _ = ph.diag_to_mesr(diagram, 1 / count)
        mean += measure - ph.float_error
    if cached:
        return np.load(destination)
    mean += ph.float_error
    np.save(destination, mean)
    return mean


def tune(points, dataset, n, power, maximum, costs, output):
    rng = np.random.RandomState(SEED + 10000 * n + power)
    prefix = f"{dataset}_n{n}_p{power}"
    counts, errors, ratios, times = [], [], [], []
    previous_mean, previous_error = None, None
    converged = False
    start = perf_counter()
    for count in range(1, maximum + 1):
        mean = compute_mean(points, n, count, CUTOFFS[dataset], rng,
                            output / f"{prefix}_B{count}_mpm.npy")
        if previous_mean is not None:
            error = ph.wass_dist(mean, previous_mean, costs, reg=1)
            ratio = np.inf
            if previous_error is not None:
                ratio = abs(error / previous_error - 1) if previous_error > 0 else (0.0 if error == 0 else np.inf)
            counts.append(count)
            errors.append(error)
            ratios.append(ratio)
            times.append(perf_counter() - start)
            converged = ratio < TOLERANCE
            np.savez_compressed(output / f"{prefix}_history.npz", counts=counts, errors=errors,
                                relative_changes=ratios, elapsed_seconds=times, converged=converged,
                                tolerance=TOLERANCE, seed=SEED + 10000 * n + power)
            print(f"{dataset}, n={n}, p={power}, B={count}: cost={error:.6g}, change={ratio:.6g}", flush=True)
            previous_error = error
        previous_mean = mean
        if converged:
            break
    print(f"{prefix}: {'stability threshold reached' if converged else 'subsample limit reached'}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", choices=list(CUTOFFS), nargs="+", default=list(CUTOFFS))
    parser.add_argument("--subsample-sizes", type=int, nargs="+", default=[100, 200, 500, 1000, 1500, 2000])
    parser.add_argument("--powers", type=int, nargs="+", default=[2, 5, 7, 10])
    parser.add_argument("--max-subsamples", type=int, default=300)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/parameter_tuning")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    ground_cost = ph.dist_mat(ph.mesh_gen(), 1)
    for dataset in args.datasets:
        points = np.load(ROOT / "data/parameter_tuning" / f"{dataset}.npy")
        for n in args.subsample_sizes:
            for power in args.powers:
                tune(points, dataset, n, power, args.max_subsamples, ground_cost ** power, args.output)


if __name__ == "__main__":
    main()
