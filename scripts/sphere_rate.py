'''
Sphere convergence rate experiment.

Run from the repository directory:
    python scripts/sphere_rate.py --max-size 3800 --b-ratio 0.1

Each mean is saved as mean_mesr_n<N>_B<B>.npy and reused on reruns.
'''

import argparse
import json
import sys
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ApproxPH as ph

SEED = 20211121
POWERS = (2, 9)
N_POINTS = 20000
DIMENSION = 4  # The 3-sphere is embedded in R^4.
RADIUS = 0.5


def compute_measures(args, output):
    points_path = output / "true-sphere-points.npy"
    reference_path = output / "true-sphere-diagram.npy"
    if not points_path.exists():
        points = np.random.default_rng(SEED).normal(size=(N_POINTS, DIMENSION))
        points *= RADIUS / np.linalg.norm(points, axis=1, keepdims=True)
        np.save(points_path, points)
    points = np.load(points_path)
    means, counts = [], []
    rng = np.random.RandomState(SEED)
    for n in args.sizes:
        count = max(1, round(args.b_ratio * n))
        cache = output / f"mean_mesr_n{n}_B{count}.npy"
        cached = cache.exists()
        mean = np.zeros(ph.mat_size)
        for _ in range(count):
            indices = rng.choice(len(points), n, replace=False)
            if cached:
                continue
            subset = points[indices]
            diag = ph.get_PD(subset, max_edge_length=0.45, min_persistence=0.01, sparse=0.3)
            diag[np.isposinf(diag[:, 1]), 1] = 0.45
            measure, _ = ph.diag_to_mesr(diag, 1 / count)
            mean += measure - ph.float_error
        if cached:
            means.append(np.load(cache))
        else:
            means.append(mean + ph.float_error)
            np.save(cache, means[-1])
            print(f"n={n}, B={count}: mean saved", flush=True)
        counts.append(count)
    if not reference_path.exists():
        print("Computing the sphere reference diagram...", flush=True)
        diagram = ph.get_PD(points, max_edge_length=0.45, min_persistence=0.01, sparse=0.3)
        np.save(reference_path, diagram)
    diagram = np.load(reference_path)
    diagram = diagram[diagram[:, 1] - diagram[:, 0] > 0.01]
    diagram[np.isposinf(diagram[:, 1]), 1] = 0.45
    reference, _ = ph.diag_to_mesr(diagram, 1)
    return np.asarray(args.sizes), np.asarray(counts), np.asarray(means), reference


def analyze(sizes, means, reference, output):
    results = {}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ground_cost = ph.dist_mat(ph.mesh_gen(), 1)
    for ax, power in zip(axes, POWERS):
        costs = ground_cost ** power
        losses = np.array([ph.wass_dist(mean, reference, costs, reg=1) for mean in means])
        params = ph.fit_power_law(sizes, losses)
        predicted = ph.func(sizes, *params)
        r2 = 1 - np.sum((losses - predicted) ** 2) / np.sum((losses - losses.mean()) ** 2)
        results[f"p{power}"] = dict(losses=losses.tolist(), a=float(params[0]),
                                    b=float(params[1]), c=float(params[2]), r_squared=float(r2))
        ax.scatter(sizes, losses, label="Empirical loss", s=24)
        smooth_n = np.linspace(min(sizes), max(sizes), 200)
        ax.plot(smooth_n, ph.func(smooth_n, *params), label=f"Fit: b={params[1]:.2f}")
        ax.set(xlabel="Subsample size n", ylabel="Transport loss", title=f"p={power}")
        ax.legend()
        print(f"p={power}: exponent=-{params[1]:.4f}, R2={r2:.4f}", flush=True)
    fig.savefig(output / "rate.png", dpi=160)
    plt.close(fig)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-size", type=int, default=3800, help="Largest subsample size, in steps of 200.")
    parser.add_argument("--b-ratio", type=float, default=0.1, help="B=round(ratio*n), at least one.")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/sphere_rate")
    args = parser.parse_args()
    args.sizes = list(range(400, args.max_size + 1, 200))
    args.output.mkdir(parents=True, exist_ok=True)
    config = dict(sizes=args.sizes, b_ratio=args.b_ratio, seed=SEED, powers=POWERS,
                  sampling="sphere_S3", population_size=N_POINTS, dimension=DIMENSION, radius=RADIUS,
                  cutoff=0.45, min_persistence=0.01, sparse=0.3, sinkhorn_regularization=1,
                  versions={name: version(name) for name in ("numpy", "scipy", "gudhi", "POT")})
    sizes, counts, means, reference = compute_measures(args, args.output)
    (args.output / "settings.json").write_text(json.dumps(config, indent=2), encoding="utf8")
    results = analyze(sizes, means, reference, args.output)
    report = {"sizes": sizes.tolist(), "counts": counts.tolist(), "fits": results}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf8")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
