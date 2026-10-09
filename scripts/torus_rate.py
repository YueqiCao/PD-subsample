"""Torus rate experiment"""

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


def load_measures(path, sizes, b_ratio):
    if path.is_file():
        with np.load(path) as saved:
            return tuple(saved[key].copy() for key in ("sizes", "counts", "means", "reference"))
    files = list(path.glob("mean_mesr_nb*.npy"))
    if len(files) != len(sizes):
        raise ValueError("--sizes must describe every saved mean_mesr_nb*.npy file in order.")
    means = np.stack([np.load(path / f"mean_mesr_nb{i}.npy") for i in range(len(sizes))])
    diagram = np.load(path / "true-torus-diagram.npy")
    diagram[np.isposinf(diagram[:, 1]), 1] = 0.9
    reference, _ = ph.diag_to_mesr(diagram, 1)
    return np.asarray(sizes), np.maximum(1, np.rint(b_ratio * np.asarray(sizes))).astype(int), means, reference


def compute_measures(args, output):
    points = np.load(args.points)
    diagram = np.load(args.reference)
    diagram[np.isposinf(diagram[:, 1]), 1] = 0.9
    reference, _ = ph.diag_to_mesr(diagram, 1)
    rng = np.random.default_rng(args.seed)
    means, counts = [], []
    for n in args.sizes:
        if n > len(points):
            raise ValueError("Subsample size exceeds the reference point cloud.")
        count = max(1, round(args.b_ratio * n))
        mean = np.zeros(ph.mat_size)
        for _ in range(count):
            subset = points[rng.choice(len(points), n, replace=False)]
            diag = ph.get_PD(subset, max_edge_length=0.9, min_persistence=0.01, sparse=0.3)
            diag[np.isposinf(diag[:, 1]), 1] = 0.9
            measure, _ = ph.diag_to_mesr(diag, 1 / count)
            mean += measure - ph.float_error
        means.append(mean + ph.float_error)
        counts.append(count)
        # Checkpoint only the means, never the sampled point clouds or simplices.
        np.savez_compressed(output / "measures.npz", sizes=args.sizes[:len(means)],
                            counts=counts, means=means, reference=reference)
        print(f"n={n}, B={count}: mean saved", flush=True)
    return np.asarray(args.sizes), np.asarray(counts), np.asarray(means), reference


def analyze(sizes, means, reference, powers, output):
    if means.shape != (len(sizes), ph.mat_size) or reference.shape != (ph.mat_size,):
        raise ValueError("Saved measures do not match the 50-unit grid.")
    if not np.all(np.isfinite(means)) or np.any(means < 0):
        raise ValueError("Saved measures must be finite and nonnegative.")
    results = {}
    fig, axes = plt.subplots(2, len(powers), figsize=(5 * len(powers), 8), squeeze=False,
                             constrained_layout=True)
    # Compute the ground costs once, then raise them to p=3 and p=8.
    ground_cost = ph.dist_mat(ph.mesh_gen(), 1)
    for column, power in enumerate(powers):
        costs = ground_cost ** power
        for row, solver in enumerate(("sinkhorn_reg1", "exact")):
            if solver == "exact":
                # Remove the legacy numerical background before exact transport.
                source = np.maximum(means - ph.float_error, 0)
                target = np.maximum(reference - ph.float_error, 0)
                reg = None
            else:
                source, target, reg = means, reference, 1
            losses = np.array([ph.wass_dist(mean, target, costs, reg=reg) for mean in source])
            params = ph.fit_power_law(sizes, losses)
            predicted = ph.func(sizes, *params)
            r2 = 1 - np.sum((losses - predicted) ** 2) / np.sum((losses - losses.mean()) ** 2)
            results[f"p{power}_{solver}"] = {
                "losses": losses.tolist(), "a": float(params[0]),
                "b": float(params[1]), "c": float(params[2]), "r_squared": float(r2),
                "exponent_at_bound": bool(params[1] < 0.0101 or params[1] > 4.9999),
            }
            ax = axes[row, column]
            ax.scatter(sizes, losses, label="Observed cost", s=24)
            smooth_n = np.linspace(min(sizes), max(sizes), 200)
            ax.plot(smooth_n, ph.func(smooth_n, *params), label=f"Fit: b={params[1]:.3f}")
            ax.set(xlabel="Subsample size n", ylabel="Powered transport cost",
                   title=f"p={power}: {'Sinkhorn, reg=1' if row == 0 else 'Exact OT'}")
            ax.legend()
            print(f"{solver}, p={power}: exponent=-{params[1]:.4f}, R2={r2:.4f}", flush=True)
    fig.savefig(output / "rate.png", dpi=160)
    plt.close(fig)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-measures", type=Path,
                        help="Replot a measures.npz file or the former notebook's outputs directory.")
    parser.add_argument("--points", type=Path, default=ROOT / "outputs/true-torus-points.npy")
    parser.add_argument("--reference", type=Path, default=ROOT / "outputs/true-torus-diagram.npy")
    parser.add_argument("--sizes", type=int, nargs="+", default=list(range(400, 3201, 200)))
    parser.add_argument("--b-ratio", type=float, default=0.1, help="B=round(ratio*n), at least one.")
    parser.add_argument("--powers", type=int, nargs="+", default=[3, 8])
    parser.add_argument("--seed", type=int, default=20211121)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/experiments/torus_rate")
    args = parser.parse_args()
    if len(set(args.sizes)) < 4 or args.sizes != sorted(set(args.sizes)) or min(args.sizes) < 2:
        parser.error("Use at least four distinct, increasing sample sizes >=2.")
    if args.b_ratio <= 0 or min(args.powers) <= 0:
        parser.error("The subsampling ratio and powers must be positive.")
    args.output.mkdir(parents=True, exist_ok=False)
    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    config.update(grid_units=ph.nb_units, grid_width=ph.grid_width, quantization="floor birth/death",
                  cutoff=0.9, min_persistence=0.01, reference_min_persistence=0.0, sparse=0.3,
                  versions={name: version(name) for name in ("numpy", "scipy", "gudhi", "POT")},
                  python=sys.version.split()[0])
    (args.output / "settings.json").write_text(json.dumps(config, indent=2), encoding="utf8")
    if args.from_measures:
        sizes, counts, means, reference = load_measures(args.from_measures, args.sizes, args.b_ratio)
    else:
        sizes, counts, means, reference = compute_measures(args, args.output)
    results = analyze(sizes, means, reference, args.powers, args.output)
    report = {"sizes": sizes.tolist(), "counts": counts.tolist(), "fits": results,
              "interpretation": "Fit a*n^(-b)+c to powered costs, not their pth roots. "
              "Sinkhorn reg=1 is the original regularized proxy, not exact OT. "
              "The theoretical upper bound does not require an empirical exponent of -0.5. "
              "These are single-run fits, not estimates of the expected loss over repeated experiments."}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf8")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
