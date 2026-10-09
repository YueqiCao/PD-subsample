"""Cluster Bearing and Motor shapes using subsampled H1 summaries."""

import argparse
import csv
import hashlib
import json
import sys
from importlib.metadata import version
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import umap
from gudhi.wasserstein import wasserstein_distance
from gudhi.wasserstein.barycenter import lagrangian_barycenter
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ApproxPH as ph


def normalize_shape(points):
    points = np.asarray(points, dtype=float)
    low, high = points.min(axis=0), points.max(axis=0)
    span = (high - low).max()
    if not np.isfinite(points).all() or span <= 0:
        raise ValueError("Expected a finite, nonconstant point cloud.")
    return (points - (low + high) / 2) * (2 / span)


def mean_measure(diagrams):
    mean = np.zeros(ph.mat_size)
    for diagram in diagrams:
        for birth, death in diagram:
            if not np.isfinite([birth, death]).all() or not 0 <= birth <= death <= ph.grid_width:
                raise ValueError("Expected finite diagrams within the persistence grid.")
            if death <= birth:
                continue
            i = min(int(np.floor(birth / ph.unit)), ph.nb_units - 1)
            j = min(int(np.ceil(death / ph.unit)), ph.nb_units)
            j = max(j, i + 1)
            mean[ph.nb_units * i + j - i * (i + 1) // 2 - 1] += 1 / len(diagrams)
    return mean


def diagram_bank(path, index, args, cache):
    settings = dict(object=f"{path.parent.name}/{path.name}",
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    seed=args.seed + 10000 * index, n=args.subsample_size, B=args.subsamples,
                    normalization="isotropic_bbox", cutoff=0.4, min_persistence=0.03, sparse=0.3)
    destination = cache / f"{path.parent.name}_{path.stem}.npz"
    if destination.exists():
        with np.load(destination) as saved:
            if json.loads(str(saved["settings"])) != settings:
                raise ValueError(f"Diagram settings or input changed: {destination}. Use a new cache.")
            return [saved[f"diagram_{i}"] for i in range(args.subsamples)]
    if args.diagrams:
        raise FileNotFoundError(f"Missing saved diagram bank: {destination}")
    points = normalize_shape(np.load(path))
    if args.subsample_size > len(points):
        raise ValueError(f"Subsample size exceeds {path.name}'s point count.")
    rng = np.random.default_rng(settings["seed"])
    diagrams = []
    for _ in range(args.subsamples):
        subset = points[rng.choice(len(points), args.subsample_size, replace=False)]
        diag = ph.get_PD(subset, max_edge_length=0.4, min_persistence=0.03, sparse=0.3)
        diag[np.isposinf(diag[:, 1]), 1] = 0.4
        diagrams.append(diag[diag[:, 1] > diag[:, 0]])
    np.savez_compressed(destination, settings=json.dumps(settings),
                        **{f"diagram_{i}": diag for i, diag in enumerate(diagrams)})
    return diagrams


def distance_matrix(summaries, method):
    count = len(summaries)
    distances = np.zeros((count, count))
    costs = ph.dist_mat(ph.mesh_gen(), 2) if method == "mpm" else None
    for i in range(count):
        for j in range(i + 1, count):
            if method == "mpm":
                value = np.sqrt(max(0, ph.wass_dist(summaries[i], summaries[j], costs, reg=None)))
            else:
                value = wasserstein_distance(summaries[i], summaries[j], order=2, internal_p=np.inf)
            distances[i, j] = distances[j, i] = value
    return distances


def cluster(distances, truth, neighbors=15, eps=1.0, min_samples=10, seed=30):
    embedding = umap.UMAP(metric="precomputed", n_neighbors=neighbors, min_dist=0.1,
                          random_state=seed, n_jobs=1).fit_transform(distances)
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(embedding)
    assigned = labels >= 0
    correct = sum(np.unique(truth[labels == k], return_counts=True)[1].max()
                  for k in np.unique(labels[assigned]))
    scores = dict(adjusted_rand_index=float(adjusted_rand_score(truth, labels)),
                  cluster_count=int(len(np.unique(labels[assigned]))),
                  noise_count=int((~assigned).sum()), coverage=float(assigned.mean()),
                  purity_assigned=float(correct / assigned.sum()) if assigned.any() else None,
                  correct_fraction_all=float(correct / len(labels)))
    return embedding, labels, scores


def plot_clusters(embedding, labels, truth, method, destination):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for label in np.unique(truth):
        mask = truth == label
        axes[0].scatter(*embedding[mask].T, s=25, label=label)
    for label in np.unique(labels):
        mask = labels == label
        options = dict(color="0.65", marker="x") if label == -1 else {}
        axes[1].scatter(*embedding[mask].T, s=25,
                        label="Noise" if label == -1 else f"Cluster {label}", **options)
    for ax, title in zip(axes, ("Known classes", "DBSCAN clusters")):
        ax.set_title(f"{method.upper()}: {title}")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(fontsize=8)
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/experiments/shape_clustering")
    parser.add_argument("--diagrams", type=Path, help="Reuse a previous run's diagrams directory.")
    parser.add_argument("--methods", choices=["mpm", "fm"], nargs="+", default=["mpm"])
    parser.add_argument("--subsample-size", type=int, default=1500)
    parser.add_argument("--subsamples", type=int, default=15)
    parser.add_argument("--seed", type=int, default=20211121)
    parser.add_argument("--neighbors", type=int, help="UMAP neighbors: MPM 15, FM 10 by default.")
    parser.add_argument("--eps", type=float, help="DBSCAN radius: MPM 1, FM 2 by default.")
    parser.add_argument("--min-samples", type=int, default=10)
    parser.add_argument("--umap-seed", type=int, help="UMAP seed: MPM 30, FM 20 by default.")
    parser.add_argument("--limit-per-class", type=int, help="Use a small subset for a smoke test.")
    args = parser.parse_args()
    files = []
    for name in ("Bearing", "Motor"):
        group = sorted((args.data / name).glob("*.npy"), key=lambda path: int(path.stem))
        files.extend(group[:args.limit_per_class])
    method_settings = {}
    for method in args.methods:
        neighbors, eps, seed = (15, 1.0, 30) if method == "mpm" else (10, 2.0, 20)
        method_settings[method] = dict(neighbors=args.neighbors if args.neighbors is not None else neighbors,
                                       eps=args.eps if args.eps is not None else eps,
                                       seed=args.umap_seed if args.umap_seed is not None else seed,
                                       min_samples=args.min_samples)
    args.output.mkdir(parents=True, exist_ok=False)
    cache = args.diagrams or args.output / "diagrams"
    if not args.diagrams:
        cache.mkdir()
    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    config.update(objects=[f"{path.parent.name}/{path.name}" for path in files],
                  clustering=method_settings,
                  normalization="isotropic_bbox", cutoff=0.4, min_persistence=0.03, sparse=0.3,
                  grid_units=ph.nb_units, grid_width=ph.grid_width, quantization="birth-death",
                  mpm_distance="OT2",
                  fm="OT2 barycenter",
                  versions={name: version(name) for name in
                            ("numpy", "scipy", "gudhi", "POT", "umap-learn", "scikit-learn")},
                  python=sys.version.split()[0])
    (args.output / "settings.json").write_text(json.dumps(config, indent=2), encoding="utf8")
    summaries = {method: [] for method in args.methods}
    start = perf_counter()
    for index, path in enumerate(files):
        diagrams = diagram_bank(path, index, args, cache)
        if "mpm" in summaries:
            summaries["mpm"].append(mean_measure(diagrams))
        if "fm" in summaries:
            # GUDHI returns a local Frechet mean; its iterative solver can be slow.
            summaries["fm"].append(lagrangian_barycenter(diagrams, init=0))
        print(f"{index + 1}/{len(files)}: {path.parent.name}/{path.name}", flush=True)
    truth = np.array([path.parent.name for path in files])
    reports = {}
    for method, values in summaries.items():
        distances = distance_matrix(values, method)
        embedding, labels, scores = cluster(distances, truth, **method_settings[method])
        np.savez_compressed(args.output / f"{method}.npz", distances=distances,
                            embedding=embedding, labels=labels, truth=truth)
        plot_clusters(embedding, labels, truth, method, args.output / f"{method}.png")
        with (args.output / f"{method}_labels.csv").open("w", newline="", encoding="utf8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["object", "class", "cluster", "umap_x", "umap_y"])
            writer.writerows((f"{p.parent.name}/{p.name}", t, int(k), float(x), float(y))
                             for p, t, k, (x, y) in zip(files, truth, labels, embedding))
        reports[method] = scores
        print(method, json.dumps(scores), flush=True)
    reports["elapsed_seconds"] = perf_counter() - start
    (args.output / "summary.json").write_text(json.dumps(reports, indent=2), encoding="utf8")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
