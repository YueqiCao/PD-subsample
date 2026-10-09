'''
Cluster the Bearing and Motor shapes using multiple subsampling.

Run from the repository directory:
    python scripts/shape_clustering.py --subsamples 15 --subsample-size 1500 --methods mpm

To compute both MPM and Frechet mean from the same subsamples:
    python scripts/shape_clustering.py --subsamples 15 --subsample-size 1500 --methods mpm fm

Diagrams are saved automatically and reused when the command is rerun.
'''

import argparse
import csv
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

SEED = 20211121
DATA = ROOT / "data"
CLUSTERING = {
    "mpm": dict(neighbors=15, eps=1.75, min_samples=10, seed=30),
    "fm": dict(neighbors=10, eps=2.0, min_samples=10, seed=20),
}


def normalize_shape(points):
    points = np.asarray(points, dtype=float)
    low, high = points.min(axis=0), points.max(axis=0)
    span = (high - low).max()
    return (points - (low + high) / 2) * (2 / span)


def mean_measure(diagrams):
    mean = np.zeros(ph.mat_size)
    for diagram in diagrams:
        for birth, death in diagram:
            i = min(int(np.floor(birth / ph.unit)), ph.nb_units - 1)
            j = min(int(np.ceil(death / ph.unit)), ph.nb_units)
            j = max(j, i + 1)
            mean[ph.nb_units * i + j - i * (i + 1) // 2 - 1] += 1 / len(diagrams)
    return mean


def diagram_bank(path, index, args, cache):
    settings = dict(object=f"{path.parent.name}/{path.name}",
                    seed=SEED + 10000 * index, n=args.subsample_size, B=args.subsamples,
                    normalization="isotropic_bbox", cutoff=0.4, min_persistence=0.03, sparse=0.3)
    destination = cache / f"{path.parent.name}_{path.stem}.npz"
    if destination.exists():
        with np.load(destination) as saved:
            return [saved[f"diagram_{i}"] for i in range(args.subsamples)]
    points = normalize_shape(np.load(path))
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


def cluster(distances, truth, neighbors=15, eps=1.75, min_samples=10, seed=30):
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
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/shape_clustering")
    parser.add_argument("--methods", choices=["mpm", "fm"], nargs="+", default=["mpm"])
    parser.add_argument("--subsample-size", type=int, default=1500)
    parser.add_argument("--subsamples", type=int, default=15)
    args = parser.parse_args()
    files = []
    for name in ("Bearing", "Motor"):
        files.extend(sorted((DATA / name).glob("*.npy"), key=lambda path: int(path.stem)))
    args.output.mkdir(parents=True, exist_ok=True)
    cache = args.output / "diagrams"
    cache.mkdir(exist_ok=True)
    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    config.update(objects=[f"{path.parent.name}/{path.name}" for path in files],
                  seed=SEED, clustering={method: CLUSTERING[method] for method in args.methods},
                  normalization="isotropic_bbox", cutoff=0.4, min_persistence=0.03, sparse=0.3,
                  grid_units=ph.nb_units, grid_width=ph.grid_width, quantization="birth-death",
                  mpm_distance="OT2",
                  fm="OT2 barycenter",
                  versions={name: version(name) for name in
                            ("numpy", "scipy", "gudhi", "POT", "umap-learn", "scikit-learn")},
                  python=sys.version.split()[0])
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
    (args.output / "settings.json").write_text(json.dumps(config, indent=2), encoding="utf8")
    reports = {}
    for method, values in summaries.items():
        distances = distance_matrix(values, method)
        embedding, labels, scores = cluster(distances, truth, **CLUSTERING[method])
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
