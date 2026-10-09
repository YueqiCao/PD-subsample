'''
Approximate the persistent homology of the Knot and Lock point clouds.

Run both datasets with their full settings from the repository directory:
    python scripts/knot_lock.py

To run each dataset separately:
    python scripts/knot_lock.py --datasets knot --subsamples 25 --subsample-size 9500
    python scripts/knot_lock.py --datasets lock --subsamples 30 --subsample-size 9000

Lock uses data/grayloc.ply. Completed diagrams are saved and reused on reruns.
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
from gudhi.representations import PersistenceImage
from gudhi.wasserstein.barycenter import lagrangian_barycenter as bary
from plyfile import PlyData

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ApproxPH as ph

SEED = 20211121
CUTOFF = 0.55
DATASETS = {
    "knot": dict(file="knot.ply", subsamples=25, subsample_size=9500, min_persistence=0.1),
    "lock": dict(file="grayloc.ply", subsamples=30, subsample_size=9000, min_persistence=0.07),
}


def compute_means(points, name, n, count, min_persistence, output):
    np.random.seed(SEED)
    subs = ph.get_subsample(points, n, count)
    prefix = f"{name}_n{n}_B{count}"
    diags = []
    for index, subset in enumerate(subs):
        destination = output / f"{prefix}_diagram_{index + 1}.npy"
        if destination.exists():
            diag = np.load(destination)
        else:
            diag = ph.get_PD(subset, max_edge_length=CUTOFF, min_persistence=min_persistence)
            diag[np.isposinf(diag[:, 1]), 1] = CUTOFF
            diag = diag[diag[:, 1] > diag[:, 0]]
            np.save(destination, diag)
        diags.append(diag)
        print(f"{name}: {index + 1}/{count} diagrams ready", flush=True)

    mean_mesr, mean_mesr_vis = ph.diag_to_mesr(np.concatenate(diags), 1 / count)
    np.save(output / f"{prefix}_mpm.npy", mean_mesr)
    np.save(output / f"{prefix}_mpm_grid.npy", mean_mesr_vis)
    image = PersistenceImage(bandwidth=0.02, weight=lambda point: point[1],
                             resolution=[50, 50], im_range=[0, CUTOFF, 0, CUTOFF])
    mean_image = image.fit_transform(diags).mean(axis=0).reshape(50, 50)
    np.save(output / f"{prefix}_mean_image.npy", mean_image)
    destination = output / f"{prefix}_fm.npy"
    if destination.exists():
        wmean = np.load(destination)
    else:
        print(f"{name}: computing the Frechet mean...", flush=True)
        wmean = bary(diags, init=0)
        np.save(destination, wmean)
    return mean_mesr_vis, wmean, mean_image


def plot_means(mean_mesr_vis, wmean, mean_image, destination):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    measure_plot = axes[0].imshow(mean_mesr_vis.T, origin="lower",
                                  extent=[-ph.unit / 2, 1 - ph.unit / 2] * 2,
                                  cmap="hot_r", interpolation="nearest")
    fig.colorbar(measure_plot, ax=axes[0], shrink=0.7, label="Mass per grid cell")
    axes[1].scatter(wmean[:, 0], wmean[:, 1], s=30, color="tab:red")
    image_plot = axes[2].imshow(mean_image, origin="lower", extent=[0, CUTOFF, 0, CUTOFF],
                                cmap="hot_r", interpolation="nearest")
    fig.colorbar(image_plot, ax=axes[2], shrink=0.7, label="Image value")
    for ax, title in zip(axes, ["Mean persistence measure", "Frechet mean", "Mean persistence image"]):
        ax.set(xlim=(0, CUTOFF + ph.unit), ylim=(0, CUTOFF + ph.unit), xlabel="Birth", title=title)
        ax.set_aspect("equal")
    for ax in axes[:2]:
        ax.plot([0, CUTOFF], [0, CUTOFF], color="0.6", linewidth=0.7)
        ax.set_ylabel("Death")
    axes[2].set_ylabel("Persistence (death - birth)")
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", choices=list(DATASETS), nargs="+", default=list(DATASETS))
    parser.add_argument("--subsamples", type=int, help="Override the dataset's default number of subsamples.")
    parser.add_argument("--subsample-size", type=int, help="Override the dataset's default subsample size.")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/knot_lock")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in args.datasets:
        settings = DATASETS[name].copy()
        count = args.subsamples if args.subsamples is not None else settings["subsamples"]
        n = args.subsample_size if args.subsample_size is not None else settings["subsample_size"]
        vertices = PlyData.read(ROOT / "data" / settings["file"])["vertex"]
        points = ph.rescale_points(np.column_stack([vertices[axis] for axis in ("x", "y", "z")]))
        prefix = f"{name}_n{n}_B{count}"
        settings.update(subsamples=count, subsample_size=n, seed=SEED, cutoff=CUTOFF, sparse=0.3,
                        normalization="per_axis_unit_box", image_bandwidth=0.02, image_resolution=[50, 50],
                        versions={package: version(package) for package in ("numpy", "gudhi", "POT", "plyfile")})
        settings_path = args.output / f"{prefix}_settings.json"
        settings_path.write_text(json.dumps(settings, indent=2), encoding="utf8")
        means = compute_means(points, name, n, count, settings["min_persistence"], args.output)
        plot_means(*means, args.output / f"{prefix}.png")
        print(f"{name}: results saved to {args.output}", flush=True)


if __name__ == "__main__":
    main()
