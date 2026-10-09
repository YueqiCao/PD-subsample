# Approximating persistent homology via multiple subsampling

This repository contains codes to implement the numerical experiments in [Approximating Persistent Homology for Large Datasets](https://arxiv.org/abs/2204.09155).

Computing persistent diagrams for extremely large data sets is prohibitive. One way to bypass this difficulty is to draw subsamples from the large data set. Then the mean of persistence diagrams of subsamples *approximates* the persistence diagram of the original data. 

In our work, we propose to use mean persistence measures and Fréchet means of persistence diagrams to estimate the persistent homology of large data sets using subsampling. In particular, let `X` be a large point cloud with a predefined probability distribution satisfying some standard assumptions. We sample `B` subsets each consisting of `n` i.i.d. samples from `X`. We then compute the mean persistence measure and Fréchet mean which can be regarded as two types of averages of persistence diagrams of subsample sets. As `B` increases, we expect the empirical means to converge to their corresponding population means, and as `n` increases, we expect the population mean converge to the true persistence diagram. It turns out that `B` controls the variance error, and `n` controls the bias error. The approximation error between mean diagram/measure and the true persistence diagram is bounded by quantities involving `B` and `n`.

## Usage

Activate the `pd-multiple-subsampling` conda environment and open `tutorial.ipynb`
from the repository directory. Select the **Python (PD multiple subsampling)**
kernel. The tutorial has two sections: three subsampling averages on an annulus,
and persistent-homology approximation of a large point cloud. It includes the
mean persistence measure (MPM), a local Fréchet mean (FM), and the mean persistence
image. Its settings are deliberately smaller than those of the paper experiments.

- `ApproxPH.py`: sampling, persistent homology, persistence measures and transport utilities.
- `tutorial.ipynb`: the two-section tutorial (formerly `multiple-subsampling-persistence.ipynb`).
- `scripts/torus_rate.py`: the torus rate experiment and fitted curves.
- `scripts/shape_clustering.py`: Bearing/Motor clustering with MPM or FM.
- `data/`: input point clouds.
- `outputs/`: saved results from the original notebook.

### Torus rate experiment

To fit the existing notebook results without recomputing persistent homology:

```bash
python scripts/torus_rate.py --from-measures outputs
```

The saved `mean_mesr_nb0.npy` through `mean_mesr_nb14.npy` must correspond to
`n=400,600,...,3200` and `B=0.1*n`, as in the former notebook. These older files
contain no settings metadata; use `--sizes` and `--b-ratio` if their generating
schedule differs. Use the regenerated measures from the corrected grid indexing,
rather than mixing them with files produced by the old indexing bug.

To recompute the subsampling experiment using the bundled 50,000-point torus and
its saved reference diagram:

```bash
python scripts/torus_rate.py --output outputs/experiments/torus_fresh
```

This retains the original radii 0.8 and 0.3, sparse Rips parameter 0.3, filtration
cutoff 0.9, subsample persistence threshold 0.01, and `B=0.1*n` schedule. Computing
all new diagrams can take several hours. The script saves only the mean measures
in a compressed `measures.npz`, checkpointed after each sample size. It does not
save subsampled clouds or recompute the expensive full-cloud reference diagram.
For the paper's longer range, pass `--sizes 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800`.
For a quick code check, use `--sizes 40 60 80 100 --b-ratio 0.025` and a new output
directory; such small settings are not a rate experiment.

Both `p=3` and `p=8` are evaluated. The fit is `a*n**(-b)+c`, with scaled inputs,
nonnegative amplitude and offset, and a free exponent (`0.01 <= b <= 5`). It does
not impose `b=0.5`. `rate.png` and `summary.json` report two distinct quantities:

- The **original Sinkhorn transport cost with regularization 1**, retained for
  comparison with the earlier experiment. This is a regularized proxy, not an
  exact Wasserstein distance; even self-comparison can have positive cost.
- The **exact powered partial-transport cost**, using the same discretization
  without the artificial numerical background mass. Its empirical rate can differ.

The fit is to powered costs, not their pth roots. An upper bound containing an
`n**(-0.5)` term does not require the observed loss to have exactly that exponent.
The curves are single realizations, not averages over independent repetitions.
Replaying the current 15 saved means gave:

| Cost | p | Fitted exponent −b | R² |
| --- | --- | --- | --- |
| Sinkhorn, regularization 1 | 3 | −0.462 | 0.9997 |
| Sinkhorn, regularization 1 | 8 | −0.447 | 0.9793 |
| Exact partial transport | 3 | −0.913 | 0.9994 |
| Exact partial transport | 8 | −3.822 | 0.99998 |

Replot a newly saved bundle with `--from-measures path/to/measures.npz` and a fresh
`--output` directory. The sizes and counts are read from that bundle.

### Shape clustering

```bash
python scripts/shape_clustering.py
```

The default computes MPM for the **74 Bearing and 52 Motor** clouds included here.
It uses 15 subsets of 1,500 points per shape, instead of a percentage of each
cloud's size. This gives equal subsampling budgets and avoids large intermediate
files. Each shape is centered and scaled by its largest coordinate span, preserving
aspect ratios. H1 uses sparse Rips parameter 0.3, cutoff 0.4 and persistence
threshold 0.03; surviving features are capped at the cutoff.

For clustering, the 50-unit persistence grid rounds births down and deaths up,
retaining positive-lifetime features without adding background mass. Exact
partial-transport costs are square-rooted to obtain 2-Wasserstein distances with
an L-infinity ground metric. UMAP receives this **precomputed distance matrix**,
with 15 neighbors, `min_dist=0.1` and seed 30. DBSCAN uses `eps=1.0` and
`min_samples=10`. The former `eps=3` merged clusters after changing the distance
pipeline; these settings were calibrated for the revised pipeline on this dataset.
The torus script retains the former floor/floor grid for comparison with its saved
means, so the two experiments explicitly use different quantization conventions.

The script writes compressed diagram banks, distance matrices, embeddings, cluster
labels, figures, scores and settings. It saves no sampled clouds or Rips complexes.
NumPy's sampling seed is recorded; sparse-Rips construction and numerical routines
may still vary across library builds. To repeat with another subsampling seed, use
`--seed 20211122 --output outputs/experiments/shapes_seed22`.

The FM comparison is available separately because the iterative barycenter can be
slow. To reuse the same diagrams after the default MPM run:

```bash
python scripts/shape_clustering.py --methods fm --diagrams outputs/experiments/shape_clustering/diagrams --output outputs/experiments/shape_fm
```

Alternatively, `--methods mpm fm` computes both from a shared bank in a new run.
GUDHI's FM is a local barycenter with an L2 ground metric; pairwise FM distances
retain the notebook's L-infinity ground metric. FM also retains the notebook's
10 UMAP neighbors, UMAP seed 20, and DBSCAN `eps=2`, with `min_samples=10`.
`--neighbors`, `--umap-seed` and `--eps` override these defaults for either method.
The FM path was checked on a small subset; a complete FM experiment was not rerun.

As a regression check, the new MPM implementation was applied to three existing
diagram banks with these settings, after verifying that all input point clouds
matched this repository. These checks recomputed distances and clustering, not
the persistent homology of all 126 shapes:

| Subsampling seed | Adjusted Rand index | Clusters | Noise points | Purity of assigned points | Correct fraction, all points |
| --- | --- | --- | --- | --- | --- |
| 20211121 | 0.5935 | 3 | 2 | 96.77% | 95.24% |
| 20211122 | 0.5825 | 3 | 3 | 96.75% | 94.44% |
| 20211123 | 0.4942 | 4 | 2 | 97.58% | 96.03% |

Purity assigns each cluster its majority class; several clusters can have the same
class. The last column counts noise as incorrect. These are clustering diagnostics,
not held-out classification accuracy or a claim of exactly two recovered clusters.
Class labels are used only for evaluation and plotting. The cached banks used for
this check are not bundled; a fresh run creates its own compact banks.

For a quick end-to-end check of both methods:

```bash
python scripts/shape_clustering.py --methods mpm fm --limit-per-class 3 --subsample-size 300 --subsamples 3 --neighbors 3 --min-samples 2 --output outputs/experiments/shape_smoke
```

Both scripts accept `--help`. They require a **new output directory** to prevent
overwriting a previous run. Generated files under `outputs/experiments/` are
ignored by Git; the original saved results are preserved.

## Required Libraries


| Package | Version |
| --- | --- |
| Python | 3.12.15 | 
| NumPy | 2.5.3 | 
| SciPy | 1.18.1 | 
| Matplotlib | 3.11.2 | 
| GUDHI | 3.12.0 | 
| POT | 0.9.7.post1 | 
| plyfile | 1.1 | 
| scikit-learn | 1.9.1 | 
| umap-learn | 0.5.12 | 


## Academic Use

Please cite

```bibtex
@misc{cao2026approximatingpersistenthomologylarge,
  title={Approximating Persistent Homology for Large Datasets},
  author={Yueqi Cao and Anthea Monod},
  year={2026},
  eprint={2204.09155},
  archivePrefix={arXiv},
  primaryClass={stat.ML},
  url={https://arxiv.org/abs/2204.09155}
}
```


 
