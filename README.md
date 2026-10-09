# Approximating persistent homology via multiple subsampling

This repository contains implementations of algorithms and experiments in [Approximating Persistent Homology for Large Datasets](https://arxiv.org/abs/2204.09155).

Computing persistent diagrams for extremely large data sets is prohibitive. One way to bypass this difficulty is to draw subsamples from the large data set. In particular, let `X` be a large point cloud. We sample `B` subsets each consisting of `n` i.i.d. samples from `X`. We then compute the "averages" of different representations of persistent homology of subsampled sets, which, in principle, should approximate the true persistent homology the original data `X`.


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
| Gensim | 4.4.0 |

(*Gensim is needed only for `scripts/train_poincare.py`.*)

## Usage

The repository is organized as follows:

```text
.
├── ApproxPH.py                 # Multiple subsampling and persistence measure utilities
├── tutorial.ipynb              # Tutorial and large point cloud illustration
├── scripts/
│   ├── torus_rate.py           # Torus convergence rate experiment
│   ├── sphere_rate.py          # Sphere convergence rate experiment
│   ├── shape_clustering.py     # Clustering of bearing and motor shapes
│   ├── knot_lock.py            # Persistent homology approximation for Knot and Lock
│   ├── parameter_tuning.py     # Subsample size, count, and transport power tuning
│   ├── poincare_embedding.py   # Permutation tests for Poincaré embeddings
│   └── train_poincare.py       # Training Poincaré embeddings from WordNet relations
├── data/
│   ├── Bearing/               # Bearing point clouds
│   ├── Motor/                 # Motor point clouds
│   ├── knot.ply               # Knot point cloud
│   ├── grayloc.ply            # Lock point cloud
│   ├── parameter_tuning/      # Torus and sphere point clouds
│   └── poincare/              # Embeddings, word labels, relations, and source licenses
├── outputs/                   # Saved experiment results (partially, due to repo size)
│   ├── torus_rate/
│   ├── sphere_rate/
│   ├── shape_clustering/
│   ├── knot_lock/
│   ├── parameter_tuning/
│   └── poincare_embedding/
├── LICENSE                    
└── README.md
```

## Academic Use

If you find this useful for your research, you can use the following BibTex entry:

```bibtex
@misc{cao2026approximating,
  title={Approximating Persistent Homology for Large Datasets},
  author={Yueqi Cao and Anthea Monod},
  year={2026},
  eprint={2204.09155},
  archivePrefix={arXiv},
  primaryClass={stat.ML},
  url={https://arxiv.org/abs/2204.09155}
}
```


 
