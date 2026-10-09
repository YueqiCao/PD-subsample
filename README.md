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

## Usage

The repository is organized as follows:




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


 
