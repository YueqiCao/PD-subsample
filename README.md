# Approximating the persistent homology of large datsets using mean persistence measures and Fréchet means

This repository contains codes to implement the numerical experiments in [Approximating Persistent Homology for Large Datasets](https://arxiv.org/abs/2204.09155).

Computing persistent diagrams for extremely large data sets is prohibitive. One way to bypass this difficulty is to draw subsamples from the large data set. Then the mean of persistence diagrams of subsamples *approximates* the persistence diagram of the original data. 

In our work, we propose to use mean persistence measures and Fréchet means of persistence diagrams to estimate the persistent homology of large data sets using subsampling. In particular, let `X` be a large point cloud with a predefined probability distribution satisfying some standard assumptions. We sample `B` subsets each consisting of `n` i.i.d. samples from `X`. We then compute the mean persistence measure and Fréchet mean which can be regarded as two types of averages of persistence diagrams of subsample sets. As `B` increases, we expect the empirical means to converge to their corresponding population means, and as `n` increases, we expect the population mean converge to the true persistence diagram. It turns out that `B` controls the variance error, and `n` controls the bias error. The approximation error between mean diagram/measure and the true persistence diagram is bounded by quantities involving `B` and `n`.

## Usage

The repository is organized as follows:

- `ApproxPH.py`: a module containing necessary functions to compute PH etc.;
- `multiple-subsampling-persistence.ipynb`: a notebook shows how to compute/use mean diagram/measure;
- `data`: a file containing input data;
- `outputs`: a file containing output results.

## Academic Use

Please cite

> @article{cao2022approximating,  
>  author = {Cao, Yueqi and Monod, Anthea},  
>  title = {Approximating Persistent Homology for Large Datasets},  
>  journal = {arXiv preprint arXiv:2204.09155},  
>  year = {2022}  
>  }


 
