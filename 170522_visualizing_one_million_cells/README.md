*Updated: February 4, 2018. First compiled: May 22, 2017.*   
Thanks to M. Lotfollahi for rerunning the computations!

## Visualizing and Clustering 1.3M neurons

This uses the 1.3M neurons [dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1M_neurons) from 10x Genomics. 

If you have at least 30 GB of memory, run [*cluster.py*](cluster.py) to produce the following result - this will take much less CPU time in the future:

<img src="figures/tsne.png" height=300px>

For simply plotting precomputed results on your laptop, use [*plot.ipynb*](plot.ipynb), which does not require much memory.

For quick tests, see [*cluster_subsampled.ipynb*](cluster_subsampled.ipynb) for a notebook that processes the subsampled dataset.

