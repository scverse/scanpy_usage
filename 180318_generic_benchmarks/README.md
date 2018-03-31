# Benchmarking nearest neighbor computations

In version 0.*, Scanpy's numpy-based neighbor-search used too much memory:

<img src="figures/mem-numpy_20K.svg" height=200>

Run [this](memory_over_time.py) to reproduce this figure. Using the approximate nearest neighbor search of UMAP [McInnes & Healy (2018)](https://arxiv.org/abs/1802.03426), this problem is resolved:

<img src="figures/mem-umap_20K.svg" height=200>

By this, also computation times have been much reduced for larger datasets:

<img src="figures/time-n_neighbors20.svg" height=200>

which can be produced with this [notebook](neighbors.ipynb).