from memory_profiler import profile
import numpy as np
import scanpy.neighbors

np.random.seed(0)
X = np.random.rand(20000, 50)

@profile
def numpy():
    scanpy.neighbors.compute_neighbors_numpy(X, n_neighbors=10)


@profile
def umap():
    scanpy.neighbors.compute_neighbors_umap(X, n_neighbors=10)


# numpy()
umap()
