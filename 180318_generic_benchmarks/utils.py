import scanpy.neighbors
from sklearn.metrics import pairwise_distances


def scanpy0(X, n_neighbors):
    scanpy.neighbors.compute_neighbors_numpy(X, n_neighbors)


def scanpy1(X, n_neighbors):
    if X.shape[0] < 8192:
        scanpy.neighbors.compute_neighbors_numpy(X, n_neighbors)
    else:
        scanpy.neighbors.compute_neighbors_umap(X, n_neighbors)
