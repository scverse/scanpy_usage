*First compiled: May 5, 2017.*   
See the [notebook](https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb).

### Benchmark: Scanpy versus Seurat

Scanpy provides a number of Seurat's features ([Satija *et al.*, Nat. Biotechnol., 2015](https://doi.org/10.1038/nbt.3192)), but at significantly higher computationally efficiency. [Here](https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb), we reproduce most of Seurat's [guided clustering tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html) as compiled on March 30, 2017. The tutorial starts with preprocessing and ends with the identification of cell types through marker genes of clusters. The data consists in *3k PBMCs from a Healthy Donor* and is freely available from 10x ([here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)). The profiling information for Seurat has been obtained within [seurat_R.ipynb](seurat_R.ipynb).

*Note:* The profiling information was obtained in June 2017 for Scanpy 0.2.1 and Seurat 1.4.0.4). In the meanwhile, both Scanpy and Seurat have become faster and the difference should not be as dramatic any more.

| | Scanpy | Seurat |
|----|-----|----|
| **preprocessing** | < 1 s | 14 s |
| **highly variable genes** | <img src="figures/scanpy_high_var_genes.png" height=200> | <img src="figures/seurat_high_var_genes.png" height=200> |
| **correction, regressing out** | 6 s | 129 s |
| **PCA** | < 1 s | 45 s |
|  | <img src="figures/scanpy_pca.png" height=200> | <img src="figures/seurat_pca.png" height=200> |
| **clustering** | 1.3 s | 65 s |
| **tSNE** | 6.5 s | 25 s |
|  | <img src="figures/scanpy_tsne.png" height=200> | <img src="figures/seurat_tsne.png" height=200> |
| **finding marker genes** | 0.8 s | 96 s |
|  | <img src="figures/scanpy_violin.png" height=200> | <img src="figures/seurat_violin.png" height=200> |
|  | <img src="figures/scanpy_diff_genes.png" height=200> | <img src="figures/seurat_diff_genes.png" height=200> |

