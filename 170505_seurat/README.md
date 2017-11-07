*First compiled: May 5, 2017.*   
See the [notebook](seurat.ipynb).

### Scanpy and Seurat

Scanpy provides a number of Seurat's features ([Satija *et al.*, Nat. Biotechnol., 2015](https://doi.org/10.1038/nbt.3192)), but at significantly higher computationally efficiency.

[Here](seurat.ipynb), we reproduce most of Seurat's [guided clustering tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html) as compiled on March 30, 2017. The tutorial starts with preprocessing and ends with the identification of cell types through marker genes in clusters. The data consists in *3k PBMCs from a Healthy Donor* and is freely available from 10x ([here](https://support.10xgenomics.com/single-cell/datasets/pbmc3k), section *Chromium Demonstration (v1 Chemistry)* in [datasets](https://support.10xgenomics.com/single-cell/datasets)).

| | Scanpy | Seurat |
|----|-----|----|
| highly variable genes | <img src="figures/scanpy_high_var_genes.png" height=200> | <img src="figures/seurat_high_var_genes.png" height=200> |
| PCA | <img src="figures/scanpy_pca.png" height=200> | <img src="figures/seurat_pca.png" height=200> |
| tSNE | <img src="figures/scanpy_tsne.png" height=200> | <img src="figures/seurat_tsne.png" height=200> |
| violin | <img src="figures/scanpy_violin.png" height=200> | <img src="figures/seurat_violin.png" height=200> |