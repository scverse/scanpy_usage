import numpy as np
import scanpy as sc
import re

filename_data = './data/pbmc3k_filtered_gene_bc_matrices/hg19/matrix.mtx'
filename_gene_names = './data/pbmc3k_filtered_gene_bc_matrices/hg19/genes.tsv'
adata = sc.read(filename_data).transpose()
adata.var_names = np.genfromtxt(filename_gene_names, dtype=str)[:, 1]

adata.smp['n_counts'] = np.sum(adata.X, axis=1).A1
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

mito_genes = np.array([name for name in adata.var_names
                       if bool(re.search("^MT-", name))])
# for each cell compute fraction of counts in mito genes vs. all genes
adata.smp['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as sample annotation to adata
adata.smp['n_counts'] = np.sum(adata.X, axis=1).A1
# some plots
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, show=True)

sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes', show=True)

adata.filter_smp(adata.smp['n_genes'] < 2500)
adata.filter_smp(adata.smp['percent_mito'] < 0.05)

sc.pp.normalize_per_cell(adata, scale_factor=1e4)
result = sc.pp.filter_genes_dispersion(adata.X, log=True, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.filter_genes_dispersion(result)

adata.filter_var(result['gene_filter'])
sc.pp.log1p(adata)

adata_corrected = sc.pp.regress_out(adata,
                                    smp_keys=['n_counts', 'percent_mito'],
                                    copy=True)
sc.pp.scale(adata_corrected, max_value=10)

sc.tl.pca(adata_corrected)
adata_corrected.smp['X_pca'] *= -1  # multiply by 1 for correspondence with R
sc.pl.pca_scatter(adata_corrected)

sc.pl.pca_ranking(adata_corrected)

sc.tl.tsne(adata_corrected, n_pcs=10)

sc.tl.dbscan(adata_corrected, basis='tsne', n_comps=2, min_samples=5)
sc.pl.tsne(adata_corrected, smp='dbscan_groups')

sc.tl.diffmap(adata)
sc.pl.diffmap(adata)
