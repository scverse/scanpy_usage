
## Cell-Cycle Scoring and Regression

This is a _Scanpy_ demo that shows how to regress cell cycle effect, following the approach showed in [Seurat's vignette](http://satijalab.org/seurat/cell_cycle_vignette.html#assign-cell-cycle-scores). As for the R example, toy dataset consists of murine hematopoietic progenitors from [Nestorowa et al., Blood 2016](https://doi.org/10.1182/blood-2016-05-716480).
Note that this notebook can only executed with scanpy version higher than 0.43


```python
import scanpy.api as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_version_and_date()
# we will soon provide an update with more recent dependencies
sc.logging.print_versions_dependencies_numerics()

```

    /Library/Frameworks/Python.framework/Versions/3.5/lib/python3.5/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
      from ._conv import register_converters as _register_converters


    Running Scanpy 0.4.3+0.gbbe5f9a.dirty on 2018-02-09 21:26.
    Dependencies: numpy==1.14.0 scipy==1.0.0 pandas==0.22.0 scikit-learn==0.19.1 statsmodels==0.8.0 python-igraph==0.7.1 louvain==0.6.1+5.g9c94269 


Load data


```python
path = "data/"
adata = sc.read_csv(path + "nestorawa_forcellcycle_expressionMatrix.txt", delimiter='\t').transpose()
```

Load cell cycle genes defined in [Tirosh et al, 2015](https://doi.org/10.1126/science.aad0501). It is a list of 97 genes, represented by their gene symbol. The list here is for humans, in case of alternate organism, a list of ortologues should be compiled.
There are major differences in the way _Scanpy_ and _Seurat_ manage data, in particular we need to filter out cell cycle genes that are not present in our dataset to avoid errors.


```python
cell_cycle_genes = [x.strip() for x in open(path + 'regev_lab_cell_cycle_genes.txt')]

```

Here we define two lists, genes associated to the S phase and genes associated to the G2M phase


```python
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes  if x in adata.var_names]
```

Standard filters applied. Note that we do not extract variable genes and work on the whole dataset, instead. This is because, for this demo, almost 70 cell cycle genes would not be scored as variable. Cell cycle scoring on ~20 genes is ineffective. 


```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.1, max_mean=8, min_disp=1)
#sc.pl.filter_genes_dispersion(filter_result)
```

    normalizing by total count per cell
        finished (0:00:00.245): normalized adata.X and added
        'n_counts', counts per cell before normalization (adata.obs)


Log-transformation of data and scaling should always be performed before scoring


```python
sc.pp.log1p(adata)
sc.pp.scale(adata)
```

We here perform cell cycle scoring. The function is actually a wrapper to `sc.tl.score_gene_list`, which is launched twice, to score separately S and G2M phases. Both `sc.tl.score_gene_list` and `sc.tl.score_cell_cycle_genes` are a port from _Seurat_ and are supposed to work in a very similar way. 
To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list.
Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.


```python
sc.tl.score_cell_cycle_genes(adata, s_genes=s_genes, g2m_genes=g2m_genes)
```

    Calculating cell cycle scores
    Adding score S_score
    Adding score G2M_score


Here comes another difference from _Seurat_. The R package stores raw data, scaled data and variable genes information in separate slots, _Scanpy_ instead keeps only one snapshot of the data. This implies that PCA is always calculated on the entire dataset. In order to calculate PCA reduction using only a subset of genes (like `cell_cycle_genes`), a trick should be used.
Basically we create a dummy object to store information of PCA projection, which is then reincorporated into original dataset.


```python
X = sc.tl.pca(adata[:, cell_cycle_genes], copy = True) #Project using a subset of genes
adata.obsm['X_pca'] = X.obsm['X_pca'] #Substitute the object
sc.pl.pca_scatter(adata, color='Phase') #Plot
```

    ... storing Phase as categorical type
        access categories as adata.obs['Phase'].cat.categories



![png](figures/output_15_1.png)


As in the original vignette, cells can be easily separated by their cell cycle status when cell cycle genes are used.
Now we can regress out both S score and G2M score.


```python
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
sc.pp.scale(adata)
```

    regressing out ['S_score', 'G2M_score']
    finished (0:00:50.822)
    --> after `sc.pp.regress_out`, consider rescaling the adata using `sc.pp.scale`


Finally, we reproject dataset using cell cycle genes again. Since we regressed the scores, no effect of cell cycle is now evident.


```python
X = sc.tl.pca(adata[:, cell_cycle_genes], copy = True)
adata.obsm['X_pca'] = X.obsm['X_pca']
sc.pl.pca_scatter(adata, color='Phase')
```


![png](figures/output_19_0.png)

