*First compiled: May 5, 2017.*   
See the [notebook](seurat.ipynb).

### Scanpy and Seurat

Scanpy provides a number of Seurat's features ([Macosko *et al.*, Cell, 2015](http://dx.doi.org/10.1016/j.cell.2015.05.002)), but at significantly higher computationally efficiency.

[Here](seurat.ipynb), we reproduce most of Seurat's [guided clustering tutorial](http://satijalab.org/seurat/pbmc-tutorial.html) as compiled on March 30, 2017. The tutorial starts with preprocessing and ends with the identification of cell types through marker genes in clusters.

The data consists in 3k PBMCs from a Healthy Donor and is freely available from 10x (section *Chromium Demonstration v1 Chemistry* in [datasets](https://support.10xgenomics.com/single-cell/datasets)).
