*First compiled: May 3, 2017*

## Scanpy computationally outperforms Cell Ranger

We compare Scanpy with *Cell Ranger* R kit, the 10x Genomics toolkit [(Zheng *et al.*, Nat. Comm. 2017)](https://dx.doi.org/10.1038/ncomms14049). The comparison shows that Scanpy requires about a factor of 10 less CPU time and less memory in crucial steps of the analysis. This enables analyzing 64000 cells without waiting times interactively on a regular laptop (MacBook Pro 13-inch, Early 2015, one 2,7 GHz Intel Core i5 processor with two cores, 16 GB RAM). Note that since recently (May 9, 2017), 10x Genomics has released a command-line based version of Cell Ranger, which shares a few of Scanpy's features, but is not available for interactive use as an API.

This Scanpy [notebook](zheng17_pbmc68k_cellranger_Py.ipynb) is benchmarked together with a [notebook](zheng17_pbmc68k_cellranger_R.ipynb) that reproduces the original Cell Ranger analysis. Both notebooks produce profiling information about CPU time and memory usage and yield exactly the same results.

The data used for this consists in 68,579 PBMC cells and is freely available [[page](https://support.10xgenomics.com/single-cell/datasets/fresh_68k_pbmc_donor_a)/[file](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz)].


#### Design of the comparison

Both the analysis of 68000 cells with [Scanpy (68000)](http://falexwolf.de/scanpy_usage/170503_zheng17/zheng17_pbmc64k_cellranger_Py_68000cells.html) and Cell Ranger [Cell Ranger (68000)](http://falexwolf.de/scanpy_usage/170503_zheng17/zheng17_pbmc64k_cellranger_R_68000cells.html) has been split into meaningful steps, in particular *preprocessing*, *further preprocessing/PCA* and *tSNE*. For Scanpy, we also profile the computation of *Diffusion Maps* and *Diffusion Pseudotime analysis*. Here, we also compile the files for lower cell numbers: [html](html).

#### Speedup

We obtain a speed up in the preprocessing routines of about a factor 10, in tSNE and PCA of about a factor 3 to 5. Particularly in standard preprocessing.

<img src="figs/cpu_time.png" height="250">
<img src="figs/cpu_time_ratio.png" height="250">

Scanpy offers *Diffusion Maps* and the *Diffusion Pseudotime Analysis* [Haghverdi *et al.*, Nat. Meth. (2016)](http://10.1038/nmeth.3971) for generating cell trajectories that are related to continuous processes as development or dose response.

<img src="figs/cpu_time_dpt.png" height="250">


#### Memory

The memory measurement here only concerns the memory at the end of a step, it is *not* the maximum memory during a process, where R tends to allocate even a lot more. Still, here, we already observe 10 times less memory use. In practice we observe a clearer memory of efficiency of Scanpy.

<img src="figs/memory.png" height="250">
<img src="figs/memory_ratio.png" height="250">
