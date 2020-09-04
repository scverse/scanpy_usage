from joblib import parallel_backend
import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import scanpy as sc

sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=400)  # set sufficiently high resolution for saving

filename = sys.argv[1]  # read filename from command line

def basic_analysis(filename):
    adata = sc.read_10x_h5(filename)
    sc.pp.recipe_zheng17(adata)
    with parallel_backend('threading', n_jobs=16):
        sc.pp.neighbors(adata)
    sc.tl.louvain(adata)
    sc.tl.umap(adata, random_state=None)
    sc.tl.rank_genes_groups(adata, 'louvain')
    adata.write('./write/result_opt.h5ad')
    # plotting
    sc.pl.umap(adata, color='louvain', save='_opt.png')
    sc.pl.rank_genes_groups(adata, save='_opt.pdf')


if __name__ == "__main__":
    basic_analysis(filename)
