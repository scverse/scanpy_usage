import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import scanpy.api as sc

sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving

filename = sys.argv[1]  # read filename from command line

def basic_analysis(filename):
    adata = sc.read_10x_h5(filename)
    sc.pp.recipe_zheng17(adata)
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata)
    sc.tl.paga(adata)
    sc.tl.umap(adata)
    sc.tl.rank_genes_groups(adata, 'louvain')
    adata.write('./write/result.h5ad')
    # plotting
    sc.pl.paga(adata)
    sc.pl.umap(adata, color='louvain')
    sc.pl.rank_genes_groups(adata, save='.pdf')


if __name__ == "__main__":
    basic_analysis(filename)
