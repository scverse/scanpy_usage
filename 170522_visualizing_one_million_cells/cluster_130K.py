import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import scanpy.api as sc

sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=400)  # set sufficiently high resolution for saving

filename = sys.argv[1]  # read filename from command line

def basic_analysis(filename):
    adata = sc.read_10x_h5(filename)
    sc.pp.subsample(adata, fraction=0.1)
    sc.pp.recipe_zheng17(adata)
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata)
    sc.tl.umap(adata)
    sc.tl.rank_genes_groups(adata, 'louvain')
    adata.write('./write/result_130K.h5ad')
    # plotting
    sc.pl.umap(adata, color='louvain', save='_130K.png')
    sc.pl.rank_genes_groups(adata, save='_130K.pdf')


if __name__ == "__main__":
    basic_analysis(filename)
