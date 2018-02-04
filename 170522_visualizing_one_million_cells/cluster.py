import sys
import scanpy.api as sc

sc.settings.verbosity = 2  # show logging output
sc.settings.autoshow = False  # do not show figures automatically

path = sys.argv[1]  # read path from command line
n_jobs = int(sys.argv[2])  # number of jobs


def tsne_cluster(path, njobs):
    adata = sc.read_10x_h5(path)
    sc.pp.recipe_zheng17(adata)
    sc.tl.tsne(adata, n_jobs=n_jobs)
    sc.tl.louvain(adata)
    sc.pl.tsne(adata, color='louvain_groups', save=True, right_margin=2)
    adata.write('one_million.h5ad')

if __name__ == "__main__":
    one_milion(path, n_jobs)
