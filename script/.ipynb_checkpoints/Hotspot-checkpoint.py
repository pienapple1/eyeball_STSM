###Author Jinpei Lin
###Times 2024.07

import os
import glob
import sys
import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import pickle
import loompy as lp
import os
import glob
from PIL import Image
from PIL import ImageDraw
import matplotlib.pyplot as plt
key = "SS200000709TR_D1_1"
path1 = "/00df_coor/"
path2 = "/01df_ex/"
outdir1 = "/result/"
pos = pd.read_csv(path1+key+"_cellpos.csv", index_col=0)
pos.index=pos.index.astype("int64")
lf = lp.connect(path2+key+"_expmtx_nvg5000.loom", mode='r', validate=False)
counts = pd.DataFrame(lf[:, :],
                      index=lf.ra["Gene"],
                      columns=lf.ca["CellID"],
                      dtype=np.int64)
lf.close()
counts.columns=counts.columns.astype("int64")
# Align the indices
outdir = outdir1+key+"/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
pos = pos.loc[counts.columns, :]
counts = counts.loc[:, pos.index]
barcodes = pos.index.values
# Swap position axes
# We swap x'=y and y'=-x to match the slides in the paper
pos = pd.DataFrame(
    {
        'X': pos.X,
        'Y': pos.Y,
    }, index=pos.index
)
num_umi = counts.sum(axis=0)
hs = hotspot.Hotspot(counts, model='normal', latent=pos, umi_counts=num_umi)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=30,
)
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=20)
#select the genes with significant spatial autocorrelation
hs_genes = hs_results.index[hs_results.FDR < 0.05]
# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=20)
modules = hs.create_modules(
    min_gene_threshold=20, core_only=False, fdr_threshold=0.05
)
plt.rcParams['figure.figsize'] = (15.0, 16.0)
hs.plot_local_correlations()
plt.savefig(''.join([outdir, "/", key, "_module_number.pdf"]), dpi = 600)
plt.close()
results = hs.results.join(hs.modules)
results.to_csv(''.join([outdir, "/", key, "_Cluster.csv"]))
module_scores = hs.calculate_module_scores()
module_scores.to_csv(''.join([outdir, "/", key, "_ModuleScore.csv"]))
for module in range(1, hs.modules.max()+1):
    scores = hs.module_scores[module]

    vmin = np.percentile(scores, 1)
    vmax = np.percentile(scores, 99)
    plt.scatter(x=hs.latent.iloc[:, 0],
                y=hs.latent.iloc[:, 1],
                s=8,
                c=scores,
                vmin=vmin,
                vmax=vmax,
                edgecolors='none'
                )
    axes = plt.gca()
    for sp in axes.spines.values():
        sp.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    plt.title('Module {}'.format(module))
    plt.savefig(f'{outdir}Module{module}.pdf')
    plt.close()