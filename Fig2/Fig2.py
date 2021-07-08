##Fig2e (Violin plots)

import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, dpi_save=300,facecolor='white',format="svg")
sc.logging.print_header()
scanpy==1.6.1 anndata==0.7.5 umap==0.4.6 numpy==1.19.4 scipy==1.5.2 pandas==1.1.4 scikit-learn==0.23.2 statsmodels==0.11.1 python-igraph==0.8.2 leidenalg==0.8.0

#load the converted h5ad (we created this object in Fig2_SuppFig2.R)
adata = sc.read_h5ad("~/Seurat.obj.h5ad")

marker_genes = ["CD4","CD8A","CD8B","CD3E","SELL","CCR7","TCF7","IL7R","PDCD1", "LAG3", "HAVCR2", "KLRG1","TIGIT", "CD244", "CD160", "BTLA","CTLA4","ENTPD1", "ID2", "GZMB", "GZMK", "GZMH", "GZMA","GZMM", "PRF1", "GNLY", "NKG7","TNF", "IFNG","TOX", "MKI67", "FOXP3"]

sc.pl.stacked_violin(adata, marker_genes, groupby='new_clustering', rotation=90,save='.svg')

##Fig2g (circa plots)

%matplotlib inline
import sys 
import numpy as np
sys.path.append("~/Software/pyCircos/")
from pycircos import * 

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, dpi_save=300,facecolor='white',format="svg")
scanpy==1.6.1 anndata==0.7.5 umap==0.4.6 numpy==1.19.4 scipy==1.5.2 pandas==1.1.4 scikit-learn==0.23.2 statsmodels==0.11.1 python-igraph==0.8.2 leidenalg==0.8.0
In [2]:

#Exec external commands with the definition of the circus plot (we created this commands in Fig2_SuppFig2.R)
exec(open("~/pycircus_commands_alpha.txt").read())

gcircle.save()

plt.savefig("circus_plot_alpha.pdf",bbox_inches="tight")

#Exec external commands with the definition of the circus plot (we created this commands in Fig2_SuppFig2.R)
exec(open("~/pycircus_commands_beta.txt").read())

gcircle.save()

plt.savefig("circus_plot_beta.pdf",bbox_inches="tight")