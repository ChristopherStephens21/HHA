#!/projects/b1217/Chris/condaenvs/SCENICPlusEnv/bin/python3
#Importing packages
print("Importing Packages")
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import SEACells
import os
import pickle
from scipy.io import mmread #Reading in mtx file 
from scipy.sparse import csr_matrix #convert from coo to csr matrix

#Directories
work_dir = "/projects/b1217/HHA/Matrix_ATAC_SEACells/"

#Loading in AnnDatas
print("Loading ATAC Object")
atac_meta_ad = sc.read(os.path.join(work_dir, "Matrix_Multiome_atac_meta_ad.h5ad"))
print(atac_meta_ad)
print("Loading RNA Object")
rna_meta_ad = sc.read(os.path.join(work_dir, "Matrix_Multiome_rna_meta_ad.h5ad"))
print(rna_meta_ad)

#-----Running Peak to Gene Linkage-----#
print("Running Peak to Gene Linkage")
gene_set = rna_meta_ad.var_names
gene_peak_cors = SEACells.genescores.get_gene_peak_correlations(atac_meta_ad, rna_meta_ad, 
                                           path_to_gtf='/projects/b1217/HHA/Dictys/Inputs/gene.gtf', 
                                           span=100000, 
                                           n_jobs=62,
                                           gene_set=gene_set)
print("Saving Results")
pickle.dump(gene_peak_cors, open(os.path.join(work_dir, 'gene_peak_cors.pkl'), 'wb'))
print("Done")
