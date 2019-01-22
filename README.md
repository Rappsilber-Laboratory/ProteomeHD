# About this repository
R scripts and data files related to the manuscript "The human proteome co-regulation map reveals functional relationships between proteins" by Kustatscher <i> et al </i> (2019). The article is available in bioRxiv {ENTER REF HERE}.

To reproduce the results of our paper download the files, set the R working directory to the download folder and execute the appropriate scripts. All scripts were tested in Ubuntu 18.04, but should generally be working on all operating systems. Below is a brief description of each file.


# R scripts
- <b> Coregulation_scores.R: </b> This script creates the co-regulation scores for protein pairs in ProteomeHD. It first calculates treeClust dissimilarities (using hyperparameters optimised with the tune_treeclust.R script) and turns them into similarities (i.e. 1 - dissimilarity). Using the WGCNA package, the script subsequently performs a sigmoid transformation to create an adjacency matrix (using hyperparameters optimised with the tune_sigmoid.R script). Finally, the script applies the topological overlap measure. The script outputs three files:
    (a) <i> treeClust_similarities.csv: </i> Contains the treeClust similarities, the adjacency matrix and the TOM values.
    (b) <i> coregulation_scores.csv: </i> Contains only the final pairwise coregulation scores (treeClust + sigmoid transformation + TOM). 
    (c) <i> ScaleFreeNess.png: </i> A plot showing that the resulting network is scale free.

- <b> tune_treeclust.R: </b> This script was used to perform a grid search to optimise treeClust / rpart hyperparameters (serule and cp) against true and false positives pairs annoated in the Reactome_TP_FP_10perc_subset_for_GS.csv file. The optimal values turned out to be cp = 0.105 and serule = 1.8, providing a ~10% improvement over default settings.

- <b> tune_sigmoid.R: </b> This script was used to perform a grid search to find the optimal parameters for WGCNA's topological overlap matrix (TOM). These parameters were mu and alpha, which relate to the sigmoid transformation taking place before calculating TOM. The optimal values turned out to be mu = 0.91 and alpha = 37. This step provided a further 10% improvement in performance. Note that this script is designed for execution on a server with 30 free cores.

# Data files
- <b> ProteomeHD_v1_1.7z: </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios. An exact copy of this file is Supplementary_Table_S1.csv {DOUBLE CHECK TABLE NAME} in the publication.

- <b> Reactome_TP_FP.7z: </b> This is a compressed csv file containing the "gold standard", i.e. a list of true and false positive protein pairs. It was generated as described in the article.

- <b> Reactome_TP_FP_10perc_subset_for_GS.7z: </b> This is a compressed csv file containing a subset of the Reactome gold standard (Reactome_TP_FP.7z). It was used to optimise treeClust hyperparameters in a grid search.





