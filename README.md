# About this repository
R scripts and data related to the manuscript "Co-regulation map of the human proteome enables identification of protein functions" by Kustatscher et al (https://www.nature.com/articles/s41587-019-0298-5)

To reproduce the results of our paper download the files, set the R working directory to the download folder and execute the appropriate scripts. All scripts were tested in Ubuntu 18.04, but should generally be working on all operating systems. Below is a brief description of each file.


# R scripts
- <b> Coregulation_scores.R: </b> This script creates the co-regulation scores for protein pairs in ProteomeHD. It first calculates treeClust dissimilarities (using hyperparameters optimised with the tune_treeclust.R script) and turns them into similarities (i.e. 1 - dissimilarity). Using the WGCNA package, the script subsequently performs a sigmoid transformation to create an adjacency matrix (using hyperparameters optimised with the tune_sigmoid.R script). Finally, the script applies the topological overlap measure, thereby creating the "co-regulation scores". The script outputs three files:
    (a) <i> treeClust_similarities.csv: </i> Contains the treeClust similarities, the adjacency matrix and the TOM-transformed values.
    (b) <i> coregulation_scores.csv: </i> Contains only the final pairwise coregulation scores (treeClust + sigmoid transformation + TOM). 
    (c) <i> ScaleFreeNess.png: </i> A plot showing that the resulting network is scale free.

- <b> Fig1_coreg_proteins.R: </b> This script handles multiple tasks related to Figure 1 of the manuscript. It produces a plot that shows the distribution of the co-regulation scores. It produces a plot showing that co-regulated protein pairs are enriched for proteins from the same complex and subcellular localisation and for enzymes catalyzing consecutive metabolic reactions. Membership of protein complexes and metabolic pathways is based on Reactome annotation (the relavant input file, homo_sapiens.interactions.txt, can be downloaded from the Reactome website). A file containing all reviewed Uniprot proteins together with their subcellular localisation was downloaded from the Uniprot website (<i>uniprot-all.tab</i>) and is available in compressed form in this repository. This script also creates the <i>coregulated_protein_pairs.csv</i> file, which contains all co-regulated protein pairs in the form “Protein X” - “Coregulated with protein X”. This file is part of the manuscript as Supplementary Table 3. Finally, the script produces a barplot showing how many co-regulation partners proteins tend to have (in bins of 0, 1-5, 6-20 etc).

- <b> topGO_enrichment_analysis_BP.R / _CC.R / _MF.R: </b> Three almost identical scripts calculating the enrichment of gene ontology terms among groups of co-regulated proteins. The analysis uses the topGO R package. The difference between the three scripts is that they analyse GO Biological Processes (BP), Cellular Component (CC) and Molecular Function (MF). Each script takes the corresponding input file (<i>GO_associations_BP.tsv, GO_associations_CC.tsv or GO_associations_MF.tsv</i>), which are available in this respository in compressed form. These files contain GO associations downloaded from QuickGO using taxon = human, qualifier = NULL and aspect = Biological process (or aspect = Molecular function, or aspect = Cellular component). 

- <b> Coregulation_score_vs_corr_part1.R / _part2.R: </b> These two scripts compare the performance of treeClust dissimilarities with three correlation metrics (Spearman's rho, Pearson correlation and biweight midcorrelation). The first script calculates the correlation coefficients for ProteomeHD and saves them in a large, temporary output file called <i> ProHD_correlations.csv </i>. The second script performs a precision - recall analysis using the Reactome gold standard, taking the correlation coefficients from <i> ProHD_correlations.csv </i> and the treeClust and treeClust + TOM similarities from <i> treeClust_similarities.csv: </i> (produced by the Coregulation_scores.R script). 

- <b> tSNE_map.R: </b> This script produces the tSNE map of ProteomeHD. It takes the coregulation_scores.csv file as input, which is itself created by the Coregulation_scores.R script (also available as download elsewhere, see see below). Annotation of the map is based on the "tSNE_map_annotation.csv" file, which contains Uniprot-based organelle annotation as well as manual annotation of zoom regions. The script outputs <i> tSNE_coordinates_and_annotation.csv </i>, which contains the tSNE coordinates and the annotations together, and which is also saved as Supplementary Table 4.

- <b> Orthogonal_Methods.R: </b> Related to Figure 2. This script compares protein - protein associations detected by co-regulation with those found in STRING, BioGRID, IntAct and BioPlex 2.0. It uses downloaded files from these databases as input, which can be found on the respective websites. These include <i> 9606.protein.links.detailed.v10.5.txt</i> and <i> BIOGRID-ALL-3.4.152.tab2.txt </i>. The IntAct_interactions.csv file is provided here, because it is a pre-processed version of the original data from IntAct (pre-processed only to reduce file size by removing unnecessary information). The file nature22366-s2_Bioplex2.xlsx is Table S1 from the BioPlex 2.0 paper (Huttlin et al, Nature, 2017). Finally, the script uses <i> SwissProt_IDs_HoSa_Nov_18.tab </i>, downloaded from Uniprot on 19.11.2018 and containing their up-to-date IDs for removing old and unused protein IDs. This file is available here.

- <b> Uncharacterised_proteins.R: </b> Related to Figure 2. This script looks at uncharacterised and disease - related proteins covered by the co-regulation map and analyses their connectivity and molecular weight distribution. It uses the following input files (see below for availability): <i> coregulation_scores.csv, SwissProt_HS_AnnotScore_Mass.tab, Cancer_gene_census_Uniprot.csv, curated_gene_disease_associations.tsv, mapa_geneid_4_uniprot_crossref.tsv, nature22366-s2_Bioplex2.xlsx</i>.

- <b> Protein_examples_zoom.R: </b> Related to Figure 2. This script creates small versions of the co-regulation map to show how it may be used to annotate uncharacterised proteins (e.g. TMEM256). It uses the input files <i> coregulation_scores.csv, tSNE_coordinates_and_annotation.csv, ProteomeHD_v1_1.csv, QuickGO_mito_IMM.tsv</i>.

- <b> Protein_vs_mRNA_part_1.R: </b> This script compares the performance of treeClust on transcriptomics and proteomics measurements for the same 59 cell lines, producing a precision - recall curve. These measurements for the same set of genes are available as <i> BattleSILAC_PickrellRPKM.csv</i>, which we have generated previously (Kustatscher et al, MSB, 2017).

- <b> Protein_vs_mRNA_part_2.R: </b> This script compares treeClust co-regulation scores based on ProteomeHD with the transcriptomics-based coexpression scores from STRING. Please note that the file "RNA_STRING_Uniprot.csv", which contains Pearson correlation coefficients based on microarrays and RNA-sequencing, is not available. This is because these data were kindly provided to us by the STRING consortium specifically for the purpose of this analysis. Note that these Pearson correlation coefficients are distinct from publically available STRING's coexpression score, which is calibrated against the KEGG database.

- <b> ProHD_data_drop_impact.R </b> This script performs PR analyses after randomly removing data points from ProteomeHD, showing that the co-regulation analysis has not reached peak performance yet. This script should be run on a multi-core server, otherwise it will take a very long time to complete.

- <b> ProHD_stats.R </b> This script outputs basic protein-centric statistics about ProteomeHD, e.g. how many proteins were detected in how many experiments etc. 

- <b> ProHD_peptide_stats.R </b> This script outputs peptide-related stats (number of peptides per protein etc) for ProteomeHD, including microproteins. It uses the MaxQuant proteinGroups and evidence files, which are available through PRIDE (see below).

- <b> MP_connectivity.R </b> This script calculates the number of interaction partners of microproteins, as opposed to larger proteins, in ProteomeHD, STRING and BioGRID. See below for a description of the input files that it uses.

- <b> Network_plot.R </b> This script plots the co-regulation matrix as a "correlation network" with six different network layouts.
    
- <b> tune_treeclust.R: </b> This script was used to perform a grid search to optimise treeClust / rpart hyperparameters (serule and cp) against true and false positives pairs annoated in the Reactome_TP_FP_10perc_subset_for_GS.csv file. The optimal values turned out to be cp = 0.105 and serule = 1.8, providing a ~10% improvement over default settings.

- <b> tune_sigmoid.R: </b> This script was used to perform a grid search to find the optimal parameters for WGCNA's topological overlap matrix (TOM). These parameters were mu and alpha, which relate to the sigmoid transformation taking place before calculating TOM. The optimal values turned out to be mu = 0.91 and alpha = 37. This step provided a further 10% improvement in performance. Note that this script is designed for execution on a server with 30 free cores.

# Data files
- <b> ProteomeHD_v1_1.7z: </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios. An exact copy of this file is Supplementary_Table_S1.csv in the publication.

- <b> Reactome_TP_FP.7z: </b> This is a compressed csv file containing the "gold standard", i.e. a list of true and false positive protein pairs. It was generated as described in the article. Note that these pairs are already sorted row-wise in alphabetical order (Protein_1 > Protein_2), for easy merging with other tables.

- <b> Reactome_TP_FP_10perc_subset_for_GS.7z: </b> This is a compressed csv file containing a subset of the Reactome gold standard (Reactome_TP_FP.7z). It was used to optimise treeClust hyperparameters in a grid search.

- <b> uniprot-all.7z </b> File containing all reviewed Uniprot proteins together with their subcellular location.

- <b> GO_associations_BP.7z / _CC.7z / _MF.7z </b> These files contain GO associations downloaded from QuickGO using taxon = human, qualifier = NULL and aspect = Biological process (or aspect = Molecular function, or aspect = Cellular component).

- <b> tSNE_map_annotation.csv: </b> Used for annotation of the tSNE map by the tSNE_map.R script. Contains Uniprot-based organelle annotation as well as manual annotation of zoom regions.

- <b> IntAct_interactions.7z </b> To reduce file size, the table downloaded from IntAct was pre-filtered to exlude interactions from non-human species, interactors that were not proteins, interactions without Uniprot ID as well as self-interactions and duplicates. Finally, annotation strings were simplified, e.g. by reducing IDs from "uniprotkb:P06748" to "P06748".

- <b> SwissProt_IDs_HoSa_Nov_18.tab </b> Up-to-date Homo Sapiens Uniprot IDs (19.11.2018) for removing old and unused protein IDs from input files.

- <b> SwissProt_HS_AnnotScore_Mass.tab </b> Contains reviewed human Uniprot IDs together with their Uniprot annotation score and molecular weight.

- <b> Cancer_gene_census_Uniprot.csv </b> A file containing annotation from the cancer gene census (COSMIC v81). Originally downloaded from http://cancer.sanger.ac.uk/census. Entrez Gene IDs were converted to reviewed Homo Sapiens Uniprot IDs, ambiguous mappings were removed.

- <b> QuickGO_mito_IMM.tsv: </b> This file was obtained from QuickGO (www.ebi.ac.uk/QuickGO/) and contains a list of UniProt IDs associated with the GO term "mitochondrion inner membrane" (GO:0005743). 

- <b> BattleSILAC_PickrellRPKM.csv: </b> This file contains mRNA and protein measurements for the same genes in 59 LCL cell lines (see Kustatscher et al, MSB, 2017 for more information about this dataset).

# Data files available elsewhere
- <b> coregulation_scores.csv: </b> This file contains the final pairwise co-regulation scores for all 12,562,578 protein pairs analysed in our study. It exceeds Github's 100 MB limit even in compressed form. It can be generated by simply running the Coregulation_scores.R script, which only takes ProteomeHD_v1_1.csv as input. Alternatively, a copy of the file can be downloaded from either www.proteomeHD.net or the PRIDE repository with identifier PXD008888.

- <b> 9606.protein.links.detailed.v10.5.txt: </b> Pairwise association scores for human proteins in STRING 10.5 (https://string-db.org).  

- <b> BIOGRID-ALL-3.4.152.tab2.txt: </b> Pairwise interaction scores from BioGRID (https://thebiogrid.org).

- <b> nature22366-s2_Bioplex2.xlsx: </b> This is Table S1 from the BioPlex 2.0 paper (Huttlin et al, Nature, 2017).

- <b> curated_gene_disease_associations.tsv </b> and <b> mapa_geneid_4_uniprot_crossref.tsv: </b> These two files can be downloaded from www.disgenet.org and contain their curated list of gene - disease associations and the mapping of their GeneIDs to Uniprot, respectively.

- <b> proteinGroups.txt </b> The complete MaxQuant search result file for ProteomeHD is available at PRIDE (project identifier PXD008888).

- <b> evidence.txt </b> The complete evidence file produced by MaxQuant for ProteomeHD is available at PRIDE (project identifier PXD008888).





