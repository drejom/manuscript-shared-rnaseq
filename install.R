bioc_pkgs <- c('knitr',
'OrganismDbi',
'Mus.musculus',
'Homo.sapiens',
'TxDb.Hsapiens.UCSC.hg38.knownGene',
'ggrepel',
'pheatmap',
'corpcor',
'Biobase',                      #Functions that are needed by many other packages or which replace R functions.
'biomaRt',                      #Interface to BioMart databases
'coop',                         #Fast implementations of the co-operations: covariance, correlation, and cosine similarity.
'doParallel',
'flashClust',
'foreach',                      #Provides Foreach Looping Construct
'ggplot2',
'goseq',
'heatmap.plus',
'GO.db',
'gplots',
'knitcitations',
'limma',
'Matrix',
'plyr',
'dplyr',
'readr',
'rtracklayer',
'parallelDist',
'preprocessCore',
'printr',
'reshape2',
'RColorBrewer',
'RCurl',
'sva',
'viridis',
'stephenturner/annotables',
'msa',
'ape',
'Biostrings',
'igraph',
'gridExtra',
'matrixStats',
'WCGNA',
'KEGGREST',
'tftargets',
'grid',
'seqinr',
'RBGL',
'Rgraphviz'
)

ip <- installed.packages()
bioc_pkgs <- bioc_pkgs[!(bioc_pkgs %in% rownames(ip))]
BiocManager::install(bioc_pkgs,  dependencies = TRUE, Ncpus = 4)
