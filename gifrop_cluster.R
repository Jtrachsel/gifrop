args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(igraph))
library(digest)
library(parallelDist)
library(purrr)

cluster_islands <- 
  function(dereplicated_island_info, dereplication_info, prefix=NULL){
    # browser()
    # clusters genomic islands into 4 levels
    #takes: 
    # 1) dereplicated genomic islands
    # 2) dereplication info to map cluster assignments to full dataset
    # 3) prefix that is added to cluster assignment
    print('constructing island_gene_presence_absence matrix')
    dat_mat <- 
      dereplicated_island_info %>% 
      # filter(island_ID %in% keep_these_islands) %>% #this should be passed a fully reduced set of islands
      select(island_ID, genes) %>% 
      separate_rows(genes, sep = '\\|') %>%
      mutate(present=1) %>% 
      pivot_wider(names_from = genes, values_from=present, values_fill=0) %>% 
      write_csv('./gifrop_out/dereplicated_island_genes_PA.csv') %>%
      column_to_rownames(var='island_ID') %>% 
      as.matrix() 
    
    
    #### IF SWITCH HERE FOR TRAD CLUSTERING VS NEW CLUSTERING?
    print('constructing graph: Nodes are islands, edges are simpson similarity index (AKA overlap coef) between joined nodes')
    
    # jaccard dist? switch to jaccard similarity for edge weights?
    print('calculating simpson similarity AKA overlap coefficient')
    sim_mat <- 
      (1 - parallelDist::parallelDist(dat_mat, method = 'simpson')) %>% 
      as.matrix() #%>%
    # Matrix::Matrix(sparse = T)
    
    g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)
    g <- delete_edges(g, E(g)[weight == 0])
    
    print('writing overlap coef graph')
    igraph::write.graph(g, file = './gifrop_out/overlap_coef_graph.dot', format = 'dot')
    # find clusters in this network
    print('Primary clustering, any islands sharing any number of genes will be in the same primary cluster')
    clust1 <- clusters(mode='weak', g)
    
    
    print('secondary clustering')
    print('pruning graph, removing edges with overlap coef of less than .5')
    print('an overlap coefficient of .5 means that at least 1/2 of the genes in the smaller island are also present in the other island')
    print('any two islands with an overlap coef of at least .5 are in the same secondary cluster')
    
    g <- delete_edges(g, E(g)[weight<.5])
    clust2 <- clusters(g)
    
    print('pruning graph, removing edges with overlap coef of less than 1')
    
    
    g <- delete_edges(g, E(g)[weight<1])
    clust3 <- clusters(g)
    
    print('constructing new graph with jaccard similarities')
    print('This will allow the identification of extremely similar genomic islands')
    
    sim_mat <- 
      (1 - parallelDist::parallelDist(dat_mat, method = 'binary')) %>% 
      as.matrix() #%>%
    # Matrix::Matrix(sparse = T)
    
    g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)
    igraph::write.graph(g, file = './gifrop_out/jaccard_coef_graph.dot', format = 'dot')
    print('removing edges representing jaccard similarities of less than .75')
    
    
    g <- delete_edges(g, E(g)[weight<.75])
    clust4 <- clusters(g)
    
    clust_info <- tibble(island_ID = names(membership(clust1)),
                         primary_cluster = membership(clust1),
                         secondary_cluster = membership(clust2),
                         tertiary_cluster = membership(clust3),
                         quat_cluster = membership(clust4))
    
    ### re-replicate and propogate cluster labels
    clust_info <- 
      dereplication_info %>%
      right_join(clust_info) %>%
      separate_rows(all_IDs, sep = '\\|') %>% 
      select(-island_ID) %>% 
      transmute(island_ID=all_IDs, 
                primary_cluster=paste0(primary_cluster, prefix), 
                secondary_cluster=paste0(secondary_cluster,prefix), 
                tertiary_cluster=paste0(tertiary_cluster, prefix), 
                quat_cluster=paste0(quat_cluster, prefix)) %>% 
      unique() # shouldn't need this...
    
    
    return(clust_info)
  }


# setwd('/90daydata/fsepru113/jtrachsel/Jan_21_mono_update/verified_i4512i/passing_QC/PLASMID_COLLECTOR/all_islands/target_islands/pan/')

island_info <- read_csv('./gifrop_out/classified_island_info.csv', 
                        col_types = c('cccddddcddlcccccccccc')) 


dereplication_info <- 
  island_info %>% 
  mutate(HASH=map_chr(.x = genes, .f = ~digest(.x))) %>%   # makes a hash of every island's "genes" column
  group_by(HASH) %>%
  summarise(all_IDs = paste(island_ID, collapse = '|'),
            island_ID=island_ID[1],
            num_elements=n()
            ) %>% 
  write_csv('./gifrop_out/dereplication_info.csv')

dereplicated_island_info <- 
  island_info %>% 
  filter(island_ID %in% dereplication_info$island_ID)


if (!exists('reduce_clustering_problem') & nrow(dereplicated_island_info) > 50000){
  reduce_clustering_problem <- T
} else{ 
  reduce_clustering_problem <- F
}



if (reduce_clustering_problem){
  print('REDUCING CLUSTERING PROBLEM BY EXCLUDING "ONLY_PHAGE" AND SMALL UNKNOWNS')
  
  phages <- 
    dereplicated_island_info %>% 
    filter(grepl('phage', island_type))
  
  only_phage <- 
    dereplicated_island_info %>% 
    filter(island_type == 'phage') %>% 
    pull(island_ID)
  
  small_unknowns <- dereplicated_island_info %>%
    filter(is.na(island_type) & num_genes < 10) %>%
    pull(island_ID)
  
  
  reduced_dereplicated_island_info <- 
    dereplicated_island_info %>%
    filter(!(island_ID %in% c(only_phage, small_unknowns)))
  
  clust_info <- cluster_islands(reduced_dereplicated_island_info, dereplication_info)
  phages_clust_info <- cluster_islands(phages, dereplication_info, prefix = 'P') 
  
  bind_rows(clust_info, phages_clust_info) %>% 
    mutate(genome=sub('(.*)_[0-9]+_[0-9]+', '\\1', island_ID)) %>%
    left_join(island_info) %>% 
    mutate(island_type=ifelse(is.na(island_type), 'unknown', island_type)) %>%
    select(island_ID, acc_frag, ends_with('cluster'), everything(), -genome)
    write_csv(file = './gifrop_out/clustered_island_info.csv')
  
  
} else {
  ### 20000x20000 clustering problem took ~2 hours? for each distance calc (so 4 hours tot)
  print('NOT REDUCING CLUSTERING PROBLEM')
  
  cluster_islands(dereplicated_island_info, dereplication_info) %>% 
    mutate(genome=sub('(.*)_[0-9]+_[0-9]+', '\\1', island_ID)) %>%
    left_join(island_info) %>% 
    mutate(island_type=ifelse(is.na(island_type), 'unknown', island_type)) %>%
    select(island_ID, acc_frag, ends_with('cluster'), everything(), -genome) %>% 
    write_csv(file = './gifrop_out/clustered_island_info.csv')
  
}
# cii <- read_csv('./gifrop_out/clustered_island_info.csv')
# 
# cii <- cii %>% left_join(dereplication_info) %>% separate_rows(all_IDs,sep='\\|')
# 
# cii %>% write_csv('TEMP_CII.csv')
