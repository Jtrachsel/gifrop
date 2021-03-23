args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(igraph))
library(digest)
library(proxy)


# setwd('./test_data/pan/gifrop_out/')

# TODO set col_types
island_info <- read_csv('./gifrop_out/classified_island_info.csv', 
                        col_types = c('cccddddcddlcccccccccc'))


island_genes_PA <- 
  island_info %>%
  select(island_ID, genes) %>% 
  separate_rows(genes, sep = '\\|') %>%
  mutate(present=1) %>% 
  spread(key = genes, value=present, fill = 0) %>% 
  write_csv('./gifrop_out/island_genes_PA.csv')

# now dereplicate identical islands


# presence/absence genes on genomic islands
dat_mat <- island_genes_PA %>% 
  column_to_rownames(var = 'island_ID') %>%
  as.matrix()


# hash the presence/absence vector for each island.
Element_hashes <- 
  tibble(ID=rownames(dat_mat),           #element ID
         hash=apply(dat_mat, 1, digest)  # calculates the hash of gene vectors
  )


# if two plasmids have identical gene vectors, they are the same plasmid
# only need to proceed with one. 
# depending on how the pangenome was constructed, plasmids with the same genes in different orders 
# will be considered the same plasmid.
print('dereplicating islands with identical gene content...')
print('generating hashes for each gene presence/absence vector')
dereplication_info <- 
  Element_hashes %>%
  group_by(hash) %>% 
  summarise(island_ID=ID[1],
            num_elements=n(), 
            all_IDs = paste(ID, collapse = '|'))


keepers <- rownames(dat_mat) %in% dereplication_info$island_ID

percent_keepers <- round(sum(keepers) / length(keepers) * 100, 1)

print(paste('Dereplication reduced islands to', percent_keepers, 'percent of the original set'))


de_rep_dat_mat <- dat_mat[keepers,]


#### IF SWITCH HERE FOR TRAD CLUSTERING VS NEW CLUSTERING?
print('constructing graph: Nodes are islands, edges are simpson similarity index (AKA overlap coef) between joined nodes')

# jaccard dist? switch to jaccard similarity for edge weights?
print('calculating simpson similarity AKA overlap coefficient')
sim_mat <- 
  proxy::simil(de_rep_dat_mat, method = 'simpson') %>% 
  as.matrix() %>%
  Matrix::Matrix(sparse = T)

g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)

 
# find clusters in this network
print('Primary clustering, any islands sharing any gene will be in the same primary cluster')
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
  proxy::simil(de_rep_dat_mat, method = 'jaccard') %>% 
  as.matrix() %>%
  Matrix::Matrix(sparse = T)

g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)

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
  left_join(clust_info) %>%
  separate_rows(all_IDs, sep = '\\|') %>% 
  select(-island_ID) %>% 
  transmute(island_ID=all_IDs, 
            primary_cluster=primary_cluster, 
            secondary_cluster=secondary_cluster, 
            tertiary_cluster=tertiary_cluster, 
            quat_cluster=quat_cluster) %>% 
  bind_rows(clust_info) %>% 
  unique()


# output clustered island info.  
clust_info <- 
  clust_info %>%
  mutate(genome=sub('(.*)_[0-9]+_[0-9]+', '\\1', island_ID)) %>%
  left_join(island_info) %>% 
  mutate(island_type=ifelse(is.na(island_type), 'unknown', island_type)) %>%
  select(island_ID, acc_frag, ends_with('cluster'), everything(), -genome) %>%
  write_csv('./gifrop_out/clustered_island_info.csv')


