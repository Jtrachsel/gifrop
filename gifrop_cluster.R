args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(purrr, quietly = TRUE, warn.conflicts = FALSE))
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

# probable sparse matricies....
# sum(dat_mat > 0)/ (nrow(dat_mat) * ncol(dat_mat)) *100 


# apply(dat_mat, 1, digest)

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
# library(vegan)
print('constructing graph: Nodes are islands, edges are simpson similarity index (AKA overlap coef) between joined nodes')

# jaccard dist? switch to jaccard similarity for edge weights?
# sim_mat <- (1 - dist(de_rep_dat_mat, method = 'binary')) %>% as.matrix()
print('calculating simpson similarity AKA overlap coefficient')
sim_mat <- 
  proxy::simil(de_rep_dat_mat, method = 'simpson') %>% 
  as.matrix() %>%
  Matrix::Matrix(sparse = T)

g <- graph_from_adjacency_matrix(adjmatrix = sim_mat,  weighted = T, mode='upper', diag = F)

############# FIX BELOW HERE #############

# 
# gpa_gath <- gpa %>%
#   gather(key='genome', value='locus_tags', -(1:14)) %>%
#   select(Gene,genome, locus_tags) #mmmk...
# 
# # this is a tibble of which locus tags are in each island (along with the island's length)
# island_ID_loc_tags <- island_info %>%
#   select(island_ID, locus_tags, island_length) %>%
#   mutate(locus_tags = strsplit(x = locus_tags, split='|', fixed=TRUE)) %>%
#   unnest(cols = locus_tags)
# 
# # 
# # # island_dat_loc <- '/home/julian.trachsel/VDL/Islands/'
# # island_dat_loc <- paste(getwd(),'./gifrop_out/my_islands/', sep = '' )
# 
# # BUILD GRAPH WITH ISLANDS AS NODES AND EDGES REPRESENTING THE NUMBER OF SHARED GENES
# # THEN RUN GRAPH BASED CLUSTERING
# #
# 
# # adds island info to the roary gene_presence_absence dataframe
# gpa_islands <- gpa_gath %>%
#   separate_rows(locus_tags, sep = '\t') %>%
#   left_join(island_ID_loc_tags) %>%
#   filter(!is.na(island_ID))%>%
#   group_by(Gene) %>%
#   summarise(num_islands=length(unique(island_ID)),
#             longest_island=island_ID[which.max(island_length)],
#             shortest_island=island_ID[which.min(island_length)],
#             len_longest=island_length[which.max(island_length)],
#             len_shortest=island_length[which.min(island_length)],
#             all_islands=list(unique(island_ID)), 
#             .groups='drop') %>%
#   
#   right_join(gpa)
# 
# 
# 
# 
# # This writes out two pangenome csvs
# # SHOULD I WAIT TO DO THIS TILL LATER?
# # I think I should write this out and then do clustering in a different R script.
# poi <- gpa_islands %>%
#   mutate(all_islands = map_chr(all_islands, paste, collapse = '|')) %>%
#   # write_csv('pan_with_island_info.csv') %>%
#   filter(num_islands > 0) #%>%
# # write_csv('pan_only_islands.csv')
# 


# ii <- read_csv('./my_islands/island_info.csv')


# poi <- read_csv('./pan_only_islands.csv')


# ipa is a matrix describing which genes are present on which islands
# rows are islands, columns are genes

# 
# ipa <- poi %>% select(Gene, all_islands) %>%
#   separate_rows(all_islands, sep = '\\|') %>%
#   group_by(Gene, all_islands) %>% tally() %>%
#   ungroup() %>%
#   spread(key=Gene, value = n, fill = 0) %>%
#   tibble::column_to_rownames(var = 'all_islands') %>% as.matrix()

#### INSERT COMMUNITY BASED CLUSTERING HERE ####
# maybe save this to do within primary clusters?
# comdist <- dist(ipa, method = 'binary')
# 
# plot(cmdscale(comdist))


###

# if you dont invert this it becomes how many times each gene occurs together in all the islands?
# USE THIS FOR SEPARATE SUBMODULE ANALYSIS
# could interpret edges to represent the number of islands that both of the linked genes co-occur within
# 
# ipa <- t(ipa)
# # ipa[1:10,1:10]
# # after these two steps the rows are genes and the columns are islands
# # TRUE means that gene is present on that island and FALSE means that gene is not present on that island
# 
# #cross product
# co_mat <- t(ipa) %*% ipa
# #co_mat is a co-occurance matrix, both rows and columns are islands
# # numbers indicate how many genes those islands share with eachother
# 
# # set diagonal to 0
# diag(co_mat) <- 0


# Create graph from adjacency matrix
# edge weights are equal to frequency of co-occurrence
# g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)

# Assign nodes weight equal to number of genes?
# g <- set.vertex.attribute(g, "v_weight", value = colSums(ipa))

# #
# plot.igraph(g, vertex.label=NA)

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



dereplication_info %>%
  left_join(clust_info) %>%
  separate_rows(all_IDs, sep = '\\|') %>% 
  select(-island_ID) %>% 
  transmute(island_ID=all_IDs, 
            primary_cluster=primary_cluster, 
            secondary_cluster=secondary_cluster, 
            tertiary_cluster=tertiary_cluster, 
            quat_cluster=quat_cluster) %>% 
  bind_rows(clust_info)






clust_info <- 
  clust_info %>%
  mutate(genome=sub('(.*)_[0-9]+_[0-9]+', '\\1', island_ID)) %>%
  left_join(island_info) %>% 
  # left_join(island_info) %>%
  # left_join(island_types) %>%
  mutate(island_type=ifelse(is.na(island_type), 'unknown', island_type)) %>%
  select(island_ID, acc_frag, ends_with('cluster'), everything(), -genome) %>%
  write_csv('./gifrop_out/clustered_island_info.csv')


# 
# # making gpa with island clustering info
# 
# gpa_clust <- gpa_gath %>%
#   left_join(island_ID_loc_tags) %>%
#   filter(!is.na(island_ID))%>%
#   left_join(clust_info, by='island_ID') %>%
#   group_by(Gene) %>%
#   summarise(num_islands=length(unique(island_ID)),
#             num_Sclusts=length(unique(secondary_cluster)),
#             all_islands=list(unique(island_ID)),
#             Pcluster = unique(primary_cluster),
#             all_Sclusters=list(unique(secondary_cluster)),
#             all_Tclusters=list(unique(tertiary_cluster)), 
#             all_Qclusters=list(unique(quat_cluster)),
#             .groups='drop') %>%
#   
#   right_join(gpa)
# 
# 
# 
# # this collapses list columns 'all_islands' and 'all_Sclusters' to pipe delimited strings
# # so that the gpa_clust dataframe can be written to csv
# 
# gpa_clust %>%
#   mutate(all_islands = map_chr(all_islands, paste, collapse = '|'),
#          all_Sclusters  = map_chr(all_Sclusters , paste, collapse = '|'), 
#          all_Tclusters = map_chr(all_Tclusters, paste, collapse = '|'), 
#          all_Qclusters = map_chr(all_Qclusters, paste, collapse = '|')) %>%
#   write_csv('./gifrop_out/pan_with_island_info.csv') %>%
#   filter(num_islands > 0) %>%
#   write_csv('./gifrop_out/pan_only_islands.csv')
# 
# 
# 
# 

# some simple stats

# what proportion of the pangenome is contained within these islands?

#total genes in pangenome
# num_tot_genes <- nrow(gpa_clust)
# 
# # genes in core genome
# num_core_genes <- sum(gpa_clust$`No. isolates` == max(gpa_clust$`No. isolates`))
# 
# # genes in accessory genome
# num_accessory_genes <- sum(gpa_clust$`No. isolates` < max(gpa_clust$`No. isolates`))
# 
# # acessory genes on islands
# num_access_on_island <- nrow(gpa_clust %>% filter(`No. isolates` < max(gpa_clust$`No. isolates`)) %>% filter(num_Sclusts > 0))
# 
# # acessory genes not on islands
# num_access_not_on_island <- nrow(gpa_clust %>% filter(`No. isolates` < max(gpa_clust$`No. isolates`)) %>% filter(is.na(num_Sclusts)))
# 
# 
# # % of the pangenome contained within these genomic islands
# 
# print('Percent of pan genome contained within genomic islands:')
# (nrow(gpa_clust %>% filter(num_Sclusts > 0)) / nrow(gpa_clust)) * 100
# 
# # 68% of the accessory genome is contained on these islands
# print('Percent of accessory genome contained within these islands:')
# (num_access_on_island / num_accessory_genes) * 100


