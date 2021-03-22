
suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(purrr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(igraph))


library(digest)
library(tibble)


# TODO set col_types
island_info <- read_csv('./gifrop_out/classified_island_info.csv')


island_genes_PA <- 
  island_info %>%
  select(island_ID, genes) %>% 
  separate_rows(genes, sep = '\\|') %>%
  mutate(present=1) %>% 
  spread(key = genes, value=present, fill = 0)

# now dereplicate identical islands





# library(tidyverse)

# make fake test data
# 
# gpa <- 
#   tibble('ID'=paste('plasmid', 1:1500, sep = '_'),
#          'Gene1'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene2'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene3'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene4'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene5'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene6'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene7'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene8'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene9'=rbinom(n = 1500, size = 1, prob = .5), 
#          'Gene10'=rbinom(n = 1500, size = 1, prob = .5))
# 
# dups <- 
#   gpa %>% 
#   mutate(ID=paste(ID, 'X', sep = '_'))
# 
# 
# dat <- bind_rows(gpa, dups)

# presence/absence genes on genomic islands
dat_mat <- island_genes_PA %>% column_to_rownames(var = 'island_ID') %>% as.matrix()

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
dereplication_info <- 
  Element_hashes %>%
  group_by(hash) %>% 
  summarise(representative=ID[1],
            num_elements=n(), 
            all_IDs = paste(ID, collapse = '|'))


keepers <- rownames(dat_mat) %in% dereplication_info$representative

de_rep_dat_mat <- dat_mat[keepers,]



# library(vegan)

# jaccard dist? switch to jaccard similarity for edge weights?
sim_mat <- (1 - dist(de_rep_dat_mat, method = 'binary')) %>% as.matrix()

g <- graph_from_adjacency_matrix(sim_mat)


E(g)
# vegdist(de_rep_dat_mat, method = 'jaccard', binary = T)




############# FIX BELOW HERE #############


gpa_gath <- gpa %>%
  gather(key='genome', value='locus_tags', -(1:14)) %>%
  select(Gene,genome, locus_tags) #mmmk...

# this is a tibble of which locus tags are in each island (along with the island's length)
island_ID_loc_tags <- island_info %>%
  select(island_ID, locus_tags, island_length) %>%
  mutate(locus_tags = strsplit(x = locus_tags, split='|', fixed=TRUE)) %>%
  unnest(cols = locus_tags)

# 
# # island_dat_loc <- '/home/julian.trachsel/VDL/Islands/'
# island_dat_loc <- paste(getwd(),'./gifrop_out/my_islands/', sep = '' )

# BUILD GRAPH WITH ISLANDS AS NODES AND EDGES REPRESENTING THE NUMBER OF SHARED GENES
# THEN RUN GRAPH BASED CLUSTERING
#

# adds island info to the roary gene_presence_absence dataframe
gpa_islands <- gpa_gath %>%
  separate_rows(locus_tags, sep = '\t') %>%
  left_join(island_ID_loc_tags) %>%
  filter(!is.na(island_ID))%>%
  group_by(Gene) %>%
  summarise(num_islands=length(unique(island_ID)),
            longest_island=island_ID[which.max(island_length)],
            shortest_island=island_ID[which.min(island_length)],
            len_longest=island_length[which.max(island_length)],
            len_shortest=island_length[which.min(island_length)],
            all_islands=list(unique(island_ID)), 
            .groups='drop') %>%
  
  right_join(gpa)




# This writes out two pangenome csvs
# SHOULD I WAIT TO DO THIS TILL LATER?
# I think I should write this out and then do clustering in a different R script.
poi <- gpa_islands %>%
  mutate(all_islands = map_chr(all_islands, paste, collapse = '|')) %>%
  # write_csv('pan_with_island_info.csv') %>%
  filter(num_islands > 0) #%>%
# write_csv('pan_only_islands.csv')



# ii <- read_csv('./my_islands/island_info.csv')


# poi <- read_csv('./pan_only_islands.csv')


# ipa is a matrix describing which genes are present on which islands
# rows are islands, columns are genes

print('constructing graph: Nodes are islands, edges are number of shared genes')

ipa <- poi %>% select(Gene, all_islands) %>%
  separate_rows(all_islands, sep = '\\|') %>%
  group_by(Gene, all_islands) %>% tally() %>%
  ungroup() %>%
  spread(key=Gene, value = n, fill = 0) %>%
  tibble::column_to_rownames(var = 'all_islands') %>% as.matrix()

#### INSERT COMMUNITY BASED CLUSTERING HERE ####
# maybe save this to do within primary clusters?
# comdist <- dist(ipa, method = 'binary')
# 
# plot(cmdscale(comdist))


###

# if you dont invert this it becomes how many times each gene occurs together in all the islands?
# USE THIS FOR SEPARATE SUBMODULE ANALYSIS
# could interpret edges to represent the number of islands that both of the linked genes co-occur within

ipa <- t(ipa)
# ipa[1:10,1:10]
# after these two steps the rows are genes and the columns are islands
# TRUE means that gene is present on that island and FALSE means that gene is not present on that island

#cross product
co_mat <- t(ipa) %*% ipa
#co_mat is a co-occurance matrix, both rows and columns are islands
# numbers indicate how many genes those islands share with eachother

# set diagonal to 0
diag(co_mat) <- 0


# Create graph from adjacency matrix
# edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)

# Assign nodes weight equal to number of genes?
g <- set.vertex.attribute(g, "v_weight", value = colSums(ipa))

# #
# plot.igraph(g, vertex.label=NA)

# find clusters in this network
print('Primary clustering, any islands sharing any gene will be in the same primary cluster')
cweak <- clusters(mode='weak', g)


# considering changing this to leiden algorithm
print('Secondary clustering: detecting densely connected communities in the graph using the Louvain method')
print('this should separate very different islands that only share a gene or two')
print('but in reality its still pretty bad, so I am adding some graph pruning steps')
clouv <- cluster_louvain(g)


# generate layout so primary and secondary clusters can be plotted on the same layout
# PRUNE LOW WEIGHT EDGES HERE #
# JUST MAKE A DAMN DESCISION AND CLEAN UP THIS MESS ALREADY #
# 
# LO <- layout_nicely(g, dim = 2)


# other clustering options (not all are appropriate)
# ceb <- cluster_edge_betweenness(graph = g)
# cfg <- cluster_fast_greedy(g)
# cim <- cluster_infomap(g, v.weights = NULL)
# clp <- cluster_label_prop(g)
# cle <- cluster_leading_eigen(g)
# copt <- cluster_optimal(g)
# cspin <- cluster_spinglass(g)
# cwt <- cluster_walktrap(g)

# these two are the same for undirected graphs
# I think this is good for very coarse clustering because all connected elements are in one cluster
# ie if islands share even 1 gene they are in the same cluster?
# cstrong <- clusters(mode='strong', g)




# partition <- membership(clouv)

# membership(ceb)
# membership(cfg)
# membership(cim)
# membership(clp)
# membership(cle)
# membership(copt)
# membership(cstrong)
# membership(cwt)
#
# names(membership(clouv)) == names(membership(cweak))

# table(partition)


# is_hierarchical(clouv)
# is_hierarchical(ceb) # TRUE
# is_hierarchical(cfg) # TRUE
# is_hierarchical(cim)
# is_hierarchical(clp)
# is_hierarchical(cle) # TRUE
# is_hierarchical(cstrong)
# is_hierarchical(cwt) # TRUE
# is_hierarchical(cweak)

# as.hclust(cfg)
# as.hclust(clouv)
# plot(x=clouv,
#      y=g,
#      col=membership(cweak),
#      mark.groups = communities(cweak),
#      edge.color = c("black",
#                     "red")[crossing(clouv, g) + 1],
#      vertex.label=NA, main='Graph based clustering, Primary clusters shown')
#
# plot(x=clouv,
#      y=g,
#      col=membership(clouv),
#      mark.groups = communities(clouv),
#      edge.color = c("black",
#                     "red")[crossing(clouv, g) + 1],
#      vertex.label=NA,
#      main='Graph based clustering, Secondary clusters shown')

### Experimental cluster zone ###

# clust_info %>% count(primary_cluster) %>% arrange(desc(n))
# 
# subme <- clust_info %>% filter(primary_cluster == 5) %>% pull(island_ID)
# subg <- induced_subgraph(g, subme)
# plot(subg)
# 




# 
# 
# E(subg)
# 
# 
# E(g)
# 
# # E(subg)$weight
# # 
# # hist(E(subg)$weight, breaks = 100)
# 
# quantile(E(subg)$weight)

# removes edges with weights 1 and 2
# subg.filt <- delete.edges(subg, which(E(subg)$weight > 3))
# #
# plot(subg.filt)
# 

# this can come in handy?
# set_edge_attr(graph, name, index = E(graph), value)
# V(g)[idx]$attr
# 
# # 
# igraph::adjacent_vertices()
# igraph::as_data_frame(g)
# # 
# # thiscanwork <- igraph::as_long_data_frame(subg)
# thiscanwork <- igraph::as_long_data_frame(g)


# special case where there are no edges in the graph...
# how to deal with this?
g_long_df <- igraph::as_long_data_frame(g)

if (nrow(g_long_df) > 0){
  
  edge_values <- 
    g_long_df %>%
    mutate(
      largestVw=ifelse(from_v_weight > to_v_weight, from_v_weight, to_v_weight), 
      LVWmEW=largestVw-weight, 
      LVWdEW=largestVw/weight, 
      absdifVW=abs(from_v_weight - to_v_weight), 
      AdifVWmEW=absdifVW - weight, 
      AdifVWdEW=absdifVW / weight, 
      maybe=LVWmEW/largestVw, 
      maybe2=LVWmEW/weight, 
      sumVWsm2xEW=(to_v_weight + from_v_weight)- (2*weight), 
      sumVWsd2xEW=(to_v_weight + from_v_weight) / (2*weight), 
      x2EWdsumVWs=(2*weight)/(to_v_weight + from_v_weight))
  
  
  
  # Largest Vertex minus edge weight divided by edge weight
  # what proportion of the 
  # hist(edge_values$weight, breaks=50)
  # hist(edge_values$largestVw, breaks=50)
  # hist(edge_values$LVWmEW, breaks=50)
  # hist(edge_values$LVWdEW, breaks=50)
  # hist(edge_values$absdifVW, breaks=50)
  # hist(edge_values$AdifVWmEW, breaks = 50)
  # hist(edge_values$AdifVWdEW, breaks = 50)
  # hist(edge_values$maybe, breaks = 50)
  # hist(edge_values$maybe2, breaks = 50)
  # 
  # hist(edge_values$sumVWsm2xEW, breaks = 50)
  # hist(edge_values$sumVWsd2xEW, breaks = 500, xlim = c(0,10))
  # 
  # hist(edge_values$x2EWdsumVWs, breaks = 100)
  # # 
  # # 
  # min(edge_values$x2EWdsumVWs)
  # max(edge_values$x2EWdsumVWs)
  # 
  # library(Rtsne)
  # 
  # # 
  # ev4t <- edge_values %>%
  #   tidyr::unite(from, to, from_name, to_name, sep='~', col=ID) %>%
  #   column_to_rownames(var='ID') %>% as.matrix() %>% scale()
  # 
  # # 
  # 
  # ev4t %>% count(ID) %>% arrange(desc(n))
  # Rtsne(unique(ev4t))
  
  # plot(cmdscale(dist(ev4t)))
  
  
  # scale(ev4t)
  # pca1 <- princomp(ev4t)
  # pca2 <- prcomp(ev4t)
  # 
  # plot(pca1$scores[,1],pca1$scores[,2] )
  # plot(pca1$scores[,1],pca1$scores[,2] )
  # 
  # tidyr::unite()
  
  # min(edge_values$sumVWsd2xEW)
  # min(edge_values$sumVWsm2xEW)
  # 
  # quantile(edge_values$sumVWsd2xEW)
  # 
  # hist(log(edge_values$sumVWsm2xEW), breaks = 50)
  # hist(log(edge_values$sumVWsd2xEW), breaks = 50)
  
  # quantile(edge_values$x2EWdsumVWs)
  # 
  
  
  # should probably just pick a cutoff
  #
  # evq <- quantile(edge_values$maybe)
  # hist(edge_values$x2EWdsumVWs)
  bad_edges <- which(edge_values$x2EWdsumVWs < .5)
  
  num_bad <- length(bad_edges)
  tot_edge <- nrow(edge_values)
  percent_bad <- round((num_bad/tot_edge) *100)
  
  
  prune_message <- paste('removing ',
                         num_bad,
                         ' edges from the graph. This is ',
                         percent_bad,
                         '% of the total edges', 
                         sep='')
  
  print(prune_message)
  
  
  g_filt <- delete_edges(g, bad_edges)
  
  # plot(g_filt)
  
  
  print('Tertiary clustering, after removing edges, any islands sharing any gene will be in the same primary cluster')
  cweak2 <- clusters(mode='weak', g_filt)
  
  
  # considering changing this to leiden algorithm
  print('Quaternary clustering: detecting densely connected communities in the pruned graph using the Louvain method')
  # print('this should separate very different islands that only share a gene or two')
  clouv2 <- cluster_louvain(g_filt)
  
  
  # generate layout so primary and secondary clusters can be plotted on the same layout
  # PRUNE LOW WEIGHT EDGES HERE #
  
  LO <- layout_nicely(g_filt, dim = 2)
  
  
  # one idea is to calculate the value of an edge based on it's weight (shared genes)
  # and the vertex weights (num genes in each island)
  
  
  
  clust_info <- tibble(island_ID = names(membership(clouv)),
                       #island_ID2 = names(membership(clouv2)),
                       primary_cluster = membership(cweak),
                       secondary_cluster = membership(clouv), 
                       tertiary_cluster = membership(cweak2), 
                       quat_cluster = membership(clouv2))
  
  
  
  
  # clust_info %>% group_by(quat_clust) %>% tally() %>% arrange(desc(n))
  # 
  # 
  # # I thnk i stopped messing with stuff here #
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # ### FOR PLOTTING NETWORKS AND COLORING BY CLUSTER ###
  # # library(RColorBrewer)
  # #
  # # node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
  # # plot(graph_object, vertex.color = node.cols)
  # 
  # ### clustering plots...
  # 
  # # TODO # Make better clustering figures.
  # # maybe make individual figures for difficult clusters
  # 
  # png('./gifrop_out/figures/Primary_clustering.png', width=600, height=600, res=120, type="cairo")
  # p <- plot(x=clouv,
  #           y=g,
  #           col=membership(cweak),
  #           mark.groups = communities(cweak),
  #           edge.color = c("black",
  #                          "red")[crossing(clouv, g) + 1],
  #           vertex.label=NA,
  #           layout=LO,
  #           main='Graph based clustering, Primary clusters shown')
  # 
  # 
  # dev.off()
  # 
  # png('./gifrop_out/figures/Secondary_clustering.png', width=600, height=600, res=120, type="cairo")
  # p <- plot(x=clouv,
  #           y=g,
  #           col=membership(clouv),
  #           mark.groups = communities(clouv),
  #           edge.color = c("black",
  #                          "red")[crossing(clouv, g) + 1],
  #           vertex.label=NA,
  #           layout=LO,
  #           main='Graph based clustering, Secondary clusters shown')
  # 
  # 
  # dev.off()
  # 
  # png('./gifrop_out/figures/Secondary_clustering.png', width=600, height=600, res=120, type="cairo")
  # p <- plot(x=clouv,
  #           y=g,
  #           col=membership(clouv),
  #           mark.groups = communities(clouv),
  #           edge.color = c("black",
  #                          "red")[crossing(clouv, g) + 1],
  #           vertex.label=NA,
  #           layout=LO,
  #           main='Graph based clustering, Secondary clusters shown')
  # 
  # 
  # dev.off()
  # 
  # 
  # 
  # png('./gifrop_out/figures/Tertiary_clustering.png', width=600, height=600, res=120, type="cairo")
  # p <- plot(x=clouv2,
  #           y=g_filt,
  #           col=membership(cweak2),
  #           mark.groups = communities(cweak2),
  #           edge.color = c("black",
  #                          "red")[crossing(clouv2, g_filt) + 1],
  #           vertex.label=NA,
  #           layout=LO,
  #           main='Graph based clustering, Secondary clusters shown')
  # 
  # 
  # dev.off()
  
  # now need to make a clust_info
  ### END INSERT
  
  
  # final clust info construction
  # joins island clustering results with island typing dataframes
  
  
  
} else {
  # if there were no edges in the clustering graph, make a fake clust_info 
  # dataframe with each island assigned to its own cluster
  print('no edges in the graph, no clusters with more than one member')
  clust_info <- tibble(island_ID = names(membership(clouv)),
                       #island_ID2 = names(membership(clouv2)),
                       primary_cluster = membership(cweak),
                       secondary_cluster = membership(clouv), 
                       tertiary_cluster = membership(cweak), 
                       quat_cluster = membership(clouv))
  
}

clust_info <- 
  clust_info %>%
  mutate(genome=sub('(.*)_[0-9]+_[0-9]+', '\\1', island_ID)) %>%
  left_join(island_info) %>%
  left_join(island_types) %>%
  mutate(island_type=ifelse(is.na(island_type), 'unknown', island_type)) %>%
  left_join(res_types) %>%
  left_join(vir_types) %>%
  left_join(plasmid_types) %>%
  left_join(viro_types) %>%
  left_join(megares_types) %>% 
  select(island_ID, acc_frag, ends_with('cluster'), everything(), -genome) %>%
  write_csv('./gifrop_out/clustered_island_info.csv')



# making gpa with island clustering info

gpa_clust <- gpa_gath %>%
  left_join(island_ID_loc_tags) %>%
  filter(!is.na(island_ID))%>%
  left_join(clust_info, by='island_ID') %>%
  group_by(Gene) %>%
  summarise(num_islands=length(unique(island_ID)),
            num_Sclusts=length(unique(secondary_cluster)),
            all_islands=list(unique(island_ID)),
            Pcluster = unique(primary_cluster),
            all_Sclusters=list(unique(secondary_cluster)),
            all_Tclusters=list(unique(tertiary_cluster)), 
            all_Qclusters=list(unique(quat_cluster)),
            .groups='drop') %>%
  
  right_join(gpa)



# this collapses list columns 'all_islands' and 'all_Sclusters' to pipe delimited strings
# so that the gpa_clust dataframe can be written to csv

gpa_clust %>%
  mutate(all_islands = map_chr(all_islands, paste, collapse = '|'),
         all_Sclusters  = map_chr(all_Sclusters , paste, collapse = '|'), 
         all_Tclusters = map_chr(all_Tclusters, paste, collapse = '|'), 
         all_Qclusters = map_chr(all_Qclusters, paste, collapse = '|')) %>%
  write_csv('./gifrop_out/pan_with_island_info.csv') %>%
  filter(num_islands > 0) %>%
  write_csv('./gifrop_out/pan_only_islands.csv')





# some simple stats

# what proportion of the pangenome is contained within these islands?

#total genes in pangenome
num_tot_genes <- nrow(gpa_clust)

# genes in core genome
num_core_genes <- sum(gpa_clust$`No. isolates` == max(gpa_clust$`No. isolates`))

# genes in accessory genome
num_accessory_genes <- sum(gpa_clust$`No. isolates` < max(gpa_clust$`No. isolates`))

# acessory genes on islands
num_access_on_island <- nrow(gpa_clust %>% filter(`No. isolates` < max(gpa_clust$`No. isolates`)) %>% filter(num_Sclusts > 0))

# acessory genes not on islands
num_access_not_on_island <- nrow(gpa_clust %>% filter(`No. isolates` < max(gpa_clust$`No. isolates`)) %>% filter(is.na(num_Sclusts)))


# % of the pangenome contained within these genomic islands

print('Percent of pan genome contained within genomic islands:')
(nrow(gpa_clust %>% filter(num_Sclusts > 0)) / nrow(gpa_clust)) * 100

# 68% of the accessory genome is contained on these islands
print('Percent of accessory genome contained within these islands:')
(num_access_on_island / num_accessory_genes) * 100


