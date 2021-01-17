test <- 
  as.matrix(comdist) %>% 
  as.data.frame() %>%
  rownames_to_column(var='island1') %>% 
  gather(key = 'island2', value = 'distance', -island1)




hist(test$distance, breaks = 100)

hclust_obj <- hclust(comdist)

test <- test %>% filter(distance != 1)

hist(edge_values$x2EWdsumVWs)


plot(test)
hist(comdist, breaks = 1000)



cutree(hclust_obj, h = .95)
cutree(hclust_obj, h = .5)
cutree(hclust_obj, h = .25)
cutree(hclust_obj, h = .00001)


# oA+oB             = oA+oB
# -----               ------
# oA+oB+AaB           
# 
# 
# 2*AaB             =  2AaB
# -------             -----
# aA + aB             oA+oB+2AaB
# 
# 
# 
# 
# 
# 
# 

###### submodule test #########



ipa <- poi %>% select(Gene, all_islands) %>%
  separate_rows(all_islands, sep = '\\|') %>%
  group_by(Gene, all_islands) %>% tally() %>%
  ungroup() %>%
  spread(key=Gene, value = n, fill = 0) %>%
  column_to_rownames(var = 'all_islands') %>% as.matrix()


nrow(ipa)
ncol(ipa)
ipa[1:10,1:10]
#### INSERT COMMUNITY BASED CLUSTERING HERE ####
# maybe save this to do within primary clusters?
# comdist <- dist(ipa, method = 'binary')
# 
# plot(cmdscale(comdist))


###

# if you dont invert this it becomes how many times each gene occurs together in all the islands?
# USE THIS FOR SEPARATE SUBMODULE ANALYSIS
# could interpret edges to represent the number of islands that both of the linked genes co-occur within

# ipa <- t(ipa)
# ipa[1:10,1:10]
# after these two steps the rows are genes and the columns are islands
# TRUE means that gene is present on that island and FALSE means that gene is not present on that island

#cross product
gene_cooccur <- t(ipa) %*% ipa
#co_mat is a co-occurance matrix, both rows and columns are islands
# numbers indicate how many genes those islands share with eachother

# set diagonal to 0
diag(gene_cooccur) <- 0

nrow(gene_cooccur) * ncol(gene_cooccur)
# Create graph from adjacency matrix
# edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(gene_cooccur, mode = "upper", weighted = TRUE)

# Assign nodes weight equal to number of genes?
g <- set.vertex.attribute(g, "v_weight", value = colSums(ipa))

# #
# plot.igraph(g, vertex.label=NA)

# find clusters in this network
print('Primary clustering, any genes that occur together in at least 1 island will be in the same primary cluster')
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




edge_values <- 
  igraph::as_long_data_frame(g) %>%
  mutate(
    largestVw=ifelse(from_v_weight > to_v_weight, from_v_weight, to_v_weight), 
    numONLY_from=from_v_weight - weight,
    numONLY_to  = to_v_weight - weight,
    NUM_EITHER  =numONLY_from + numONLY_to + weight,
    JACARD=weight / NUM_EITHER,
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

hist(edge_values$JACARD, breaks = 200)

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

