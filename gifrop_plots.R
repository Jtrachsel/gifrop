### PLOTS! ###

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

### FOR TESTING ONLY ###
# setwd('~/Documents/IslandR/test_data/pan')
# setwd('~/Documents/gifrop_examples/test3/pan/')



suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(purrr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(igraph))





####### PLOTS #########
# should move this to its own script?
# check for cairo?

clust_info <- read_csv('./gifrop_out/clustered_island_info.csv',
                       col_types = c('ccccccccddddcddlccccccccc'))

### USE ME
p <- clust_info %>%
  ggplot(aes(x=island_length, fill=island_type)) +
  geom_histogram(bins = 50, color='black') + theme_bw() +scale_fill_brewer(palette = 'Set1')+
  ggtitle('Histogram of the lengths of all detected genomic islands')



png('./gifrop_out/figures/island_length_histogram.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()




### cluster types ###
# probably adjust the height and width based on num clusters and num types

p <- clust_info %>%
  ggplot(aes(x=num_genes, y=primary_cluster, fill=island_type)) + 
  geom_point(position = position_jitter(width = .2), color='white')

png('./gifrop_out/figures/primary_cluster_types.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()



p <- clust_info %>% ggplot(aes(x=num_genes, y=secondary_cluster, fill=island_type)) + 
  geom_point(position = position_jitter(width = .2), color='white')

png('./gifrop_out/figures/secondary_cluster_types.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()



p <- clust_info %>% ggplot(aes(x=num_genes, y=tertiary_cluster, fill=island_type)) + 
  geom_point(position = position_jitter(width = .2), color='white')

png('./gifrop_out/figures/tertiary_cluster_types.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()



p <- clust_info %>% ggplot(aes(x=num_genes, y=quat_cluster, fill=island_type)) + 
  geom_point(position = position_jitter(width = .2), color='white')

png('./gifrop_out/figures/quat_cluster_types.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()

# 
# ####  USE ME? ####
# 
# islands_per_isolate <- clust_info %>%
#   group_by(genome_name) %>%
#   tally()
# 
# width <- (nrow(islands_per_isolate)/10) *120
# 
# if (width < 600){
#   width <- 600
# }
# 
# 
# 
# p <- clust_info %>%
#   group_by(genome_name, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('number of genomic islands per isolate (All Islands)')
# 
# 
# png('./gifrop_out/figures/islands_per_isolate.png', width=width, height=600, res=120, type="cairo")
# print(p)
# dev.off()
# 
# ####  Excluding unknowns ###
# 
# p <- clust_info %>% filter(island_type != 'unknown') %>%
#   group_by(genome_name, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('number of genomic islands per isolate (All Islands)')
# 
# 
# png('./gifrop_out/figures/islands_per_isolate_no_unknowns.png', width=width, height=600, res=120, type="cairo")
# print(p)
# dev.off()
# 
# 
# 
# # #### Now limiting islands by length for the next few ##
# # p <- clust_info %>% filter(island_length > 5000) %>%
# #   group_by(genome_name, island_type) %>%
# #   tally() %>%
# #   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
# #   geom_col(color='black') +
# #   scale_fill_brewer(palette = 'Set1') + theme_bw()+
# #   theme(axis.text.x = element_text(angle = 90)) +
# #   ggtitle('number of genomic islands per isolate (>5kb)')
# #
# #
# # width <- (length(p$data$genome_name)/10) *120
# #
# # if (width < 600){
# #   width <- 600
# # }
# #
# # png('./figures/islands_per_isolate_5kb.png', width=width, height=600, res=120, type="cairo")
# # print(p)
# # dev.off()
# #
# #
# # p <- clust_info %>% filter(island_length > 10000) %>%
# #   group_by(genome_name, island_type) %>%
# #   tally() %>%
# #   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
# #   geom_col(color='black') +
# #   scale_fill_brewer(palette = 'Set1') + theme_bw()+
# #   theme(axis.text.x = element_text(angle = 90)) +
# #   ggtitle('number of genomic islands per isolate (>10kb)')
# #
# # width <- (length(p$data$genome_name)/10)*120
# # if (width < 600){
# #   width <- 600
# # }
# #
# #
# # png('./figures/islands_per_isolate_10kb.png', width=width, height=600, res=120, type="cairo")
# # print(p)
# # dev.off()
# #
# #
# #
# #
# #
# # p <- clust_info %>% filter(island_length > 20000) %>%
# #   group_by(genome_name, island_type) %>%
# #   tally() %>%
# #   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
# #   geom_col(color='black') +
# #   scale_fill_brewer(palette = 'Set1') + theme_bw()+
# #   theme(axis.text.x = element_text(angle = 90)) +
# #   ggtitle('number of genomic islands per isolate (>20kb)')
# #
# # width <- (length(p$data$genome_name)/9)*120
# # if (width < 600){
# #   width <- 600
# # }
# #
# #
# # png('./figures/islands_per_isolate_20kb.png', width=width, height=600, res=120, type="cairo")
# # print(p)
# # dev.off()
# #
# #
# #
# #
# #
# # p <- clust_info %>% filter(island_length > 40000) %>%
# #   group_by(genome_name, island_type) %>%
# #   tally() %>%
# #   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
# #   geom_col(color='black') +
# #   scale_fill_brewer(palette = 'Set1') + theme_bw()+
# #   theme(axis.text.x = element_text(angle = 90)) +
# #   ggtitle('number of genomic islands per isolate (>40kb)')
# #
# # width <- (length(p$data$genome_name)/7)*120
# # if (width < 600){
# #   width <- 600
# # }
# #
# #
# # png('./figures/islands_per_isolate_40kb.png', width=width, height=600, res=120, type="cairo")
# # print(p)
# # dev.off()
# #
# 
# 
# 
# ####
# p <- clust_info %>%
#   mutate(primary_cluster = factor(primary_cluster)) %>%
#   group_by(primary_cluster, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=primary_cluster, y=n, fill=island_type)) +
#   geom_col() +
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('Number of occurances of each genomic island cluster in the pangenome')
# 
# 
# png('./gifrop_out/figures/Number_of_occurances.png', width=840, height=600, res=120, type="cairo")
# print(p)
# dev.off()
# 
# 
# ###############
# 
# p <- clust_info %>%
#   mutate(secondary_cluster = factor(secondary_cluster)) %>%
#   group_by(secondary_cluster, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=secondary_cluster, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   theme_bw() +
#   # theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('Number of occurances of each genomic island cluster in the pangenome') +
#   scale_fill_brewer(palette = 'Set1')
# 
# 
# png('./gifrop_out/figures/Number_of_occurances_secondary.png', width=840, height=600, res=120, type="cairo")
# print(p)
# dev.off()
# 
# ################
# 
# # 
# # gpa_clust <- read_csv('gifrop_out/pan_only_islands.csv')
# # 
# # 
# # 
# # 
# # ### TEST ZONE # COMMENT OUT
# # # Island PA heatmap? #
# # 
# # 
# # # pseudo #
# # 
# # # for all secondary clusters that have variability in the number of genes
# # # make
# # # heatmap where rows are genes and columns are individual islands
# # # histogram of number of genes
# # ## this dataframe tries to assess how consistent all islands within a cluster are
# # # looks at how variable the number of genes is for all islands in a cluster
# # # it seems like ((max_genes - min_genes) / min_genes) gives a pretty good indication
# # # of variability in the cluster
# # cluster_qual <-
# #   clust_info %>% group_by(primary_cluster, secondary_cluster) %>%
# #   summarise(mean_genes = mean(num_genes),
# #             med_genes  = median(num_genes),
# #             var_genes  = var(num_genes),
# #             sd_genes   = sd(num_genes),
# #             min_genes  = min(num_genes),
# #             max_genes  = max(num_genes),
# #             num_occur  = length(num_genes),
# #             maxmin_divmin   = (max_genes - min_genes)/min_genes)
# # 
# # 
# # # variable_clusters <- cluster_qual %>% filter(maxmin_divmin != 0)
# # 
# # 
# # 
# # filt_helper <- function(test_vec, int_vec){
# #   res <- any(int_vec %in% test_vec)
# #   return(res)
# # }
# # 
# # 
# # 
# # 
# # 
# # #######
# # # this makes a heatmap of all the genes in a specified cluster (or vector of clusters) in each ISLAND
# # 
# # gene_by_island_heatmap <- function(gpa_clust, SECONDARY_CLUSTER){
# #   # gpa clust need to be the gene_presence_absence.csv file with the added clustering information
# #   # in addition the 'secondary_cluster' column needs to be a list formatted column
# #   # secondary cluster can be a vector of secondary clusters you want to see together
# #   # ISSUE!! currently genes from other sclusts are pulled in, if a gene is in two different sclusts
# #   # it will show up in these heatmaps but will not have annotation info with it.
# #   # need to include another filtering step to remove islands not belonging to Sclust at hand
# #   
# #   
# #   anno <- clust_info %>%
# #     filter(secondary_cluster %in% SECONDARY_CLUSTER) %>%
# #     select(island_ID, island_type) %>%
# #     column_to_rownames(var = 'island_ID')
# #   
# #   PA <- gpa_clust %>%
# #     filter(map_lgl(.x=all_Sclusters, .f=filt_helper, SECONDARY_CLUSTER)) %>%
# #     unnest(cols = all_islands) %>%
# #     select(Gene, all_islands) %>%
# #     group_by(all_islands, Gene) %>%
# #     tally() %>%
# #     ungroup() %>%
# #     spread(key = Gene, value = n, fill = 0) %>%
# #     column_to_rownames(var = 'all_islands') %>%
# #     as.matrix() %>%
# #     t()   # to get genes as rows and cols as islands
# #   
# #   SECONDARY_CLUSTER <- paste(SECONDARY_CLUSTER, collapse = '_', sep = '_')
# #   MAIN=paste('Presence/absence of genes among islands in secondary cluster',SECONDARY_CLUSTER)
# #   FILENAME=paste('./gifrop_out/Sclust_', SECONDARY_CLUSTER,'_gene_heatmap.jpeg', sep = '')
# #   
# #   width=ncol(PA)/5
# #   if (width < 6){
# #     width <- 6
# #   }
# #   
# #   height=nrow(PA)/10
# #   
# #   if (height < 6){
# #     height <- 6
# #   }
# #   
# #   
# #   pheatmap(PA, filename = FILENAME,
# #            height = height,
# #            width = width,
# #            main=MAIN,
# #            annotation_col = anno)
# #   
# # }
# # 
# # imperfect_clusters <- cluster_qual %>%
# #   filter(maxmin_divmin > 0) %>%
# #   pull(secondary_cluster)
# # 
# # lapply(imperfect_clusters, gene_by_island_heatmap, gpa_clust = gpa_clust)
# # # dev.off()
# # 
# # 
# # ########
# # # gene_by_island_heatmap(gpa_clust = gpa_clust, SECONDARY_CLUSTER = 120)
# # 
# # 
# # ## secondary cluster by genome heatmaps here
# # anno <- clust_info %>%
# #   select(genome_name, secondary_cluster, island_type) %>%
# #   group_by(secondary_cluster) %>%
# #   summarise(consensus_type=paste(unique(island_type), sep = '~', collapse = '~')) %>%
# #   left_join(cluster_qual) %>%
# #   transmute(secondary_cluster = secondary_cluster,
# #             # primary_cluster = factor(primary_cluster),
# #             cluster_variability = maxmin_divmin,
# #             consensus_type = consensus_type) %>%
# #   
# #   column_to_rownames(var = 'secondary_cluster')
# # 
# # 
# # anno_length <- max(nchar(anno$consensus_type))/4
# # 
# # # this produces a heat map of the number of times islands from each secondary cluster show up in each genome
# # sec_clust_by_genome <- clust_info %>%
# #   select(genome_name, secondary_cluster) %>%
# #   group_by(genome_name, secondary_cluster) %>%
# #   tally() %>%
# #   spread(key=genome_name, value = n, fill = 0) %>%
# #   column_to_rownames(var = 'secondary_cluster') %>%
# #   as.matrix()
# # 
# # 
# # width=(ncol(sec_clust_by_genome)+anno_length)/4
# # if (width < 6){
# #   width <- 6
# # }
# # 
# # height=nrow(sec_clust_by_genome)/10
# # 
# # if (height < 6){
# #   height <- 6
# # }
# # 
# # 
# # pheatmap(sec_clust_by_genome,
# #          annotation_row = anno,
# #          height=height,
# #          width=width,
# #          filename = './gifrop_out/secondary_clusters_by_genome.jpeg',
# #          main='Presence of genomic island clusters in each genome',
# #          sub='highly variable clusters are indicated in green',
# #          fontsize = 5)
# # 
# # # dev.off()
# # 
# # 
