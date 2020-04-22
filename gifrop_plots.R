### PLOTS! ###

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

### FOR TESTING ONLY ###
# setwd('~/Documents/IslandR/test_data/pan')

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

clust_info <- read_csv('clustered_island_info.csv',
                       col_types = c('ccccccddddcddlcccccc'))


### USE ME
p <- clust_info %>%
  ggplot(aes(x=island_length, fill=island_type)) +
  geom_histogram(bins = 50, color='black') + theme_bw() +scale_fill_brewer(palette = 'Set1')+
  ggtitle('Histogram of the lengths of all detected genomic islands')


png('./gifrop_out/figures/island_length_histogram.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()



####  USE ME? ####

islands_per_isolate <- clust_info %>%
  group_by(genome_name) %>%
  tally()

width <- (nrow(islands_per_isolate)/10) *120

if (width < 600){
  width <- 600
}



p <- clust_info %>%
  group_by(genome_name, island_type) %>%
  tally() %>%
  ggplot(aes(x=genome_name, y=n, fill=island_type)) +
  geom_col(color='black') +
  scale_fill_brewer(palette = 'Set1') + theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('number of genomic islands per isolate (All Islands)')


png('./gifrop_out/figures/islands_per_isolate.png', width=width, height=600, res=120, type="cairo")
print(p)
dev.off()

####  Excluding unknowns ###

p <- clust_info %>% filter(island_type != 'unknown') %>%
  group_by(genome_name, island_type) %>%
  tally() %>%
  ggplot(aes(x=genome_name, y=n, fill=island_type)) +
  geom_col(color='black') +
  scale_fill_brewer(palette = 'Set1') + theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('number of genomic islands per isolate (All Islands)')


png('./gifrop_out/figures/islands_per_isolate_no_unknowns.png', width=width, height=600, res=120, type="cairo")
print(p)
dev.off()



# #### Now limiting islands by length for the next few ##
# p <- clust_info %>% filter(island_length > 5000) %>%
#   group_by(genome_name, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('number of genomic islands per isolate (>5kb)')
#
#
# width <- (length(p$data$genome_name)/10) *120
#
# if (width < 600){
#   width <- 600
# }
#
# png('./figures/islands_per_isolate_5kb.png', width=width, height=600, res=120, type="cairo")
# print(p)
# dev.off()
#
#
# p <- clust_info %>% filter(island_length > 10000) %>%
#   group_by(genome_name, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('number of genomic islands per isolate (>10kb)')
#
# width <- (length(p$data$genome_name)/10)*120
# if (width < 600){
#   width <- 600
# }
#
#
# png('./figures/islands_per_isolate_10kb.png', width=width, height=600, res=120, type="cairo")
# print(p)
# dev.off()
#
#
#
#
#
# p <- clust_info %>% filter(island_length > 20000) %>%
#   group_by(genome_name, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('number of genomic islands per isolate (>20kb)')
#
# width <- (length(p$data$genome_name)/9)*120
# if (width < 600){
#   width <- 600
# }
#
#
# png('./figures/islands_per_isolate_20kb.png', width=width, height=600, res=120, type="cairo")
# print(p)
# dev.off()
#
#
#
#
#
# p <- clust_info %>% filter(island_length > 40000) %>%
#   group_by(genome_name, island_type) %>%
#   tally() %>%
#   ggplot(aes(x=genome_name, y=n, fill=island_type)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90)) +
#   ggtitle('number of genomic islands per isolate (>40kb)')
#
# width <- (length(p$data$genome_name)/7)*120
# if (width < 600){
#   width <- 600
# }
#
#
# png('./figures/islands_per_isolate_40kb.png', width=width, height=600, res=120, type="cairo")
# print(p)
# dev.off()
#



####
p <- clust_info %>%
  mutate(primary_cluster = factor(primary_cluster)) %>%
  group_by(primary_cluster, island_type) %>%
  tally() %>%
  ggplot(aes(x=primary_cluster, y=n, fill=island_type)) +
  geom_col() +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Number of occurances of each genomic island cluster in the pangenome')


png('./gifrop_out/figures/Number_of_occurances.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()


###############

p <- clust_info %>%
  mutate(secondary_cluster = factor(secondary_cluster)) %>%
  group_by(secondary_cluster, island_type) %>%
  tally() %>%
  ggplot(aes(x=secondary_cluster, y=n, fill=island_type)) +
  geom_col(color='black') +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Number of occurances of each genomic island cluster in the pangenome') +
  scale_fill_brewer(palette = 'Set1')


png('./gifrop_out/figures/Number_of_occurances_secondary.png', width=840, height=600, res=120, type="cairo")
print(p)
dev.off()

