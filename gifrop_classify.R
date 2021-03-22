args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

# CHANGE THIS TO GIFROP_CLASSIFY.R
suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(purrr, quietly = TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(igraph))



### ONLY FOR HERE FOR TESTING ###
# setwd('/home/julian/Documents/gifrop_examples/test5')
# setwd('/project/fsep_004/jtrachsel/klima/assembly/both/second_flye_polish/pananal/plasmids/pan/')
# setwd('/home/julian/Documents/gifrop/test_data/pan')
#getwd()

## read in island info data ##
island_info <- read_csv('./gifrop_out/my_islands/island_info.csv', col_types = c('cccddddcddlccc'))


# current_directory <- getwd()
# seq_dat_path <- paste0(current_directory, '/gifrop_out/sequence_data/')
gff_files <- list.files(path = './gifrop_out/sequence_data/', pattern = 'short.gff', full.names = TRUE)

# This creates a vector of column specifications to be passed to the read_csv function
# I had trouble with some of the locus tag column types being guessed as logical
pan_cols <- c('ccciidiiiiciii')
locus_tag_cols <- rep_len('c', length(gff_files)) %>% paste(sep = '', collapse = '')
all_cols <- paste(pan_cols, locus_tag_cols, sep = '', collapse = '')

# this is where the gene, presence/absense is read in
gpa <- read_csv('./gene_presence_absence.csv', col_types = all_cols)

# READ IN ABRICATE RESULTS #
pwd <- getwd()
print('reading in abricate files...')

# read in plasmidfinder results and bind together
plasfiles <- list.files(path = './gifrop_out/my_islands/abricate/', pattern = 'plasmidfinder', full.names = TRUE)
plasfinders <- lapply(plasfiles, read_tsv, col_types = c('ccddcccccddcccc')) #
plasfinders <- bind_rows(plasfinders)

# concatenate all plasmid genes and produce plasmid type per island
plasmid_types <- plasfinders %>%
  mutate(gene_percent=paste(GENE,'<',  `%COVERAGE`,'%', '>', sep = '')) %>%
  group_by(SEQUENCE) %>%
  summarise(plasmid_type=paste(gene_percent, sep = '~', collapse = '~'), 
            .groups='drop') %>%
  transmute(island_ID=SEQUENCE, plasmid_type=plasmid_type)

# read in vfdb results and bind together
vfdbfiles <- list.files(path = './gifrop_out/my_islands/abricate/', pattern = 'vfdb', full.names = TRUE)
vfdbs <- lapply(vfdbfiles, read_tsv, col_types = c('ccddcccccddcccc'))
vfdbs <- bind_rows(vfdbs)

# concatenate all vfdb genes and produce plasmid type per island
vir_types <- vfdbs%>%
  mutate(gene_percent=paste(GENE,'<',  `%IDENTITY`,'%', '>', sep = '')) %>%
  group_by(SEQUENCE) %>%
  summarise(vir_type=paste(gene_percent, sep = '~', collapse = '~'), 
            .groups='drop') %>%
  transmute(island_ID=SEQUENCE, vir_type=vir_type)

# read in ncbi resistance results 
resfiles <- list.files(path = './gifrop_out/my_islands/abricate/', pattern = 'ncbi', full.names = TRUE)
resfinders <- lapply(resfiles, read_tsv, col_types = c('ccddcccccddcccc'))
resfinders <- bind_rows(resfinders) %>% filter(`%COVERAGE` > 66)

# concatenate all resfinder genes and produce res_type for all islands
res_types <- resfinders%>%
  mutate(gene_percent=paste(GENE,'<',  `%IDENTITY`,'%', '>', sep = '')) %>%
  group_by(SEQUENCE) %>%
  summarise(res_type=paste(gene_percent, sep = '~', collapse = '~'), 
            .groups='drop') %>%
  transmute(island_ID=SEQUENCE, res_type=res_type)

# virotypes


virofiles <- list.files(path = './gifrop_out/my_islands/abricate/', pattern = 'viroseqs', full.names = TRUE)
virofinders <- lapply(virofiles, read_tsv, col_types = c('ccddcccccddcccc')) # check this coltypes
virofinders <- bind_rows(virofinders) %>% filter(`%COVERAGE` > 66)

# concatenate all resfinder genes and produce viro_type for all islands
viro_types <- virofinders%>%
  mutate(gene_percent=paste(GENE,'<',  `%IDENTITY`,'%', '>', sep = '')) %>%
  group_by(SEQUENCE) %>%
  summarise(viro_type=paste(gene_percent, sep = '~', collapse = '~'), 
            .groups='drop') %>%
  transmute(island_ID=SEQUENCE, viro_type=viro_type)


# megares/bacmet : metals and biocides

megares_files <- list.files(path = './gifrop_out/my_islands/abricate/', pattern = 'megares', full.names = TRUE)

megares <- bind_rows(lapply(megares_files, read_tsv, col_types = c('ccddcccccddcccc'))) %>%
      filter(`%COVERAGE` > 66)

megares_types <- megares%>%
      mutate(gene_percent=paste(GENE,'<',  `%IDENTITY`,'%', '>', sep = '')) %>%
      group_by(SEQUENCE) %>%
      summarise(megares_type=paste(gene_percent, sep = '~', collapse = '~')) %>%
      transmute(island_ID=SEQUENCE, megares_type=megares_type)



#
# might need to check in on things here # in case island_IDs got garbled somehow
allbricates <- bind_rows(plasfinders, vfdbs, resfinders, virofinders, megares) %>%
  mutate(island_ID=SEQUENCE) %>%
  select(island_ID, everything(), -SEQUENCE)





print('Done reading in abricate files')


# this block creates a 'resistance type' by concatenating all the detected resistances into a string.
res_info <- allbricates %>%
  filter(!is.na(RESISTANCE)) %>% group_by(island_ID) %>%
  mutate(RESISTANCE=paste(unique(RESISTANCE), collapse = '|', sep = '|')) %>%
  select(island_ID, RESISTANCE) %>% unique()



# this block creates a broad 'island type'
# if an island has a hit to one of the five database types it gets assigned that type
# all types are then concatenated to form the final island type
island_types <- allbricates %>% 
  select(island_ID, DATABASE) %>% 
  unique() %>%
  mutate(DATABASE=case_when(
    DATABASE == 'plasmidfinder' ~ 'plasmid',
    DATABASE == 'ncbi'          ~ 'AMR',
    DATABASE == 'vfdb'          ~ 'virulence',
    DATABASE == 'PHAGE'         ~ 'phage',
    DATABASE == 'BacMet'        ~ 'metals/biocides')) %>%
  group_by(island_ID) %>%
  mutate(island_type = paste(DATABASE, sep = '_', collapse = '_')) %>%
  ungroup() %>% select(island_ID, island_type) %>%
  unique() %>%
  left_join(res_info)

island_info %>%
  left_join(island_types) %>% 
  left_join(res_types) %>%
  left_join(vir_types) %>%
  left_join(plasmid_types) %>%
  left_join(viro_types) %>%
  left_join(megares_types) %>%
  write_csv('./gifrop_out/classified_island_info.csv')

print('Done with island classification')

### END READ IN ABRICATE STUFF ###




### TEST ZONE # COMMENT OUT
# Island PA heatmap? #


# pseudo #

# for all secondary clusters that have variability in the number of genes
    # make
      # heatmap where rows are genes and columns are individual islands
      # histogram of number of genes
## this dataframe tries to assess how consistent all islands within a cluster are
# looks at how variable the number of genes is for all islands in a cluster
# it seems like ((max_genes - min_genes) / min_genes) gives a pretty good indication
# of variability in the cluster
# THESE ARE THE CLUSTER QUALITIES OF THE quat_clusters NOW
# cluster_qual <-
#   clust_info %>% group_by(primary_cluster, secondary_cluster, tertiary_cluster, quat_cluster) %>%
#   summarise(mean_genes = mean(num_genes),
#             med_genes  = median(num_genes),
#             var_genes  = var(num_genes),
#             sd_genes   = sd(num_genes),
#             min_genes  = min(num_genes),
#             max_genes  = max(num_genes),
#             num_occur  = length(num_genes),
#             maxmin_divmin   = (max_genes - min_genes)/min_genes)
# 
# 
# 
# # variable_clusters <- cluster_qual %>% filter(maxmin_divmin != 0)
# 
# 
# filt_helper <- function(test_vec, int_vec){
#     res <- any(int_vec %in% test_vec)
#     return(res)
# }





#######
# This might be nice but is broken when there are some islands with only 1 gene on them.
# this makes a heatmap of all the genes in a specified cluster (or vector of clusters) in each ISLAND
# 
# gene_by_island_heatmap <- function(gpa_clust, QUAT_CLUST){
#   # gpa clust need to be the gene_presence_absence.csv file with the added clustering information
#   # in addition the 'secondary_cluster' column needs to be a list formatted column
#   # secondary cluster can be a vector of secondary clusters you want to see together
#   # ISSUE!! currently genes from other sclusts are pulled in, if a gene is in two different sclusts
#   # it will show up in these heatmaps but will not have annotation info with it.
#   # need to include another filtering step to remove islands not belonging to Sclust at hand
# 
# 
#   anno <- clust_info %>%
#     filter(quat_cluster %in% QUAT_CLUST) %>%
#     select(island_ID, island_type) %>%
#     column_to_rownames(var = 'island_ID')
# 
#   PA <- gpa_clust %>%
#     filter(map_lgl(.x=all_Qclusters, .f=filt_helper, QUAT_CLUST)) %>%
#     unnest(cols = all_islands) %>%
#     select(Gene, all_islands) %>%
#     group_by(all_islands, Gene) %>%
#     tally() %>%
#     ungroup() %>%
#     spread(key = Gene, value = n, fill = 0) %>%
#     column_to_rownames(var = 'all_islands') %>%
#     as.matrix() %>%
#     t()   # to get genes as rows and cols as islands
# 
#   QUAT_CLUST <- paste(QUAT_CLUST, collapse = '_', sep = '_')
#   MAIN=paste('Presence/absence of genes among islands in quaternary cluster',QUAT_CLUST)
#   FILENAME=paste('./gifrop_out/Qclust_', QUAT_CLUST,'_gene_heatmap.jpeg', sep = '')
# 
#   width=ncol(PA)/5
#     if (width < 6){
#     width <- 6
#   }
# 
#   height=nrow(PA)/10
# 
#   if (height < 6){
#     height <- 6
#   }
# 
# 
#   pheatmap(PA, filename = FILENAME,
#            height = height,
#            width = width,
#            main=MAIN,
#            annotation_col = anno)
# 
# }
# 
# imperfect_clusters <- cluster_qual %>%
#   filter(maxmin_divmin > 0) %>%
#   pull(quat_cluster)
# 
# lapply(imperfect_clusters, gene_by_island_heatmap, gpa_clust = gpa_clust)
# dev.off()

# gene_by_island_heatmap(gpa_clust = gpa_clust, SECONDARY_CLUSTER = 9)


########
# gene_by_island_heatmap(gpa_clust = gpa_clust, SECONDARY_CLUSTER = 120)

# THIS NEEDS TO BE UPDATED TO QUAT CLUSTERS
## secondary cluster by genome heatmaps here
# anno <- clust_info %>%
#   select(genome_name, quat_cluster, island_type) %>%
#   group_by(quat_cluster) %>%
#   summarise(consensus_type=paste(unique(island_type), sep = '~', collapse = '~')) %>%
#   left_join(cluster_qual) %>%
#   transmute(quat_cluster = quat_cluster,
#            # primary_cluster = factor(primary_cluster),
#             cluster_variability = maxmin_divmin,
#             consensus_type = consensus_type) %>%
# 
#   column_to_rownames(var = 'quat_cluster')
# 
# 
# anno_length <- max(nchar(anno$consensus_type))/4
# 
# # this produces a heat map of the number of times islands from each secondary cluster show up in each genome
# quat_clust_by_genome <- clust_info %>%
#   select(genome_name, quat_cluster) %>%
#   group_by(genome_name, quat_cluster) %>%
#   tally() %>%
#   spread(key=genome_name, value = n, fill = 0) %>%
#   column_to_rownames(var = 'quat_cluster') %>%
#   as.matrix()
# 
# 
# width=(ncol(quat_clust_by_genome)+anno_length)/4
# if (width < 6){
#   width <- 6
# }
# 
# height=nrow(quat_clust_by_genome)/10
# 
# if (height < 6){
#   height <- 6
# }
# 
# 
# pheatmap(quat_clust_by_genome,
#          annotation_row = anno,
#          height=height,
#          width=width,
#          filename = './gifrop_out/quat_clusters_by_genome.jpeg',
#          main='Presence of genomic island clusters in each genome',
#          sub='highly variable clusters are indicated in green',
#          fontsize = 5)

# dev.off()
# 
# # COPY FOR TERTIARY
# # COPTY FOR QUATERNARY
# 
# 
