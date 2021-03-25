args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# adds island info to the roary pangenome files

gff_files <- list.files(path = './gifrop_out/sequence_data/', pattern = 'short.gff', full.names = TRUE)

# This creates a vector of column specifications to be passed to the read_csv function
# I had trouble with some of the locus tag column types being guessed as logical
pan_cols <- c('ccciidiiiiciii')
locus_tag_cols <- rep_len('c', length(gff_files)) %>% paste(sep = '', collapse = '')
all_cols <- paste(pan_cols, locus_tag_cols, sep = '', collapse = '')

# this is where the gene, presence/absense is read in
gpa <- read_csv('./gene_presence_absence.csv', col_types = all_cols)

island_info <- read_csv('./gifrop_out/clustered_island_info.csv', col_types = c('ccccccccddddcddlccccccccc'))


gpa_gath <- gpa %>%
  gather(key='genome', value='locus_tags', -(1:14)) %>%
  select(Gene,genome, locus_tags) #mmmk...

# this is a tibble of which locus tags are in each island (along with the island's length)
island_ID_loc_tags <- island_info %>%
  select(island_ID, locus_tags) %>%
  mutate(locus_tags = strsplit(x = locus_tags, split='|', fixed=TRUE)) %>%
  unnest(cols = locus_tags)
 
# making gpa with island clustering info
# 
gpa_clust <- gpa_gath %>%
  left_join(island_ID_loc_tags) %>%
  filter(!is.na(island_ID))%>%
  left_join(island_info, by='island_ID') %>%
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

# 
# 
# 
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
#







# some simple stats

# what proportion of the pangenome is contained within these islands?

# total genes in pangenome
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

# % of the accessory genome contained on these islands
print('Percent of accessory genome contained within these islands:')
(num_access_on_island / num_accessory_genes) * 100


