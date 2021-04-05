library(ggplot2)
library(cowplot)
library(tidyverse)

cii <- read_csv('/home/Julian.Trachsel/Documents/gifrop/test_data4/gifrop_out/clustered_island_info.csv')
pgff <- read_csv('/home/Julian.Trachsel/Documents/gifrop/test_data4/gifrop_out/islands_pangenome_gff.csv')


cii %>%
  ggplot(aes(x=log(num_genes), y=factor(primary_cluster), fill=island_type)) +
  geom_point(shape=21, color='white', size=3) + 
  theme_cowplot() + 
  theme(panel.grid.major = element_line(color='grey'))



cii %>% 
  filter(!is.na(RESISTANCE)) %>% 
  mutate(num_res=map_dbl(.x=strsplit(RESISTANCE, split='|', fixed = T), .f=length)) %>%# pull(num_res)
  ggplot(aes(x=num_res, y=num_genes, color=island_type)) + geom_point() + xlim(0,7)


cii8 <- cii %>% filter(primary_cluster == 8) 


# what's the 'core genome' of these islands?
cii8 %>% 
  ggplot(aes(x=log(num_genes), y=factor(secondary_cluster), fill=island_type)) +
  geom_point(shape=21, color='white', size=3) + 
  theme_cowplot() + 
  theme(panel.grid.major = element_line(color='grey'))


cii8_PA <- 
  cii8 %>%
  select(island_ID, genes) %>% 
  separate_rows(genes, sep = '\\|') %>% 
  mutate(present=1) %>% 
  pivot_wider(names_from = genes, values_from=present, values_fill=0)

#

cii8_matrix <- cii8_PA %>% column_to_rownames(var='island_ID') %>% as.matrix()

#

#

test <- cii8 %>%
  select(island_ID, genes) %>%
  separate_rows(genes, sep='\\|') %>% 
  group_by(island_ID) %>%
  nest() %>% 
  mutate(data=map(data, pull))

EDGE_LIST <- expand_grid(test, test,.name_repair ='universal') %>%
  filter(island_ID...1 != island_ID...3) %>% 
  mutate(in_common=map2_chr(.x = data...2, .y =data...4 , .f = ~ paste(.x[.x %in% .y] , collapse = '|')), 
         num_in_common=map2_dbl(.x = data...2, .y =data...4 , .f = ~ length(.x[.x %in% .y]))) %>% 
  select(-c(data...2, data...4)) %>% 
  transmute(from=island_ID...1, 
            to=island_ID...3, 
            in_common=in_common, 
            num_in_common=num_in_common) %>% 
  filter(num_in_common > 0)

library(igraph)


graph_from_edgelist(EDGE_LIST)

#


###

expand_grid(x = 1:3, y = 1:2)
expand_grid(l1 = letters, l2 = LETTERS)

# Can also expand data frames
expand_grid(df = data.frame(x = 1:2, y = c(2, 1)), z = 1:3)
# And matrices
expand_grid(x1 = matrix(1:4, nrow = 2), x2 = matrix(5:8, nrow = 2))


max(colSums(cii8_matrix))

cii8_matrix[,which(colSums(cii8_matrix) == 6)]



pgff %>% filter(Gene == 'group_5474')
