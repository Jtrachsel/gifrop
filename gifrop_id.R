# deal with command line args
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
min_genes <- as.numeric(args[2])
flankingDNA <- as.numeric(args[3])

# type check of args
stopifnot(exprs =
            {is.numeric(min_genes)
             is.numeric(flankingDNA)
             min_genes >= 1}
          )

### FOR TESTING ONLY ####
# setwd('./test_data3/')
# min_genes <- 4
# flankingDNA <-1000
# getwd()

###
# load packages
print('loading packages')
suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(purrr, quietly = TRUE, warn.conflicts = FALSE))
print('done loading packages')
# 
# # Functions #


### return orders for locus tags along seqids ###
# filters gffs to only contain CDS type features (just like roary does)
# then determines the order of these remaining locus tags for each seqid
# takes a gff produced by prokka and parsed by parse_gff3()
get_loc_tag_order <- function(gff){
  loc_tag_orders <- 
    gff %>%
    filter(type == 'CDS') %>% 
    mutate(num_loc_tag=as.numeric(sub('(.*)_([0-9]+)','\\2',locus_tag))) %>%
    arrange(seqid, num_loc_tag) %>%
    mutate(loc_tag_order=seq_along(num_loc_tag)) 
  return(loc_tag_orders)
}

### remove_core_genome
# takes a roary gene_presence_absence.csv file
# removes the core genome so that genomic islands can be determined from the remaining
# locus tags.
# keeps genes that dont occur in every isolate
# keeps genes that occur more than once within any genome
remove_core_genome <- function(roary_gpa){
  
  tot_isolates <- ncol(roary_gpa) - 14
  not_present_in_all <- roary_gpa$`No. isolates` < tot_isolates
  
  keepers <- roary_gpa$`Avg sequences per isolate` > 1 | not_present_in_all
  access_frags <- roary_gpa[keepers,]
  ####
  # this step will remove genomes from the analysis which do not
  # contain any non-core genes
  ISNA <- apply(access_frags, 1, is.na)
  ISNA <- rowSums(!ISNA) < 0      # this remove genomes with no non-core genes.
  ISNA[1:15] <- FALSE  # hacky way around this issue of removing columns that are not genomes....
  access_frags <- access_frags[,!ISNA]
  return(access_frags)
}


# ID_islands should be called in a mutate(map()) situation on the nested dataframe
ID_islands <- function(datfrm, min_genes){
  # island_names <- paste(datfrm)
  step0 <- datfrm$loc_tag_order
  names(step0) <- datfrm$locus_tag
  step1 <- sort(step0) # vector
  step2 <- split(step1, cumsum(c(1, diff(step1) != 1))) # list
  step3 <- step2[lapply(step2, length) >= min_genes] # list
  if (length(step3) == 0){
    # mess <- paste('NO ISLANDS LONGER THAN', min_genes, 'GENES FOUND IN seqid', seqid)
    # print(mess)
    return(NULL)
  } else {
    #need to swap orders back for locus tags
    #locus tags are names of the orders
    # browser()
    RESULTS <- lapply(step3, names)
  }
  return(RESULTS)
}



get_islands <- function(island_info, genome){
  # using the island info dataframe, this function
  # extracts the sequence data for each island from the corresponding
  # fasta.
  locs <- island_info %>% select(seqid, Istart, Iend) %>%
    transmute(chrom=seqid, start=Istart, end=Iend)

  islands <- BSgenome::getSeq(genome, as(locs, "GRanges")) # this needs Biostrings and BSgenome loaded
  names(islands) <- island_info$island_ID

  return(islands)
}

#### another gff parsing option ####

gff_parse3 <- function(path){
  # only excludes 'gene' type annotations, might be better than to only allow CDS and trna...
  # added mutate statement to remove extra stuff that sometimes comes along with the locus tag
  gff <- read_delim(path,
                    delim = '\t',
                    col_names = c("seqid", "source", "type", "start", "end", "score", "strand","phase","attributes"),
                    comment = '#', progress = FALSE, col_types = c('cccddcccc')) %>%
    filter(type != 'gene') %>%
    tidyr::extract(attributes,
                   into = c('ID', 'locus_tag', 'product'),
                   regex ='ID=(.*);.*locus_tag=(.*_[0-9]+);.*product=(.*)',
                   remove = FALSE) %>%
    mutate(locus_tag = sub('([A-Za-z]_[0-9]+).*', '\\1', locus_tag))
  return(gff)
}


# helper function to get island IDs associated with locus_tags
enframe_island_list <- function(island_list){
  tibble::enframe(island_list,name='island_id', value='locus_tag') %>% 
    unnest(cols = locus_tag)
}


# change this to select 10? flanking genes on either side, then can cluster spots of insertion
# using same clustering alg as islands themselves
ID_flanking_genes <- function(datfrm, seqid_loc_tags){
  extract_these <- 
    datfrm %>% 
    summarise(flank_low_ord=min(loc_tag_order)-1, 
              flank_hig_ord=max(loc_tag_order)+1)
  # browser()
  loc_tag_low <- seqid_loc_tags$locus_tag[which(seqid_loc_tags$loc_tag_order == extract_these$flank_low_ord)]
  loc_tag_hig <- seqid_loc_tags$locus_tag[which(seqid_loc_tags$loc_tag_order == extract_these$flank_hig_ord)]
  
  # sometimes there are not flanking genes (contig boundaries, plasmids etc)
  if (identical(loc_tag_low, character(0))){
    loc_tag_low <- 'none'
  }
  if (identical(loc_tag_hig, character(0))){
    loc_tag_hig <- 'none'
  }
  
  return(paste(loc_tag_low, loc_tag_hig, sep = '|'))
  
}


# done with functions #



##### read in files #

current_directory <- getwd()
gff_files <- list.files(path = './gifrop_out/sequence_data/', pattern = 'short.gff', full.names = TRUE)



# This creates a vector of column specifications to be passed to the read_csv function
# I had trouble with some of the locus tag column types being guessed as logical

pan_cols <- c('ccciidiiiiciii')
locus_tag_cols <- rep_len('c', length(gff_files)) %>%
  paste(sep = '', collapse = '')
all_cols <- paste(pan_cols, locus_tag_cols, sep = '', collapse = '')

gpa <- read_csv('./gene_presence_absence.csv', col_types = all_cols)


# read in gffs
print('reading in gffs...')
# parallel this? meh.
gffs <- lapply(gff_files, gff_parse3)

gff_names <- sub('./gifrop_out/sequence_data/(.*)_short.gff','\\1',gff_files)
gff_names <- gsub('/?','',gff_names)
names(gffs) <- gff_names


# this is an important line
# helps solve issue #1, where some islands are inappropriately interrupted
# because raw locus tags were being converted to numeric.
# now, locus tags are sorted according to the converted numeric, then 
# a new order is assigned, this way there are no gaps caused by 
# tRNA rRNA type locus tags being filtered out.
loc_tag_orders <- bind_rows(lapply(gffs, get_loc_tag_order))



# this tibble is to help determine if an island occupies the whole contig
# also used to extract flanking locus tags ( core genes that border islands)
loctags_on_seqids <- 
  tibble::enframe(gffs, name = 'genome', value='gff_df') %>% 
  unnest(cols = gff_df) %>% 
  select(genome, seqid, locus_tag) %>% 
  left_join(loc_tag_orders) %>% 
  select(genome, seqid, locus_tag,loc_tag_order) %>% 
  nest(seqid_loc_tags=c(locus_tag, loc_tag_order)) %>% 
  ungroup() %>% 
  select(genome, seqid, seqid_loc_tags)


###

genome_filenames <- list.files(path = './gifrop_out/sequence_data', pattern = '.fna', full.names = TRUE)

print('reading in fastas...')
genome_seqs_list <- sapply(genome_filenames, Biostrings::readDNAStringSet)
print('done reading in fastas')

# makes a dataframe of seqid and seqid length
seq_lens <- tibble(seqid=unlist(lapply(genome_seqs_list, names)),
                   seqid_len=unlist(lapply(genome_seqs_list, Biostrings::width)))

# identify islands in the context of the pangenome
islands_pangenome_gff <- 
  remove_core_genome(gpa) %>%    # removes genes found in all genomes
  pivot_longer(cols=-c(1:14),
               names_to = 'genome',
               values_to='locus_tag', 
               values_drop_na = TRUE) %>% 
  separate_rows(locus_tag, sep = '\t') %>% # helps solve single pangenome genes with tab separated locus tags
  left_join(loc_tag_orders) %>%            # brings in locus_tag orders, no longer converting locus tags to numeric
  group_by(genome, seqid) %>%              
  nest() %>% 
  mutate(ISLANDS=map(.x = data, .f = ID_islands, min_genes = min_genes)) %>% 
  filter(!map_lgl(ISLANDS, is.null)) %>%              # NULL values generated when min_genes removes small islands
  mutate(island_IDs=map(.x = ISLANDS, enframe_island_list)) %>% # converts the detected island list to a tibble
  select(-ISLANDS) %>% 
  mutate(newdat=map2(.x = data, .y=island_IDs, .f = left_join)) %>% 
  select(genome, seqid, newdat) %>% 
  unnest(cols = newdat) %>% 
  # mutate(names(loc_tag_order)=locus_tag) %>% pull(loc_tag_order)
  filter(!is.na(island_id)) %>% 
  ungroup() %>% 
  mutate(island_ID=paste(seqid, island_id, sep = '_')) %>% 
  select(island_ID, everything(), -island_id) %>% 
  group_by(island_ID, genome, seqid) %>%
  nest() %>% 
  mutate(island_loc_tags=map(.x = data, .f = pull, locus_tag)) %>% 
  left_join(loctags_on_seqids) %>% 
  left_join(seq_lens) %>% 
  mutate(only_island=map2_lgl(.x = island_loc_tags, .y = seqid_loc_tags, .f= ~ all(.y$locus_tag%in%.x)), 
         flanking_genes=map2_chr(.x=data, .y=seqid_loc_tags, .f=ID_flanking_genes)) %>% 
  select(island_ID, genome,seqid, data, only_island, seqid_len, flanking_genes) %>% 
  unnest(cols = data) %>% 
  ungroup() 

islands_pangenome_gff %>% write_csv('./gifrop_out/islands_pangenome_gff.csv')

# creates the island_info summary of the islands
island_info <- 
  islands_pangenome_gff %>% 
  group_by(island_ID) %>%
  summarise(seqid=unique(seqid),
            genome_name=unique(genome),
            start=min(start),
            end=max(end),
            island_length=max(end)-min(start),
            num_genes=length(locus_tag),
            locus_tags=paste(locus_tag, collapse = '|'),
            seqid_len = unique(seqid_len),
            percent_island = island_length/unique(seqid_len),
            only_island = unique(only_island),
            Istart=case_when(
              only_island ~ 1, 
              !only_island ~ min(start)),
            Iend=case_when(
              only_island ~ as.numeric(unique(seqid_len)),   # was having int vs double conflicts here
              !only_island ~ max(end)),
            acc_frag = paste(unique(`Accessory Fragment`), collapse = '|'), 
            genes=paste(unique(Gene), collapse = '|'), 
            flank_loc_tags=paste(unique(flanking_genes), collapse = '|')) 


island_info %>% 
  select(-Istart, -Iend) %>% 
  write_csv('./gifrop_out/my_islands/island_info.csv')

# writing out the fastas for the islands
# creates a tibble containing fasta files for each genome name
genome_seqs <- tibble(genome_name=sub('.*/gifrop_out/sequence_data/(.*).fna','\\1',names(genome_seqs_list)),
                      genome=genome_seqs_list)

nesty_res <- island_info %>%
  group_by(genome_name) %>%
  nest() %>%
  ungroup() %>%
  left_join(genome_seqs)

print('Extracting island fastas... this can take a second...')
nesty_res <- nesty_res %>%
  mutate(islands_mon = map2(.x = data, .y = genome, .f = get_islands))

all_islands <- Biostrings::DNAStringSetList(nesty_res$islands_mon)
all_islands <- unlist(all_islands)

# I could probably just have abricate run on this fasta
# instead of splitting everything out into their own fastas and running abricate 
# on all those files independently...
Biostrings::writeXStringSet(all_islands, './gifrop_out/my_islands/All_islands.fasta')



#IF BLOCK HERE?
# #### OUTPUT GFFS HERE ####
# Gffs can now be re-constructed from the islands_pangenome_gff.csv file


# get_gffs <- function(island_info, gff){
#   # using the island info dataframe, this function
#   # extracts the annotation data for each island from the corresponding
#   # gff
#   
#   # uses helper function, get_island_gff
#   island_gffs <- map(island_info$locus_tags, .f = get_island_gff, gff)
#   names(island_gffs) <- island_info$island_ID
#   return(island_gffs)
# }
# 
# 
# get_island_gff <- function(locus_tag_vec, gff){
#   #this little guy should take the locus tag column from the island info dataframe
#   # and return the corresponding annotations from the gff file for the genome,
#   # basically a little mini gff that only covers the island
#   gff[gff$locus_tag %in% locus_tag_vec,]
#   
# }

# # now need way to incorporate ID_non_island_loc_tags
# print('collecting island gffs')
# 
# 
# # # fix extra info on locus tags...
# #moved this up
# # res_4_real <- bind_rows(gffs) %>%
# #   group_by(seqid) %>%
# #   select(seqid, locus_tag) %>%
# #   nest() %>%
# #   mutate(seqid_loc_tags=list(unlist(data, use.names = FALSE))) %>% 
# #   select(-data) %>%
# #   right_join(res_4_real) %>%
# #   mutate(only_island=map2_lgl(.x = seqid_loc_tags, .y = locus_tags, .f = ~ all(.x %in% .y)))
# 
# 
# gff_table <- island_info %>%
#   group_by(genome_name) %>%
#   nest() %>% ungroup() %>%
#   left_join(tibble(genome_name=names(gffs), gffs=gffs)) # this tibble is just genome names and gffs
# 
# ###### THIS ISNT WORKING ####
# gff_table <- gff_table %>%
#   mutate(island_gffs = map2(.x = data, .y = gffs, .f = get_gffs))
# 
# ## I SHOULD JUST WRITE OUT PANGENOME_GFF???
# 
# 
# all_island_gffs <- unlist(gff_table$island_gffs, recursive = FALSE)
# 
# # names(all_island_gffs)
# 
# 
# invisible(sapply(names(all_island_gffs),
#        function (x) write_tsv(all_island_gffs[[x]], path = paste('./gifrop_out/my_islands/', x, ".gff", sep=""))))
# 
# 
# print('done writing island gffs')
# END IF BLOCK



### ADJUST START AND END COORDS HERE ###
# might want to take leading and trailing DNA sequences before and after coding sequences
# could look at insertion sites etc
if(flankingDNA > 0){

  print('ALSO OUTPUTTING ISLANDS WITH FLANKING DNA')
  print(paste('taking', flankingDNA, 'bp on either side of the island'))

  #EXTEND CHROMOSOMAL ISLANDS STARTS AND ENDS TO flankingDNA_len
  island_info <- island_info %>%
    mutate(Istart=Istart-flankingDNA,
           Iend=Iend + flankingDNA,
           Istart=ifelse(Istart < 1, 1, Istart),
           Iend= ifelse(Iend > seqid_len, seqid_len, Iend))


  system('mkdir ./gifrop_out/my_islands/with_flanking')

  nesty_res <- island_info %>%
    group_by(genome_name) %>%
    nest() %>%
    ungroup() %>%
    left_join(genome_seqs)

  nesty_res <- nesty_res %>%
    mutate(islands_mon = map2(.x = data, .y = genome, .f = get_islands))
  all_islands <- Biostrings::DNAStringSetList(nesty_res$islands_mon)
  all_islands <- unlist(all_islands)
  Biostrings::writeXStringSet(all_islands, './gifrop_out/my_islands/with_flanking/All_islands_with_flanking.fasta')

}
