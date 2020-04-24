
# deal with command line args
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
min_genes <- as.numeric(args[2])
flankingDNA <- as.numeric(args[3])

# type check of args
stopifnot(exprs =
            {is.numeric(min_genes)
             is.numeric(flankingDNA)
             min_genes > 0}
          )


### FOR TESTING ONLY ####
# setwd('../IslandR/test_data/pan/')
# min_genes <- 4
# flankingDNA <- 0
# getwd()




###
# load packages
print('loading packages')
suppressPackageStartupMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(readr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(purrr, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(BSgenome, quietly = TRUE, warn.conflicts = FALSE))
print('done loading packages')

# Functions #

find_islands <- function(roary_gpa, min_genes = 4){
  # This function identifies genomic islands from roary output gene_presence_absence.csv
  # dataframes.

  # This step filters to only non core genes of the pangenome
  # I am suspicious that this is the best way of doing things....
  # I am allowing genes that occur more than 1 time in the same genome to be kept
  keepers <- gpa$`Avg sequences per isolate` > 1 | gpa$`No. isolates`< max(gpa$`No. isolates`)
  access_frags <- gpa[keepers,]
  ####
  # this step will remove genomes from the analysis which do not
  # contain any non-core genes
  ISNA <- apply(access_frags, 1, is.na)
  ISNA <- rowSums(!ISNA) < 0      # this remove genomes with no non-core genes.
  ISNA[1:15] <- FALSE  # hacky way around this issue of removing columns that are not genomes....
  access_frags <- access_frags[,!ISNA]


  # trying a big ol loop
  genomes <- colnames(access_frags)[15:length(colnames(access_frags))]
  master_res <- list()
  for (genome in genomes){
    # browser()
    # this should be a function?
    # for each column of locus tags (genome) in a dataframe,
    # change the locus tags to numeric vectors
    # sort them,
    # looks for stretches of consecutive numbers
    print(paste("finding islands for", genome))

    # If a gene is detected more than once in a genome roary will
    # put that locus tag in a tab delimited list within one of the csv cells
    # this block selects only the locus tag column from the genome at hand
    # and then separates the tab delimited values
    # this vector of locus tags is then used to identify genomic islands
    
    # this is inefficient because it does this select for each genome, 
    # should split into list of each genome's loc_tags and then apply subsequent function to it.
    step0 <- access_frags %>%
      select(all_of(genome)) %>%
      separate_rows(genome, sep = '\t') %>%
      na.omit()

    locus_tag_prefix <- unique(na.omit(sub('(.*_)(.*)', '\\1', step0[[genome]])))
    #
    # Maybe add if statement, if locus_tag_prefix
    #
    if (length(locus_tag_prefix) > 1){
      print(locus_tag_prefix)
      print(length(locus_tag_prefix))
      print('YOU SHOULD NEVER SEE THIS MESSAGE')
      print('MULTIPLE LOCUS TAG PREFIXES DETECTED')
      print('YOU NEED TO RE-ANNOTATE YOUR GENOMES')
      stop()
    }

    step1 <- sub('(.*_)(.*)', '\\2', step0[[genome]])        # converts loc.tags to numbers
    step2 <- as.numeric(step1)                               # could do in step 1...
    step3 <- sort(step2)                                     # orders loc tags for this genome
    step4 <- split(step3, cumsum(c(1, diff(step3) != 1)))    # identifies streches of sequential loc.tags
    #browser()
    step5 <- step4[lapply(step4, length) >= min_genes]  # this sets the minimum size of the islands to 4 genes
    if (length(step5) == 0){
      mess <- paste('NO ISLANDS LONGER THAN', min_genes, 'GENES FOUND IN GENOME', genome)
      print(mess)
    } else {
      names(step5) <- paste(genome, names(step5), sep = '_')
      island_list <- step5
      res <- lapply(island_list, restore_locus_tags, locus_tag_prefix)  # restores locus tags to original state
      master_res[[genome]] <- res

    }


  }
  return(master_res)
}

restore_locus_tags <- function(island, locus_tag_prefix){
  # This function is used to restore the original state of the locus tags.
  # The locus tags are stripped of their contig identifier during
  # the search for strings of consequtive locus tag numbers
  # and this puts them back the way they were
  tmp_res <- c()
  for (index in 1:length(island)){
    locus_tag <- island[[index]]
    z <- as.character(rep(0, 5-nchar(locus_tag)))
    zz <- paste(c(locus_tag_prefix, z, locus_tag), sep = '', collapse = '')
    tmp_res[index] <- zz
  }
  return(tmp_res)
}

get_islands <- function(island_info, genome){
  # using the island info dataframe, this function
  # extracts the sequence data for each island from the corresponding
  # fasta.
  locs <- island_info %>% select(seqid, Istart, Iend) %>%
    transmute(chrom=seqid, start=Istart, end=Iend)

  islands <- getSeq(genome, as(locs, "GRanges")) # this needs Biostrings and BSgenome loaded
  names(islands) <- island_info$island_ID

  return(islands)
}

# id_subcols <- function(attributes_column){
#   # was trying to make a better way of parsing gff files and
#   # extract subcolumns from the attribute column?
#
#   # should return a list of 2.
#   #[[1]] = text vector of column names to split attributes column into, fed to extract's "into"
#   #[[2]] = regex to extract these columns
# }


# this one only extracts ID locus_tag and product
# needs more testing.
# assumes that order is ID, locus_tag, product, with product at the very end.

# gff_parse2 <- function(path){
#   # used to parse the shortened prokka gff files.
#   # can take full gff (with sequence data appended)
#   # but is much slower
#   # ADDING tRNA TO TYPES ALLOWED BECAUSE
#   gff <- read_delim(path,
#                     delim = '\t',
#                     col_names = c("seqid", "source", "type", "start", "end", "score", "strand","phase","attributes"),
#                     comment = '#', progress = FALSE, col_types = c('cccddcccc')) %>%
#     filter(type %in% c('CDS', 'tRNA')) %>%
#     tidyr::extract(attributes,
#                    into = c('ID', 'locus_tag', 'product'),
#                    regex ='ID=(.*);.*locus_tag=(.*_[0-9]+);.*product=(.*)',
#                    remove = FALSE)
#   return(gff)
# }



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



######## this is where the gene, presence/absence is read in
gpa <- read_csv('./gene_presence_absence.csv', col_types = all_cols)


# ran into a weird thng where one of the isolates locus tag columns is given all NAs...
# really need to QC the inputs...



##### call to find islands function
print('Identifying genomic islands... (streches of consecutive locus tags')
master_res <- find_islands(roary_gpa = gpa)


if ( length(master_res) == 0){
  stop_message <- paste('NO ISLANDS LONGER THAN', min_genes, 'GENES FOUND IN ANY OF YOUR GENOMES')
  # does this throw an error?
  # should cause the main script to exit if this happens.
  stop(stop_message)

}

# read in gffs
print('reading in gffs...')
# parallel this? meh.
gffs <- lapply(gff_files, gff_parse3)

gff_names <- sub('./gifrop_out/sequence_data/(.*)_short.gff','\\1',gff_files)
gff_names <- gsub('/?','',gff_names)

names(gffs) <- gff_names


###### main loop ###
# should this be a function?
# it could be applied to each genome column of the gpa dataframe

res <- list()
for (genome_index in 1:length(master_res)){
  genome_name <- names(master_res[genome_index])
  genome <- master_res[[genome_index]]
  genome_res <- list()

  print(paste('working_on_genome_', genome_name))
  for (island_index in 1:length(genome)){
    island <- genome[[island_index]]
    island_name <- names(genome[island_index])
    tmp <- gffs[[genome_name]] %>% filter(locus_tag %in% island)
    if (nrow(tmp) == 0){
      messa <- paste('locus tags for island "', island_name, ' " not found in gffs...')

    } else {
      tmp2 <- tmp %>%
        group_by(seqid) %>%
        summarise(genome_name=genome_name,
                  start=min(start),
                  end=max(end),
                  island_number=unique(group_indices()), # this helps discriminate when two islands look like one but are on different contigs
                  island_name=island_name,
                  island_ID=paste(island_name, island_number, sep = '_'),
                  island_length=(end-start)+1,
                  locus_tags=list(locus_tag), # list of all the locus tags in the island
                  num_genes = length(locus_tag))

      genome_res[[island_index]] <- tmp2

    }


  }
  print(paste('Finished',genome_name))
  tmp_res <- bind_rows(genome_res)
  res[[genome_index]] <- tmp_res
}

##### bind all the islands together and select appropriate cols

res_4_real <- bind_rows(res) %>%
  select(seqid, genome_name, start, end, island_ID, island_length, num_genes, locus_tags) %>%
  filter(num_genes > min_genes)

##### read in fastas
genome_filenames <- list.files(path = './gifrop_out/sequence_data', pattern = '.fna', full.names = TRUE)


print('reading in fastas...')
genome_seqs_list <- sapply(genome_filenames, readDNAStringSet)
print('done reading in fastas')

# makes a dataframe of seqid and seqid length
seq_lens <- tibble(seqid=unlist(lapply(genome_seqs_list, names)),
                   seqid_len=unlist(lapply(genome_seqs_list, width)))


## IMPORTANT ##
# this adds in the seqid lenghts and then if the island occupies greater than
# 90 % of the whole contig it selects the whole thing.  Islands will be output with Istart and Iend
res_4_real <- res_4_real %>%
  left_join(seq_lens) %>%
  mutate(percent_island=round(island_length/seqid_len*100, 3),
         Istart=ifelse(percent_island > 90, 1, start),
         Iend= ifelse(percent_island >90, seqid_len, end))


# this block associates each locus tag with an accessory fragment value from roary
acc_frag_tmp <- gpa %>%
  gather(key='genome', value='locus_tags', -c(1:14)) %>%
  select(locus_tags, `Accessory Fragment`) %>%
  transmute(acc_frag = `Accessory Fragment`,
            locus_tags = locus_tags) %>%
  filter(!is.na(acc_frag) & !is.na(locus_tags))




# creates a tibble containing fasta files for each genome name
genome_seqs <- tibble(genome_name=sub('.*/gifrop_out/sequence_data/(.*).fna','\\1',names(genome_seqs_list)),
                      genome=genome_seqs_list)

nesty_res <- res_4_real %>%
  group_by(genome_name) %>%
  nest() %>%
  ungroup() %>%
  left_join(genome_seqs)


nesty_res <- nesty_res %>%
  mutate(islands_mon = map2(.x = data, .y = genome, .f = get_islands))

all_islands <- DNAStringSetList(nesty_res$islands_mon)
all_islands <- unlist(all_islands)

Biostrings::writeXStringSet(all_islands, './gifrop_out/my_islands/All_islands.fasta')


#### OUTPUT GFFS HERE ####

# now need way to incorporate ID_non_island_loc_tags
print('collecting island gffs')


# fix extra info on locus tags...
res_4_real <- bind_rows(gffs) %>%
  group_by(seqid) %>%
  select(locus_tag) %>%
  nest() %>%
  mutate(seqid_loc_tags=list(unlist(data, use.names = FALSE))) %>% 
  select(-data) %>%
  right_join(res_4_real) %>%
  mutate(only_island=map2_lgl(.x = seqid_loc_tags, .y = locus_tags, .f = ~ all(.x %in% .y)))


get_island_gff <- function(locus_tag_vec, gff){
  #this little guy should take the locus tag column from the island info dataframe
  # and return the corresponding annotations from the gff file for the genome,
  # basically a little mini gff that only covers the island
  gff[gff$locus_tag %in% locus_tag_vec,]

}

get_gffs <- function(island_info, gff){
  # using the island info dataframe, this function
  # extracts the annotation data for each island from the corresponding
  # gff

  # uses helper function, get_island_gff
  island_gffs <- map(island_info$locus_tags, .f = get_island_gff, gff)
  names(island_gffs) <- island_info$island_ID
  return(island_gffs)
}


gff_table <- res_4_real %>%
  group_by(genome_name) %>%
  nest() %>% ungroup() %>%
  left_join(tibble(genome_name=names(gffs), gffs=gffs)) # this tibble is just genome names and gffs


gff_table <- gff_table %>%
  mutate(island_gffs = map2(.x = data, .y = gffs, .f = get_gffs))




all_island_gffs <- unlist(gff_table$island_gffs, recursive = FALSE)

# names(all_island_gffs)


invisible(sapply(names(all_island_gffs),
       function (x) write_tsv(all_island_gffs[[x]], path = paste('./gifrop_out/my_islands/', x, ".gff", sep=""))))


print('done writing island gffs')

############## island info df ###########

# this builds the final island info dataframe
island_info <- res_4_real %>%
  ungroup() %>%
  tidyr::unnest(cols = 'locus_tags') %>%
  left_join(acc_frag_tmp) %>% group_by(island_ID) %>%
  summarise(seqid=unique(seqid),
            genome_name=unique(genome_name),
            start=unique(start),
            end=unique(end),
            island_length = unique(island_length),
            num_genes=unique(num_genes),
            locus_tags=paste(locus_tags, collapse = '|'),
            seqid_len = unique(seqid_len),
            percent_island = unique(percent_island),
            only_island = unique(only_island),
            Istart=unique(Istart),
            Iend=unique(Iend),
            acc_frag = paste(unique(acc_frag), collapse = '|')) %>%
  select(-Istart, -Iend)


write_csv(island_info, './gifrop_out/my_islands/island_info.csv')



### ADJUST START AND END COORDS HERE ###
# might want to take leading and trailing DNA sequences before and after coding sequences
# could look at insertion sites etc

if(flankingDNA > 0){

  print('ALSO OUTPUTTING ISLANDS WITH FLANKING DNA')
  print(paste('taking', flankingDNA, 'bp on either side of the island'))

  #EXTEND CHROMOSOMAL ISLANDS STARTS AND ENDS TO flankingDNA_len
  res_4_real <- res_4_real %>%
    mutate(Istart=Istart-flankingDNA,
           Iend=Iend + flankingDNA,
           Istart=ifelse(Istart < 1, 1, Istart),
           Iend= ifelse(Iend > seqid_len, seqid_len, Iend))


  system('mkdir ./gifrop_out/my_islands/with_flanking')

  nesty_res <- res_4_real %>%
    group_by(genome_name) %>%
    nest() %>%
    ungroup() %>%
    left_join(genome_seqs)

  nesty_res <- nesty_res %>%
    mutate(islands_mon = map2(.x = data, .y = genome, .f = get_islands))
  all_islands <- DNAStringSetList(nesty_res$islands_mon)
  all_islands <- unlist(all_islands)
  Biostrings::writeXStringSet(all_islands, './gifrop_out/my_islands/with_flanking/All_islands_with_flanking.fasta')

}
