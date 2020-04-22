# Need to check that these libraries are installed

reqs <- c('dplyr','tidyr','readr','tibble','ggplot2','purrr','Biostrings','BSgenome', 'igraph', 'pheatmap')


unsatisfied_reqs <- reqs[!(reqs %in% installed.packages()[,"Package"])]

if(length(unsatisfied_reqs) > 0){
  print('The following R libraries were not detected:')
  print(unsatisfied_reqs)
  print('please install these packages and try again')

} else {

  print('All required R packages were detected')

}





