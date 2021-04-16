
[![install with
bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gifrop/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gifrop/badges/downloads.svg)](https://anaconda.org/bioconda/gifrop)  

# gifrop

**G**enomic **I**slands **f**rom **Ro**ary **P**angenomes.
## This readme is a bit out of date right now.  I'll try and update soon.  

This program is supposed to identify ‘genomic islands’ from roary
pangenomes. The intent was to try and identify units of horizontal gene
transfer within very closely related strains using a pangenome
framework.

In general, this type of analysis works best for very closely related
strains, for example, a collection of *Salmonella enterica serovar
typhimurium*.

In theory, this program will work on pangenomes of any size, though
things get very messy with too many genomes or genomes that are too
different from one another.

**Prerequisites:**

1.  A pangenome calculated by roary: Any pangenome should work, even
    those that do not split paralogues (the -s option of roary).

2.  the \*.gff files used by roary to create the pangenome. I have only
    used those produced by prokka but others may work. 

If all you have are nucleotide fastas then you can use the included
`pan_pipe` script to take care of the whole pipeline from annotation
with prokka to pangenome calculation with roary to genomic island
extraction with gifrop.

## conda installation

**create new conda environment**  
`conda create -n gifrop python=3.7`

**activate new environment**  
`conda activate gifrop`

**install gifrop**  
`conda install -c conda-forge -c bioconda -c defaults gifrop=0.0.9`

or

## manual installation

Install the dependencies and make sure they are in your path:

1)  [GNU parallel](https://www.gnu.org/software/parallel/)  
2)  [abricate](https://github.com/tseemann/abricate)  
    \- You will also need these custom abricate databases
      - <https://github.com/Jtrachsel/megares_db_4_abricate>  
      - <https://github.com/Jtrachsel/viroseqs>  
3)  [R 3.6.\*](https://www.r-project.org/) and these R packages:  
    \-
    ‘dplyr’,‘tidyr’,‘readr’,‘tibble’,‘ggplot2’,‘purrr’,‘Biostrings’,‘BSgenome’,
    ‘igraph’, ‘parallelDist’, ‘digest’  
4)  [roary](https://sanger-pathogens.github.io/Roary/) (if using
    pan\_pipe script)  
5)  [prokka](https://github.com/tseemann/prokka) (if using pan\_pipe
    script)  
6)  make the scripts from this repo available in your path

<br>

## Running gifrop

  - Navigate to a directory containing prokka annotated gff files and
    the output of roary, specifically the `gene_presence_absence.csv`
    file.  
  - run gifrop:
      - `gifrop --get_islands --threads 8`  
  - all outputs go to the `gifrop_out` directory.

## pan\_pipe script

The included script `pan_pipe` is a complete pipeline for extracting
genomic islands from a collection of nucleotide fastas.

These are the basic steps:  
1\) Annotate with prokka  
2\) Calculate pangenome with roary  
3\) Extract genomic islands with gifrop

You are able to pass arguments to prokka, roary and gifrop using quoted
strings of arguments, see `pan_pipe --help` for more information.

<br>

## How this program identifies pangenomic islands

1)  **remove the core genome from the pangenome.**
    
      - This removes any gene that occurs in all isolates
        exactly 1 time. This means that genes that only occur in some
        genomes and genes that occur twice or more in any genome are
        currently kept for consideration.  
      - Alternatively, you can find genomic islands relative to a reference with the '--reference / -r' flag.
            - In this mode, any genes present in the reference are removed from the pangenome.  

2)  **Identify strings of consecutive genes** (locus tag orders)  

3)  **Remove strings that contain fewer than the minimum number of
    genes.**
    
      - These remaining strings of consecutive genes are ‘genomic
        islands’ ‘pan-genomic islands?’  

4)  **Classify pangenomic islands** by running abricate with the
    following databases:
    1)  The NCBI amr database packaged with abricate.  
    2)  [Megares2.0](https://megares.meglab.org/) database (bc it has
        metal tolerance)  
    3)  [plasmid finder](https://cge.cbs.dtu.dk/services/PlasmidFinder/)
        (replicon genes)  
    4)  [vfdb](http://www.mgc.ac.cn/VFs/main.htm) (virulence genes)  
    5)  [ProphET](https://github.com/jaumlrc/ProphET) phage db  

5)  **Cluster pangenomic islands**, this can get messy.
    
      - De-replicate islands with identical gene content
      - calculate pairwise distances between all unique islands using two separate distance metrics
            - Simpson (AKA overlap distance), and Jaccard distance.  
            - Convert distances to similarities.  
      - make two graphs where pangenomic islands are nodes and they are
        joined by edges representing their similarities to other islands.
        One graph for Overlap Coeficient and one for Jaccard.  Edges with 0 weight are 
        removed from these graphs.
      - four levels of clustering:  
        \- primary clustering: Uses the Overlap Coeficient graph any island connected by non-zero weight edges are in the same primary cluster.  
        \- secondary clustering: First The Overlap Coeficient graph is pruned to remove edges with weights less than 0.5.
        an overlap similarity of 0.5 linking two islands means that at least half of the genes in the smaller island are present in the larger island.  
        The Louvain clustering on algorithm is then applied to detect densly connected communities.  
        \- tertiary clustering: The Overlap Coeficient graph is pruned to remove edges with weights less than 0.75.
        an overlap similarity of 0.75 linking two islands means that at least 3/4 of the genes in the smaller island are present in the larger island.
        The Louvain clustering on algorithm is then applied to detect densly connected communities.
        \- quaternary clustering: The Jaccard Coeficient graph is pruned to remove edges with weights less than 0.75.
        a jaccard similarity of 0.75 linking two islands means that both connected islands contain at least 75% of the total gene content of both islands combined.  
        The Louvain clustering on algorithm is then applied to detect densly connected communities.


6)  **output pangenomic island info.**  
    All outputs go to the `gifrop_out` folder.  
    main outputs:
    
    ``` 
      1. clustered_island_info.csv  
      2. pan_with_island_info.csv  
      3. pan_only_islands.csv  
      4. figures directory  
         - some basic exploratory figs  
      5. All the pangenomic islands  
         - as a single fasta file 'All_islands.fasta'  
         - can be split into individual fastas if requested.   
      6. logfiles  
    ```

<br>

| file                        | description                                                                                                                           |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| clustered\_island\_info.csv | This is a csv file that contains a detailed description of every genomic island detected in this pangenome.                           |
| pan\_with\_island\_info.csv | The roary ‘gene\_presence\_absence.csv’ file annotated with information about which genomic islands and clusters the genes belong to. |
| pan\_only\_islands.csv      | The same as the previous file but filtered to only include genes that occur on at least 1 genomic island                              |

## graphical overview:

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->

## Known issues

**Draft assemblies**  
gifrop can handle incomplete, draft assemblies, but the genomic islands
it identifies **will never span contig breaks**. That is, if there is a
genomic island in your assembly that is divided across two contigs, it
will be split into 2 different genomic islands.

**Complete genomes**  
If you input complete genomes, (with closed circular chromosomes) and a
genomic island spans the arbitrarily chosen start/end coordinates, it
will be split into two genomic islands.

**Combining of adjacent/nested genomic islands**  
If there is a situation where two distinct genomic islands are inserted
adjacent to one another or one nested within another etc., they will be
treated as one genomic island.

**Error when identical contig names in different assemblies**  


## TODO

1)  link abricate results to pangenome, use coords of hits with gffs to
    get loc\_tags?  
2)  Include example data?

Suggestions welcome.
