# Pangenome_complex_loci
visualizaion of gene locations and paths of complex loci in human pangenome graph

## bandage plot
- Select reference ranges of interested loci:

  Based on the loci/genes we are interested in, we added ~10k flanking sequence to each side of the genes to ensure all bubbles of the region in the pangenome graph are included. Take RHD/RHCE loci as an example, positions of RHD, TMEM50A and RHCE genes in GRCh38 are:
  ```
  grch38#chr1	25272393	25330445	RHD
  grch38#chr1	25362249	25430192	RHCE
  grch38#chr1	25338317	25362361	TMEM50A
  ```
  The reference range is: 
  ```
  grch38#chr1	25240000	25460000	RH_locus
  ```

- Subgraph use gfabase / odgi:  
  Upload gfa to bandage, run; randomly, if overlap, rerun.  
  Command line tool 


## color viz to reveal gene position  
- tsv generation
- loaded to bandage plot 
- potential automatic command line to do this 


## arrows/lines to show positions of genes and gene elements
- Software: adobe illustrator  
Load bandage plot into adobe illustrator,  
Exon position: Nodes of exons and genes are found by align sequences to pangenome graph by GraphAligner. Mark those nodes in bandage plot.   
Lines are drawn by curvature tool along with the gene. Arrows at the end of each line are added by changing stroke.  
To make lines match the curve of the bandage plot, we adjusted the curve and position of each anchor by direct selection tool.   
Exons are drawn by rectangle tool. Widths and positions are based on alignment results.  



