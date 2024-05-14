# Pangenome_complex_loci
visualize and genotype structural variants in the human pangenome graphs (Figure 5 of [A draft human pangenome reference](https://doi.org/10.1038/s41586-023-05896-x))

## Visualize structures of complex loci using [bandage plot](https://github.com/rrwick/Bandage) 
- Select reference ranges of specific loci:

  The region we used for the loci is ~10kb upstream downstream to the start and end of the genes to ensure all bubbles at the gene region are included in the pangenome subgraph. Take RHD/RHCE loci as an example, positions of RHD, TMEM50A and RHCE genes in GRCh38 are:
  ```
  grch38#chr1	25272393	25330445	RHD
  grch38#chr1	25362249	25430192	RHCE
  grch38#chr1	25338317	25362361	TMEM50A
  ```
  The region of RHD&RHCE loci we used is: 
  ```
  grch38#chr1	25240000	25460000	RH_locus
  ```

- extract subgraph using gfabase / odgi:  
  ```
  # minigraph-cactus graph (index graph first to get .gfab index file)
  gfabase sub GRCh38-f1g-90-mc-aug11.gfab GRCh38.chr1:25240000-25460000 --range --connected --view --cutpoints 1 --guess-ranges -o RH_locus.walk.gfa
  ## remove --connected parameter to include nodes and edges only in the output subgraph without paths, to make the file small and easy to load to bandage
  
  # PGGB graph (put region to a bed file)
  odgi extract -i chr1.pan.smooth.og -o chr1.pan.RH_locus.og -b chr1.RH_locus.bed -E -P
  # conver og to gfa
  odgi view -i chr1.pan.RH_locus.og -g > chr1.pan.RH_locus.gfa
  ```

  then upload sugraph gfa file to bandage, and visualize the structure using bandage 

## Identify gene position in the subgraph and visualize in bandage plot
- We aligned Ensembl (release 106) GRCh38 version gene sequences to the MC graph and PGGB graph using GraphAligner (v.1.0.13) with parameter settings -x vg --try-all-seeds --multimap-score-fraction 0.1 to identify the gene positions within the graphs.
  ```
  GraphAligner -g chr1.pan.RH_locus.no_pline.gfa -f RHD.RHCE.fa -a RHD.RHCE.aln.more.gaf -x vg --try-all-seeds --multimap-score-fraction 0.3
  # --multimap-score-fraction could be change to include more or less alignment results based on your need
  ```
- To show locations of genes on Bandage plots, we applied colour gradients from green to blue to the nodes of each gene. 
  
  Customized color on bandage can be added by uploada a csv file to bandage, detail can be found here: https://github.com/rrwick/Bandage/wiki/Colour-schemes
  (command line of bandage might also accomplish this, but I have never tried)

- Lines alongside the Bandage plots showing approximate gene positions, exons and transcription start sites based on Ensembl Canonical transcripts were drawn by hand.
  (to be edited)        
  Software: adobe illustrator  
  Load bandage plot into adobe illustrator,  
  Exon position: Nodes of exons and genes are found by align sequences to pangenome graph by GraphAligner. Mark those nodes in bandage plot.   
  Lines are drawn by curvature tool along with the gene. Arrows at the end of each line are added by changing stroke.  
  To make lines match the curve of the bandage plot, we adjusted the curve and position of each anchor by direct selection tool.   
  Exons are drawn by rectangle tool. Widths and positions are based on alignment results.
      <img width="500" alt="Screenshot 2024-05-10 at 1 44 42 AM" src="https://github.com/lushjia/Pangenome_complex_loci/assets/38059727/3c11f344-188e-4fdf-b196-c8ba5e20bb6c">

  
  
## Genotype and visualize structural alleles 
- Sequences of each assembly are represented by paths in a GFA file, and variants are represented by bubbles in the graph.
- Copy number of genes in the graph could be told from where the gene sequence has been aligned to the grpah based on GraphAligner results. Copy number of genes in each assembly can then be told by tracing the path of assembly through the regions that the genes have been aligned. 
- Big insertions, deletions of each assembly are identified by tracing paths of each assembly (which represent sequences) through big bubbles.
- Gene conversions are more complex to identiry, because they are not shown as bubbles in the graph.         
  we identified nodes that were different between a gene and its homologous gene (for example, RHD and RHCE) based on the GraphAligner alignments described above. We refer to these as paralogous sequence variants. A gene conversion event was detected if a path of a gene goes through more than four paralogous sequence variants of its homologous gene in a row.
- Other small insertion, deletion and SNPs are not covered by this analysis. 
           
  <img width="500" alt="Screenshot 2024-05-10 at 1 48 33 AM" src="https://github.com/lushjia/Pangenome_complex_loci/assets/38059727/cc4d03ca-3e99-4289-a800-e2546745cad4">


