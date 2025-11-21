# Pangenome complex loci
genotype and visualize structural variants in the human pangenome graphs (Figure 5 of [A draft human pangenome reference](https://doi.org/10.1038/s41586-023-05896-x))

## Visualize structures of complex loci using [bandage plot](https://github.com/rrwick/Bandage) 
- Select reference ranges of targeted loci:

  For each locus, we defined a window extending ~10 kb upstream and downstream of the gene boundaries to ensure all gene related bubbles are captured in the pangenome subgraph. As an example for the RHD/RHCE locus, the GRCh38 positions of RHD, TMEM50A, and RHCE are:
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
  ```bash
  # minigraph-cactus graph (index graph first to get .gfab index file)
  gfabase sub GRCh38-f1g-90-mc-aug11.gfab GRCh38.chr1:25240000-25460000 --range --connected --view --cutpoints 1 --guess-ranges -o RH_locus.walk.gfa
  ## remove --connected parameter to include nodes and edges only in the output subgraph without paths, to make the file small and easy to load to bandage
  
  # PGGB graph (put region to a bed file)
  odgi extract -i chr1.pan.smooth.og -o chr1.pan.RH_locus.og -b chr1.RH_locus.bed -E -P
  # conver og to gfa
  odgi view -i chr1.pan.RH_locus.og -g > chr1.pan.RH_locus.gfa
  ```

  then upload sugraph gfa file to bandage, and visualize the structure using bandage 

## Identify gene position in the subgraph and visualize it in bandage plot
- We aligned Ensembl (release 106) GRCh38 version gene sequences to the MC graph and PGGB graph using GraphAligner (v.1.0.13) with parameter settings -x vg --try-all-seeds --multimap-score-fraction 0.1 to identify the gene positions within the graphs.
  ```bash
  GraphAligner -g chr1.pan.RH_locus.no_pline.gfa -f RHD.RHCE.fa -a RHD.RHCE.aln.more.gaf -x vg --try-all-seeds --multimap-score-fraction 0.3
  # --multimap-score-fraction could be change to include more or less alignment results based on your need
  ```
- To display gene positions on Bandage plots, we applied a green-to-blue color gradient to the nodes belonging to each gene.
  
  Customized color on bandage can be added by uploada a csv file to bandage, detail can be found here: https://github.com/rrwick/Bandage/wiki/Colour-schemes
  (command line of bandage might also accomplish this, but I have never tried)

- We manually drew guides next to each Bandage plot to indicate approximate gene positions, exon structures, and transcription start sites from Ensembl canonical transcripts.
  (to be edited)        
  Software: Adobe Illustrator        
  Export the Bandage plot and open it in Adobe Illustrator. Lock this layer.        
  Exon position: Use GraphAligner alignments to locate nodes corresponding to each exon and gene; mark these nodes on a new annotation layer.        
  With the Curvature Tool, trace lines along each gene’s path. Add direction arrowheads via the Stroke panel.        
  Use the Direct Selection Tool to fine-tune anchor points and handles so each line follows the underlying Bandage curve precisely.        
  Draw exons with the Rectangle Tool; set their widths and positions according to the alignment results.        
  <p align="center">
      <img width="500" alt="Screenshot 2024-05-10 at 1 44 42 AM" src="https://github.com/lushjia/Pangenome_complex_loci/assets/38059727/3c11f344-188e-4fdf-b196-c8ba5e20bb6c">
  </p>
  
  
## Genotype and visualize structural alleles 
- Sequences from each assembly are represented as paths in the GFA graph, and variants appear as bubbles.
- We infer gene copy number in the graph by aligning gene sequences with GraphAligner and examining where they map on the graph; copy number in each assembly is then obtained by tracing each assembly’s path through the aligned gene regions.
- We identify big insertion nd deletion structure in each assembly by tracing each assembly path through large bubbles.
- Gene conversions are more complex to identify, because they do not show as bubbles in the graph. we identified nodes that were different between a gene and its homologous gene (for example, RHD and RHCE) based on the GraphAligner alignments described above. We refer to these as paralogous sequence variants (PSVs). A gene conversion event was detected when a gene's path traverses more than four PSVs of its homologous gene consecutively.
- Small insertions, deletions, and SNPs are out of scope for this analysis.        
- **Example analysis pipeline:**       
>_CY2D6/7_: `CYP2D6_7_locs.analysis.sh`           
_C4A/B_:       
_RHD/RHCE_:       
_LPA_:       
_HLA-A_:

  <p align="center">      
  <img width="500" alt="Screenshot 2024-05-10 at 1 48 33 AM" src="https://github.com/lushjia/Pangenome_complex_loci/assets/38059727/cc4d03ca-3e99-4289-a800-e2546745cad4">
  </p>

