# pangenome graphs
# minigraph-cactus
GRCh38-f1g-90-mc-aug11.gfa
# pggb 
$chr.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz

#======================================================
# CY2D6 gene, chr22:42,115,000-42,155,000
mc_path=""
pggb_path=""
local_path=""


#===================
# bandage to visualize minigraph-cactus graph 

wkdir = $mc_path/cyp2d_mini_2021_08_11

# index gfa 
gfabase load -o GRCh38-f1g-90-mc-aug11.gfab GRCh38-f1g-90-mc-aug11.gfa

# extract cyp2d6/7 locus 
gfabase sub GRCh38-f1g-90-mc-aug11.gfab GRCh38.chr22:42115000-42155000 --range --view --cutpoints 1 --guess-ranges -o cyp2d6_locus.gfa 

# csv for coloring GRCh39 node in bandage plot  
echo "Name,Colour" > hg38.cyp2d6_locus.csv
less cyp2d6_locus.gfa | grep "^S" |grep "GRCh38"| awk '{print $2",red"}' >> hg38.cyp2d6_locus.csv

# extract path + walk in this region 
gfabase sub GRCh38-f1g-90-mc-aug11.gfab GRCh38.chr22:42115000-42155000 --range --connected --view --cutpoints 1 --guess-ranges -o cyp2d6_locus.walk.gfa 
less cyp2d6_locus.walk.gfa | grep -E "^W|^P" > path_walk_cyp2d.gfa 

# trasnfer to my mac to visualize in bandage 

rsync -arvzp $user@server:$mc_path/cyp2d6_locus.gfa $local_path
rsync -arvzp $user@server:$mc_path/hg38.cyp2d6_locus.csv $local_path


#===================
# bandage to visualize PGGB graph  
wkdir = $pggb_path
cd $wkdir 
chr=chr22

# gfa to og 
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/$chr.hprc-v1.0-pggb.gfa.gz
gunzip $chr.hprc-v1.0-pggb.gfa.gz
odgi build -t 48 -P -g $chr.hprc-v1.0-pggb.gfa -o $chr.hprc-v1.0-pggb.og

# extract cyp2d6/7 locus from chr22 og
odgi extract -i ./$chr.hprc-v1.0-pggb.og -o ./$chr.pan.CYP2D6_locus.og -b ./chr22.CYP2D6_locus.bed -E --threads 2 -P

# transfer to gfa
odgi view -i chr22.pan.CYP2D6_locus.og -g > chr22.pan.CYP2D6_locus.gfa

# extract p lines for assembly paths in pangenome 
zless chr22.pan.CYP2D6_locus.gfa| grep "^P" > chr22.pan.CYP2D6_locus.pline.txt

# csv (bandage color)
less chr22.pan.CYP2D6_locus.pline.txt| sed -n '4,6p'| awk '{print $3}'| sed 's/,/\n/g' | sed 's/+/,red/g'| sed 's/-/,red/g' > hg38.chr22.pan.CYP2D6_locus.csv
# vi hg38.chr22.pan.CYP2D6_locus.csv to add header Name,Colour

# try pggb png visualization
odgi sort -i chr22.pan.CYP2D6_locus.og -o - -O | \
    odgi viz -i - -o chr22.pan.CYP2D6_locus.png -s '#' -P 20

# transfer to my mac
rsync -arvzP $user@server:$pggb_path/chr22.pan.CYP2D6_locus.gfa $local_path
rsync -arvzP $user@server:$pggb_path/hg38.chr22.pan.CYP2D6_locus.csv $local_path
# on my mac to view gfa by bandage 

# why two components in pggb graph 
$pggb_path/try
# conclusion - not a bubble in a larger region at graph

# seperate two components in og file 
odgi explode -i chr22.pan.CYP2D6_locus.og -p chr22.pan.CYP2D6_locus
# transfer to gfa
odgi view -i chr22.pan.CYP2D6_locus.0.og -g > chr22.pan.CYP2D6_locus.0.gfa
odgi view -i chr22.pan.CYP2D6_locus.1.og -g > chr22.pan.CYP2D6_locus.1.gfa



#===================
# gene conversion analysis 
#=========
# pggb 
# graph align to understand the graph pggb
wkdir = $pggb_path
cd $wkdir 

# fasta of cyp2d6, cpy2d7, spacer
cyp2d6_cyp2d7_spacer.fasta
# gene reigon based on gencode gtf https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz 
# download 
<gencode_folder>
# extract only gene
less gencode.v39.annotation.gtf| awk '$3=="gene"{ print $0}' > gencode.v39.annotation.geneonly.gtf
# to bed
less gencode.v39.annotation.geneonly.gtf| awk -F "\t" '{OFS="\t"; print $1,$4-1,$5,$7,$9}' > gencode.v39.annotation.geneonly.bed

# gene sequence get from ucsc genome browser
# spacer from curius papcer
cyp2d6: chr22:42,126,499-42,130,865
cyp2d7: chr22:42,140,203-42,149,455
cyp2d7 spacer: chr22:42,138,124-42,139,676

GraphAligner -g chr22.pan.CYP2D6_locus.gfa -f cyp2d6_cyp2d7_spacer.fasta -a cyp2d6_cyp2d7_spacer.aln.gaf -x vg




=========
# distinguish cyp2d6 cyp2d7 cyp2d6-7 pggb (gene conversion analysis)

# pggb path - with all diff&uniq cyp2d6,cyp2d7 nodes
wkdir = $pggb_path
cd $wkdir 

# cyp2d6 nodes
less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d'> cyp2d6.pan.nodes.txt
# cyp2d7 nodes
less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '2p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d'> cyp2d7.pan.nodes.txt
# spacer nodes 
less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '3p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d'> spacer.pan.nodes.txt

# cyp2d7 longer than cyp2d6 
# remove both cyp2d6 & cyp2d7 nodes after 345128 because graph branch after this node

# find cyp2d6 uniq / cyp2d7 uniq nodes
comm -23 <(sort cyp2d6.pan.nodes.txt) <(sort cyp2d7.pan.nodes.txt)| awk '{print $1}' > cyp2d6.pan.uniq.nodes.txt
comm -13 <(sort cyp2d6.pan.nodes.txt) <(sort cyp2d7.pan.nodes.txt)| awk '{print $1}' > cyp2d7.pan.uniq.nodes.txt

# visualize cyp2d6 uniq / cyp2d7 uniq nodes in bandage 
# csv
# cyp2d6, cyp2d7, spacer: '#A6B692', '#FBF9D0','#535341'
cat <(less cyp2d6.pan.uniq.nodes.txt| awk '{print $1",\"#A6B692\""}') <(less cyp2d7.pan.uniq.nodes.txt| awk '{print $1",\"#FBF9D0\""}') <(less spacer.pan.nodes.txt| awk '{print $1",\"#535341\""}') > CYP2D6_uniq_nodes.csv 
cat <(less cyp2d6.pan.uniq.nodes.txt| awk '{print $1",red"}') <(less cyp2d7.pan.uniq.nodes.txt| awk '{print $1",blue"}') <(less spacer.pan.nodes.txt| awk '{print $1",\"#535341\""}') > CYP2D6_uniq_nodes.csv 
# vi CYP2D6_uniq_nodes.csv to add header Name,Colour
rsync -arvzP $user@server:$pggb_path/CYP2D6_uniq_nodes.csv $local_path


# grep header and improtant nodes
# make bash file 
echo -n "less chr22.pan.CYP2D6_locus.pline.txt | sed 's/\t/\n/g'|sed 's/,/\n/g'| grep -E \"grch38|chm13|HG0|NA|h1tg|h2tg" > pggb.find.cyp2d67.path.sh
for node in `less cyp2d6.pan.uniq.nodes.txt`;
do
	echo -n "|" >> pggb.find.cyp2d67.path.sh
	echo -n $node >> pggb.find.cyp2d67.path.sh
done
for node in `less cyp2d7.pan.uniq.nodes.txt`;
do
	echo -n "|" >> pggb.find.cyp2d67.path.sh
	echo -n $node >> pggb.find.cyp2d67.path.sh
done
for node in `less spacer.pan.nodes.txt`;
do
	echo -n "|" >> pggb.find.cyp2d67.path.sh
	echo -n $node >> pggb.find.cyp2d67.path.sh
done
# end node
echo -n "|" >> pggb.find.cyp2d67.path.sh
echo -n 345128 >> pggb.find.cyp2d67.path.sh
echo "\" > path_walk.cyp2d.pan.alluniq.txt" >> pggb.find.cyp2d67.path.sh

sbatch pggb.find.cyp2d67.path.sh


# clean not useful rows (assemblies not have these nodes) and turn rows to columns
# 91 lines (3 HG01071)
less path_walk.cyp2d.pan.alluniq.txt| tr '\n' '\t'|sed 's/\tchm13/\nchm13/g'|sed 's/\tgrch38/\ngrch38/g'|sed 's/\tNA/\nNA/g' | sed 's/\tHG/\nHG/g'|awk 'NF>1{print $0}' > path_walk.CYP2D6_locus.end.alluniq.pan.clean.txt

# turn positive to negative 
"""
convert positive strand to negative strand to make all path direction consistent
"""
python3

input_file=open('path_walk.CYP2D6_locus.end.alluniq.pan.clean.txt','r')
output_file='path_walk.CYP2D6_locus.end.alluniq.pan.clean.direction.txt'

with open(output_file,'w') as output_f:
	for line in input_file:
		line_split=line.strip().split('\t')
		node_list=line_split[1:]
		if node_list[0][:-1] != "345128": #end node
			output_f.write(line)
		else:
			# need to invert order
			output_f.write(line_split[0] + '\t')
			# deal with last node
			if node_list[-1][-1] == "+":
				node_1 = node_list[-1][0:len(node_list[-1])-1] + '-'
			else:
				node_1 = node_list[-1][0:len(node_list[-1])-1] + '+'
			output_f.write(node_1)
			# other nodes
			for node in node_list[len(node_list)-2::-1]:
				if node[-1] == "+":
					node = node[0:len(node)-1] + '-'
				else:
					node = node[0:len(node)-1] + '+'
				output_f.write("\t" + node )
			output_f.write("\n")

# 58 uniq path
less -S path_walk.CYP2D6_locus.end.alluniq.pan.clean.direction.txt| cut -f2-| sort| uniq | wc -l

# turn nodes into cyp2d6-6; cyp2d7-7; spacer-8; end-0
python3

input_file=open('path_walk.CYP2D6_locus.end.alluniq.pan.clean.direction.txt','r')
output_file='path_walk.CYP2D6_locus.end.alluniq.pan.clean.direction.6780.txt'

with open('cyp2d6.pan.uniq.nodes.txt') as f:
	cyp2d6_nodes = f.read().splitlines()

with open('cyp2d7.pan.uniq.nodes.txt') as f:
	cyp2d7_nodes = f.read().splitlines()

with open('spacer.pan.nodes.txt') as f:
	spacer_nodes = f.read().splitlines()

with open(output_file,'w') as output_f:
	for line in input_file:
		line_split=line.strip().split('\t')
		output_f.write(line_split[0])
		for node in line_split[1:]:
			node = node[0:len(node)-1]
			if node in cyp2d6_nodes:
				output_f.write('\t6')
			elif node in cyp2d7_nodes:
				output_f.write('\t7')
			elif node in spacer_nodes:
				output_f.write('\t8')
			elif node == '345128':
				output_f.write('\t0'*20) # changed to add more cells to seperate 2 gene
			else:
				output_f.write('\t' + node)
		output_f.write('\n') 

# 57 uniq paths, 2 missing
less -S path_walk.CYP2D6_locus.end.alluniq.pan.clean.direction.6780.txt | cut -f2-| sort | uniq | wc -l

rsync -arvzP $user@server:$pggb_path/path_walk.CYP2D6_locus.end.alluniq.pan.clean.direction.6780.txt $local_path


# add exon1 sequence back 
? where is the exon1 - based on duplication region?
#cyp2d6: chr22:42,126,499-42,130,865
#cyp2d7: chr22:42,140,203-42,149,455
#cyp2d7 spacer: chr22:42,138,124-42,139,676

chr22:42,130,417-42,131,001
chr22:42,144,078-42,144,724


GraphAligner -g chr22.pan.CYP2D6_locus.gfa -f $pggb_path/exon1.fa -a exon1.aln.gaf -x vg --try-all-seeds

# cyp2d6 exon1  
less exon1.aln.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d' > cyp2d6.exon1.nodes.txt
# cyp2d7 exon1 nodes
less exon1.aln.gaf| sed -n '2p;3p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d'> cyp2d7.exon1.nodes.txt

# find cyp2d6 uniq / cyp2d7 uniq nodes in exon1 
comm -23 <(sort cyp2d6.exon1.nodes.txt) <(sort cyp2d7.exon1.nodes.txt|uniq)| awk '{print $1}' > cyp2d6.pan.exon1.uniq.nodes.txt
comm -13 <(sort cyp2d6.exon1.nodes.txt) <(sort cyp2d7.exon1.nodes.txt|uniq)| awk '{print $1}' | sed '1d'> cyp2d7.pan.exon1.uniq.nodes.txt

# combine with uniq nodes before
cat cyp2d6.pan.uniq.nodes.txt cyp2d6.pan.exon1.uniq.nodes.txt > cyp2d6.pan.wholegene.uniq.nodes.txt
cat cyp2d7.pan.uniq.nodes.txt cyp2d7.pan.exon1.uniq.nodes.txt > cyp2d7.pan.wholegene.uniq.nodes.txt


# grep header and improtant nodes
# make bash file 
echo -n "less chr22.pan.CYP2D6_locus.pline.txt | sed 's/\t/\n/g'|sed 's/,/\n/g'| grep -E \"grch38|chm13|HG0|NA|h1tg|h2tg" > pggb.find.cyp2d67.wholegene.path.sh
for node in `less cyp2d6.pan.wholegene.uniq.nodes.txt`;
do
	echo -n "|" >> pggb.find.cyp2d67.wholegene.path.sh
	echo -n $node >> pggb.find.cyp2d67.wholegene.path.sh
done
for node in `less cyp2d7.pan.wholegene.uniq.nodes.txt`;
do
	echo -n "|" >> pggb.find.cyp2d67.wholegene.path.sh
	echo -n $node >> pggb.find.cyp2d67.wholegene.path.sh
done
for node in `less spacer.pan.nodes.txt`;
do
	echo -n "|" >> pggb.find.cyp2d67.wholegene.path.sh
	echo -n $node >> pggb.find.cyp2d67.wholegene.path.sh
done
# end node
echo -n "|" >> pggb.find.cyp2d67.wholegene.path.sh
echo -n 346111 >> pggb.find.cyp2d67.wholegene.path.sh
echo -n "|" >> pggb.find.cyp2d67.wholegene.path.sh
echo -n 343868 >> pggb.find.cyp2d67.wholegene.path.sh
echo "\" > path_walk.cyp2d.pan.alluniq.wholegene.txt" >> pggb.find.cyp2d67.wholegene.path.sh
# add sbatch header 
#!/bin/bash
#SBATCH --job-name=grep_node
#SBATCH --out="slurm-%j.out"
#SBATCH --time=05:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=None
#SBATCH --partition=pi_hall
sbatch pggb.find.cyp2d67.wholegene.path.sh

# clean not useful rows (assemblies not have these nodes) and turn rows to columns
# still 91 lines
less path_walk.cyp2d.pan.alluniq.wholegene.txt| tr '\n' '\t'|sed 's/\tgrch38/\ngrch38/g'|sed 's/\tchm13/\nchm13/g'|sed 's/\tNA/\nNA/g' | sed 's/\tHG/\nHG/g'|awk 'NF>1{print $0}' > path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.txt

# turn positive to negative 
"""
convert positive strand to negative strand
"""
python3
input_file=open('path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.txt','r')
output_file='path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.direction.txt'

with open(output_file,'w') as output_f:
	for line in input_file:
		line_split=line.strip().split('\t')
		node_list=line_split[1:]
		if (node_list[0][:-1] != "346111") and (node_list[0][:-1] != "343868"):
			output_f.write(line)
		else:
			# need to invert order
			output_f.write(line_split[0] + '\t')
			# deal with last node
			if node_list[-1][-1] == "+":
				node_1 = node_list[-1][0:len(node_list[-1])-1] + '-'
			else:
				node_1 = node_list[-1][0:len(node_list[-1])-1] + '+'
			output_f.write(node_1)
			# other nodes
			for node in node_list[len(node_list)-2::-1]:
				if node[-1] == "+":
					node = node[0:len(node)-1] + '-'
				else:
					node = node[0:len(node)-1] + '+'
				output_f.write("\t" + node )
			output_f.write("\n")

# 60 uniq path
less -S path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.direction.txt| cut -f2-| sort| uniq | wc -l

# turn nodes into cyp2d6-6; cyp2d7-7; spacer-8; end-0
python3

input_file=open('path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.direction.txt','r')
output_file='path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.direction.6780.txt'

with open('cyp2d6.pan.wholegene.uniq.nodes.txt') as f:
	cyp2d6_nodes = f.read().splitlines()

with open('cyp2d7.pan.wholegene.uniq.nodes.txt') as f:
	cyp2d7_nodes = f.read().splitlines()

with open('spacer.pan.nodes.txt') as f:
	spacer_nodes = f.read().splitlines()

with open(output_file,'w') as output_f:
	for line in input_file:
		line_split=line.strip().split('\t')
		output_f.write(line_split[0])
		for node in line_split[1:]:
			node = node[0:len(node)-1]
			if node in cyp2d6_nodes:
				output_f.write('\t6')
			elif node in cyp2d7_nodes:
				output_f.write('\t7')
			elif node in spacer_nodes:
				output_f.write('\t8')
			elif node in ['346111','343868']:
				output_f.write('\t0'*20) # changed to add more cells to seperate 2 gene
			else:
				output_f.write('\t' + node)
		output_f.write('\n') 

# 61 uniq paths, 1 missing
less -S path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.direction.6780.txt | cut -f2-| sort | uniq | wc -l

rsync -arvzP $user@server:$pggb_path/path_walk.CYP2D6_locus.end.alluniq.pan.wholegene.clean.direction.6780.txt $local_path




=========
# mini_cactus path with all diff&uniq cyp2d6,cyp2d7 nodes
wkdir=$mc_path
cd $wkdir

# align
GraphAligner -g cyp2d6_locus.gfa -f ./cyp2d6_cyp2d7_spacer.fasta -a cyp2d6_cyp2d7_spacer.aln.gaf -x vg --try-all-seeds

# no cyp2d7 at c2, so extend alignment cretieria 
GraphAligner -g cyp2d6_locus.gfa -f ./cyp2d6_cyp2d7_spacer.fasta -a cyp2d6_cyp2d7_spacer.aln.more.gaf -x vg --try-all-seeds --multimap-score-fraction 0.5


# mini_cactus path with all diff&uniq cyp2d6,cyp2d7 nodes

c2 cyp2d6(3) 21/ cyp2d7(59) 155
c1 cyp2d6(2) 5
s  cyp2d6(1) 0/ cyp2d7(12) 1
c3 cyp2d6(7) 237/ cyp2d7(11) 0

spacer 
 
c1 (61) 0 76912741 - 76913157 * need to modify 
c2l (63) 5 76913155 - 76913157
c2r (64) 8 76912893 - 76912900 * need to modify
c3 (62) 0 76913992 - 76914035 


# combine cyp2d6/cyp2d7/spacer nodes together
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '1p;2p;3p;7p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' > cyp2d6.all.nodes.txt
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '11p;12p;59p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' > cyp2d7.all.nodes.txt

less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '61,64p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sort| uniq| sed '1d' > spacer.all.nodes.txt
# remove 76912741 - 76912752 (not only space, but also on stem); also 76912893
# remove 76912900 from spacer due to overlap with gene 

>76913183 >76913673 s
<76912776 <76912830 c1
<76912900 <76913137 c2 
 >76914055 >76914314 c3

 >76914055 >76914445
 >76913183 >76914445
 <76912900 <76913137 

# since cyp2d7 longer than cyp2d6, remove following nodes from cyp2d7.all.nodes.txt
c3: >76914315 >76914445
s: >76913675 >76914445

# uniq nodes
comm -23 <(sort cyp2d6.all.nodes.txt) <(sort cyp2d7.all.nodes.txt)| awk '{print $1}' > cyp2d6.uniq.nodes.mini_cactus.txt
comm -13  <(sort cyp2d6.all.nodes.txt) <(sort cyp2d7.all.nodes.txt)| awk '{print $1}' > cyp2d7.uniq.nodes.mini_cactus.txt


# visualize cyp2d6 uniq / cyp2d7 uniq nodes in bandage 
# csv
# cyp2d6, cyp2d7, spacer: '#A6B692', '#FBF9D0','#535341'
cat <(less cyp2d6.uniq.nodes.mini_cactus.txt| awk '{print $1",red"}') <(less cyp2d7.uniq.nodes.mini_cactus.txt| awk '{print $1",blue"}') <(less spacer.all.nodes.txt| awk '{print $1",\"#535341\""}') > CYP2D6_uniq_nodes.csv 
# vi CYP2D6_uniq_nodes.csv to add header Name,Colour
rsync -arvzP $user@server:$mc_path/CYP2D6_uniq_nodes.csv $local_path



# modify path_walk_cyp2d.gfa to combine sample and contig name (diff p line and w line)
less path_walk_cyp2d.gfa | sed -n '1p' > path_walk_cyp2d.mod.gfa
less path_walk_cyp2d.gfa | sed '1d' | awk '{print $1"\t"$2"#"$3"#"$4"\t"$5"\t"$6"\t"$7}' >> path_walk_cyp2d.mod.gfa

# grep header and important nodes
# make bash file 
echo -n "less path_walk_cyp2d.mod.gfa | sed 's/\t/\n/g'|sed 's/,/\n/g'|sed 's/>/\n>/g' | sed 's/</\n</g'| grep -E \"GRCh38|CHM13|HG0|NA|JA|chr" > find.cyp2d67.modified.path.sh
for node in `less cyp2d6.uniq.nodes.mini_cactus.txt`;
do
	echo -n "|" >> find.cyp2d67.modified.path.sh
	echo -n $node >> find.cyp2d67.modified.path.sh
done
for node in `less cyp2d7.uniq.nodes.mini_cactus.txt`;
do
	echo -n "|" >> find.cyp2d67.modified.path.sh
	echo -n $node >> find.cyp2d67.modified.path.sh
done
for node in `less spacer.all.nodes.txt`;
do
	echo -n "|" >> find.cyp2d67.modified.path.sh
	echo -n $node >> find.cyp2d67.modified.path.sh
done
for node in 76913673 76912833 76913137 76914317; # end s c1, c2, c3 
do
	echo -n "|" >> find.cyp2d67.modified.path.sh
	echo -n $node >> find.cyp2d67.modified.path.sh
done

echo "\" > path_walk.cyp2d.modified.end.mini_cactus.alluniq.txt" >> find.cyp2d67.modified.path.sh
# add sbatch header 
#!/bin/bash
#SBATCH --job-name=grep_node
#SBATCH --out="slurm-%j.out"
#SBATCH --time=05:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=None
#SBATCH --partition=general
sbatch find.cyp2d67.modified.path.sh 


less path_walk.cyp2d.modified.end.mini_cactus.alluniq.txt| tr '\n' '\t'| sed 's/CHM13/\nCHM13/g'| sed 's/NA/\nNA/g' | sed 's/HG/\nHG/g'| awk 'NF>2{print $0}' > path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.pre.txt

less path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.pre.txt | sed 's/+\t/\t>/g' > path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.txt
# manully move the end of the fitst line > to head;



# turn positive to negative 
"""
convert positive stran to negative strand
"""
python3
input_file=open('path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.txt','r')
output_file='path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.direction.txt'

with open(output_file,'w') as output_f:
	for line in input_file:
		line_split=line.strip().split('\t')
		node_list=line_split[1:] # change
		if node_list[0][1:] != "76913673" and node_list[0][1:] != "76914317" :
			output_f.write(line)
		else:
			# need to invert order
			output_f.write(line_split[0] + '\t' )
			# deal with last node
			if node_list[-1][0] == ">":
				node_1 = '<' + node_list[-1][1:len(node_list[-1])] 
			else:
				node_1 = '>' + node_list[-1][1:len(node_list[-1])] 
			output_f.write(node_1)
			# other nodes
			for node in node_list[len(node_list)-2::-1]:
				if node[0] == ">":
					node = '<' + node[1:len(node)] 
				else:
					node = '>' + node[1:len(node)] 
				output_f.write("\t" + node )
			output_f.write("\n")

# 91 line, 68 uniq path
less -S path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.direction.txt| cut -f2-| sort| uniq | wc -l


# turn nodes into cyp2d6-6; cyp2d7-7; spacer-8; end-11,12,13,14
# use diff number to represent diff end
# 76913673 76912833 76913137 76914317; # end s c1, c2, c3 

python3

input_file=open('path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.direction.txt','r')
output_file='path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.direction.6780.txt'

with open('cyp2d6.uniq.nodes.mini_cactus.txt') as f:
	cyp2d6_nodes = f.read().splitlines()

with open('cyp2d7.uniq.nodes.mini_cactus.txt') as f:
	cyp2d7_nodes = f.read().splitlines()

with open('spacer.all.nodes.txt') as f:
	spacer_nodes = f.read().splitlines()

with open(output_file,'w') as output_f:
	for line in input_file:
		line_split=line.strip().split('\t')
		output_f.write(line_split[0] )
		for node in line_split[1:]:
			node = node[1:len(node)]
			if node in cyp2d6_nodes:
				output_f.write('\t6')
			elif node in cyp2d7_nodes:
				output_f.write('\t7')
			elif node in spacer_nodes:
				output_f.write('\t8') # spacer
			elif node == '76912833': # c1 end
				output_f.write('\t11'*20)
			elif node == '76913137': # c2 end 
				output_f.write('\t12'*20)
			elif node == '76913673': # s end
				output_f.write('\t13'*20)
			elif node == '76914317': # c3 end
				output_f.write('\t14'*20)
			else:
				output_f.write('\t' + node)
		output_f.write('\n') 

# 58 uniq paths, 2 missing
less -S path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.direction.6780.txt | cut -f2-| sort | uniq | wc -l

rsync -arvzP $user@server:$mc_path/path_walk.cyp2d.modified.end.mini_cactus.alluniq.clean.direction.6780.txt $local_path






===================
# odgi 2d graph of cyp2d6 mini - show diff haplotype 
odgi build -P -g cyp2d6_locus.walk.gfa -o cyp2d6_locus.og
# sort, layout, draw 
odgi sort -P -t 16 -i cyp2d6_locus.og -o cyp2d6_locus.sort.og -O; odgi sort -P -p Ygs -t 16 -i cyp2d6_locus.sort.og -o cyp2d6_locus.sort2.og; odgi layout -i cyp2d6_locus.sort2.og -o cyp2d6_locus.sort2.og.lay -t 16 -P ; odgi draw -i cyp2d6_locus.sort2.og -c cyp2d6_locus.sort2.og.lay -p CYP2D6_locus.sort2.og.lay.draw.png -H 1000; odgi draw -i cyp2d6_locus.sort2.og -c cyp2d6_locus.sort2.og.lay -p CYP2D6_locus.sort2.og.lay.draw_hap.png -H 1000 -C -w 20 -H 1000

module load dSQ
dsq --job-file odgi.sh --mem-per-cpu 5G --cpus-per-task=8 -t 12:00:00 --mail-type None --partition=pi_hall -J odgi_cyp2d6 --max-jobs 1
dsqa -j 81058502  

rsync -arvzP $user@server:$mc_path/CYP2D6_locus.sort2.og.lay.draw.png $local_path

rsync -arvzP $user@server:$mc_path/CYP2D6_locus.sort2.og.lay.draw_hap.png $local_path

# odgi 1D # show structure: psv/ins/del $user@server
odgi viz -i cyp2d6_locus.sort.og -o cyp2d6_locus.odgi_viz.png
rsync -arvzP $user@server:$mc_path/cyp2d6_locus.odgi_viz.png $local_path

# doesn't work, need to change the w line into p line, or some other things, will do later..

===================
# pggb odgi viz 
odgi sort -P -t 16 -i chr22.pan.CYP2D6_locus.0.og -o chr22.pan.CYP2D6_locus.0.sort.og -O;
odgi sort -P -p Ygs -t 16 -i chr22.pan.CYP2D6_locus.0.sort.og -o chr22.pan.CYP2D6_locus.0.sort2.og

odgi viz -i chr22.pan.CYP2D6_locus.0.sort2.og -o chr22.pan.CYP2D6_locus.0.odgi_viz.png -m 

rsync -arvzP $user@server:$pggb_path/chr22.pan.CYP2D6_locus.0.odgi_viz.png $local_path

odgi viz -i chr22.pan.CYP2D6_locus.0.sort2.og -o chr22.pan.CYP2D6_locus.0.path_position.odgi_viz.png -u -d
rsync -arvzP $user@server:$pggb_path/chr22.pan.CYP2D6_locus.0.path_position.odgi_viz.png $local_path




odgi layout -i chr22.pan.more.CYP2D6_locus.sort2.og -o chr22.pan.more.CYP2D6_locus.sort2.og.lay -t 16 -P
odgi draw -i cyp2d6_locus.sort2.og -c cyp2d6_locus.sort2.og.lay -p CYP2D6_locus.sort2.og.lay.draw.png -H 1000; 
odgi draw -i cyp2d6_locus.sort2.og -c cyp2d6_locus.sort2.og.lay -p CYP2D6_locus.sort2.og.lay.draw_hap.png -H 1000 -C -w 20 -H 1000

odgi viz -i chr22.pan.CYP2D6_locus.0.og -o chr22.pan.CYP2D6_locus.0.odgi_viz.png

rsync -arvzP $user@server:$pggb_path/chr22.pan.CYP2D6_locus.odgi_viz.png $local_path





===================
# bandage color by path direction start to end 
(mark the PSV nodes)


# hg38 path 
# 200 color from blue to red
module load R/3.6.1-foss-2018b
# rainbow(8)
# colorRampPalette(c("blue", "red"))(1242)
R
write(colorRampPalette(c("red", "blue"))(1241), file = "color1241_red_blue.txt",
      ncolumns = 1, append = FALSE)
# write(rainbow(1242), file = "color1242_rainbow.txt",
#       ncolumns = 1, append = FALSE)

# change csv file from red to color changing 
less path_walk_cyp2d.gfa| sed -n '1p'| grep -o '76912359.*76914652'| sed 's/+,/\n/g'| sed 's/-,/\n/g'| wc -l  # 1241

paste -d "," <(less path_walk_cyp2d.gfa| sed -n '1p'| grep -o '76912359.*76914652'| sed 's/+,/\n/g'| sed 's/-,/\n/g') color1241_red_blue.txt > hg38.cyp2d6_locus.change.csv 
# vi change the frst line to Name,Colour

# change csv file from red to color rainbow 
# paste -d "," hg38.cyp2d6_locus.node.txt color1242_rainbow.txt|awk -F, '{print $1",\""$2"\""}' > hg38.cyp2d6_locus.rainbow.csv 
# vi change the frst line to Name,Colour

rsync -arvzP $user@server:$mc_path/hg38.cyp2d6_locus.change.csv $local_path

# rsync -arvzP $user@server:$mc_path/hg38.cyp2d6_locus.rainbow.csv $local_path





# chm13 pathes & other tipical pathes
# csv (bandage color)
less path_walk_cyp2d.gfa| sed -n '2p'| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1236
R 
write(colorRampPalette(c("red", "blue"))(1236), file = "color1236_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| sed -n '2p'| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" ) color1236_red_blue.txt > chm13.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour
rsync -arvzP $user@server:$mc_path/chm13.cyp2d6_locus.change.csv $local_path

# HG01891#1#JAGYVO010000117.1 path
# csv (bandage color)
less path_walk_cyp2d.gfa| grep JAGYVO010000117| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" | wc -l  # 887
R 
write(colorRampPalette(c("red", "blue"))(887), file = "color887_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAGYVO010000117| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" ) color887_red_blue.txt > HG01891_1_JAGYVO010000117.1.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour

rsync -arvzP $user@server:$mc_path/HG01891_1_JAGYVO010000117.1.cyp2d6_locus.change.csv $local_path


# HG02080#2#JAHEOV010000194.1 path
# csv (bandage color)
less path_walk_cyp2d.gfa| grep JAHEOV010000194| grep -o '76914652.*76912359'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1269
R 
write(colorRampPalette(c("red", "blue"))(1269), file = "color1269_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAHEOV010000194| grep -o '76914652.*76912359'| tr ">" "\n" | tr "<" "\n" ) color1269_red_blue.txt > HG02080_2_JAHEOV010000194.1.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour
rsync -arvzP $user@server:$mc_path/HG02080_2_JAHEOV010000194.1.cyp2d6_locus.change.csv $local_path


# HG03540#1#JAGYVY010000096.1 path
# csv (bandage color)
less path_walk_cyp2d.gfa| grep JAGYVY010000096| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1325
R 
write(colorRampPalette(c("red", "blue"))(1325), file = "color1325_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAGYVY010000096| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" ) color1325_red_blue.txt > HG03540_1_JAGYVY010000096.1.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour
rsync -arvzP $user@server:$mc_path/HG03540_1_JAGYVY010000096.1.cyp2d6_locus.change.csv $local_path


# HG00733#1#JAHEPQ010000046.1 path
# csv (bandage color)
less path_walk_cyp2d.gfa| grep JAHEPQ010000046| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1134
R 
write(colorRampPalette(c("red", "blue"))(1134), file = "color1134_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAHEPQ010000046| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" ) color1134_red_blue.txt > HG00733_1_JAHEPQ010000046.1.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour

rsync -arvzP $user@server:$mc_path/HG00733_1_JAHEPQ010000046.1.cyp2d6_locus.change.csv $local_path


# HG00438#2#JAHBCA010000050.1 path 
# csv (bandage color)
less path_walk_cyp2d.gfa| grep JAHBCA010000050| grep -o '76914652.*76912359'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1127
R 
write(colorRampPalette(c("red", "blue"))(1127), file = "color1127_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAHBCA010000050| grep -o '76914652.*76912359'| tr ">" "\n" | tr "<" "\n") color1127_red_blue.txt > HG00438_2_JAHBCA010000050.1.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour
rsync -arvzP $user@server:$mc_path/HG00438_2_JAHBCA010000050.1.cyp2d6_locus.change.csv $local_path

# HG00438#2#JAHBCA010000050.1 path - color reverse
less path_walk_cyp2d.gfa| grep JAHBCA010000050| grep -o '76914652.*76912359'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1127
R 
write(colorRampPalette(c("blue", "red"))(1127), file = "color1127_blue_red.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAHBCA010000050| grep -o '76914652.*76912359'| tr ">" "\n" | tr "<" "\n") color1127_blue_red.txt > HG00438_2_JAHBCA010000050.1.cyp2d6_locus.change.rev.csv
# vi change the frst line to Name,Colour
rsync -arvzP $user@server:$mc_path/HG00438_2_JAHBCA010000050.1.cyp2d6_locus.change.rev.csv $local_path




# HG01258#1#JAGYYV010000053.1
# csv (bandage color)
less path_walk_cyp2d.gfa| grep JAGYYV010000053| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n" | wc -l  # 1493
R 
write(colorRampPalette(c("red", "blue"))(1493), file = "color1493_red_blue.txt",
      ncolumns = 1, append = FALSE)
paste -d "," <(less path_walk_cyp2d.gfa| grep JAGYYV010000053| grep -o '76912359.*76914652'| tr ">" "\n" | tr "<" "\n") color1493_red_blue.txt > HG01258_1_JAGYYV010000053.1.cyp2d6_locus.change.csv
# vi change the frst line to Name,Colour
rsync -arvzP $user@server:$mc_path/HG01258_1_JAGYYV010000053.1.cyp2d6_locus.change.csv $local_path




===================
# bandage color by cyp2d6 cyp2d7 graph alignment  mini-cactus

# mini_cactus path with all diff&uniq cyp2d6,cyp2d7 nodes

c2 cyp2d6(3) 21/ cyp2d7(59) 155
c1 cyp2d6(2) 5
s  cyp2d6(1) 0/ cyp2d7(12) 1
c3 cyp2d6(7) 237/ cyp2d7(11) 0

# combine cyp2d6/cyp2d7/spacer nodes together
# cyp2d6 s
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' > bandage.gene.nodes.txt
# 200 color from blue to red
module load R/3.6.1-foss-2018b
R
write(rainbow(324), file = "color324_rainbow.txt",ncolumns = 1, append = FALSE)
less color324_rainbow.txt| cut -c1-7 > color324_rainbow.mod.txt
paste -d "," bandage.gene.nodes.txt color324_rainbow.mod.txt|awk -F, '{print $1","$2}' > hg38.cyp2d6_gene.rainbow.csv 
# cyp2d6 c1
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '2p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' > bandage.gene.nodes.txt
R
write(rainbow(37), file = "color37_rainbow.txt",ncolumns = 1, append = FALSE)
less color37_rainbow.txt| cut -c1-7 > color37_rainbow.mod.txt
cat hg38.cyp2d6_gene.rainbow.csv  <(paste -d "," bandage.gene.nodes.txt color37_rainbow.mod.txt) > hg38.cyp2d6_gene.2.rainbow.csv
# cyp2d6 c2
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '3p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' > bandage.gene.nodes.txt
R
write(rainbow(155), file = "color155_rainbow.txt",ncolumns = 1, append = FALSE)
less color155_rainbow.txt| cut -c1-7 > color155_rainbow.mod.txt
cat hg38.cyp2d6_gene.2.rainbow.csv  <(paste -d "," bandage.gene.nodes.txt color155_rainbow.mod.txt) > hg38.cyp2d6_gene.3.rainbow.csv
# cyp2d7 c3
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '7p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' > bandage.gene.nodes.txt
less bandage.gene.nodes.txt| wc -l
R
write(rainbow(172), file = "color172_rainbow.txt",ncolumns = 1, append = FALSE)
less color172_rainbow.txt| cut -c1-7 > color172_rainbow.mod.txt
cat hg38.cyp2d6_gene.3.rainbow.csv  <(paste -d "," bandage.gene.nodes.txt color172_rainbow.mod.txt) > hg38.cyp2d6_gene.4.rainbow.csv
# vi hg38.cyp2d6_gene.4.rainbow.csv change the frst line to Name,Colour

rsync -arvzP $user@server:$mc_path/hg38.cyp2d6_gene.4.rainbow.csv $local_path



===================
# bandage color by cyp2d6 cyp2d7 graph alignment pggb 
$pggb_path

less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d'| wc -l  # 424
R 
write(rainbow(424), file = "color424_rainbow.txt",ncolumns = 1, append = FALSE)
paste -d "," <(less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d' ) <(less color424_rainbow.txt| cut -c1-7) > hg38.cyp2d6_gene.rainbow.csv
# vi change the frst line to Name,Colour

rsync -arvzP $user@server:$pggb_path/hg38.cyp2d6_gene.rainbow.csv $local_path


===================
# bandage color by hg38 path pggb 
$pggb_path

less -S chr22.pan.CYP2D6_locus.pline.txt|grep grch38#chr22| grep 345128| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'| wc -l # 1604
module load R/3.6.1-foss-2018b
R
write(colorRampPalette(c("red", "blue"))(1604), file = "color1604_red_blue.txt",
      ncolumns = 1, append = FALSE)

paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep grch38#chr22| grep 345128| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g') color1604_red_blue.txt > hg38.cyp2d6_locus.change.csv 
# vi change the frst line to Name,Colour
# vi remove - in last row
rsync -arvzP $user@server:$pggb_path/hg38.cyp2d6_locus.change.csv $local_path





# hg38 path 
# 200 color from blue to red
module load R/3.6.1-foss-2018b
# rainbow(8)
# colorRampPalette(c("blue", "red"))(1242)
R
write(colorRampPalette(c("red", "blue"))(1241), file = "color1241_red_blue.txt",
      ncolumns = 1, append = FALSE)
# write(rainbow(1242), file = "color1242_rainbow.txt",
#       ncolumns = 1, append = FALSE)

# change csv file from red to color changing 
less path_walk_cyp2d.gfa| sed -n '1p'| grep -o '76912359.*76914652'| sed 's/+,/\n/g'| sed 's/-,/\n/g'| wc -l  # 1241

paste -d "," <(less path_walk_cyp2d.gfa| sed -n '1p'| grep -o '76912359.*76914652'| sed 's/+,/\n/g'| sed 's/-,/\n/g') color1241_red_blue.txt > hg38.cyp2d6_locus.change.csv 
# vi change the frst line to Name,Colour

# change csv file from red to color rainbow 
# paste -d "," hg38.cyp2d6_locus.node.txt color1242_rainbow.txt|awk -F, '{print $1",\""$2"\""}' > hg38.cyp2d6_locus.rainbow.csv 
# vi change the frst line to Name,Colour

rsync -arvzP $user@server:$mc_path/hg38.cyp2d6_locus.change.csv $local_path

# rsync -arvzP $user@server:$mc_path/hg38.cyp2d6_locus.rainbow.csv $local_path



===================
# linear visualization mini
# align hg38 cyp2d6/cyp2d7 to assemblies 
$mc_path/align2assembly


module load SAMtools/1.13-GCCcore-10.2.0

# HG01891#1#JAGYVO010000117.1
#select chromosomes from genbank assemblies
samtools faidx ./data/94assembly_v2/assemblies/HG01891.paternal.f1_assembly_v2_genbank.fa HG01891#1#JAGYVO010000117.1 > HG01891_1_JAGYVO010000117.1.fa

# cyp2d6/ cyp2d7 / spacer fa
$pggb_path/cyp2d6_cyp2d7_spacer.fasta

#index gene fasta
samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta
# cyp2d6 alignment 
minimap2 -ax sr HG01891_1_JAGYVO010000117.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) | samtools sort | samtools view -S -b > cyp2d6_HG01891_1_JAGYVO010000117.1.sorted.bam
samtools index cyp2d6_HG01891_1_JAGYVO010000117.1.sorted.bam
# cyp2d7 alignment 
minimap2 -ax sr HG01891_1_JAGYVO010000117.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_chr22:42,140,203-42,149,455 ) | samtools sort | samtools view -S -b > cyp2d7_HG01891_1_JAGYVO010000117.1.sorted.bam
samtools index cyp2d7_HG01891_1_JAGYVO010000117.1.sorted.bam
# cyp2d7 space 
minimap2 -ax sr HG01891_1_JAGYVO010000117.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_spacer_chr22:42,138,124-42,139,676 ) | samtools sort | samtools view -S -b > cyp2d7_spacer_HG01891_1_JAGYVO010000117.1.sorted.bam
samtools index cyp2d7_spacer_HG01891_1_JAGYVO010000117.1.sorted.bam


rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6_HG01891_1_JAGYVO010000117.1.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6_HG01891_1_JAGYVO010000117.1.sorted.bam.bai $local_path/align2assembly

rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7_HG01891_1_JAGYVO010000117.1.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7_HG01891_1_JAGYVO010000117.1.sorted.bam.bai $local_path/align2assembly

rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7_spacer_HG01891_1_JAGYVO010000117.1.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7_spacer_HG01891_1_JAGYVO010000117.1.sorted.bam.bai $local_path/align2assembly






# HG00438#2#JAHBCA010000050.1

#select chromosomes from genbank assemblies
samtools faidx ./data/94assembly_v2/assemblies/HG00438.maternal.f1_assembly_v2_genbank.fa HG00438#2#JAHBCA010000050.1 > HG00438_2_JAHBCA010000050.1.fa

# cyp2d6/ cyp2d7 / spacer fa
$pggb_path/cyp2d6_cyp2d7_spacer.fasta

# change parameter 
# cyp2d6 alignment 
minimap2 -ax asm5 HG00438_2_JAHBCA010000050.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG00438_2_JAHBCA010000050.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d6.HG00438_2_JAHBCA010000050.1.asm5.sam| samtools view -S -b > cyp2d6.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam
samtools index cyp2d6.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam.bai $local_path/align2assembly

# cyp2d7 alignment 
minimap2 -ax asm5 HG00438_2_JAHBCA010000050.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_chr22:42,140,203-42,149,455 ) > cyp2d7.HG00438_2_JAHBCA010000050.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d7.HG00438_2_JAHBCA010000050.1.asm5.sam| samtools view -S -b > cyp2d7.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam
samtools index cyp2d7.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam.bai $local_path/align2assembly

# spacer alignment 
minimap2 -ax asm5 HG00438_2_JAHBCA010000050.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_spacer_chr22:42,138,124-42,139,676 ) > spacer.HG00438_2_JAHBCA010000050.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort spacer.HG00438_2_JAHBCA010000050.1.asm5.sam| samtools view -S -b > spacer.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam
samtools index spacer.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG00438_2_JAHBCA010000050.1.asm5.sorted.bam.bai $local_path/align2assembly



# minimap2 -ax asm5 HG00438_2_JAHBCA010000050.1.fa $pggb_path/cyp2d6_cyp2d7_spacer.fasta >  cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm5.sorted.sam
# samtools sort cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm5.sorted.sam| samtools view -S -b > cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm5.sorted.2.bam
# samtools index cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm5.sorted.2.bam


# rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm5.sorted.2.bam $local_path/align2assembly
# rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm5.sorted.2.bam.bai $local_path/align2assembly


# minimap2 -ax asm10 HG00438_2_JAHBCA010000050.1.fa $pggb_path/cyp2d6_cyp2d7_spacer.fasta | samtools sort | samtools view -S -b > cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm10.sorted.bam

# minimap2 -ax asm20 HG00438_2_JAHBCA010000050.1.fa $pggb_path/cyp2d6_cyp2d7_spacer.fasta | samtools sort | samtools view -S -b > cyp2d6_7_s.HG00438_2_JAHBCA010000050.1.asm20.sorted.bam




# HG01258#1#JAGYYV010000053.1

#select chromosomes from genbank assemblies
samtools faidx ./data/94assembly_v2/assemblies/HG01258.paternal.f1_assembly_v2_genbank.fa HG01258#1#JAGYYV010000053.1 > HG01258_1_JAGYYV010000053.1.fa
rsync -arvzP $user@server:$mc_path/align2assembly/HG01258_1_JAGYYV010000053.1.fa $local_path/align2assembly

# cyp2d6/ cyp2d7 / spacer fa
$pggb_path/cyp2d6_cyp2d7_spacer.fasta

# change parameter 
# cyp2d6 alignment 
# minimap2 -ax asm20 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sam

minimap2 -ax sr --secondary=yes -p 0.1 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sam| samtools view -S -b > cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam
samtools index cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam.bai $local_path/align2assembly

# cyp2d7 alignment 
minimap2 -ax sr --secondary=yes -p 0.1 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_chr22:42,140,203-42,149,455 ) > cyp2d7.HG01258_1_JAGYYV010000053.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d7.HG01258_1_JAGYYV010000053.1.asm5.sam| samtools view -S -b > cyp2d7.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam
samtools index cyp2d7.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam.bai $local_path/align2assembly

# spacer alignment 
minimap2 -ax sr --secondary=yes -p 0.1 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_spacer_chr22:42,138,124-42,139,676 ) > spacer.HG01258_1_JAGYYV010000053.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort spacer.HG01258_1_JAGYYV010000053.1.asm5.sam| samtools view -S -b > spacer.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam
samtools index spacer.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG01258_1_JAGYYV010000053.1.asm5.sorted.bam.bai $local_path/align2assembly




# HG00733#1#JAHEPQ010000046.1

#select chromosomes from genbank assemblies
samtools faidx ./data/94assembly_v2/assemblies/HG00733.paternal.f1_assembly_v2_genbank.fa HG00733#1#JAHEPQ010000046.1 > HG00733_1_JAHEPQ010000046.1.fa
rsync -arvzP $user@server:$mc_path/align2assembly/HG00733_1_JAHEPQ010000046.1.fa $local_path/align2assembly

# cyp2d6/ cyp2d7 / spacer fa
$pggb_path/cyp2d6_cyp2d7_spacer.fasta

# change parameter 
# cyp2d6 alignment 
# minimap2 -ax asm20 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sam

minimap2 -ax sr --secondary=yes -p 0.1 HG00733_1_JAHEPQ010000046.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG00733_1_JAHEPQ010000046.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d6.HG00733_1_JAHEPQ010000046.1.asm5.sam| samtools view -S -b > cyp2d6.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam
samtools index cyp2d6.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam.bai $local_path/align2assembly

# cyp2d7 alignment 
minimap2 -ax sr --secondary=yes -p 0.1 HG00733_1_JAHEPQ010000046.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_chr22:42,140,203-42,149,455 ) > cyp2d7.HG00733_1_JAHEPQ010000046.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d7.HG00733_1_JAHEPQ010000046.1.asm5.sam| samtools view -S -b > cyp2d7.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam
samtools index cyp2d7.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam.bai $local_path/align2assembly

# spacer alignment 
minimap2 -ax sr --secondary=yes -p 0.1 HG00733_1_JAHEPQ010000046.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_spacer_chr22:42,138,124-42,139,676 ) > spacer.HG00733_1_JAHEPQ010000046.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort spacer.HG00733_1_JAHEPQ010000046.1.asm5.sam| samtools view -S -b > spacer.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam
samtools index spacer.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG00733_1_JAHEPQ010000046.1.asm5.sorted.bam.bai $local_path/align2assembly



# HG03540#1#JAGYVY010000096.1 (2)


#select chromosomes from genbank assemblies
samtools faidx ./data/94assembly_v2/assemblies/HG03540.paternal.f1_assembly_v2_genbank.fa HG03540#1#JAGYVY010000096.1 > HG03540_1_JAGYVY010000096.1.fa
rsync -arvzP $user@server:$mc_path/align2assembly/HG03540_1_JAGYVY010000096.1.fa $local_path/align2assembly

# cyp2d6/ cyp2d7 / spacer fa
$pggb_path/cyp2d6_cyp2d7_spacer.fasta

# change parameter 
# cyp2d6 alignment 
# minimap2 -ax asm20 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sam

minimap2 -ax sr --secondary=yes -p 0.1 HG03540_1_JAGYVY010000096.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG03540_1_JAGYVY010000096.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d6.HG03540_1_JAGYVY010000096.1.asm5.sam| samtools view -S -b > cyp2d6.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam
samtools index cyp2d6.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam.bai $local_path/align2assembly

# cyp2d7 alignment 
minimap2 -ax sr --secondary=yes -p 0.1 HG03540_1_JAGYVY010000096.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_chr22:42,140,203-42,149,455 ) > cyp2d7.HG03540_1_JAGYVY010000096.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d7.HG03540_1_JAGYVY010000096.1.asm5.sam| samtools view -S -b > cyp2d7.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam
samtools index cyp2d7.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam.bai $local_path/align2assembly

# spacer alignment 
minimap2 -ax sr --secondary=yes -p 0.1 HG03540_1_JAGYVY010000096.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_spacer_chr22:42,138,124-42,139,676 ) > spacer.HG03540_1_JAGYVY010000096.1.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort spacer.HG03540_1_JAGYVY010000096.1.asm5.sam| samtools view -S -b > spacer.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam
samtools index spacer.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.HG03540_1_JAGYVY010000096.1.asm5.sorted.bam.bai $local_path/align2assembly





# hg38 path (80 asm)

#select chromosomes from genbank assemblies
samtools faidx ./reference/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta chr22 > chr22.fa
rsync -arvzP $user@server:$mc_path/align2assembly/chr22.fa $local_path/align2assembly

# cyp2d6/ cyp2d7 / spacer fa
$pggb_path/cyp2d6_cyp2d7_spacer.fasta

# change parameter 
# cyp2d6 alignment 
# minimap2 -ax asm20 HG01258_1_JAGYYV010000053.1.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.HG01258_1_JAGYYV010000053.1.asm5.sam

minimap2 -ax sr --secondary=yes -p 0.1 chr22.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d6_chr22:42,126,499-42,130,865 ) > cyp2d6.chr22.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d6.chr22.asm5.sam| samtools view -S -b > cyp2d6.chr22.asm5.sorted.bam
samtools index cyp2d6.chr22.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.chr22.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d6.chr22.asm5.sorted.bam.bai $local_path/align2assembly

# cyp2d7 alignment 
minimap2 -ax sr --secondary=yes -p 0.1 chr22.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_chr22:42,140,203-42,149,455 ) > cyp2d7.chr22.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort cyp2d7.chr22.asm5.sam| samtools view -S -b > cyp2d7.chr22.asm5.sorted.bam
samtools index cyp2d7.chr22.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.chr22.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/cyp2d7.chr22.asm5.sorted.bam.bai $local_path/align2assembly

# spacer alignment 
minimap2 -ax sr --secondary=yes -p 0.1 chr22.fa <(samtools faidx $pggb_path/cyp2d6_cyp2d7_spacer.fasta cyp2d7_spacer_chr22:42,138,124-42,139,676 ) > spacer.chr22.asm5.sam
# vi 
# change seuqence from * to original sequence
samtools sort spacer.chr22.asm5.sam| samtools view -S -b > spacer.chr22.asm5.sorted.bam
samtools index spacer.chr22.asm5.sorted.bam
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.chr22.asm5.sorted.bam $local_path/align2assembly
rsync -arvzP $user@server:$mc_path/align2assembly/spacer.chr22.asm5.sorted.bam.bai $local_path/align2assembly





===================
# color by Green and blue 
$mc_path

# mini_cactus path with all diff&uniq cyp2d6,cyp2d7 nodes

c2 cyp2d6(3) 21/ cyp2d7(59) 155
c1 cyp2d6(2) 5
s  cyp2d6(1) 0/ cyp2d7(12) 1
c3 cyp2d6(7) 237/ cyp2d7(11) 0

# combine cyp2d6/cyp2d7/spacer nodes together
# cyp2d6 s
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' | wc -l #324
# cyp2d6 c1
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '2p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' | wc -l #37
# cyp2d6 c2
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '3p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' | wc -l #155
# cyp2d7 c3
less cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '7p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d' | wc -l #172

$mc_path/green_blue
# color
python3
from colour import Color
green = Color("green")
colors = list(green.range_to(Color("blue"),324))
output_file = 'color324_green_blue.txt'
with open(output_file,'w') as f:
	for i in colors:
		f.write(str(i) + "\n")

colors = list(green.range_to(Color("blue"),37))
output_file = 'color37_green_blue.txt'
with open(output_file,'w') as f:
	for i in colors:
		f.write(str(i) + "\n")

colors = list(green.range_to(Color("blue"),155))
output_file = 'color155_green_blue.txt'
with open(output_file,'w') as f:
	for i in colors:
		f.write(str(i) + "\n")

colors = list(green.range_to(Color("blue"),172))
output_file = 'color172_green_blue.txt'
with open(output_file,'w') as f:
	for i in colors:
		f.write(str(i) + "\n")

echo Name,Colour > hg38.CYP2D_gene.green_blue.csv
paste -d "," <(less ../cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d') <(less color324_green_blue.txt) >> hg38.CYP2D_gene.green_blue.csv
paste -d "," <(less ../cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '2p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d') <(less color37_green_blue.txt) >> hg38.CYP2D_gene.green_blue.csv
paste -d "," <(less ../cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '3p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d') <(less color155_green_blue.txt) >> hg38.CYP2D_gene.green_blue.csv
paste -d "," <(less ../cyp2d6_cyp2d7_spacer.aln.more.gaf| sed -n '7p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' | sed '1d') <(less color172_green_blue.txt) >> hg38.CYP2D_gene.green_blue.csv

rsync -arvzP $user@server:$mc_path/green_blue/hg38.CYP2D_gene.green_blue.csv $local_path



===================
# change to green-blue gradient color - PGGB
$pggb_path
# GraphAligner -g chr22.pan.CYP2D6_locus.gfa -f cyp2d6_cyp2d7_spacer.fasta -a cyp2d6_cyp2d7_spacer.aln.gaf -x vg

less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d'| wc -l  # 424

# color
python3
from colour import Color
green = Color("green")
colors = list(green.range_to(Color("blue"),424))
output_file = 'color424_green_blue.txt'
with open(output_file,'w') as f:
	for i in colors:
		f.write(str(i) + "\n")

echo Name,Colour > hg38.CYP2D_gene.green_blue.csv
paste -d "," <(less cyp2d6_cyp2d7_spacer.aln.gaf| sed -n '1p' | cut -f6 | sed 's/>/\n/g' |sed 's/</\n/g' |sed '1d') <(less color424_green_blue.txt) >> hg38.CYP2D_gene.green_blue.csv

rsync -arvzP $user@server:$pggb_path/hg38.CYP2D_gene.green_blue.csv $local_path


===================
# bandage color path PGGB 
$pggb_path

# GRCh38
less -S chr22.pan.CYP2D6_locus.pline.txt|grep grch38| cut -f3| grep -o "346438.*343650"| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 1540

module load R/4.1.0-foss-2020b
R
write(colorRampPalette(c("red", "blue"))(1540), file = "color1540_red_blue.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > hg38.cyp2d_locus.change.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep grch38| cut -f3| grep -o "346438.*343650"| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color1540_red_blue.txt >> hg38.cyp2d_locus.change.csv 

rsync -arvzP $user@server:$pggb_path/hg38.cyp2d_locus.change.csv $local_path

# HG01891#1#JAGYVO010000117.1
less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAGYVO010000117| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 1771
R 
write(colorRampPalette(c("red", "blue"))(1771), file = "color1771_red_blue.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > HG01891_1_JAGYVO010000117.1.cyp2d_locus.path.change.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAGYVO010000117| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color1771_red_blue.txt >> HG01891_1_JAGYVO010000117.1.cyp2d_locus.path.change.csv 

rsync -arvzP $user@server:$pggb_path/HG01891_1_JAGYVO010000117.1.cyp2d_locus.path.change.csv $local_path

# HG03540#1#JAGYVY010000096.1
less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAGYVY010000096| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 3060
R 
write(colorRampPalette(c("red", "blue"))(3060), file = "color3060_red_blue.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > HG03540_1_JAGYVY010000096.1.cyp2d_locus.path.change.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAGYVY010000096| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color3060_red_blue.txt >> HG03540_1_JAGYVY010000096.1.cyp2d_locus.path.change.csv 

rsync -arvzP $user@server:$pggb_path/HG03540_1_JAGYVY010000096.1.cyp2d_locus.path.change.csv $local_path

# HG00733#1#JAHEPQ010000046.1
less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAHEPQ010000046| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 3110
R 
write(colorRampPalette(c("red", "blue"))(3110), file = "color3110_red_blue.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > HG00733_1_JAHEPQ010000046.1.cyp2d_locus.path.change.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAHEPQ010000046| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color3110_red_blue.txt >> HG00733_1_JAHEPQ010000046.1.cyp2d_locus.path.change.csv 

rsync -arvzP $user@server:$pggb_path/HG00733_1_JAHEPQ010000046.1.cyp2d_locus.path.change.csv $local_path

# HG00438#2#JAHBCA010000050.1
less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAHBCA010000050| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 3099
R 
write(colorRampPalette(c("red", "blue"))(3099), file = "color3099_red_blue.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAHBCA010000050| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color3099_red_blue.txt >> HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.csv 

rsync -arvzP $user@server:$pggb_path/HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.csv $local_path


# HG00438#2#JAHBCA010000050.1 - color reverse
less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAHBCA010000050| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 3099
R 
write(colorRampPalette(c("blue", "red"))(3099), file = "color3099_blue_red.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAHBCA010000050| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color3099_blue_red.txt >> HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.csv 

rsync -arvzP $user@server:$pggb_path/HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.csv $local_path



# HG01258#1#JAGYYV010000053.1
less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAGYYV010000053| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g'| wc -l # 4599
R 
write(colorRampPalette(c("red", "blue"))(4599), file = "color4599_red_blue.txt",
      ncolumns = 1, append = FALSE)

echo Name,Colour > HG01258_1_JAGYYV010000053.1.cyp2d_locus.path.change.csv 
paste -d "," <(less -S chr22.pan.CYP2D6_locus.pline.txt|grep JAGYYV010000053| cut -f3| sed 's/+,/\n/g'| sed 's/-,/\n/g'|sed 's/+//g'|sed 's/-//g') color4599_red_blue.txt >> HG01258_1_JAGYYV010000053.1.cyp2d_locus.path.change.csv 

rsync -arvzP $user@server:$pggb_path/HG01258_1_JAGYYV010000053.1.cyp2d_locus.path.change.csv $local_path



# change order head to tail on my mac 
$local_path
HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.csv
echo Name,Colour > HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.order.csv
less HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.csv| sed '1d' |tac >> HG00438_2_JAHBCA010000050.1.cyp2d_locus.path.change.rev.order.csv








