#!/usr/bin/env bash


cd /home/projects/beyondMARS/animals10/
cat /home/projects/beyondMARS/data/10_animals.fasta | awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1==0){file=sprintf("seq_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }'


name=$(grep -oE '^>[a-zA-Z0-9_.]+' seq_3.fa) # The accession number
name="${name:1}"


# Construct vg graph from the chosen file and circularize it
/home/ctools/vg/bin/vg construct -M seq_3.fa > graph.vg
/home/ctools/vg/bin/vg circularize -p $name graph.vg > graph_circ.vg  
/home/ctools/vg/bin/vg stats -z graph_circ.vg # Nodes and edges
rm graph.vg
rm seq_3.fa


for d in {2,1,4,6,9,5,8,7,0}  
do

	# Index graph
	/home/ctools/vg/bin/vg index -x graph_circ.xg graph_circ.vg 
	/home/ctools/vg/bin/vg prune -k 48 graph_circ.vg > graph_circ_pruned.vg 
	/home/ctools/vg/bin/vg index -g graph_circ.gcsa -Z 400 graph_circ_pruned.vg 
	rm -f graph_circ_pruned.vg


	name=$(grep -oE '^>[a-zA-Z0-9_.]+' seq_${d}.fa) 
	name="${name:1}"

	# Make it a string and align to graph
	tail -n +2 seq_${d}.fa | tr -d '\n' > mitogenome_str  
	/home/ctools/vg/bin/vg map -s $(< mitogenome_str) -V $name -g graph_circ.gcsa -x graph_circ.xg > seq_${d}.gam
	/home/ctools/vg/bin/vg view -a seq_${d}.gam | jq ".identity"
	/home/ctools/vg/bin/vg augment -i -S graph_circ.vg seq_${d}.gam > graph_aug.vg

	rm mitogenome_str
	mv graph_aug.vg graph_circ.vg

	rm seq_${d}.gam
	rm seq_${d}.fa

	/home/ctools/vg/bin/vg stats -z graph_circ.vg

done


# Convert to ODGI and GFA formats
/home/ctools/vg/bin/vg convert -o graph_circ.vg > graph_circ.odgi
/home/ctools/vg/bin/vg convert -f graph_circ.vg > graph_circ.gfa

rm graph_circ.xg
rm graph_circ.gcsa.lcp
rm graph_circ.gcsa 
