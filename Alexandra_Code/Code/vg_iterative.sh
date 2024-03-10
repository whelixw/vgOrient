#!/usr/bin/env bash


echo $(date)
fasta_filename=$1;
output_graph_name=$2;

# Make temporary directory
tmp_dir=$(mktemp -d -t vg-XXXXXXXXXX)
echo "Temporary directory: $tmp_dir"

cd $tmp_dir
mkdir fasta_files
mkdir alignments


# Split fasta file into individual sequences and save to the fasta_files/ dir
# Count number of sequences in fasta file
if [[ $fastafilename == *.gz ]]; then
    nr=$(zcat $fasta_filename | grep '>' -c)
    zcat $fasta_filename | awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1==0){file=sprintf("fasta_files/seq_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }'
else
    nr=$(grep '>' $fasta_filename -c)
    cat $fasta_filename | awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1==0){file=sprintf("fasta_files/seq_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }'
fi

echo "Number of sequences in total: " $nr


cd fasta_files

# Choose a random file/sequence
random_file=$(ls | shuf -n 1)
name=$(grep -oE '^>[a-zA-Z0-9_.]+' ${random_file}) # The accession number
name="${name:1}"
echo "Choosing file $random_file ($name) to build graph from"

cd ..

# Construct vg graph from the chosen file and circularize it
/home/ctools/vg/bin/vg construct -M fasta_files/$random_file > graph.vg
/home/ctools/vg/bin/vg circularize -p $name graph.vg > graph_circ.vg  
/home/ctools/vg/bin/vg stats -z graph_circ.vg 


rm graph.vg
rm fasta_files/$random_file 
cd fasta_files


for ((i = 1; i <= $nr-1; i++))
do

	identity=0

	# Index graph
	/home/ctools/vg/bin/vg index -x $tmp_dir/graph_circ.xg $tmp_dir/graph_circ.vg # XG index
	/home/ctools/vg/bin/vg prune -k 48 $tmp_dir/graph_circ.vg > $tmp_dir/graph_circ_pruned.vg 
	/home/ctools/vg/bin/vg index -g $tmp_dir/graph_circ.gcsa -Z 400 $tmp_dir/graph_circ_pruned.vg # GCSA index
	rm -f $tmp_dir/graph_circ_pruned.vg

	for FILE in *; do

		# Extract path name
		name=$(grep -oE '^>[a-zA-Z0-9_.]+' ${FILE})
		name="${name:1}"

		# convert mitogenome to a string and align to graph
		tail -n +2 ${FILE} | tr -d '\n' > mitogenome_str  
		/home/ctools/vg/bin/vg map -s $(< mitogenome_str) -V $name -g $tmp_dir/graph_circ.gcsa -x $tmp_dir/graph_circ.xg > $tmp_dir/alignments/${FILE}.gam

		identity_new=$(/home/ctools/vg/bin/vg view -a $tmp_dir/alignments/${FILE}.gam | jq ".identity" -e)
		
		# Set identity to 0.0 if it is "null"
		if [[ $? -ne 0 ]]; then
			echo "it is null!"
			identity_new=0.0
			#echo "null set to 0"
		fi

		# Save the file with the highest %identity
		if [[ $identity_new > $identity ]]; then
			identity=$identity_new
			best_file=${FILE}
			best_name=$name
		fi

		echo "${FILE} ($name) aligned with identity score $identity_new"

	done

	echo "$best_file has the highest %identity ($identity)" 


	# Augment graph with that sequence
	/home/ctools/vg/bin/vg augment -i -S $tmp_dir/graph_circ.vg $tmp_dir/alignments/${best_file}.gam  > $tmp_dir/graph_aug.vg

	rm $tmp_dir/fasta_files/$best_file
	rm $tmp_dir/fasta_files/mitogenome_str
	mv $tmp_dir/graph_aug.vg $tmp_dir/graph_circ.vg

	/home/ctools/vg/bin/vg stats -z $tmp_dir/graph_circ.vg

	# Clear the alignments dir for the next round
	rm $tmp_dir/alignments/*


done


mv $tmp_dir/graph_circ.vg /home/projects/beyondMARS/iterative/${output_graph_name}.vg

cd /home/projects/beyondMARS/iterative/

# Convert to ODGI and GFA formats
/home/ctools/vg/bin/vg convert -o $output_graph > ${output_graph_name}.odgi
/home/ctools/vg/bin/vg convert -f $output_graph > ${output_graph_name}.gfa


rm -rf $tmp_dir
echo $(date)


