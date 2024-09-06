#!/bin/bash

# Activate the .bashrc source file so that we can use the conda environments in the script.
source ~/.bashrc

input_dir="/home/kevin.muirhead/phytoplasma_hybrid_genes/gene_assemblies"
#input_dir="/Users/kevin.muirhead/Desktop/phytoplasma_hybrid_genes/gene_assemblies"

output_dir="/home/kevin.muirhead/phytoplasma_hybrid_genes"
#output_dir="/Users/kevin.muirhead/Desktop/phytoplasma_hybrid_genes"

mkdir -p $output_dir

gene_output_dir="${output_dir}/gene_sequences"
mkdir -p $gene_output_dir

# The minimum sequence length for the ORFs.
min_size=500

# The maximum sequence length for the ORFs.
max_size=4000

# Find the top three longest ORFs.
orf_find=3

# The strain name
strain_name="BbSP"

# The organism name.
organism_name="Candidatus Phytoplasma asteris"

# List of gene names to match in the hybrid assembly files and order you want the sequences to be concatenated.
gene_list=("cpn60" "secY" "secA" "nusA" "tuf")

# Join strings by delimeter character.
join_by_string() {
  local separator="$1"
  shift
  local first="$1"
  shift
  printf "%s" "$first" "${@/#/$separator}"
}

# Get the gene regular expression string. i.e. "cpn60\|secY\|secA\|nusA\|tuf"
gene_regex=$(join_by_string "\|" "${gene_list[@]}")

# The list of hybrid assembly file paths.
assembly_list_file="${input_dir}/hybrid_assembly_files.txt"

# Find all the hybrid assembly files and write the file paths to the assembly_list_file.
#find $input_dir -name "*\.txt" | grep "cpn60\|secY\|secA\|nusA\|tuf" > $assembly_list_file
find $input_dir -name "*\.txt" -type f | grep "${gene_regex}" > $assembly_list_file

# Iterate over the hybrid assembly files.
for fasta_file in $(cat $assembly_list_file); 
do 

	# Get the gene name. Use basename to get the filename. Remove the ".txt" extension using sed. Grab 	
	gene_name=$(basename $fasta_file | sed 's/\.txt//g' | sed 's/.*\(cpn60\|secY\|secA\|nusA\|tuf\).*/\1/g');
	
	echo $gene_name;

	#exit 0;
	
	# The path of the gene output file.
	gene_fasta_outfile="${gene_output_dir}/${gene_name}.fasta"
	
	echo $gene_fasta_outfile

	# Activate the EMBOSS conda environment.
	conda activate emboss_env

	# Use the getorf command to get ORFs from each assembly.
	echo -e "getorf -sequence $fasta_file -minsize $min_size -maxsize $max_size -find $orf_find -methionine N -outseq $gene_fasta_outfile"
	getorf -sequence $fasta_file -minsize $min_size -maxsize $max_size -find $orf_find -methionine N -outseq $gene_fasta_outfile

done

gene_fasta_list_file="${gene_output_dir}/hybrid_genes_files.txt"

# The list of hybrid gene file paths.
find $gene_output_dir -name "*\.fasta" -type f | grep "${gene_regex}" > $gene_fasta_list_file

# Get a concatenated ordered list of gene names separated by commas.
gene_name_list=$(join_by_string "," "${gene_list[@]}")

# Activate the biopython conda environment.
conda activate biopython_env

echo "Concatinating ${gene_name_list} files..."

# Concatenate hybrid gene panel sequences in order of the gene_name_list.
echo -e "python concat_seqs_order_by_sample.py --fasta_file_list_infile $gene_fasta_list_file --gene_name_list $gene_name_list --strain_name $strain_name --organism_name  $organism_name --output_dir $output_dir"
python concat_seqs_order_by_sample.py --fasta_file_list_infile $gene_fasta_list_file --gene_name_list $gene_name_list --strain_name "\"${strain_name}\"" --organism_name "\"${organism_name}\"" --output_dir $output_dir



