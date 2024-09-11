#!/bin/bash
#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G
#SBATCH --output=run_install_phy_hyb_databases.%A.out
#SBATCH --error=run_install_phy_hyb_databases.%A.err

# Get the source ~/.bashrc file so that we can use conda.
source ~/.bashrc

# Database path.
database_dir="${HOME}/phytoplasma_hybrid_genes/databases"

# The reference protein sequence fasta file for building the bowtie2 reference database index.
ref_fasta_infile="${database_dir}/ref_db/ref_protein_seqs.fasta"

# Activate the bowtie2 conda environment.
conda activate bowtie2_env

# Build the bowtie2 reference database.
echo "bowtie2-build -f ${ref_fasta_infile} ${ref_fasta_infile}"
bowtie2-build -f ${ref_fasta_infile} ${ref_fasta_infile}

# The header line content from the reference protein sequence fasta file in a list file for retaining fastq reads aligned to the reference subgroup protein sequences.
subgroup_gene_list_file="${database_dir}/ref_db/ref_protein_seqs_header_list.txt"

# Get the contents of each fasta header in the reference fasta database and write to a file so we can extract mapped reads to a fastq file for assembly.
grep ">" ${ref_fasta_infile} | sed 's/>//g' > ${subgroup_gene_list_file}


# Activate the blast conda environment.
conda activate blast_env

for gene_name in $(cat $subgroup_gene_list_file | cut -d '_' -f2 | uniq);
do

	# The blast database fasta file.
	blast_db_fasta_file="${database_dir}/blast_dbs/${gene_name}/${gene_name}.fasta"

	# Make the blast databases for each gene.
	echo "makeblastdb -in ${blast_db_fasta_file} -dbtype nucl -out ${blast_db_fasta_file}"
	makeblastdb -in ${blast_db_fasta_file} -dbtype nucl -out ${blast_db_fasta_file}

done

