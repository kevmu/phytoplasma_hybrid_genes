#!/bin/bash
#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G
#SBATCH --output=run_phyto_hyb_analysis.%A.out
#SBATCH --error=run_phyto_hyb_analysis.%A.err

#Phytoplasma multi-gene hyb panel
#new algorithm for mapping and assembly 220122
#step 1: trimmomatic/flash2
#step 2: mapping vs protein coding genes (all types)
#step 3: assemble protein coding genes and determine phytoplasma type - eg cpnclassiphyr
#step 4: trim all assemblies to retain only contigs >=500 bp
#final product: set of fasta files containing long assemblies for each of 6 protein-coding genes/loci plus 16S

# Get the source ~/.bashrc file so that we can use conda.
source ~/.bashrc

# The list of fastq filenames.
fastq_list_file="/home/AGR.GC.CA/muirheadk/phytoplasma_hybrid_genes/fastq_files_list.txt"

# The fastq input file directory.
fastq_input_dir="/archive/dumonceauxt/230801_M01666_0203_000000000-KTT65_IMPACTTlilacBBSPTWMQhybsEPL/Fastq"

# The read1 fastq suffix.
#read1_suffix="_R1.fastq"
read1_suffix="_L001_R1_001.fastq.gz"

# The read2 fastq suffix.
#read2_suffix="_R2.fastq"
read2_suffix="_L001_R2_001.fastq.gz"

# The fastq file extension suffix.
#fastq_file_ext=".fastq.gz"

# The number of cpu threads to use.
num_threads=40

# The base output directory.
output_dir="/home/AGR.GC.CA/muirheadk/phytoplasma_hybrid_genes/output"
mkdir -p $output_dir

# The preprocessing output directory.
preprocessing_output_dir="${output_dir}/preprocessing"
mkdir -p $preprocessing_output_dir

#fastq_list_file="${preprocessing_output_dir}/fastq_files_list.txt"

### Trimmomatic program parameters.

# Cut bases off the start of a read, if below a threshold quality < N.
cut_leading_bases=3

# Cut bases off the end of a read, if below a threshold quality < N.
cut_trailing_bases=20

## Use a sliding window of size sliding_window_size that will remove bases if their phred score is below sliding_window_phred_score.

# Sliding window of size N.
sliding_window_size=4

# Remove bases if their phred score is below 20
sliding_window_phred_score=15

# Minimum length for fastq reads.
min_read_length=36


### Flash2 program parameters.

## Flash2 command configuration data.

# The fragment length.
fragment_length=550

# The fragment length standard deviation.
fragment_length_stddev=200

# The read length.
read_length=300

### Databases.

# Database path.
database_dir="/home/AGR.GC.CA/muirheadk/phytoplasma_hybrid_genes/databases"
### Bowtie2 database.

# The reference protein sequence fasta file for building the bowtie2 reference database index.
ref_fasta_infile="${database_dir}/ref_db/ref_protein_seqs.fasta"

# The header line content from the reference protein sequence fasta file in a list file for retaining fastq reads aligned to the reference subgroup protein sequences.
subgroup_gene_list_file="${database_dir}/ref_db/ref_protein_seqs_header_list.txt"

# Get the contents of each fasta header in the reference fasta database and write to a file so we can extract mapped reads to a fastq file for assembly.
grep ">" ${ref_fasta_infile} | sed 's/>//g' > ${subgroup_gene_list_file}

### Transabyss program parameters.

# The transabyss kmer size.
transabyss_kmers=32

# Length Filtering program parameters.
min_seq_length=500

### EMBOSS getorf program parameters.

# The minimum sequence length for the ORFs.
min_size=500

# The maximum sequence length for the ORFs.
max_size=4000

# Find the top three longest ORFs.
orf_find=3

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


# Find all the fastq files in the dataset and write the full path to the file.
#echo "find ${fastq_input_dir} -name \"*${fastq_file_ext}\" -type f | sed \"s/${read1_suffix}\|${read2_suffix}//g\" | rev | cut -d '/' -f1 | rev | sort -V | uniq > ${fastq_list_file}"
#find ${fastq_input_dir} -name "*${fastq_file_ext}" -type f | sed "s/${read1_suffix}\|${read2_suffix}//g" | rev | cut -d '/' -f1 | rev | sort -V | uniq > ${fastq_list_file}


#step 1A: trimmomatic

# Activate the trimmomatic conda environment.
conda activate trimmomatic_env

# Make the trimmomatic trim output directory.
trim_fastq_dir="${preprocessing_output_dir}/trim"
mkdir -p $trim_fastq_dir

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";
    
    fastq_read1_file="${fastq_input_dir}/${fastq_filename}${read1_suffix}";
    fastq_read2_file="${fastq_input_dir}/${fastq_filename}${read2_suffix}";

    trim_paired_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.fq";
    trim_unpaired_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.fq";
    trim_paired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.reverse.fq";
    trim_unpaired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.reverse.fq";
    
    echo "trimmomatic PE ${fastq_read1_file} ${fastq_read2_file} ${trim_paired_fastq_file} ${trim_unpaired_fastq_file} ${trim_paired_rev_fastq_file} ${trim_unpaired_rev_fastq_file} LEADING:${cut_leading_bases} TRAILING:${cut_trailing_bases} SLIDINGWINDOW:${sliding_window_size}:${sliding_window_phred_score} MINLEN:${min_read_length}";    
    trimmomatic PE ${fastq_read1_file} ${fastq_read2_file} ${trim_paired_fastq_file} ${trim_unpaired_fastq_file} ${trim_paired_rev_fastq_file} ${trim_unpaired_rev_fastq_file} LEADING:${cut_leading_bases} TRAILING:${cut_trailing_bases} SLIDINGWINDOW:${sliding_window_size}:${sliding_window_phred_score} MINLEN:${min_read_length};

done


### Start merging trimmed fastq files using flash2.

# Activate the flash2 conda environment.
conda activate flash2_env

# Make the flash merge output directory.
flash2_merge_dir="${preprocessing_output_dir}/merge"
mkdir -p $flash2_merge_dir

for fastq_filename in $(cat $fastq_list_file);
do

    trim_paired_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.fq";
    trim_paired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.reverse.fq";
    merged_fastq_outfile="${fastq_filename}.merged";
    
    echo "flash2 -f ${fragment_length} -s ${fragment_length_stddev} -r ${read_length} ${trim_paired_fastq_file} ${trim_paired_rev_fastq_file} -d ${flash2_merge_dir} -o ${merged_fastq_outfile}";
    flash2 -f ${fragment_length} -s ${fragment_length_stddev} -r ${read_length} ${trim_paired_fastq_file} ${trim_paired_rev_fastq_file} -d ${flash2_merge_dir} -o ${merged_fastq_outfile};

done


##step 3: mapping vs protein-coding genes
##Install bowtie2
#
#conda install -c conda-forge tbb
#
#conda install -c bioconda bowtie2
#
##make bowtie reference files
##Ensure you are using a compute node (qlogin)
#
##build one ref file for all protein-coding seqs (cpn60; nusA, secY, secA, rp, tuf - call this file ref_protein_seqs.fasta)
##then, build separate ref files containing type-specific 16S genes with junk; eg typeI_16S_junk.fasta
#

conda activate bowtie2_env

# Build the bowtie2 reference database.
echo "bowtie2-build -f ${ref_fasta_infile} ${ref_fasta_infile}"
bowtie2-build -f ${ref_fasta_infile} ${ref_fasta_infile}

# The mapped reference output directory.
mapped_reference_output_dir="${output_dir}/mapped_reference"
mkdir -p $mapped_reference_output_dir

# Map protein-coding genes.
for fastq_filename in $(cat $fastq_list_file);
do

    merged_fastq_outfile="${flash2_merge_dir}/${fastq_filename}.merged.extendedFrags.fastq";
    sam_file="${mapped_reference_output_dir}/${fastq_filename}.ref_genes.sam";
    
    echo "bowtie2 --local -p ${num_threads} -x ${ref_fasta_infile} -U ${merged_fastq_outfile} -S ${sam_file}";
    bowtie2 --local -p ${num_threads} -x ${ref_fasta_infile} -U ${merged_fastq_outfile} -S ${sam_file};

done

conda activate samtools_env

for fastq_filename in $(cat $fastq_list_file);
do
    sam_file="${mapped_reference_output_dir}/${fastq_filename}.ref_genes.sam";
    bam_file="${mapped_reference_output_dir}/${fastq_filename}.ref_genes.bam";
    sorted_bam_file="${mapped_reference_output_dir}/${fastq_filename}.ref_genes.bam.sorted";

    echo "samtools view -bS --threads ${num_threads} ${sam_file} -o ${bam_file}"
    samtools view -bS --threads ${num_threads} ${sam_file} -o ${bam_file};
    
    echo "samtools sort ${bam_file} -o ${sorted_bam_file}"
    samtools sort ${bam_file} -o ${sorted_bam_file};

    echo "samtools index ${sorted_bam_file}"
    samtools index ${sorted_bam_file};

done


## The mapped bam file output directory.
mapped_bam_output_dir="${output_dir}/mapped_bam"
mkdir -p $mapped_bam_output_dir

## Extract mapped reads and corresponding fastq for each protein-coding gene.
for fastq_filename in $(cat $fastq_list_file);
do
    sorted_bam_file="${mapped_reference_output_dir}/${fastq_filename}.ref_genes.bam.sorted";
    
    for subgroup_gene_name in $(cat $subgroup_gene_list_file);
    do
        
        # The mapped bam file for each gene to write all the mapped reads in bam format.
        mapped_bam_file="${mapped_bam_output_dir}/${fastq_filename}.mapped.${subgroup_gene_name}.bam";
        
        # Obtain all the reads aligned to each gene and write to the mapped bam file.
        echo "samtools view -bh --threads ${num_threads} ${sorted_bam_file} ${subgroup_gene_name} -o ${mapped_bam_file}";
        samtools view -bh --threads ${num_threads} ${sorted_bam_file} ${subgroup_gene_name} -o ${mapped_bam_file};

    done
done


## The mapped bam output directory.
mapped_fastq_output_dir="${output_dir}/mapped_fastq"
mkdir -p $mapped_fastq_output_dir

for fastq_filename in $(cat $fastq_list_file);
do

    for subgroup_gene_name in $(cat $subgroup_gene_list_file);
    do
        mapped_bam_file="${mapped_bam_output_dir}/${fastq_filename}.mapped.${subgroup_gene_name}.bam";
        mapped_subgroup_gene_fastq_file="${mapped_fastq_output_dir}/${fastq_filename}.mapped.${subgroup_gene_name}.fastq";
	singleton_fastq_file="${mapped_fastq_output_dir}/${fastq_filename}_singletons.fastq" 
        
	echo "samtools fastq ${mapped_bam_file} -s ${singleton_fastq_file} --threads ${num_threads} > ${mapped_subgroup_gene_fastq_file}";
        samtools fastq ${mapped_bam_file} -s ${singleton_fastq_file} --threads ${num_threads} > ${mapped_subgroup_gene_fastq_file};

        # Get the gene name so that we can append each subgroup gene into the same file.
        gene_name=$(echo $subgroup_gene_name | cut -d '_' -f2);

        # The gene fastq file.
        mapped_gene_fastq_file="${mapped_fastq_output_dir}/${fastq_filename}.mapped.${gene_name}.fastq";
        
        echo "cat ${mapped_subgroup_gene_fastq_file} >> ${mapped_gene_fastq_file}";
        cat ${mapped_subgroup_gene_fastq_file} >> ${mapped_gene_fastq_file};
        
    done
    
done

## Calculate the number of fastq reads mapped to each gene for each sample.
# The mapped fastq read count summary per gene.
num_fastq_reads_file="${output_dir}/read_count_summary.tsv"

echo "echo -e \"sample_name\tgene_name\tnum_fastq_reads\" > ${num_fastq_reads_file}"
echo -e "sample_name\tgene_name\tnum_fastq_reads" > ${num_fastq_reads_file}

# Calculate the number of fastq reads and report to the summary file.
for fastq_filename in $(cat $fastq_list_file);
do

    for gene_name in $(cat $subgroup_gene_list_file | cut -d '_' -f2 | uniq);
    do
        # The gene fastq file.
        mapped_gene_fastq_file="${mapped_fastq_output_dir}/${fastq_filename}.mapped.${gene_name}.fastq";
        
        # Calculate the number of fastq reads.
	num_fastq_reads=$(echo $(cat ${mapped_gene_fastq_file} | wc -l) / 4 | bc)

        echo "echo -e \"${fastq_filename}\t${gene_name}\t${num_fastq_reads}\" >> ${num_fastq_reads_file}"
        echo -e "${fastq_filename}\t${gene_name}\t${num_fastq_reads}" >> ${num_fastq_reads_file}
        
    done
done

##Assembly

## The assembly output directory.
assembly_output_dir="${output_dir}/gene_assemblies"
mkdir -p $assemblies_output_dir

# Activate the transabyss conda environment.
conda activate transabyss_env

# Run the assembly step on each mapped gene fastq file.
for fastq_filename in $(cat $fastq_list_file);
do

    for gene_name in $(cat $subgroup_gene_list_file | cut -d '_' -f2 | uniq);
    do
        # The gene fastq file.
        mapped_gene_fastq_file="${mapped_fastq_output_dir}/${fastq_filename}.mapped.${gene_name}.fastq";
        gene_assembly_output_dir="${assembly_output_dir}/${fastq_filename}/${gene_name}"
        
        echo "transabyss --se ${mapped_gene_fastq_file} -k ${transabyss_kmers} --threads ${num_threads} --outdir ${gene_assembly_output_dir}"
        transabyss --se ${mapped_gene_fastq_file} -k ${transabyss_kmers} --threads ${num_threads} --outdir ${gene_assembly_output_dir}

    done
done

##step6 - trim assemblies to retain only those >500 bp 

# Activate the biopython conda environment.
conda activate biopython_env

for fastq_filename in $(cat $fastq_list_file);
do

    for gene_name in $(cat $subgroup_gene_list_file | cut -d '_' -f2 | uniq);
    do
        
	gene_assembly_output_dir="${assembly_output_dir}/${fastq_filename}/${gene_name}"
        
	assembly_file="${gene_assembly_output_dir}/transabyss-final.fa"

        filtered_assembly_file="${gene_assembly_output_dir}/${fastq_filename}__${gene_name}_assembly.fasta"
        
        echo "python filter_sequences_by_length.py -i ${assembly_file} -l ${min_seq_length} -o ${filtered_assembly_file}"
        python filter_sequences_by_length.py -i ${assembly_file} -l ${min_seq_length} -o ${filtered_assembly_file}
    done
done

### BLAST

for fastq_filename in $(cat $fastq_list_file);
do

    for gene_name in $(cat $subgroup_gene_list_file | cut -d '_' -f2 | uniq);
    do

        gene_assembly_output_dir="${assembly_output_dir}/${fastq_filename}/${gene_name}"

        filtered_assembly_file="${gene_assembly_output_dir}/${fastq_filename}__${gene_name}_assembly.fasta"
	
	blast_db_fasta_file="${database_dir}/blast_dbs/${gene_name}/${gene_name}.fasta"

	blast_results_file="${gene_assembly_output_dir}/${fastq_filename}__${gene_name}_blastn_results.tsv"

	# Activate the blast conda environment.
	conda activate blast_env

	echo "makeblastdb -in ${blast_db_fasta_file} -dbtype nucl -out ${blast_db_fasta_file}"
	makeblastdb -in ${blast_db_fasta_file} -dbtype nucl -out ${blast_db_fasta_file}

	echo "echo \"qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs qcovhsp nident positive mismatch gaps qframe sframe\" | sed 's/ /\t/g' > ${blast_results_file}"
	echo "qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs qcovhsp nident positive mismatch gaps qframe sframe" | sed 's/ /\t/g' > ${blast_results_file}

	echo "blastn -query ${filtered_assembly_file} -db ${blast_db_fasta_file} -out "-" -evalue 1e-05 -max_target_seqs 10 -num_threads ${num_threads} -outfmt '6 qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs qcovhsp nident positive mismatch gaps qframe sframe' >> ${blast_results_file}"
	blastn -query ${filtered_assembly_file} -db ${blast_db_fasta_file} -out "-" -evalue 1e-05 -max_target_seqs 10 -num_threads ${num_threads} -outfmt '6 qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs qcovhsp nident positive mismatch gaps qframe sframe' >> ${blast_results_file}

	# Activate the biopython conda environment.
	conda activate biopython_env 
	
	echo "python parse_best_hit_fasta.py --fasta_infile ${filtered_assembly_file} --blast_results_infile ${blast_results_file} --output_dir ${gene_assembly_output_dir}"
	python parse_best_hit_fasta.py --fasta_infile ${filtered_assembly_file} --blast_results_infile ${blast_results_file} --output_dir ${gene_assembly_output_dir}

    done
done

### EMBOSS getorf.

gene_seqs_output_dir="${output_dir}/gene_sequences"
mkdir -p $gene_seqs_output_dir

for fastq_filename in $(cat $fastq_list_file);
do

    for gene_name in $(cat $subgroup_gene_list_file | cut -d '_' -f2 | uniq);
    do

	# The gene assembly output directory.
        gene_assembly_output_dir="${assembly_output_dir}/${fastq_filename}/${gene_name}"
	
	# The gene assembly best hit fasta file with contig with best match to the database.
	gene_assembly_fasta_file="${gene_assembly_output_dir}/${fastq_filename}__${gene_name}_assembly_best_hit.fasta"	

	# The gene output directory.
	gene_output_dir="${gene_seqs_output_dir}/${fastq_filename}"
	mkdir -p $gene_output_dir

	# The gene fasta output file that has the gene in correct orientation and ORF.
	gene_fasta_outfile="${gene_output_dir}/${fastq_filename}__${gene_name}.fasta"

	# Activate the emboss conda environment.
	conda activate emboss_env	

        # Use the getorf command to get ORFs from each assembly.
	echo -e "getorf -sequence ${gene_assembly_fasta_file} -minsize ${min_size} -maxsize ${max_size} -find ${orf_find} -methionine N -outseq ${gene_fasta_outfile}"
        getorf -sequence ${gene_assembly_fasta_file} -minsize ${min_size} -maxsize ${max_size} -find ${orf_find} -methionine N -outseq ${gene_fasta_outfile}

    done
done

### Concatenate genes into one fasta file per sample.

# Get a concatenated ordered list of gene names separated by commas.
gene_name_list=$(join_by_string "," "${gene_list[@]}")

for fastq_filename in $(cat $fastq_list_file);
do
        # The gene output directory.
        gene_output_dir="${gene_seqs_output_dir}/${fastq_filename}"
        mkdir -p $gene_output_dir

	# Concatenate the genes into one sequence so that we can generate a phylogenetic tree.
	gene_fasta_list_file="${gene_output_dir}/${fastq_filename}_hybrid_genes_files.txt"

	# The list of hybrid gene file paths.
	find $gene_output_dir -name "*\.fasta" -type f | grep "${gene_regex}" > $gene_fasta_list_file

	# Activate the biopython conda environment.
	conda activate biopython_env

	echo "Concatinating ${gene_name_list} files..."

	# Concatenate hybrid gene panel sequences in order of the gene_name_list.
	echo -e "python concat_seqs_order.py --fasta_file_list_infile ${gene_fasta_list_file} --gene_name_list ${gene_name_list} --sample_name ${fastq_filename} --output_dir ${gene_output_dir}"
	python concat_seqs_order.py --fasta_file_list_infile ${gene_fasta_list_file} --gene_name_list ${gene_name_list} --sample_name ${fastq_filename} --output_dir ${gene_output_dir}

done

