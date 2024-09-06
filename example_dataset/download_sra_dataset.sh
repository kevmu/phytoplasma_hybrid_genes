#!/bin/bash
#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=5G
#SBATCH --output=run_download_sra_subsampled_dataset.%A.out
#SBATCH --error=run_download_sra_subsampled_dataset.%A.err

# Get the conda environment path at the start of the shell script.
source ~/.bashrc

# Activate the NCBI SRA tools conda environment.
conda activate sra_tools_env

# The SRA run information summary file.
sra_run_info_infile="phyto_hyb_SraRunInfo.csv"

# The output directory to download SRA files.
fastq_output_dir="sra_fastq_files"

# The output directory.
output_dir="fixed_sra_fastq_files"

mkdir -p $output_dir

# The temporary file directory for fasterq-dump.
tmp_dir="$fastq_output_dir/tmp_dir"

# The number of threads to use for fasterq-dump.
num_threads=40

# Internal Field Separator (IFS) used when using cat in a for loop so that lines are separated by newlines instead of separated based on spaces.
IFS=$'\n'

# The fastq file list of sample names to use for running the script.
fastq_list_file="${output_dir}/fastq_files_list.txt"

# Iterate over the SraRunInfo.csv file 
for sra_run_info_entry in $(tail -n+2 $sra_run_info_infile); 
do 
	echo $sra_run_info_entry

	# Get the SRA Run ID from the sra_run_info_infile
	sra_run_id=$(echo $sra_run_info_entry | cut -d ',' -f1); 
	echo $sra_run_id

	# Get the sample_id from the sra_run_info_infile row.
	sample_id=$(echo $sra_run_info_entry | cut -d ',' -f30 | sed 's/ /_/g')
	echo $sample_id

	# Download the SRA fastq files for each sample_id using the sra_run_id SRA Run ID. Split them into read1 and read2 fastqs. Use num_threads for the number of threads for fasterq-dump. 
	echo "fasterq-dump --seq-defline '@$ac:$sn:$ri:$rl/$ri' --threads $num_threads --outdir ${fastq_output_dir} --temp $tmp_dir --split-files $sra_run_id"
	fasterq-dump --seq-defline '@$ac:$sn:$ri:$rl/$ri' --threads $num_threads --outdir ${fastq_output_dir} --temp $tmp_dir --split-files $sra_run_id

	# The sra run id fastq files.
	fastq_read1_infile="${fastq_output_dir}/${sra_run_id}_1.fastq"
	fastq_read2_infile="${fastq_output_dir}/${sra_run_id}_2.fastq"

	# Rename the fastq files to the sample name in the sra run file.
	new_fastq_read1_infile="${output_dir}/${sample_id}_R1.fastq"
	new_fastq_read2_infile="${output_dir}/${sample_id}_R2.fastq"

	# Remove the fastq comment line from the fastq files.
	sed 's/^+SRR.*/+/g' ${fastq_read1_infile} > ${new_fastq_read1_infile}
	sed 's/^+SRR.*/+/g' ${fastq_read2_infile} > ${new_fastq_read2_infile}

	# Write the sample names to the fastq list file.
	echo "${sample_id}" >> ${fastq_list_file}

done

echo "Use the following for the fastq_list_file and fastq_input_dir variables in phytoplasma_hybrid_analysis.sh"
echo "fastq_list_file=\"${fastq_list_file}\""
echo "fastq_input_dir=\"${output_dir}\""


