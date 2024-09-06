
# Get the source ~/.bashrc file so that we can use conda.
source ~/.bashrc

# The concatenated phytoplasma hybridization gene fasta file.
fasta_file="/home/AGR.GC.CA/muirheadk/phytoplasma_hybrid_genes/output/all_concat_hybrid_genes.fasta"

# The output directory.
output_dir="/home/AGR.GC.CA/muirheadk/phytoplasma_hybrid_genes/tree"
mkdir -p $output_dir

# The number of threads to use.
num_threads=40

# The number of bootstrap replicates.
num_bootstraps=1000

# The model to use for generating the phylogenetic tree.
model_test="TEST"

# The fasta file basename.
fasta_basename=$(basename ${fasta_file})

# The fasta file filename.
fasta_filename=$(echo $fasta_basename | cut -d '.' -f1)

# The alignment file.
align_file="${output_dir}/${fasta_filename}.align.txt"

# Activate the mafft conda environment.
conda activate mafft_env

# The mafft command.
echo "mafft --maxiterate 1000 --localpair --thread ${num_threads} ${fasta_file} > ${align_file}"
mafft --maxiterate 1000 --localpair --thread ${num_threads} ${fasta_file} > ${align_file}

# Activate the iqtree conda environment.
conda activate iqtree_env

# The iqtree command.
echo "iqtree -s ${align_file} -m ${model_test} -mtree -B ${num_bootstraps}"
iqtree -s ${align_file} -m ${model_test} -mtree -B ${num_bootstraps}


