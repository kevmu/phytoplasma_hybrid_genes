### Instructions on how to install the conda environments for the shell script and python script.

## Install EMBOSS conda environment to get getorf command.

# Create the EMBOSS conda environment.
conda create --name emboss_env

# Activate the EMBOSS conda environment
conda activate emboss_env

# Install the EMBOSS conda environment.
conda install -c bioconda emboss

## Install biopython conda environment to get getorf command.

# Create the biopython conda environment.
conda create --name biopython_env

# Activate the biopython conda environment
conda activate biopython_env

# Install the biopython conda environment.
conda install -c bioconda biopython


### How to use the extract_hybrid_genes.sh script.

Use a text editor to edit the input and output directories in the extract_hybrid_genes.sh file.

### Run using the following command.
sh extract_hybrid_genes.sh


