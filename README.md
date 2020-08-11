# ribosomal_snakemake

# Ribosomal protein tree generation

A workflow to generate ribosomal protein phylogenetic trees, using 15 ribosomal protein sequences.  Either protein or DNA sequences can be used to build the tree, it may be of benefit to use DNA sequences for more closely related genomes, and protein sequences for those that are more divergent.

# Installation via conda:
Use conda to install to a linux or mac osx environment:

# Requirements:
Snakemake>=3.5<br>
Python>=3.3<br>
Biopython>=1.77<br>
Diamond BLAST version>=0.9.32<br>
Mafft>=7.429<br>
Fasttree>=2.1.10 OR:<br>
Iqtree>=1.6.11<br>

Python scripts from this github folder

An environment.yaml is included for conda in env_yaml directory if required.

# Install conda:
Install conda here https://conda.io/en/latest/miniconda.html
Install snakemake:

conda install -c conda-forge mamba<br>
mamba create -c conda-forge -c bioconda -n snakemake snakemake<br>
<br>
conda activate snakemake<br>
snakemake --help<br>

# Install the ribosomal protein tree workflow from github:
git clone XXX<br>
cd XXX<br>
cd testfiles<br>
snakemake -n<br>

# Usage
Configure workflow:
Configure the workflow according to your needs by editing the file config.yaml. For use without any sequence database downloads, you need to provide the files and pathnames in cleannames.txt and atccs.txt.

1.	Gather your ribosomal protein sequences as protein sequence fasta files in 15 separate files.  These will be used for either protein or DNA sequence trees as per your choice laid out in the config.yaml.  They should be named within the file as L14_rplN, L16_rplP, L18_rplR, L2_rplB, L22_rplV, L24_rplX, L3_rplC, L4_rplD, L5_rplE, L6_rplF, S10_rpsJ, S17_rpsQ, S19_rpsS, S3_rpsC and S8_rpsH, respectively.  Staphylococcus sequences are included in the github folder.
2.	Gather all of your genome sequences as annotated genbank files in the same directory. 
3.	Create a file, “genus.txt” containing the name of the genus and species you want to build a tree for, one per line.  Note that the genus or species name should be enclosed with quotes.
1.	 If you do not want to download any genomes from the sequence databases, create a file, “cleannames.txt” containing the name of each genome file you want to include in your tree, one per line, and place “no” in the config.yaml under download_genbank options.  To use a mixture of download and provided genomes place “yes” in the config.yaml and do not provide a cleannames.txt file.
4.	To update any previous trees, collect any files generated as concatenated deduplicated fasta files that you want to update, and edit appropriately the config.yaml file with the filenames.  
Add Genbank files, ribosomal protein sequence files and a list of the files:
If you do not want to download any genomes, add the names of all provided genbank format files to cleannames.txt on a one-line per file basis.  The suffix of the genbank files should end in .gbff (following the genbank download nomenclature). Add the names of your ribosomal sequence files on a one-line per file basis to atccs.txt.  An easy way to generate these files is with: 
ls *gbff > cleannames.txt

# Execute workflow:
Test your configuration by performing a dry-run via
snakemake -n<br>
Execute the workflow locally via<br>
snakemake --cores $N<br>
using $N cores<br>

or run it in a cluster environment such as<br>
snakemake --use-conda --cluster qsub --jobs 100<br>

Further information on running snakemake in a cluster environment can be found on the snakemake website<br>
