#!/usr/bin/bash

# Run multiGenomeLoader on all Gammaproteobacteria
# Download assembly files
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

# Unzip taxonomy names and tree files
tar xz taxonomyNames.pkl.xz
tar xz taxonomyParsedTree.pkl.xz

# Run Loader
python3 multigenomeLoader.py Gammaproteobacteria --download
