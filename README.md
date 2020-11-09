# scripts_for_the_people
Commonly used scripts for basic computational biology tasks

# Scripts description
1. Headerparser.py

    Headerparser.py is a script useful to modify headers in protein fasta files
    (.faa) in 3 ways. Using the --long, --short (DEFAULT) and --veryshort
    arguments the user is allowed to get the header in the 3 forms. 

2. Intergenic_finder.py

    This script is useful for the extraction of Intergenic regions (or custumized intergenic  spots) from a
    genome, giving a fasta (.fna) file and file reporting gene and region features (.gff, .gff3) at the command line.
    The script finds the starting and ending point of genes (coding regions) and allows the extraction of the sequences
    nearby the starting point. Default is --intergenic mode, that will allow to output the sequences between two
    adjacent genes. If --custom mode is selected, the region that will be considered will include the first 50 bp from
    the starting point of each gene and an arbitrary lenght specified at the command line (max lenght 400 bp). The
    output will be a fasta file, and the header of each sequence will contain infos like the ChromosomeID, gene locus
    tag, strand and, if --custom mode is indicated, if the starting point of the region selected is internal or external
    to a nearby gene ( useful to find internal promoters ).

3. multigenomeLoader.py

This script allows to massively download the Genome, Proteome, and Features of all the organisms belonging to the same taxonomic unit.
The user only have to specify at the command line the taxonomic unit (e.g. Burkholderiaceae) and the script will automatically download
all the related species and subspecies data available (that at least have data for "Complete Genomes" or "Chromosome") from NCBI database (using Refseq accessions). The files will be found into Genomes, Proteomes, Features folders, in the same directory. The program only need in the same directory the parsed taxonomic tree files taxonomyNames.pkl and taxonomyParsedTree.pkl, available xzipped (xz -d to decompress) in the same subfolder, and the assembly_summary_refseq.txt file, to be downloaded at ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/ .
The script has the below options:
The flag --num allows to choose the number of organisms the user want to get/download; the organisms data will be downloaded randomly from the total pool of the available ones
The flag --getstrain allows to get the strain name for each downloaded organism data: it will be shown in the organism file names.
The flag --exclude_sp allows to discard these organisms when species level is not identified (sp.)
The flag --download is necessary to download all the genomes of the organisms. If this flag is not specified, the script will only show
the number of genomes available for each taxonomic units, and the related ftps ncbi link for the download of the data.

Example on how to run the script:
python3 multigenomeLoader.py --num 10 --exclude_sp --getstrain Escherichia

A future release could:
- Include the possibility to update the Parsed tree files, generating new ones downloading new names.dmp and nodes.dmp from ncbi ftp taxdump folder (needs to be done every once in a while)
- automatically download assembly_summary_refseq.txt.
- Increase download speed.

