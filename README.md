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
