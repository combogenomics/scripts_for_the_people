#!/usr/bin/env python3

#Import libraries
import os
import argparse
#Making of Help menu
from argparse import _MutuallyExclusiveGroup

parser = argparse.ArgumentParser(
    description=
    ''' This script is useful for the extraction of Intergenic regions (or custumized intergenic  spots) from a
    genome, giving a fasta (.fna) file and file reporting gene and region features (.gff, .gff3) at the command line.
    The script finds the starting and ending point of genes (coding regions) and allows the extraction of the sequences
    nearby the starting point. Default is --intergenic mode, that will allow to output the sequences between two
    adjacent genes. If --custom mode is selected, the region that will be considered will include the first 50 bp from
    the starting point of each gene and an arbitrary lenght specified at the command line (max lenght 400 bp). The
    output will be a fasta file, and the header of each sequence will contain infos like the ChromosomeID, gene locus
    tag, strand and, if --custom mode is indicated, if the starting point of the region selected is internal or external
    to a nearby gene ( useful to find internal promoters ). ''', epilog = """*********************************""")
parser.add_argument('fasta', type=str, help='a genome file')
parser.add_argument('gff', type=str, help='a gff gene feature file')
group: _MutuallyExclusiveGroup = parser.add_mutually_exclusive_group(required=False)
group.add_argument('--custom', type=int, help="Custom mode: choose the region size", choices=[100,200,300,400])
group.add_argument('--intergenic', '-i', help="Intergenic mode [DEFAULT]", action="store_true")
args = parser.parse_args()




########################## FUNCTION SECTION ###############################
#The first Function Extracts the desired fields from each gff line, in which
#n is the field number ("\t" is the separator)

def locus_extract(line, n):
    field = line.split("\t")
    locus = field[n]
    if n > 8 or n < 0:
        print("Locus out of range")
    if n == 8:
        line_list = locus.split(";")
        locus_ID = line_list[-1]
        ID = ""
        #Here we get the locus ID from the last column (8)
        for n in range(10,len(locus_ID)-1):
            ID += locus_ID[n]
        return ID
    else:
        return locus

#This function takes a fasta file converted in a list (with readlines() method)
#and  returns a list with n arguments reporting each sequence of each chromosome
#(n=number of chromosomes): this allow to work with different chromosomes
#separately as different arguments of the same list.
def define_chromosomes(list):
    outer_list = []
    inner_list = []
    n = 0
    for elem in list:
        inner_list.append(elem)
        n += 1
        if ">" in elem:
            if n == 1:
                inner_list.remove(inner_list[0])
            if n > 1:
                inner_list.remove(inner_list[-1])
                outer_list.append(inner_list)
                inner_list = []
    outer_list.append(inner_list)
    global allchr
    allchr = []
    for elem in outer_list:
        joined = ''.join(elem)
        nonewline = joined.replace("\n", "")
        allchr.append(nonewline)
    return allchr

#---------------------------------CUSTOM MODE--------------------------------#
#This function is needed to iterate through each chromosome and each gene start,
#depending on the strand, and finding the intergenic region starting point
#and to determinate if this starting point is located internally or externally
#to the nearest gene, depending on the strand. Then the intergenic sequence
#and occasionally a part of nearby genic sequence (if int status) is extracted.
def intergenic_custom(strand, basepairs):
    global status
    global seq
    if strand == "+":
        gene_start = (int(start) - 1) + 50
        sequence_start = gene_start - basepairs
        if sequence_start < 0:
            newstart = len(allchr[a]) + sequence_start
            near = locus_extract(outer[a][-1], 4)
            seq = allchr[a][newstart:]+ allchr[a][:gene_start]+"\n"
            if len(allchr[a][:newstart]) <= int(near)-1:
                status = "int"
            else:
                status = "ext"
        else:
            near = locus_extract(outer[a][num-1], 4)
            seq = allchr[a][sequence_start:gene_start]+"\n"
            if sequence_start <= int(near)-1:
                status = "int"
            else:
                status = "ext"
    if strand == "-":
        gene_start = (int(start) - 1) - 50
        sequence_start = gene_start + basepairs
        if sequence_start > len(allchr[a]):
            newstart = sequence_start - len(allchr[a])
            near = locus_extract(outer[a][0], 3)
            seq = allchr[a][gene_start:]+allchr[a][:newstart]+"\n"
            if len(allchr[a][:newstart]) >= int(near)-1:
                status = "int"
            else:
                status = "ext"
        else:
            seq = allchr[a][gene_start:sequence_start]+"\n"
            if outer[a][num] != outer[a][-1]:
                near = locus_extract(outer[a][num+1], 3)
                if sequence_start >= int(near)-1:
                    status = "int"
                else:
                    status = "ext"
            else:
                status = "ext"
    return status,seq
#This function creates a fasta output in which are reported all the genes
#associated nearby sequences, according with the option specified in --custom
#argument. The header of each sequence includes the Chromosome,the locus tag,
#the strand and the internal or external status of the sequence.
def output_custom(chr, tag, strand, status, seq):
    header = [chr, tag, strand, status]
    header = "_".join(header)
    header = ">"+header+"\n"
    total_file = header+seq
    newfile.write(total_file)

#------------------------------INTERGENIC MODE-------------------------------#
#This funcion is selected by DEFAULT, in the --intergenic mode, which recalls
#the chromosome/region name, the gene locus tag and the strand of each feature,
#and appends it to the headers an new output file.
def output_intergenic(chr,tag,strand):
    header = [chr, tag, strand]
    header = "_".join(header)
    header = ">"+header+"\n"
    newfile.write(header)
#This function extracts the intergenic sequences of each gene, of each
#chromosome/region with the DEFAULT --intergenic mode.
def second_line_intergenic(strand, start):
    gene_start = (int(start)-1)
    #For the forward strand
    if strand == "+":
        if outer[a][num] == outer[a][0]:
            near = locus_extract(outer[a][-1], 4)
            seq = allchr[a][int(near)-1:]+allchr[a][:gene_start]+"\n"
        else:
            near = locus_extract(outer[a][num-1], 4)
            seq = allchr[a][int(near)-1:gene_start]+"\n"
    else:
        if outer[a][num] == outer[a][-1]:
            near = locus_extract(outer[a][0], 3)
            seq = allchr[a][gene_start:]+allchr[a][:int(near)-1]+"\n"
        else:
            near = locus_extract(outer[a][num+1], 3)
            seq = allchr[a][gene_start:int(near)-1]+"\n"
    newfile.write(seq)
    newfile.close()
    return



###############################################################################
###############################################################################
#Reading of gff file. If not ".gff" operation fails
GFF = args.gff
print(GFF)
print(os.path.splitext(GFF)[1])
if ".gff" in os.path.splitext(GFF)[1]:
    GFF = open(GFF, "rt")
    list_lines = GFF.readlines()
    start_position = []
    int_intergenic = []
    outer = []
    #Extracting fields: the region for the entire chromosome, the gene for each
    #coding sequence, and generating a list of lists: each region is a list,
    #and its genes(features) represent the arguments of each list.
    for line in list_lines:
        if "#" in line:
            pass
        else:
            Field = locus_extract(line,2)
            if "region" in Field:
                outer.append(int_intergenic)
                int_intergenic = []
            if "gene" in Field:
                int_intergenic.append(line)
    GFF.close()
else:
    print("ERROR: Only gff (.gff) files are allowed")
outer.append(int_intergenic)
outer = outer[1:]
num = 0
a = 0
#Reading input fasta file
fasta = args.fasta
name = os.path.splitext(fasta)[0]
#Define output file name
if args.custom:
    newname = name+"_Intergenic"+str(args.custom)+".fa"
else:
    newname = name+"_Intergenic.fa"
first = True
#Checking if file is a fasta (.fna)
if  os.path.splitext(fasta)[1] == ".fna":
    for list in outer:
        #For each chromosome, extract the needed locus (chr, strand, start, end)
        for line in outer[a]:
            chr = locus_extract(line, 0)
            print("The current chromosome is:", chr)
            locus = locus_extract(line, 8)
            print("The current locus is:", locus)
            strand = locus_extract(line, 6)
            print("The current strand is:", strand)
            if strand == "+":
                start = locus_extract(line, 3)
                print("The gene start is:", start)
            else:
                start = locus_extract(line, 4)
                print("The gene start is:", start)
            #Check if the file to write already exists
            if first:
                with open(newname, "wt") as newfile:
                    #Open the fna file given at the command line as an argument
                    #if output already exists, open it in append mode.
                    with open(fasta, "rt") as oldfile:
                        line_list = oldfile.readlines()
                        #Divide chromosomes with the function define_chromosomes
                        define_chromosomes(line_list)
                        #if --custom option is indicated at the command line:
                        if args.custom:
                            intergenic_custom(strand, args.custom)
                            output_custom(chr, locus, strand, status, seq)
                        #if not --custom, --intergenic is the default
                        else:
                            output_intergenic(chr, locus, strand)
                            second_line_intergenic(strand, start)
            else:
                with open(newname, "a") as newfile:
                    if args.custom:
                        intergenic_custom(strand, args.custom)
                        output_custom(chr, locus, strand, status, seq)
                    else:
                        output_intergenic(chr, locus, strand)
                        second_line_intergenic(strand, start)
            first = False
            num += 1
        num = 0
        a += 1
else:
    print("ERROR: Only genome fasta files (.fna) are allowed")
