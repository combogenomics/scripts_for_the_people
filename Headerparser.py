#!/usr/bin/env python3

#Load libraries
import os
import argparse

#Making of Help menu
parser=argparse.ArgumentParser(
    description=
    ''' Headerparser.py is a script useful to modify headers in protein fasta
    files (.faa) in 3 ways. Using the --long, --short (DEFAULT) and --veryshort arguments the user is
    allowed to get the header in this exemple format | Long mode: >AOJ17836.1_Burkholderia_cenocepacia
    | Short mode: >AOJ17836.1_B_cenocepacia | Veryshort mode: >AOJ17836.1_B_cen          *************************************************************************''',
    epilog="""*************************************************************************""")
parser.add_argument('file', type=str, help='file to be parsed')
group = parser.add_mutually_exclusive_group(required=False)
group.add_argument('--long', '-L', help="Long header mode", action="store_true")
group.add_argument('--short', '-S', help="Short header mode", action="store_true")
group.add_argument('--veryshort', '-v', help="VeryShort header mode", action="store_true")
args = parser.parse_args()

#This function allows to convert header in Short Mode
def shortparser(filename):
    #Checking if the file is a .faa file
    faa = os.path.splitext(filename)
    if faa[1] == ".faa":
        file =  faa[0] + "_PARSED_short.faa"
        #Copy old file to new file, checking where the header (">") is and modifying it
        with open(file, "wt") as newfile:
            with open(filename, "rt") as oldfile:
                for line in oldfile.readlines():
                    if ">" in line:
                        list = line.split()
                        genera = list[-2]
                        onlyfirst = ''.join([c for c in genera if c.isupper()])
                        str2= list[-1]
                        str3=list[0]
                        str = str3+"_"+onlyfirst+"_"+str2
                        finalname = str.replace("]", "\n")
                        newfile.write(line.replace(line, finalname))
                    if ">" not in line:
                        newfile.write(line)
    #If not a .faa file, print the message
    else:
        print("ERROR: Only Protein fasta files (.faa) are allowed")

#This function allows to convert header in Long Mode
def longparser(filename):
    #Checking if the file is a .faa file
    faa = os.path.splitext(filename)
    if faa[1] == ".faa":
        file =  faa[0] + "_PARSED_long.faa"
        #Copy old file to new file, checking where the header (">") is and modifying it
        with open(file, "wt") as newfile:
            with open(filename, "rt") as oldfile:
                for line in oldfile.readlines():
                    if ">" in line:
                        list = line.split()
                        str1= list[-2]
                        str2= list[-1]
                        str3=list[0]
                        str = str3+"_"+str1+"_"+str2
                        str = str.replace("[", "")
                        finalname = str.replace("]", "\n")
                        newfile.write(line.replace(line, finalname))
                    if ">" not in line:
                        newfile.write(line)
    #If not a .faa file, print the message
    else:
        print("ERROR: Only Protein fasta files (.faa) are allowed")

#This function allows to convert header in Veryshort Mode
def vshortparser(filename):
    faa = os.path.splitext(filename)
    if faa[1] == ".faa":
        file =  faa[0] + "_PARSED_vshort.faa"
        #Copy old file to new file, checking where the header (">") is and modifying it
        with open(file, "wt") as newfile:
            with open(filename, "rt") as oldfile:
                for line in oldfile.readlines():
                    if ">" in line:
                        #Modifying the header through iteration
                        list = line.split()
                        str1= list[-1]
                        str = ""
                        for n in range(3):
                            str += str1[n]
                        genera = list[-2]
                        str2 = ''.join([c for c in genera if c.isupper()])
                        str3=list[0]
                        finalstr = str3+"_"+str2+"_"+str+"\n"
                        newfile.write(line.replace(line, finalstr))
                    if ">" not in line:
                        newfile.write(line)

    #If not a .faa file, print the message
    else:
        print("ERROR: Only Protein fasta files (.faa) are allowed")
#Run the tool using mode short, long or vshort
if args.long:
    longparser(args.file)
if args.veryshort:
    vshortparser(args.file)
else:
    shortparser(args.file)
