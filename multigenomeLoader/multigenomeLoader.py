#!/usr/bin/env python3.6
# encoding: utf-8

import argparse
import pickle
import random
import subprocess
import wget
import sys

parser = argparse.ArgumentParser(
    description=''' This script is useful to write down a new table with a custom number of rows and columns and
    a boxplot with the stats of the experiment
''', epilog="""*********************************""")
parser.add_argument('taxa', type=str, help="Taxonomic unit you want to get all genomes")
parser.add_argument('--num', type=int, help="Choose the number of genomes you want to download. A representative,"
                                            "randomic pool of #num genomes will be downloaded")
parser.add_argument('--getstrain', help="From given species, get all strains on database", action="store_true")
parser.add_argument('--exclude_sp', help="Exclude organism if species level is not identified (sp.)", action="store_true")
parser.add_argument('--download', help="Download the selected genomes, features and proteomes", action="store_true")
parser.add_argument('--quiet', help="File downloading is not showed on shell")
args = parser.parse_args()

""" Parse the NCBI taxonomy """

import logging
import random
from typing import Dict, List, Tuple

LOGGER = logging.getLogger()
LOGGER.setLevel(0)


class Node:
    """ Definition of the class Node """

    def __init__(self):
        self.tax_id = "0"  # Number of the tax id.
        self.parent = "0"  # Number of the parent of this node
        self.children = []  # List of the children of this node
        self.division = None  # Division.
        self.is_tip = True  # Tip = True if it's a terminal node, False if not.
        self.name = ""  # Name of the node: taxa if it's a terminal node, number if not.


def get_genealogy(name_object, leaf_node: str) -> List[str]:
    """ Trace genealogy from root to leaf """
    ancestors = []  # Initialise the list of all nodes from root to leaf.
    gen_tax_id = leaf_node  # Define leaf
    while 1:
        if gen_tax_id in name_object:
            ancestors.append(gen_tax_id)
            gen_tax_id = name_object[gen_tax_id].parent  # Move up to parents
        else:
            break
        if gen_tax_id == "1":
            # If it is the root, we reached the end.
            # Add it to the list and break the loop
            ancestors.append(gen_tax_id)
            break
    return ancestors  # Return the list


def _get_all_descendant_nodes(name_object, taxid: str) -> List[str]:
    """ Get all descendant of a node """
    descendant_nodes: List[str] = [taxid]
    if len(name_object[taxid].children) > 0:
        for child in name_object[taxid].children:
            descendant_nodes = descendant_nodes + _get_all_descendant_nodes(name_object, child)
    return descendant_nodes


def _keep_terminal(name_object, nodes_list) -> List[str]:
    """ Keep only terminal nodes """
    terminal_nodes = [x for x in nodes_list if name_object[x].is_tip]
    return terminal_nodes


def _keep_division(name_object, nodes_list, target_division) -> List[str]:
    """ Keep only division nodes """
    division_nodes = [x for x in nodes_list if name_object[x].division == target_division]
    return division_nodes


def get_all_descendants(name_object, target_division: str, taxid: str) -> List[str]:
    """ Get all taxa of a node """

    terminal_nodes = _get_all_descendant_nodes(name_object, taxid)
    terminal_nodes = _keep_division(name_object, terminal_nodes, target_division)

    return terminal_nodes  # Return a list


def load_ncbi_names(task, filename: str = "names.dmp") -> Tuple[Dict, Dict]:
    """Load NCBI names definition ("names.dmp")
    Args:
        filename (str): filename of NCBI names
    Returns:
        name_dict, name_dict_reverse
        :param filename:
        :param task:
    """

    name_dict = {}  # Initialise dictionary with TAX_ID:NAME
    name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

    LOGGER.warning(f"Load {filename}")
    name_file = open(filename, "r")
    global lista
    lista = []
    while 1:
        line = name_file.readline()
        if line == "":
            break
        line = line.rstrip()
        line = line.replace("\t", "")
        tab = line.split("|")
        if tab[3] == "scientific name":
            tax_id, name = tab[0], tab[1]  # Assign tax_id and name ...
            name_dict[tax_id] = name  # ... and load them
            name_dict_reverse[name] = str(tax_id)  # ... into dictionaries
    name_file.close()
    b_file = open("dict_names.pkl", "wb")
    pickle.dump(name_dict_reverse, b_file)
    b_file.close()
    return name_dict, name_dict_reverse


def load_ncbi_taxonomy(name_dict, filename: str = "nodes.dmp"):
    """Load taxonomy NCBI file ("nodes.dmp")
    Args:
        filename (str): filename of ncbi taxonomy
        name_dict (dict): name_dict
    Returns:
    """

    # Define taxonomy variable
    # global name_object
    name_object: Dict = {}

    LOGGER.warning(f"Load {filename}")
    taxonomy_file = open(filename, "r")
    while 1:
        line = taxonomy_file.readline()
        if line == "":
            break
        line = line.replace("\t", "")
        tab = line.split("|")

        tax_id = str(tab[0])
        tax_id_parent = str(tab[1])
        division = str(tab[2])

        # Define name of the taxonomy id
        name = "unknown"
        if tax_id in name_dict:
            name = name_dict[tax_id]

        if tax_id not in name_object:
            name_object[tax_id] = Node()
        name_object[tax_id].tax_id = tax_id  # Assign tax_id
        name_object[tax_id].parent = tax_id_parent  # Assign tax_id parent
        name_object[tax_id].name = name  # Assign name
        name_object[tax_id].division = division  # Assign name

        # Add it has children to parents
        children_list = []
        if tax_id_parent in name_object:
            children_list = name_object[
                tax_id_parent
            ].children  # If parent is in the object
        else:
            name_object[tax_id_parent] = Node()
            name_object[tax_id_parent].tax_id = tax_id_parent  # Assign tax_id
        children_list.append(tax_id)  # ... we found its children.
        name_object[
            tax_id_parent
        ].children = children_list  # ... so add them to the parent

        # As the parent node is found, it is not a terminal node then
        name_object[tax_id_parent].is_tip = False

    a_file = open("dict_species.pkl", "wb")
    pickle.dump(name_object, a_file)
    a_file.close()
    return name_object

def getOrganisms(descendants, assembly, getStrain, excludesp, final_list= [""]):
    for line in assembly:
        if "#" in line:
            continue
        if line == "":
            break
        else:
            line = line.split("\t")
            if line[11] == "Chromosome" or line[11] == "Complete Genome":
                if str(line[5]) in descendants or str(line[6]) in descendants:
                    if excludesp:
                        if "sp." in line[7]:
                            continue
                    if getStrain:
                        if "strain=" in line[8]:
                            strain = line[8].split('=')[1]
                            organism_name = line[7].replace(" ", "_") + "_str." + strain.replace(' ', '_')
                            final_list.append(organism_name + '\t' + line[19] + '\n')
                        else:
                            final_list.append(line[7].replace(" ", "_") + '\t' + line[19] + '\n')
                    else:
                        final_list.append(line[7].replace(" ", "_") + '\t' + line[19] + '\n')
    return final_list


def ftpAvaiability(listOfdescendants, out={}):
    for elem in listOfdescendants:
        organismName = elem.split('\t')
        out[organismName[0]] = organismName[1]
    return out


def main():
    """ Main function """
    # task = False
    # multikeyes = []
    # Load name_dict, name_dict_reverse and taxonomy
    # name_dict, name_dict_reverse = load_ncbi_names(task, filename="names.dmp")  # Load names
    # ncbi_taxonomy = load_ncbi_taxonomy(filename="nodes.dmp", name_dict=name_dict)
    a_file = open("taxonomyNames.pkl", "rb")
    output_names = pickle.load(a_file)
    b_file = open("taxonomyParsedTree.pkl", "rb")
    output_species = pickle.load(b_file)

#  Get all species of a node
    organismName = args.taxa
    try:
        tax_id = output_names[organismName]
    # If taxa doesn't exist raise an exception and exit the execution.
    except KeyError:
        print("ERROR:", args.taxa, "is not a valid Taxa.")
        raise
    target_division = "species"
    print(organismName, "taxa id:", tax_id)

    all_descendants = _get_all_descendant_nodes(output_species, tax_id)
    all_descendants = _keep_division(output_species, all_descendants, target_division)

    print(f"Number of all species and sub-species of {organismName}: {len(all_descendants)}")
    return all_descendants

#  Execute the taxonomy tree parsing in order to find all species under indicated Taxa
if __name__ == "__main__":
    descendants = main()

#  Check the number of descendants
print("The total number of species found in", args.taxa, ":", len(descendants))

#  Open assembly summary from refseq
#  Download it at ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
assembly = open("assembly_summary_refseq.txt", "rt", encoding='utf8')  # This will be taken from server in the future!

#  Set a control over the arguments, that will be used by the function getOrganisms
if args.getstrain:
    getStrain = True
else:
    getStrain = False
if args.exclude_sp:
    exclude_sp = True
else:
    exclude_sp = False

#  Call the function getOrganism and assign it to variable final_list
final_list = getOrganisms(descendants, assembly, getStrain, exclude_sp)

#  Getting the Genome, Feature, Genbank and Proteome data
final_list = final_list[1:]
print("Number of genomes available:", len(final_list))
if len(final_list) == 0:
    print("No Genomes available")
    print("Exiting program . . .")
    sys.exit(0)

# Return a dictionary coupling name of organism and its ftp path on ncbi database
ftps = ftpAvaiability(final_list)
print("ftps disponibility:")
for k, v in ftps.items():
    print(f"For organism {k}, ftps url at: {v}")
#  Check the num argument at command line: numbers of genomes to be shown
if args.num:
    if args.num > len(final_list):
        print("The number of genomes to select is above the total available number of genomes for this taxa")
        print("Exiting program . . .")
        sys.exit(0)
    argument = random.sample(final_list, args.num)
    for elem in argument:
        field1 = elem.split('\t')[0]
        print(f"Species list: {field1}")

#  If num argument not specified, show all genomes found
if not args.num:
    argument = final_list
    for elem in argument:
        field1 = elem.split('\t')[0]
        print(f"Species list: {field1}")

#  If specified at command line, download all genomes previously found
if args.download:
    # System command: make directory for Genome, Proteome and Feature storage
    bashCommand = ["mkdir -p Features", "mkdir -p Genomes", "mkdir -p Proteomes"]
    for command in bashCommand:
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
    # Download data using the ftp url in argument list
    n = 0
    for elem in argument:
        n += 1
        url1 = elem.split('\t')[1].replace("\n", "")
        name = elem.split('\t')[0]
        getname = url1.split('/')[9]
        finalURLfeature = url1 + '/' + getname + '_genomic.gff.gz'
        finalURLgenome = url1 + '/' + getname + '_genomic.fna.gz'
        finalURLproteome = url1 + '/' + getname + '_protein.faa.gz'

        #print(f"\nfiles to be dowloaded: {args.num + 3 - n}\nStarting download of {name} files")
        outputFeature = "Features/" + name.replace('/','-') + '_genomic.gff.gz'
        outputGenome = "Genomes/" + name.replace('/', '-') + '_genomic.fna.gz'
        outputProteome = 'Proteomes/' + name.replace('/', '-') + '_protein.faa.gz'
        final = wget.download(finalURLgenome, outputGenome)
        final2 = wget.download(finalURLfeature, outputFeature)
        final3 = wget.download(finalURLproteome, outputProteome)
    print("\nDONE")
