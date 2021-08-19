from collections import defaultdict
import argparse
import fastq2matrix as fm
from uuid import uuid4
from Bio import Phylo 
import ete3
from tqdm import tqdm

def nexus2newick(intree,outtree):
    trees = Phylo.parse(intree, 'nexus')
    for tree in trees:
        pass

    for x in tree.find_clades():
        x.comment = None
        if x.confidence!=None:
            x.name = x.confidence
            x.confidence = None

    Phylo.write([tree],outtree,"newick")

def main(args):
    id = "7ae737a4-27b7-42f7-89c4-ce3e3e3fc970" 
    # id = str(uuid4())
    # fm.run_cmd(f"treetime ancestral --aln {args.vcf} --vcf-reference {args.ref} --tree {args.tree}  --outdir {id}")


    # nexus2newick(f"{id}/annotated_tree.nexus", f"{args.out}.annotated_tree.nwk")

    tree = ete3.Tree(f"{args.out}.annotated_tree.nwk",format=1)
    node_names = set([tree.name] + [n.name.split("/")[0] for n in tree.get_descendants()])
    leaf_names = set(tree.get_leaf_names())
    internal_node_names = node_names - leaf_names

    states = defaultdict(dict)
    for l in tqdm(fm.cmd_out(f"bcftools query -f '[%POS %REF %ALT %SAMPLE %GT\\n]' {id}/ancestral_sequences.vcf")):
        row = l.strip().split()

        alleles = [row[1]]+row[2].split(",")
        if row[4]==".":
            a = row[1]
        else:
            a = alleles[int(row[4][0])]
        states[row[0]][row[3]] = a


    acgt = set(["A", "C", "G", "T", "a", "c", "g", "t"])
    convergent_sites = []
    with open(f"{args.out}.homoplasies.txt","w") as O:
        for site in tqdm(list(states)):
            nucleotides = set([states[site][n] for n in node_names])
            if len(nucleotides)==1: continue

            # Set up storage objects
            origins = []

            tree.add_feature("state", states[site][tree.name])
            for n in tree.traverse():
                if n == tree: continue
                node_state = states[site][n.name]
                if node_state != n.get_ancestors()[0].state and node_state in acgt and n.get_ancestors()[0].state in acgt:
                    origins.append(n)
                    O.write("%s\t%s\t%s\n" % (site,n.name,len(n.get_leaf_names())))
                    
                n.add_feature("state", node_state)
            if len(origins) > 1:
                convergent_sites.append((site,  origins))

    with open(f"{args.out}.num_changes.txt","w") as O:
        for site in convergent_sites:
            O.write("%s\t%s\n" % (site[0],len(site[1])))

    fm.run_cmd(f"rm -r {id}")



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='VCF file',required=True)
parser.add_argument('--ref',type=str,help='VCF file',required=True)
parser.add_argument('--tree',type=str,help='VCF file',required=True)
parser.add_argument('--out',type=str,help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)