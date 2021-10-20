from collections import defaultdict
import argparse
from uuid import uuid4
from Bio import Phylo 
import ete3
from tqdm import tqdm
import subprocess as sp
import os

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
    # id = "fc22cdcd-380b-49a3-9548-6f0403ce1d12" 
    id = str(uuid4())
    if args.filter_low_af:
        sp.call(f"bcftools view -c 3 {args.vcf} -Oz -o {id}.reduced.vcf.gz",shell=True)
        sp.call(f"treetime ancestral --aln {id}.reduced.vcf.gz --vcf-reference {args.ref} --tree {args.tree}  --outdir {id}",shell=True)
    else:
        sp.call(f"treetime ancestral --aln {args.vcf} --vcf-reference {args.ref} --tree {args.tree}  --outdir {id}",shell=True)


    nexus2newick(f"{id}/annotated_tree.nexus", f"{args.out}.temp.annotated_tree.nwk")

    tree = ete3.Tree(f"{args.out}.temp.annotated_tree.nwk",format=1)
    node_names = set([tree.name] + [n.name.split("/")[0] for n in tree.get_descendants()])
    leaf_names = set(tree.get_leaf_names())

    states = defaultdict(dict)
    for l in tqdm(sp.Popen(f"bcftools query -f '[%POS %REF %ALT %SAMPLE %GT\\n]' {id}/ancestral_sequences.vcf",shell=True,stdout=sp.PIPE).stdout):
        row = l.decode().strip().split()

        alleles = [row[1]]+row[2].split(",")
        if row[4]==".":
            a = row[1]
        else:
            a = alleles[int(row[4][0])]
        states[row[0]][row[3]] = a

    name_conversion = {}
    i=0
    tree_renamed = tree.copy()
    if args.map_to_tree:
        import re
        tmp = open(args.tree).read()
        tmp = re.sub("\[.+?\]","",tmp)
        guide_tree = ete3.Tree(tmp,format=1)

    for n in tree_renamed.traverse("preorder"):
        if not n.is_leaf():
            i+=1
            if args.map_to_tree:
                new_id = guide_tree.get_common_ancestor(n.get_leaf_names()).name
            else:
                new_id = "N%s" % str(i)
            name_conversion[n.name] = new_id
            n.name =  new_id
        else:
            name_conversion[n.name] = n.name
        
    tree_renamed.write(format=1,outfile = f"{args.out}.annotated_tree.nwk")

    acgt = set(["A", "C", "G", "T", "a", "c", "g", "t"])
    convergent_sites = []
    with open(f"{args.out}.branch_changes.txt","w") as O:
        O.write("Position\tAncestor_Node\tDerived_Node\tAncestor_Call\tDerived_Call\tClade_Size\n")
        for site in tqdm(list(states)):
            nucleotides = set([states[site][n] for n in node_names])
            if len(nucleotides)==1: continue

            # Set up storage objects
            origins = []

            tree.add_feature("state", states[site][tree.name])
            for n in tree.traverse("preorder"):
                if n == tree: continue
                node_state = states[site][n.name]
                parent_state = n.get_ancestors()[0].state
                if node_state != parent_state and node_state in acgt and parent_state in acgt:
                    origins.append(n)
                    O.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (site,name_conversion[n.get_ancestors()[0].name], name_conversion[n.name],parent_state,node_state,len(n.get_leaf_names())))
                    
                n.add_feature("state", node_state)
            if len(origins) > 1:
                convergent_sites.append((site,  origins))



    with open(f"{args.out}.num_changes.txt","w") as O:
        O.write("Position\tNum_Homoplasies")
        for site in convergent_sites:
            O.write("%s\t%s\n" % (site[0],len(site[1])))

    sp.call(f"rm -r {id}*",shell=True)



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='VCF file',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta file',required=True)
parser.add_argument('--tree',type=str,help='Tree file',required=True)
parser.add_argument('--out',type=str,help='Ouput prefix for results',required=True)
parser.add_argument('--map-to-tree',action="store_true",help='Use the internal node names of the input tree')
parser.add_argument('--filter-low-af',action="store_true",help="Don't reconstruct variants which appear in just one sample")
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
