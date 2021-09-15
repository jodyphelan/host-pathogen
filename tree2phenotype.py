import sys
import argparse
import ete3

def main(args):
    tree = ete3.Tree(args.tree,format=1)
    if len(tree.children)!=2:
        R = tree.get_midpoint_outgroup()
        tree.set_outgroup(R)

    all_samples = set(tree.get_leaf_names())
    phenos = {}
    i = 0

    for node in tree.traverse("preorder"):

        if node.is_leaf():
            continue
        i+=1
        # if node.is_root():
        #     continue
        if node.name=="":
            node.name="100"

        bootstrap = node.name
        node_id = "N%s" % (i)
        node.name = "%s" % (node_id)
        if float(bootstrap)>args.bootstrap:
            tmp = {}
            node_samples = set(node.get_leaf_names())
            for s in node_samples:
                tmp[s] = "2"
            for s in all_samples - node_samples:
                tmp[s] = "1"
            phenos[node_id] = tmp

    
    tree.write(format=1, outfile=args.out+".annotated_tree.nwk")
    with open(args.out+".phenos.txt","w") as O:
        O.write("#IID\t%s\n" % "\t".join(list(phenos)))
        for s in tree.get_leaf_names():
            O.write("%s\t%s\n" % (s,"\t".join([phenos[p][s] for p in phenos])))

    print("Written %s phenotypes to %s.phenos.txt" % (len(phenos),args.out))
    print("Written annotated tree to %s.annotated_tree.nwk" % (args.out))
    

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--tree',type=str,help='Newick formatted tree with bootstraps')
parser.add_argument('--out',type=str,help='Prefix for output files')
parser.add_argument('--bootstrap',type=int,default=95,help='Minimum bootstrap to use clade as phenotype')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)