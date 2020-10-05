import sys
import argparse
import fastq2matrix
from tqdm import tqdm


def main(args):
    vcf = fastq2matrix.vcf_class(args.vcf)
    for l in tqdm(fastq2matrix.cmd_out("bcftools query -f '%%POS[\\t%%GT]\\n' %(vcf)s" % vars(args))):
        row = l.strip().split()
        with open("pheno_%s.txt" % row[0],"w") as O:
            if args.format=="plink2":
                O.write("#IID\tvariant\n")

            for i in range(1,len(row)):
                if args.format=="plink1":
                    O.write("0\t%s" % vcf.samples[i-1])
                elif args.format=="plink2":
                    O.write("%s" % vcf.samples[i-1])

                if row[i]=="0/0":
                    O.write("\t1\n")
                elif row[i]=="1/1":
                    O.write("\t2\n")
                elif row[i]=="2/2":
                    O.write("\t2\n")
                elif row[i]=="3/3":
                    O.write("\t2\n")
                else:
                    O.write("\t-9\n")

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--format',choices=["plink1","plink2"],help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
