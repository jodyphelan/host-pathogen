import sys
import argparse
import fastq2matrix as fm
from tqdm import tqdm
from collections import defaultdict

def main(args):
    vcf = fm.vcf_class(args.vcf)
    if args.type=="dna" or args.type=="both":
        sys.stderr.write("Loading nucleotide variants\n")
        nt_variants = defaultdict(dict)
        for l in tqdm(fm.cmd_out("bcftools query -i 'GT!=\"ref\"' -f '%%POS[\\t%%SAMPLE:%%GT]\\n' %(vcf)s" % vars(args))):
            if args.debug:
                print(l)
            row = l.strip().split()
            for d in row[1:]:
                s,gt = d.split(":")
                if gt.replace("|","/")!="./.":
                    nt_variants[row[0]][s] = "2"
                else:
                    nt_variants[row[0]][s] = "-9"

        if args.maf:
            filter_out = []
            for pos in nt_variants:
                a = sum([1 for x in nt_variants[pos].values() if x=="2"])
                m = sum([1 for x in nt_variants[pos].values() if x=="-9"])
                alt_freq = a/(len(vcf.samples) - m)
                if alt_freq<args.maf or alt_freq>(1-args.maf):
                    filter_out.append(pos)
            sys.stderr.write("Filtering out %s variants which don't pass maf criteria\n" % len(filter_out))
            for pos in filter_out:
                del nt_variants[pos]

        sys.stderr.write("Writing output\n")
        with open(args.out+".dna.pheno.txt","w") as O:
            if args.format=="plink2":
                O.write("#IID\t%s\n" % ("\t".join(["var_%s" % p for p in nt_variants])))

            for s in tqdm(vcf.samples):
                if args.format=="plink1":
                    O.write("0\t%s" % s)
                else:
                    O.write(s)
                O.write("\t%s\n" % ("\t".join([nt_variants[pos].get(s,"1") for pos in nt_variants])))

        if args.format=="plink1":
            with open(args.out+".dna.phenotype_map.txt","w") as O:
                for i,pos in enumerate(nt_variants):
                    O.write("%s\t%s\n" % (i+1,pos))


    if args.type=="aa" or args.type=="both":

        bcsq_found = False
        for l in fm.cmd_out("bcftools view %(vcf)s -h | grep BCSQ" % vars(args)):
            if "BCSQ" in l: bcsq_found = True
        if not bcsq_found:
            quit("\n########## ERROR ##########\n\nBCSQ tag not found in vcf...Quitting!\nHave you annotated it with bcftools csq yet?")
        sys.stderr.write("Loading amino acid variants\n")
        variants = defaultdict(set)
        missing = defaultdict(set)
        csq2pos = {}
        vcf = fm.vcf_class(args.vcf)
        for l in tqdm(fm.cmd_out("bcftools query -f '[%%POS\\t%%SAMPLE\\t%%TBCSQ\\n]' %(vcf)s" % vars(args))):
            if args.debug:
                print(l)
            row = l.strip().split()
            tbcsq = row[2].split("|")

            if row[2][0]=="@": continue
            if "synonymous" in tbcsq[0]: continue
            if "non_coding" in tbcsq[0]: continue
            if "start_lost" in tbcsq[0]: tbcsq.append("1M>1*")
            variants[(tbcsq[2],tbcsq[5])].add(row[1])
            csq2pos[(tbcsq[2],tbcsq[5])] = row[0]


        for l in tqdm(fm.cmd_out("bcftools query -i 'GT=\"mis\"' -f '%%POS\\t%%BCSQ[\\t%%SAMPLE]\\n' %(vcf)s" % vars(args))):
            if args.debug:
                print(l)
            row = l.strip().split()
            if row[1]==".": continue
            for csq in row[1].split(","):
                tmp = csq.split("|")

                if csq[0]=="@": continue
                if "synonymous" in csq: continue
                if "non_coding" in csq: continue
                if "start_lost" in csq: tmp.append("1M>1*")
                mut = (tmp[2],tmp[5])

                missing[mut] = set(row[2:])

        if args.maf:
            filter_out = []
            for key in variants:
                a = len(variants[key])
                m = len(missing[key])
                alt_freq = a/(len(vcf.samples) - m)
                if alt_freq<args.maf or alt_freq>(1-args.maf):
                    filter_out.append(key)
            sys.stderr.write("Filtering out %s variants which don't pass maf criteria\n" % len(filter_out))
            for key in filter_out:
                del variants[key]

        with open(args.out+".aa.pheno.txt","w") as O:
            if args.format=="plink2":
                O.write("#IID\t%s\n" % ("\t".join(["%s_%s" % (mut[0],mut[1].replace(">","_")) for mut in variants])))
            for s in tqdm(vcf.samples):
                sample_vector = []
                for mut in variants:
                    if s in missing[mut]:
                        sample_vector.append("-9")
                    elif s in variants[mut]:
                        sample_vector.append("2")
                    else:
                        sample_vector.append("1")

                if args.format=="plink1":
                    O.write("0\t%s\t%s\n" % (s,"\t".join(sample_vector)))
                else:
                    O.write("%s\t%s\n" % (s,"\t".join(sample_vector)))

        with open(args.out+".aa.phenotype_map.txt","w") as O:
            for i,mut in enumerate(variants):
                O.write("%s\t%s\t%s\t%s\t%s_%s\n" % (i+1,mut[0],mut[1],csq2pos[(mut[0],mut[1])], mut[0],mut[1].replace(">","_")))



        dna_patterns = defaultdict(list)
        for i,l in enumerate(fm.cmd_out("cat %s.dna.pheno.txt | datamash transpose" % args.out)):
            row = l.strip().split()
            if args.format=="plink2":
                if i<1: continue
                dna_patterns[tuple(row[1:])].append(row[0])
            elif args.format=="plink1":
                if i<2: continue
                dna_patterns[tuple(row)].append(str(i-1))

        with open(args.out+".dna.pheno.compressed.txt","w") as O:
            if args.format=="plink2":
                O.write("#IID\t%s\n" % ("\t".join(["variant_"+str(i).zfill(4) for i in range(len(dna_patterns))])))
            for i,s in enumerate(vcf.samples):
                if args.format=="plink1":
                    O.write("0\t%s\t%s\n" % (s,"\t".join([key[i] for key in list(dna_patterns)])))
                else:
                    O.write("%s\t%s\n" % (s,"\t".join([key[i] for key in list(dna_patterns)])))

        with open(args.out+".dna.pheno.compressed.map.txt","w") as O:
            for i in range(len(dna_patterns)):
                O.write("%s\t%s\n" % ("variant_"+str(i).zfill(4), ",".join(list(dna_patterns.values())[i])))

        patterns = defaultdict(list)
        for i,l in enumerate(fm.cmd_out("cat %s.aa.pheno.txt | datamash transpose" % args.out)):
            row = l.strip().split()
            if args.format=="plink2":
                if i<1: continue
                patterns[tuple(row[1:])].append(row[0])
            elif args.format=="plink1":
                if i<2: continue
                patterns[tuple(row)].append(str(i-1))

        with open(args.out+".aa.pheno.compressed.txt","w") as O:
            if args.format=="plink2":
                O.write("#IID\t%s\n" % ("\t".join(["variant_"+str(i).zfill(4) for i in range(len(patterns))])))
            for i,s in enumerate(vcf.samples):
                if args.format=="plink1":
                    O.write("0\t%s\t%s\n" % (s,"\t".join([key[i] for key in list(patterns)])))
                else:
                    O.write("%s\t%s\n" % (s,"\t".join([key[i] for key in list(patterns)])))

        with open(args.out+".aa.pheno.compressed.map.txt","w") as O:
            for i in range(len(patterns)):
                O.write("%s\t%s\n" % ("variant_"+str(i).zfill(4), ",".join(list(patterns.values())[i])))

parser = argparse.ArgumentParser(description='Convert vcf to phenotype files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--out',help='Prefix for the output files',required=True)
parser.add_argument('--format',choices=["plink1","plink2"],help='VCF file',required=True)
parser.add_argument('--type',choices=["dna","aa","both"],help='The type of output file generated',required=True)
parser.add_argument('--maf',type=float,help='Minor allele frequency cutoff')
parser.add_argument('--debug',action="store_true",help='Print out lines to debug')

parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
