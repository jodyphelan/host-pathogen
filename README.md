# host-pathogen

Scripts for host-pathogen project

## vcf2pheno
```
usage: vcf2pheno.py [-h] --vcf VCF --out OUT --format {plink1,plink2} --type
                    {dna,aa,both} [--maf MAF]

Convert vcf to phenotype files

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             VCF file (default: None)
  --out OUT             Prefix for the output files (default: None)
  --format {plink1,plink2}
                        VCF file (default: None)
  --type {dna,aa,both}  The type of output file generated (default: None)
  --maf MAF             Minor allele frequency cutoff (default: None)
```

## hp-vcf2fasta.py
```
usage: hp-vcf2fasta.py [-h] --vcf VCF --ref REF [--threads THREADS]

Convert vcf to fasta

optional arguments:
  -h, --help         show this help message and exit
  --vcf VCF          VCF file (default: None)
  --ref REF          Reference fasta file (default: None)
  --threads THREADS  Number of threads for parallel operations (default: 4)
```

# Detecting homoplasy 

First you will need to run `vcf_find_homoplasy.py` from the fastq2matrix package. If you can't find it please update the package on your computer. The inputs are a VCF file, the reference fasta and a rooted newick formatted tree. 

```
vcf_find_homoplasy.py  --vcf test_data.vcf.gz --ref MTB.fa --tree rooted.nwk  --out output
```

Three files will be produced, but we will only use the output.homoplasies.txt file in the next step which will use `select_homoplasy.R` either interactively or through using `Rscript`. Three parameters are set to classify mutations are homoplastic: 1) the minimum size for the biggest clade which aquired the mutation, 2) the minimum numer of other acquisitions on the tree and 3) the minimum number of samples in other (non-main) clades. When all these parameters pass a SNP will be called homoplastic.
