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
