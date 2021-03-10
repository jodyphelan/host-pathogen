#! /usr/bin/env python
import sys
import subprocess
import argparse
import random
import os
rand_generator = random.SystemRandom()

def which(program):
	import os
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

software_ok = True
for p in ["bedtools","parallel","bcftools","datamash"]:
	sys.stderr.write("Looking for {}...".format(p))
	if not which(p):
		software_ok = False
		sys.stderr.write("Not found\n")
	else:
		sys.stderr.write("Ok\n")

if not software_ok:
	sys.exit("\nError: some dependencies not found, please install before continuing\n")

def get_random_file(prefix = None,extension=None):
	randint = rand_generator.randint(1,999999)
	if prefix:
		if extension:
			return "%s.%s.%s" % (prefix,randint,extension)
		else:
			return "%s.%s.txt" % (prefix,randint)
	else:
		if extension:
			return "%s.tmp.%s" % (randint,extension)
		else:
			return "%s.tmp.txt" % (randint)

def log(msg,ext=False):
	sys.stderr.write("\n"+str(msg)+"\n")
	if ext:
		exit(1)

def run_cmd(cmd,verbose=1,target=None):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	if target and filecheck(target): return True
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/stdout","w")
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")
	else:
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		print("Command Failed! Please Check!")
		exit(1)

def nofile(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isfile(filename):
		return True
	else:
		return False

def get_vcf_prefix(filename):
	if filename[-4:]==".bcf":
		return filename[:-4]
	elif filename[-5:]==".gbcf":
		return filename[:-5]
	elif filename[-7:]==".vcf.gz":
		return filename[:-7]
	elif filename[-4:]==".vcf":
		return filename[:-4]
	else:
		return filename

class vcf_class:
	def __init__(self,filename,threads=4):
		self.samples = []
		self.filename = filename
		self.threads = threads
		self.prefix = get_vcf_prefix(filename)
		if nofile(filename+".csi"):
			run_cmd("bcftools index  %(filename)s" % vars(self))
		self.temp_file = get_random_file()
		run_cmd("bcftools query -l %(filename)s > %(temp_file)s" % vars(self))
		for l in open(self.temp_file):
			self.samples.append(l.rstrip())
		os.remove(self.temp_file)

	def vcf_to_fasta(self,ref_file,threads=4,chunk_size = 50000):
		self.ref_file = ref_file
		self.chunk_size = chunk_size
		self.cmd_split_chr = "bedtools makewindows -g %(ref_file)s.fai -w %(chunk_size)s -s  %(chunk_size)s | awk '{print $1\":\"$2\"-\"$3}'" % vars(self)
		self.tmp_file = "%s.tmp.txt" % self.prefix
		self.threads = threads
		cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"bcftools view  %(filename)s -r {} | bcftools view -V indels | bcftools query -f '%%POS[\\t%%IUPACGT]\\n' | sed 's/\.[\/|]\./N/g' |  datamash transpose > %(prefix)s.{}.tmp.txt\"" % vars(self)
		run_cmd(cmd)
		cmd = "paste `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'` > %(tmp_file)s" % vars(self)
		run_cmd(cmd)
		cmd = "rm `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'`" % vars(self)
		run_cmd(cmd)
		with open(self.prefix+".snps.fa","w") as O:
			for i,l in enumerate(open(self.tmp_file)):
				row = l.rstrip().split()
				if i==0: continue
				s = self.samples[i-1]
				seq = "".join(row)
				O.write(">%s\n%s\n" % ( s,seq))
		run_cmd("rm %s" % self.tmp_file)

def main(args):
	if nofile(args.vcf): quit("Can't find %s... Exiting!" % args.vcf)
	vcf = vcf_class(args.vcf)
	vcf.vcf_to_fasta(args.ref)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--ref',help='VCF file',required=True)
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
