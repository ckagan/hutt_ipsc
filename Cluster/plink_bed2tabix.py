#!/usr/bin/env python
import subprocess
import time
import glob
import sys
sys.path.append('/mnt/lustre/home/cusanovich/Programs/')
sys.path.append('/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/')
from DarrenTools import ifier

#step 1 - chrom beds
#step 2 - chrom raws
#step 3 - chrom bim (snpid + position)
#step 4 - load raw
#step 5 - transpose raw
#step 6 - add chrom, start, end, snpid, "0", "."
#step 7 - save table
#step 8 - compress with bgzip
#step 9 - index with tabix

genodir = '/mnt/gluster/data/internal_supp/hutt_ipsc/Genotypes/'
outname = 'hutt.all.imputed'
print 'Creating raw files...'
for j in range(1,23):
	#plinker = 'echo "plink --noweb --nonfounders --maf 0.05 --geno 0.05 --bfile ' + genodir + 'imputed_cgi --chr ' + str(j) + ' --make-bed --out ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + '; plink --bfile ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + ' --recodeA --out ' + genodir + 'ByChr/hutt.imputed.chr' + str(j) + '; touch ' + genodir + 'ByChr/chr' + str(j) + '.done" | qsub -l h_vmem=2g -o ~/dump/ -e ~/dump/'
	plinker = 'echo "plink --noweb --nonfounders --bfile ' + genodir + 'hutt.imputed.73subset --chr ' + str(j) + ' --make-bed --out ' + genodir + 'ByChr/' + outname + '.chr' + str(j) + '; plink --bfile ' + genodir + 'ByChr/' + outname + '.chr' + str(j) + ' --recodeA --out ' + genodir + 'ByChr/' + outname + '.chr' + str(j) + '; touch ' + genodir + 'ByChr/chr' + str(j) + '.done" | qsub -l h_vmem=2g -o ~/dump/ -e ~/dump/'
	ifier(plinker)

while len(glob.glob(genodir + 'ByChr/*.done')) < 22:
	time.sleep(5)

cleanup = "rm " + genodir + "ByChr/*.done"
ifier(cleanup)

print 'Creating bed files...'
for j in range(1,23):
	converter = 'echo "python /mnt/gluster/data/internal_supp/hutt_ipsc/Scripts/raw2txt.py ' + str(j) + ' ' + genodir + ' ' + outname + '" | qsub -l h_vmem=8g -o ~/dump/ -e ~/dump/'
	ifier(converter)

while len(glob.glob(genodir + 'ByChr/*.done')) < 22:
	time.sleep(5)

ifier(cleanup)

print 'Creating tabix files...'
for j in range(1,23):
	tabixer = 'echo "bgzip ' + genodir + 'ByChr/' + outname + '.chr' + str(j) + '.txt; tabix -p bed ' + genodir + 'ByChr/' + outname + '.chr' + str(j) + '.txt.gz; touch ' + genodir + 'ByChr/' + str(j) + '.done" | qsub -l h_vmem=2g -o ~/dump/ -e ~/dump/ -wd `pwd`'
	ifier(tabixer)

while len(glob.glob(genodir + 'ByChr/*.done')) < 22:
	time.sleep(5)

print 'Cleaning up a bit...'
ifier(cleanup)
for j in range(1,23):
	for k in ['.bed','.bim','.fam','.log','.nof','.raw']:
		cleaner = 'rm ' + genodir + 'ByChr/' + outname + '.chr' + str(j) + k
		ifier(cleaner)

print 'Work complete.'
