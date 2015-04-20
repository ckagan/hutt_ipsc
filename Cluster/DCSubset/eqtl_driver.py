#!/usr/bin/env python

import subprocess
import glob
import time
import sys

if len(sys.argv) != 2:
	print "Usage = python eqtl_driver.py pc"
	sys.exit()

pcs = sys.argv[1]

def ifier(commander):
	ify = subprocess.Popen(commander,shell=True)
	ify.wait()

for i in range(22,0,-1):
	eqtler = 'echo "python /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/Scripts/alt_gemma_noplink_eqtl_mapper.py chr' + str(i) + ' ' + str(pcs) + '" | qsub -l h_vmem=2g,bigio=4 -V -o ~/dump/ -e ~/dump/ -N "eQTLs.chr' + str(i) + '.PC' + str(pcs) + '"'
	ifier(eqtler)

while len(glob.glob('/mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/ByChr/*.PC' + str(pcs) + '.bonferroni.done')) < 22:
	time.sleep(300)

cleanup = "rm /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/ByChr/*.PC" + str(pcs) + ".bonferroni.done"
ifier(cleanup)

masterer = 'cat /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/ByChr/*.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt | sort -k1,1 -k2,2 > /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/eQTLs/master.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt; cat /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/ByChr/*.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt | sort -k1,1 >  /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/eQTLs/master.PC' + str(pcs) + '.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt'
ifier(masterer)
