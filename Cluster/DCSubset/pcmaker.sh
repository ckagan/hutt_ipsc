#!/bin/bash
baselayer='/mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/PCBaseLayer.txt'

for i in `seq 1 43`
do
   paste <(cat $baselayer) <(cut -f1-$i /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/qqnorm.newck.gccor.newcovcor.pcs) > /mnt/gluster/data/internal_supp/hutt_ipsc/DCsubset/qqnorm.newck.gccor.newcovcor.pc$i
done
