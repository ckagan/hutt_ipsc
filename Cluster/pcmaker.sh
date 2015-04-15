#!/bin/bash
baselayer='/mnt/gluster/data/internal_supp/hutt_ipsc/GEMMAFiles/PCBaseLayer.txt'

for i in `seq 1 23`
do
   paste <(cat $baselayer) <(cut -f1-$i /mnt/gluster/data/internal_supp/hutt_ipsc/GEMMAFiles/hutt.DConly.ordered.pcs) > /mnt/gluster/data/internal_supp/hutt_ipsc/GEMMAFiles/hutt.DConly.ordered.pc$i
done
