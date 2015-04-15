#!/bin/bash
for i in $(seq 1 1 23)
do
   nohup python /mnt/gluster/data/internal_supp/hutt_ipsc/Scripts/eqtl_driver.py $i &
done
