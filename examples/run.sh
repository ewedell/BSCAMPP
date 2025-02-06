#!/bin/bash

time=/usr/bin/time

model=data/aln_dna.fa.raxml.bestModel
queryalnpath=data/query100.fa
backbonealnpath=data/ref10000.fa
backbonetreepath=data/ref10000.tre
#bin=../BSCAMPP_code/EPA-ng-BSCAMPP.py
bin=../run_bscampp.py

outdir=bscampp_output
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

placement_method=pplacer
{ $time -v python $bin -i ${model} -t ${backbonetreepath} \
    -d ${outdir} -o placement -a ${backbonealnpath} \
    -q ${queryalnpath} -b 2000 \
    --placement-method ${placement_method} ; } 2> $outdir/runtime.txt
