#!/bin/bash

time=/usr/bin/time

export BSCAMPP_LOGGING_LEVEL=debug

queryalnpath=data/query100.fa.gz
backbonealnpath=data/ref10000.fa.gz
backbonetreepath=data/ref10000.tre
bin=../run_scampp.py
t=16

outdir=scampp_output
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

placement_method=epa-ng
model=data/aln_dna.fa.raxml.bestModel
if [[ $1 == "pplacer" ]]; then
    placement_method=pplacer
    model=data/ref10000.fasttree.log
fi

$time -v python3 $bin -i ${model} -t ${backbonetreepath} \
    -d ${outdir} -o placement -a ${backbonealnpath} \
    -q ${queryalnpath} -b 2000 \
    --placement-method ${placement_method} \
    --num-cpus $t \
    #--keeptemp True \
