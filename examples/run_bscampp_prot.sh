#!/bin/bash

time=/usr/bin/time

export BSCAMPP_LOGGING_LEVEL=debug

indir=data/protein_example
queryalnpath=$indir/query.fasta
backbonealnpath=$indir/backbone.fasta
backbonetreepath=$indir/backbone.tre
bin=../run_bscampp.py
t=16

outdir=bscampp_output_prot
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

# to demonstrate that subtree placements are working properly
subtreesize=200

placement_method=epa-ng
model=$indir/true.fasta.raxml.bestModel
if [[ $1 == "pplacer" ]]; then
    placement_method=pplacer
    model=$indir/fasttree.log
fi

$time -v python $bin -i ${model} -t ${backbonetreepath} \
    -d ${outdir} -o placement -a ${backbonealnpath} \
    -q ${queryalnpath} -b ${subtreesize} \
    --placement-method ${placement_method} \
    --num-cpus $t \
    --keeptemp True
