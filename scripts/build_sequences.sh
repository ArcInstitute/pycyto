#!/bin/bash

NUM_LANES=8;
NUM_SAMPLES=10;
NUM_RUNS=4;
NUM_READS=10000;
COMPRESSION=zst
OUTDIR="./sequences"

for lane in $(seq 1 $NUM_LANES); do
    for sample in $(seq 1 $NUM_SAMPLES); do
        for run in $(seq 1 $NUM_RUNS); do
            for mode in GEX CRISPR; do
                subname=${OUTDIR}/S${sample}_Lane${lane}_${mode}_S${run}
                echo "Generating reads for ${subname}"
                nucgen -n $NUM_READS ${subname}_R1.fq.${COMPRESSION}
                nucgen -n $NUM_READS ${subname}_R2.fq.${COMPRESSION}
            done
        done
    done
done
