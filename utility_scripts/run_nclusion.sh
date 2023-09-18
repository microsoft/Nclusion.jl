#!/bin/bash
INPUTFILE=$1
JULIAENV=$2
DATASET=$3
OUTPUTDIR=$4
NUM_THREADS=auto

a=1.0
b=1.0
k=25
elbo_ep="1.0"
julia --project=${JULIAENV} --thread=${NUM_THREADS} run_nclusion.jl  "${INPUTFILE}" "${k}" "${a}" "${b}"  "12345" "${elbo_ep}" "${DATASET}" "${OUTPUTDIR}" 