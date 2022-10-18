#!/bin/bash

# INPUT
VCF=$1
PLOT=$2

echo "[INFO] VCF: "$(basename $VCF)

vcf2circos --input=$VCF --output=$PLOT.html --options=options/options.default.json

