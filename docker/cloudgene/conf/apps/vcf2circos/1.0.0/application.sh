#!/bin/bash

# INPUT
VCF=$1
PLOT=$2

echo "# VCF file is "$(basename $VCF)"."

#vcf2circos --input=$VCF --output=$PLOT.html --options=options/options.example.json
vcf2circos --input=$VCF --output=$PLOT.html --options=/app/demo_data/options.example.json
