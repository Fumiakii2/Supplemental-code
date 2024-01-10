#!/bin/sh

#Set derectory
Workdir=/Users/fumiakii2/Desktop/temp/POPSICLE_TUTORIAL
cd ${Workdir}

# Pipeline adapted from: https://popsicle-admixture.sourceforge.io/AnalyticalPipelinePopulStr.html

# run name
RUN="POPSICLE_TUTORIAL"
IDS="${RUN}.txt"

# tool directory
BIN="${Workdir}"

# directories
WD="${Workdir}"
RUND="${WD}/${RUN}_analysis"
REF="${Workdir}/ref"
BAM="${Workdir}/bam"
block_size=10000

##Only BAM files are needed, do not put bai files together
SOMIES="${RUND}/somies"
ALLELE="${RUND}/alleles"
LIST="${RUND}/listfiles"
SNPS="${RUND}/snps"

# Before getting circos plot, review http://circos.ca/documentation/tutorials/configuration/
echo "Circos"
cd ${RUND}/circos_${block_size}
circos -conf circos.conf -outputfile ${block_size}.png
