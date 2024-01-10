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

# step 10 find local ancestries using Popsicle, n is window size.  start with 500
#Here, -i is the POPSICLE baseline file generated in Step 7. -j is the clusters file generated in Step 9. -o is the output ancestry file with local ancestry information. -k is the baseline ARFF file generated in Step 8. -l is the directory where temporary files are placed.
#Mix_6_5K: 6 is coming from number of cluster and 1k is the -n number which is the block size in bp. To change the block size, open the FilteredPopsicleBaselineSampled.clusters (Step9) file and change the cluster value and rerun the bottom 3 commands.
echo "step 10: find local ancestries"

#Have to make directory temp_ancestry
java -Xmx730G -jar $BIN/LPDtools.jar POPSICLEIntermediate \
    -i ${RUND}/SampledBaselineFile.txt \
    -j ${RUND}/Popsicle_clusters \
    -o ${RUND}/${RUN}_popsicle_ancestry_win_${block_size}.txt \
    -n ${block_size} \
    -k ${RUND}/${RUN}_popsicle_baseline.arff \
    -l ${RUND}/temp_ancestry

#step 11 polish local ancestries. start with 5. This may not be necessary for sars genome size
#Here, -i is the ancestry file generated in Step 10. -n is the number of blocks that are searched for consistency. -n can be any odd number. If -n 5, the current block, two blocks before the current block and two blocks after the current block are searched and ancestry that is present in maximum blocks is reported. If no maximum ancestry is found, the current ancestry is retained.
echo "step 11: polish local ancestries"
block=5
java -Xmx730G -jar $BIN/LPDtools.jar PolishPosicle -i ${RUND}/${RUN}_popsicle_ancestry_win_${block_size}.txt -n $block

# step 12 arrange ancestries
# Here, -i is the ancestry polished ancestry file generated in Step 11. -o is the file with samples arranged by their cluster membership as specified in the clusters file generated in Step 9
echo "step 12: arrange ancestries"
java -Xmx730G -jar $BIN/LPDtools.jar ArrangePopsicleByClustersFormed -i ${RUND}/${RUN}_popsicle_ancestry_win_${block_size}_polished.txt \
    -o ${RUND}/${RUN}_popsicle_ancestry_win_${block_size}_polished_arranged.txt \
    -j ${RUND}/Popsicle_clusters

echo "To generate circos plot, create an output folder called circos"
mkdir ${RUND}/circos_${block_size}
java -Xmx730G -jar $BIN/LPDtools.jar ConvertPopsicle2CircosHighlights \
    -i ${RUND}/${RUN}_popsicle_ancestry_win_${block_size}_polished_arranged.txt \
    -j ${RUND}/${RUN}_popsicle_baseline.txt \
    -k ${RUND}/Popsicle_clusters \
    -o ${RUND}/circos_${block_size}
