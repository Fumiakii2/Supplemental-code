#!/bin/sh

#Set derectory
Workdir=/Users/fumiakii2/Desktop/temp/POPSICLE_TUTORIAL
cd ${Workdir}

# This script adapted from: https://popsicle-admixture.sourceforge.io/AnalyticalPipelinePopulStr.html

# run name
RUN="POPSICLE_TUTORIAL"
IDS="${RUN}.txt"

# tool directory
BIN="${Workdir}"

# directories
WD="${Workdir}"
RUND="${WD}/${RUN}_analysis"
REF="${WD}/ref"
BAM="${WD}/bam"
block_size=10000

#move your reference file to "ref"

# Only BAM files are needed, do not put bai files together
SOMIES="${RUND}/somies"
ALLELE="${RUND}/alleles"
LIST="${RUND}/listfiles"
SNPS="${RUND}/snps"

# create run directories
mkdir ${RUN}_analysis
mkdir -p $BAM $SOMIES $ALLELE $SNPS $LIST
mkdir ${RUND}/temp_ancestry
mkdir ${RUND}/plot

# step 1 generate somies files
#Here, -i is the directory where sorted aligned bam files are placed. -o is the output file containing somies of chromosomes in each of the samples. -m is the ploidy of the organism and -k is the file with 2 columns. The first column is the chromosome names and the second column is size of the chromosome.
echo "step 1: generate somies file..."
java -Xmx730G -jar $BIN/LPDtools.jar FindSomies -i $BAM -o $SOMIES/${RUN}_somies.txt -m 1 -k $REF/Genome.txt

# step 2 find alleles
# Here, -i is the directory where sorted aligned Binary Alignment map files are placed. -o is the output directory containing the allele files. -n is the bp window. -n 1 finds allele composition at each base
echo "step 2: find alleles..."
java -Xmx730G -jar $BIN/LPDtools.jar findAlleles -i $BAM -n 1 -o $ALLELE

echo "Step 2.5:Listalleles and findsnps"
java -jar $BIN/LPDtools.jar listAlleles -i $ALLELE/ -m 0.4 -n 5 -o $LIST/

java -jar $BIN/LPDtools.jar findSNPs -i $LIST/ -j $REF/ToxoDB-57_TgondiiME49_Genome.fasta -k "vcf" -m 0.25 -o $SNPS/

echo "step 3: generate popsicle input..."
# step 3 generate popsicle input
#Here, -i is the directory where allele files are present. These are the files generated in Step 2 of this pipeline. -j is the directory where the Single nucleotide Polymorphism (SNP files) generated using a utility such as samtools are placed. Alternately, one may use their own markers of choice. -k is the format of the snp files (VCF is the preferred format. Markers can be submitted in tab delimited format using "tab" tag). -l is the somies file generated using Step 1 of the POPSICLE pipeline

java -Xmx730G -jar $BIN/LPDtools.jar GenerateInputFromAlleleFiles -i $ALLELE -j $SNPS -o ${RUND}/${RUN}_popsicle_in.txt -k "vcf" -l $SOMIES/${RUN}_somies.txt

# step 4 remove non-variant sites - not going to do
#/Users/fumiakii2/Desktop/temp/popsicle_gifu/GUY-2004-JAG1_AF.vcf
#Here, -i is the POPSICLE input file generated using the utility GenerateInputFromAlleleFiles utility of POPSICLE pipeline (see Step 3). -o is the filtered output file after removing markers that dont pass the filter. -m is the factor that determines which markers are filtered out. If -m is set to 0.8 and if a marker has 80% of the samples with identical allele, it is filtered out.
#Not popsicleInputFile, use output file from step3

echo "step4"
java -Xmx730G -jar $BIN/LPDtools.jar RemoveInsignificantLoci -i ${RUND}/${RUN}_popsicle_in.txt -o ${RUND}/OutputFilteredPopsicleFile1.txt -m 0.9

# step 5 remove loci with missing
#Here -i is the popsicle input file generated in Step 3 or filtered popsicle file generated in Step 4. -o is the output file generated after removing loci with lots of missing data. -m is the maximum missing data tolerable. Eg. if -m is set to 0.7, markers with 70% or more missing data are ignored.
 echo "step 5: remove loci with missing..."
java -Xmx730G -jar $BIN/LPDtools.jar RemoveLociWithLotsOfMissingData -i ${RUND}/OutputFilteredPopsicleFile1.txt -o ${RUND}/${RUN}_popsicle_filtered.txt -m 0.7

# step 6 find divergence from baseline
#Here -i is the input POPSICLE file generated in Step 3 or any of the filtered versions generated in Steps 5 or 6. -o is the output file.
 echo "step 6: find divergence from baseline..."
java -Xmx730G -jar $BIN/LPDtools.jar FindDivergenceFromBaseline -i ${RUND}/${RUN}_popsicle_filtered.txt -o ${RUND}/${RUN}_popsicle_baseline.txt

# step 7 sample baseline to generate a smaller file - not necessary, probably not recommended with low variation This sampled file is used for faster processing such as for clustering.
#Here, -i is the baseline file generated in Step 6. -o is the sampled version of the file. -m is the factor that determines the number of markers that are to be retained per kb. If -m 2, then 2 markers every kb are retained.
echo "Step7"
java -Xmx730G -jar $BIN/LPDtools.jar SampleBaseLineFile -i ${RUND}/${RUN}_popsicle_baseline.txt -o ${RUND}/SampledBaselineFile.txt -m 2

# step 8 convert files from 6/7 to .ARFF format
#Here, -i is the baseline files generated using Steps 6 and 7. -o are the output files in .arff format (See WEKA for details on arff file format)
echo "step 8:  convert files to arff"
echo "input from Step6"
java -Xmx730G -jar $BIN/LPDtools.jar Convert2ARFFformat -i ${RUND}/${RUN}_popsicle_baseline.txt -o ${RUND}/${RUN}_popsicle_baseline.arff
echo "input from Step7"
java -Xmx730G -jar $BIN/LPDtools.jar Convert2ARFFformat -i ${RUND}/SampledBaselineFile.txt -o ${RUND}/${RUN}_popsicle_baselineSampled.arff

# step 9 cluster samples to find major sample groups. - for dev253 number of clusters recommended was 5

#Here, -i is the sampled ARFF file generated in Step 8. -o is the output clusters file that indicates which samples are assigned to which clusters and a  score associated with such clustering. -n is the minimum number of clusters to which the samples are assigned. -m is the maximum number of clusters to which the samples are assigned.
echo "step 9: cluster samples to find major groups"
java -Xmx120G -jar $BIN/LPDtools.jar PerformKmeansClustering -i ${RUND}/${RUN}_popsicle_baselineSampled.arff -o ${RUND}/Popsicle_clusters -n 2 -m 4
