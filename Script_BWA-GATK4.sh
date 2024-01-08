#!/bin/sh

##Set Working directory
Workdir='''Set Working directory'''
cd ${Workdir}

##Set directories
Tgindex='''Set Path to Reference fasta'''
fastqdir='''Set Path to fastq files'''
PROJECT_PATH=${Workdir}

mkdir -p ${PROJECT_PATH}
cd ${PROJECT_PATH}

mkdir fastq
mkdir cleaned_fastq
mkdir bam
mkdir bqsr
mkdir vcf

cd $fastqdir
for fpath in `ls *_R1.fastq.gz`
do
    fname=${fpath%_R1.fastq.gz}
    trimmomatic -XX:MaxRAM=730G -XX:MaxRAMPercentage=75 \
    PE -threads 24 \
    $fastqdir/${fname}_R1.fastq.gz $fastqdir/${fname}_R2.fastq.gz \
    ${fname}_cleaned_1.fastq.gz ${fname}_unpaired_1.fastq.gz \
    ${fname}_cleaned_2.fastq.gz ${fname}_unpaired_2.fastq.gz \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20
done

cd ${PROJECT_PATH}/cleaned_fastq
for fpath in `ls *_cleaned_1.fastq.gz`
do
    fname=${fpath%_cleaned_1.fastq.gz}
    bamRG="@RG\tID:L\tSM:"${fname}"\tPL:illumina\tLB:lib1\tPU:unit1"
    bwa mem -R ${bamRG} $Tgindex \
            ${fname}_cleaned_1.fastq.gz ${fname}_cleaned_2.fastq.gz > ${fname}.sam
done

mv *.sam ${PROJECT_PATH}/bam/

cd ${PROJECT_PATH}/bam
for fpath in `ls *.sam`
do
    samtools sort -@ 24 ${fpath} > ${fpath%.sam}.bam
done

for fpath in `ls *.bam`
do
    picard MarkDuplicates \
        -I ${fpath%.bam}.bam \
        -O ${fpath%.bam}.markdup.bam \
        -M ${fpath%.bam}.markdup.metrics.txt \
        --CREATE_INDEX
done

cd ${PROJECT_PATH}/bam
for fpath in `ls *.markdup.bam`
do
    fname=${fpath%.markdup.bam}

    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
         HaplotypeCaller \
         -R $Tgindex --emit-ref-confidence GVCF \
         -I ${fname}.markdup.bam \
         -O ${fname}.g.vcf --sample-ploidy 1
done

mv *.g.vcf ../bqsr/
mv *.g.vcf.idx ../bqsr/

cd ${PROJECT_PATH}/bqsr
gvcf_files=""
for gvcf_file in `ls *.g.vcf`
do
    gvcf_files=${gvcf_files}"-V ${gvcf_file} "
done

picard VcfToIntervalList I=ME49.g.vcf O=sample.interval_list

gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"  \
    GenomicsDBImport -R ${Tgindex} ${gvcf_files} -L sample.interval_list \
                      --genomicsdb-workspace-path gvcfs_db
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs -R ${Tgindex} -V gendb://gvcfs_db -O merged.vcf

cd ${PROJECT_PATH}/bqsr
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"\
    SelectVariants -R ${Tgindex} -V merged.vcf --select-type-to-include SNP -O merged_snps.vcf

gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
   VariantFiltration -R ${Tgindex} -V merged_snps.vcf -O merged_snps_filtered.vcf \
       -filter "MQ < 20.0" --filter-name "MQ20"     \
       -filter "QD < 2.0" --filter-name "QD2"       \
       -filter "FS > 60.0" --filter-name "FS60"     \
       -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
       -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

cd ${PROJECT_PATH}/bqsr
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    SelectVariants -R ${Tgindex} -V merged.vcf --select-type-to-include INDEL -O merged_indels.vcf

gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    VariantFiltration -R ${Tgindex} -V merged_indels.vcf -O merged_indels_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2"       \
    -filter "FS > 200.0" --filter-name "FS200"   \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

cd ${PROJECT_PATH}/bam
for fpath in `ls *.markdup.bam`
do
    fname=${fpath%.markdup.bam}
    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        BaseRecalibrator -R ${Tgindex} -I ${fname}.markdup.bam \
            --known-sites ${PROJECT_PATH}/bqsr/merged_snps_filtered.vcf \
            --known-sites ${PROJECT_PATH}/bqsr/merged_indels_filtered.vcf \
            -O ${fname}_recal_data.table

    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    ApplyBQSR -R ${Tgindex} -I ${fname}.markdup.bam -bqsr ${fname}_recal_data.table -O ${fname}_bqsr.bam
done

cd ${PROJECT_PATH}/bam
for fpath in `ls *_bqsr.bam`
do
    fname=${fpath%_bqsr.bam}
    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    BaseRecalibrator -R ${Tgindex} -I ${fname}_bqsr.bam \
        --known-sites ${PROJECT_PATH}/bqsr/merged_snps_filtered.vcf \
        --known-sites ${PROJECT_PATH}/bqsr/merged_indels_filtered.vcf \
        -O ${fname}_recal_data.table.2

    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        AnalyzeCovariates -before ${fname}_recal_data.table -after ${fname}_recal_data.table.2 \
        -plots ${fname}_recalibration_plots.pdf
done

cd ${PROJECT_PATH}/bam
for fpath in `ls *_bqsr.bam`
do
    fname=${fpath%_bqsr.bam}
    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    HaplotypeCaller -R ${Tgindex} --emit-ref-confidence GVCF \
                        -I ${fname}_bqsr.bam -O ${fname}.g.vcf --sample-ploidy 1
done

mv *.g.vcf ../vcf/
mv *.g.vcf.idx ../vcf/

cd ${PROJECT_PATH}/vcf
gvcf_files=""
for gvcf_file in `ls *.g.vcf`
do
    gvcf_files=${gvcf_files}"-V ${gvcf_file} "
done

picard VcfToIntervalList I=ME49.g.vcf O=sample.interval_list
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport -R ${Tgindex} ${gvcf_files} -L sample.interval_list \
    --genomicsdb-workspace-path gvcfs_db
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs -R ${Tgindex} -V gendb://gvcfs_db -O merged.vcf

cd ${PROJECT_PATH}/vcf
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    SelectVariants -R ${Tgindex} -V merged.vcf --select-type-to-include SNP -O merged_snps.vcf

gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    VariantFiltration -R ${Tgindex} -V merged_snps.vcf -O merged_snps_filtered.vcf \
        -filter "MQ < 20.0" --filter-name "MQ20"     \
        -filter "QD < 2.0" --filter-name "QD2"       \
        -filter "FS > 60.0" --filter-name "FS60"     \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

cd ${PROJECT_PATH}/vcf
gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
SelectVariants -R ${Tgindex} -V merged.vcf --select-type-to-include INDEL -O merged_indels.vcf

gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
VariantFiltration -R ${Tgindex} -V merged_indels.vcf -O merged_indels_filtered.vcf \
                       -filter "QD < 2.0" --filter-name "QD2"       \
                       -filter "FS > 200.0" --filter-name "FS200"   \
                       -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

cd ${PROJECT_PATH}/vcf
items=3
for x in "${items[@]}";do
    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        VariantFiltration -R ${Tgindex} -V merged_snps.vcf -O merged_DP${x}_filtered.vcf \
        --genotype-filter-expression "DP<=${x}" \
        --genotype-filter-name "DP${x}" \
        -filter "MQ < 20.0" --filter-name "MQ20"     \
        -filter "QD < 2.0" --filter-name "QD2"       \
        -filter "FS > 60.0" --filter-name "FS60"     \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -filter "AC < 2.0" --filter-name "AC1"


    gatk --java-options "-Xmx730G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        SelectVariants -R ${Tgindex}  \
        -V merged_DP${x}_filtered.vcf \
        -O merged_snps_DP${x}_filtered.vcf \
        --set-filtered-gt-to-nocall \
        --restrict-alleles-to BIALLELIC \
        --select-type-to-include SNP \
        --exclude-filtered  \
        --selectExpressions "AC > 1"
done
