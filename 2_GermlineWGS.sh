#!/bin/bash
#PBS -l walltime=999:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline WGS Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, Harriet Jackson All Wales Medical Genetics Lab
#Mode: BY_COHORT

# Script 2 runs in panel folder, requires final Bams, gVCFs and a PED file
# Variant filtering assumes non-related samples. If familiy structures are known they MUST be provided in the PED file

#load run & pipeline variables
. variables

annotateVCF(){
    #annotate VCF
    perl /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl \
    --verbose \
    --no_progress \
    --everything \
    --fork 12 \
    --species homo_sapiens \
    --assembly GRCh38 \
    --input_file "$1" \
    --format vcf \
    --output_file "$2" \
    --force_overwrite \
    --no_stats \
    --cache \
    --dir /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --fasta /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --no_intergenic \
    --offline \
    --cache_version 86 \
    --allele_number \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --minimal \
    --refseq
    
    #check VEP has produced annotated VCF
    if [ ! -e "$2" ]; then
        cp "$1" "$2"
    fi

    #index annotated VCF
    /share/apps/igvtools-distros/igvtools_2.3.75/igvtools index "$2"
}

### Joint variant calling and filtering ###

#Joint genotyping
/usr/java/jdk1.7.0_51/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-V GVCFs.list \
-o "$seqId"_variants.vcf \
-ped "$seqId"_pedigree.ped \
-nt 12

#select only canonical chromosomes
/usr/java/jdk1.7.0_51/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX \
-V "$seqId"_variants.vcf \
-o "$seqId"_variants_canonical.vcf

#Build the SNP recalibration model
/usr/java/jdk1.7.0_51/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-input "$seqId"_variants_canonical.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/db/human/gatk/v0/hg38/hapmap_3.3.hg38.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 /data/db/human/gatk/v0/hg38/1000G_omni2.5.hg38.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/db/human/gatk/v0/hg38/1000G_phase1.snps.high_confidence.hg38.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile "$seqId"_SNP.recal \
-tranchesFile "$seqId"_SNP.tranches \
-rscriptFile "$seqId"_SNP_plots.R \
-ped "$seqId"_pedigree.ped \
-nt 12

#Apply the desired level of recalibration to the SNPs in the call set
/usr/java/jdk1.7.0_51/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-input "$seqId"_variants_canonical.vcf \
-mode SNP \
--ts_filter_level 99.5 \
-recalFile "$seqId"_SNP.recal \
-tranchesFile "$seqId"_SNP.tranches \
-o "$seqId"_recalibrated_snps_raw_indels.vcf \
-ped "$seqId"_pedigree.ped \
-nt 12

#Build the Indel recalibration model
/usr/java/jdk1.7.0_51/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-input "$seqId"_recalibrated_snps_raw_indels.vcf \
-resource:mills,known=false,training=true,truth=true,prior=12.0 /data/db/human/gatk/v0/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--maxGaussians 4 \
-recalFile "$seqId"_INDEL.recal \
-tranchesFile "$seqId"_INDEL.tranches \
-rscriptFile "$seqId"_INDEL_plots.R \
-ped "$seqId"_pedigree.ped \
-nt 12

#Apply the desired level of recalibration to the Indels in the call set
/usr/java/jdk1.7.0_51/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-input "$seqId"_recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level 99.5 \
-recalFile "$seqId"_INDEL.recal \
-tranchesFile "$seqId"_INDEL.tranches \
-o "$seqId"_recalibrated_variants.vcf \
-ped "$seqId"_pedigree.ped \
-nt 12

#filter genotypes
/usr/java/jdk1.7.0_51/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.4.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-V SIGNAL_WGS_recalibrated_variants.vcf \
-ped SIGNAL_WGS_pedigree.ped \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
-o SIGNAL_WGS_recalibrated_variants_genotype_filtered.vcf

#annotate with VEP
annotateVCF SIGNAL_WGS_recalibrated_variants_genotype_filtered.vcf SIGNAL_WGS_recalibrated_variants_genotype_filtered_vep.vcf

#add gnomad allele frequencies
/share/apps/jre-distros/jre1.8.0_101/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-V SIGNAL_WGS_recalibrated_variants_genotype_filtered_vep.vcf \
--resource:GNOMAD_2.0.1_Genome_chr1 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.1.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr2 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.2.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr3 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.3.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr4 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.4.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr5 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.5.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr6 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.6.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr7 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.7.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr8 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.8.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr9 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.9.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr10 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.10.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr11 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.11.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr12 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.12.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr13 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.13.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr14 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.14.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr15 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.15.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr16 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.16.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr17 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.17.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr18 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.18.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr19 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.19.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr20 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.20.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr21 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.21.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr22 /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.22.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chrX /data/db/human/gnomad/hg38/gnomad.genomes.r2.0.1.sites.X.hg38.vcf.gz \
--resource:GNOMAD_2.0.1_Exome /data/db/human/gnomad/hg38/gnomad.exomes.r2.0.1.sites.hg38.vcf.gz \
-E GNOMAD_2.0.1_Genome_chr1.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr2.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr3.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr4.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr5.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr6.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr7.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr8.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr9.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr10.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr11.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr12.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr13.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr14.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr15.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr16.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr17.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr18.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr19.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr20.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr21.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr22.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chrX.AF_POPMAX \
-E GNOMAD_2.0.1_Exome.AF_POPMAX \
-ped SIGNAL_WGS_pedigree.ped \
--resourceAlleleConcordance \
-o SIGNAL_WGS_recalibrated_variants_genotype_filtered_vep_gnomad.vcf

#validate final VCF
/share/apps/jre-distros/jre1.8.0_101/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
--reference_window_stop 300 \
-V SIGNAL_WGS_recalibrated_variants_genotype_filtered_vep_gnomad.vcf

#report variants to text
/share/apps/jre-distros/jre1.8.0_131/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx48g -jar /data/diagnostics/apps/VariantReporterSpark/VariantReporterSpark-1.3.2/VariantReporterSpark.jar \
-V SIGNAL_WGS_recalibrated_variants_genotype_filtered_vep_gnomad.vcf \
-P SIGNAL_WGS_pedigree.ped \
-T 8 \
-N

### QC ###

#relatedness test
/share/apps/vcftools-distros/vcftools-0.1.14/build/bin/vcftools \
--relatedness2 \
--out SIGNAL_WGS_recalibrated_variants_gcp_phased_gtfiltered_vep_gnomad_relatedness \
--vcf SIGNAL_WGS_recalibrated_variants_gcp_phased_gtfiltered_vep_gnomad.vcf

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3//picard.jar CollectVariantCallingMetrics \
INPUT=SIGNAL_WGS_recalibrated_variants_gcp_phased_gtfiltered_vep_gnomad.vcf.gz \
OUTPUT=CollectVariantCallingMetrics.txt \
DBSNP=/data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
THREAD_COUNT=12

#compare test NA12878 with gold standard
/share/apps/rtg-distros/rtg-tools-3.8.4/rtg vcfeval \
-b /data/db/human/giab/3.3.2-hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
-c "$2" \
-e /data/db/human/giab/3.3.2-hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
-o "RTG" \
-t /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38-SDF \
--sample "HG001,$1"

#Validate VCF file


