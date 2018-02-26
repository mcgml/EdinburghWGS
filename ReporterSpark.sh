#!/bin/bash
#PBS -l walltime=999:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Report cohort VCF to text file
#Author: Matt Lyon, Harriet Jackson All Wales Medical Genetics Lab
#Mode: BY_COHORT

#load run & pipeline variables
. variables

#report variants to text
/share/apps/jre-distros/jre1.8.0_131/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx48g -jar /data/diagnostics/apps/VariantReporterSpark/VariantReporterSpark-1.3.2/VariantReporterSpark.jar \
-V "$seqId"_recalibrated_variants_genotype_filtered_vep_gnomad.vcf \
-P "$seqId"_pedigree.ped \
-T 8 \
-N