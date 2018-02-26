#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

sampleId=$(basename $PWD)

#break point analysis
/share/apps/manta-distros/manta-1.3.0.centos6_x86_64/bin/configManta.py \
--bam "$sampleId".bam \
--referenceFasta /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
--runDir manta

manta/runWorkflow.py \
--quiet \
-m local \
-j 12

#rename sv vcf
mv manta/results/variants/diploidSV.vcf.gz "$sampleId"_sv_filtered.vcf.gz
mv manta/results/variants/diploidSV.vcf.gz.tbi "$sampleId"_sv_filtered.vcf.gz.tbi
