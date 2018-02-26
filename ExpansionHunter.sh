#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
set -euo pipefail
#PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
#cd $PBS_O_WORKDIR

#set variables
sampleId=$(basename $PWD)
sex=$(grep $(basename $PWD) ../*ped | awk '$6==2 {print $5}' | sed 's/2/female/g' |sed 's/1/male/g')

#print variables
echo "$sampleId $sex"

#determine repeat sizes
/share/apps/expansionhunter-distros/ExpansionHunter-v2.5.5-linux_x86_64/bin/ExpansionHunter \
--bam "$sampleId".bam \
--ref-fasta /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
--repeat-specs /share/apps/expansionhunter-distros/ExpansionHunter-v2.5.5-linux_x86_64/data/repeat-specs/hg38 \
--vcf "$sampleId".vcf \
--json "$sampleId".json \
--log "$sampleId".log \
--sex "$sex"
