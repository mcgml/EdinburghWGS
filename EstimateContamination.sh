#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

set -euo pipefail
file=$(ls *bam)

#Calculate dna contamination: sample-to-sample contamination
/share/apps/verifyBamID-distros/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID \
--vcf ../1kg_highconfidence_autosomal_monoallelic_snps.vcf \
--bam "$file" \
--out "$file"_Contamination \
--verbose \
--ignoreRG \
--chip-none \
--minMapQ 20 \
--maxDepth 65
