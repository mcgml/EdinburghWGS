set -euo pipefail

#Extract 1kg autosomal snps for contamination analysis
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38.fasta \
-T SelectVariants \
--variant /data/db/human/gatk/v0/hg38/1000G_phase1.snps.high_confidence.hg38.vcf \
-o 1kg_highconfidence_autosomal_monoallelic_snps.vcf \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-env \
-ef \
-L chr1 \
-L chr2 \
-L chr3 \
-L chr4 \
-L chr5 \
-L chr6 \
-L chr7 \
-L chr8 \
-L chr9 \
-L chr10 \
-L chr11 \
-L chr12 \
-L chr13 \
-L chr14 \
-L chr15 \
-L chr16 \
-L chr17 \
-L chr18 \
-L chr19 \
-L chr20 \
-L chr21 \
-L chr22
