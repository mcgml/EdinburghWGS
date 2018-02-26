#compare test NA12878 with gold standard
/share/apps/rtg-distros/rtg-tools-3.8.4/rtg vcfeval \
-b /data/db/human/giab/3.3.2-hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
-c "$2" \
-e /data/db/human/giab/3.3.2-hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
-o "RTG" \
-t /data/db/human/gatk/v0/hg38/Homo_sapiens_assembly38-SDF \
--sample "HG001,NA12878"