gunzip CosmicCodingMuts.vcf.gz
gunzip CosmicNonCodingVariants.vcf.gz
cat CosmicCodingMuts.vcf CosmicNonCodingVariants.vcf | gzip -c > cosmic.vcf.gz
rm CosmicCodingMuts.vcf
rm CosmicNonCodingVariants.vcf