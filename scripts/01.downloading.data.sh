#!/bin/bash



# Download the  HPRC_PLUS data from the S3 repository
# release desciption https://humanpangenome.org/hprc-data-release-2/
# sample data description https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/README.md

# 1. The VCF (.dip.vcf.gz) shows the variants found in the descendant that differ from the parents, each line containing specific positions (SNPs and Indels).
# Example
# ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
# #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT syndip
# chr1 89177 . A G 30 . . GT:AD 1|1:0,2

# 2. The BED (.dip.bed) shows regions and confidence (regions of perfect assembly compared to the reference) where the VCF variants are contained. Each line is a start and end interval.
# Example
# chr1 89153 180104 # the sequencing is correct in this region.



# 1. HG002, HG002.f1_assembly_v2.dip.bed, HG002.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG002.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG002.f1_assembly_v2.dip.vcf.gz

# 2. HG005/, HG005.f1_assembly_v2.dip.bed, HG005.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG005.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG005.f1_assembly_v2.dip.vcf.gz

# 3. HG00733/, HG00733.f1_assembly_v2.dip.bed, HG00733.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG00733.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG00733.f1_assembly_v2.dip.vcf.gz

# 4. HG01109/, HG01109.f1_assembly_v2.dip.bed, HG01109.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG01109/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG01109.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG01109/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG01109.f1_assembly_v2.dip.vcf.gz

# 5. HG01243/, HG01243.f1_assembly_v2.dip.bed, HG01243.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG01243/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG01243.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG01243/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG01243.f1_assembly_v2.dip.vcf.gz

# 6. HG02055/, HG02055.f1_assembly_v2.dip.bed, HG02055.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02055/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02055.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02055/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02055.f1_assembly_v2.dip.vcf.gz

# 7. HG02080/, HG02080.f1_assembly_v2.dip.bed, HG02080.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02080/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02080.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02080/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02080.f1_assembly_v2.dip.vcf.gz

# 8. HG02109/, HG02109.f1_assembly_v2.dip.bed, HG02109.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02109/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02109.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02109/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02109.f1_assembly_v2.dip.vcf.gz

# 9. HG02145/, HG02145.f1_assembly_v2.dip.bed, HG02145.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02145/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02145.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02145/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02145.f1_assembly_v2.dip.vcf.gz

# 10. HG02723/, HG02723.f1_assembly_v2.dip.bed, HG02723.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02723/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02723.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02723/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02723.f1_assembly_v2.dip.vcf.gz

# 11. HG02818/, HG02818.f1_assembly_v2.dip.bed, HG02818.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02818/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02818.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG02818/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG02818.f1_assembly_v2.dip.vcf.gz

# 12. HG03098/, HG03098.f1_assembly_v2.dip.bed, HG03098.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG03098/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG03098.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG03098/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG03098.f1_assembly_v2.dip.vcf.gz

# 13. HG03486/, HG03486.f1_assembly_v2.dip.bed, HG03486.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG03486/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG03486.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG03486/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG03486.f1_assembly_v2.dip.vcf.gz

# 14. HG03492/, HG03492.f1_assembly_v2.dip.bed, HG03492.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG03492/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG03492.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG03492/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/HG03492.f1_assembly_v2.dip.vcf.gz

# 15. NA18906/, NA18906.f1_assembly_v2.dip.bed, NA18906.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA18906/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA18906.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA18906/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA18906.f1_assembly_v2.dip.vcf.gz

# 16. NA19240/, NA19240.f1_assembly_v2.dip.bed, NA19240.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA19240/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA19240.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA19240/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA19240.f1_assembly_v2.dip.vcf.gz

# 17. NA20129/, NA20129.f1_assembly_v2.dip.bed, NA20129.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA20129/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA20129.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA20129/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA20129.f1_assembly_v2.dip.vcf.gz

# 18. NA21309/, NA21309.f1_assembly_v2.dip.bed, NA21309.f1_assembly_v2.dip.vcf.gz,.
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA21309/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA21309.f1_assembly_v2.dip.bed
wget -P /home/DATA/HPRC_PLUS/  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA21309/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/NA21309.f1_assembly_v2.dip.vcf.gz

# 19. HG01442/, não tem vcf faseado,.
# 20. HG02970/, não tem vcf faseado,.
# 21. HG06807/, não tem vcf faseado,.
# 22. NA19030/, não tem vcf faseado,.
# 23. NA20300/, não tem vcf faseado,.



# end