#!/bin/bash



# MHC region from hla-mapper
# zgrep -v "^#" /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz | cut -f1-10 | head -n 1
# chr6	29700000	rs76597150;rs1773811957	TTGG	CTGG,T	.	.	AC=72;AF=0.00672646,0.000280269;CM=0,0	GT	0|0
# zgrep -v "^#" /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz | cut -f1-10 | tail -n 1
# chr6	33149972	rs546881640	T	G	.	.	AC=1;AF=9.3423e-05;CM=3.44997	GT	0|0
#    HLA coordinates = chr6:29700000-33149972

# MHC region from Kulski 2022
#    Human leukocyte antigen super-locus: nexus of genomic supergenes, SNPs, indels, transcripts, and haplotypes.
#    HLA genomic region that corresponds to the genomic coordinates of 29602228 (GABBR1) to 33410226 (KIFC1).



# paths
path_int="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/"



# sample names
n1="HG002"
n2="HG005"
n3="HG00733"
n4="HG01109"
n5="HG01243"
n6="HG02055"
n7="HG02080"
n8="HG02109"
n9="HG02145"
n10="HG02723"
n11="HG02818"
n12="HG03098"
n13="HG03486"
n14="HG03492"
n15="NA18906"
n16="NA19240"
n17="NA20129"
n18="NA21309"
names=(${n1} ${n2} ${n3} ${n4} ${n5} ${n6} ${n7} ${n8} ${n9} ${n10} ${n11} ${n12} ${n13} ${n14} ${n15} ${n16} ${n17} ${n18})

# checking name

bcftools query -l "/dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz" > \
    "${path_int}/mapper.samples.in.hprc.txt"

for name in "${names[@]}"; do
    
    #grep "${name}"

done



# HPRC vcfs

# create folder to VCFs with the MHC region, one by sample (HPRC)
mkdir -p "${path_int}/mhc.hprc"

for name in "${names[@]}"; do
    
    # selecting chr6 from the phased VCFs
    bcftools view \
        -r chr6:29700000-33149972 \
        "/home/DATA/HPRC_PLUS/${name}.f1_assembly_v2.dip.vcf.gz" \
        -Oz \
        -o "${path_int}/${name}.dip.vcf.gz"

done



# HLA-MAPPER vcfs

# create folder to VCFs with the MHC region, one by sample (HPRC)
mkdir -p "${path_int}/mhc.hlamapper"

# selecting MHC from the populational (SABE_1KGEN_HGDP) phased VCFs (hla-mapper), one VCF by sample.

for name in "${names[@]}"; do
    
    # checking name
    bcftools query \
        -l "/dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz" | \
    grep "${name}"

    # selecting chr6 from the phased VCFs
    bcftools view \
        -s "${name}" \
        -r chr6:29700000-33149972 \
        "/dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz" \
        -Oz \
        -o "${path_int}/${name}.mapper.vcf.gz"

done




# end