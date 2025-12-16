#!/bin/bash



# SAMPLES

# sample names
#n1="HG002"
#n2="HG005"
n3="HG00733"
n4="HG01109"
n5="HG01243"
n6="HG02055"
n7="HG02080"
#n8="HG02109"
n9="HG02145"
n10="HG02723"
n11="HG02818"
n12="HG03098"
n13="HG03486"
n14="HG03492"
n15="NA18906"
n16="NA19240"
n17="NA20129"
#n18="NA21309"
names=(${n3} ${n4} ${n5} ${n6} ${n7} ${n9} ${n10} ${n11} ${n12} ${n13} ${n14} ${n15} ${n16} ${n17})

mkdir -p "/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/graphic"
pathgraph="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/graphic"

#log="${pathgraph}/graph.log"
#> "${log}"

for name in "${names[@]}"; do

    echo -e "\nSAMPLE: ${name} ------------------------------------------------" >> "${log}"

    hrpcvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hprc.mhc/${name}.dip.reheaded.vcf.gz"
    hlamvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/${name}.mapper.vcf.gz"
    pathisec="${pathgraph}/isec${name}"

    mkdir -p "${pathisec}"


    
    # HRPC

    bcftools \
        query -f "%CHROM:%POS:%REF:%ALT\n" \
        "${hrpcvcf}" \
        -o s"${pathisec}/${name}.hrpc.idcomp.vcf"
    
    hrpcvcfcomp="${pathisec}/${name}.hrpc.idcomp.vcf"
    
    bcftools index "${hrpcvcfcomp}"


    
    # HLA-mapper
    
    bcftools \
        query -f "%CHROM:%POS:%REF:%ALT\n" \
        "${hlamvcf}" \
        -o s"${pathisec}/${name}.hlamapper.idcomp.vcf"
    
    hlamvcfcomp="${pathisec}/${name}.hlamapper.idcomp.vcf"
    
    bcftools index "${hlamvcfcomp}"


    # isec; equivalent expressions -n=2, -n~11, -n+2
    bcftools \
        isec -n+2  \
        "${hrpcvcfcomp}" \
        "${hlamvcfcomp}" \
        -Oz \
        -p ${pathisec}
    
    # &>> "${log}"

done



# end