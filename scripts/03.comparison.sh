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

#conda activate whatshap_env

log="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/all.samples.log"
> "${log}"

for name in "${names[@]}"; do

    echo -e "\nSAMPLE: ${name} ------------------------------------------------" >> "${log}"

    tsvpairwise="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/${name}.tsv.pairwse.txt"
    truthvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hprc.mhc/${name}.dip.reheaded.vcf.gz"
    phasedvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/${name}.mapper.vcf.gz"

    echo -e "${name}\tchr6\ttruth\twhatshap\t${name}.dip.vcf.gz\t${name}.mapper.vcf.gz" > "${tsvpairwise}"

    # comparison hla.mapper vs trios
    whatshap compare \
        --names truth,phased \
        --tsv-pairwise "${tsvpairwise}" \
        "${truthvcf}" \
        "${phasedvcf}" &>> "${log}"

done



# end