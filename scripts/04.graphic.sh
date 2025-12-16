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

log="${pathgraph}/graph.log"
> "${log}"

for name in "${names[@]}"; do

    echo -e "SAMPLE: ${name} ------------------------------------------------"

    hrpcvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hprc.mhc/${name}.dip.reheaded.vcf.gz"
    hlamvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/${name}.mapper.vcf.gz"
    pathisec="${pathgraph}/isec${name}"

    mkdir -p "${pathisec}"


    
    # HRPC

    bcftools \
        annotate --set-id '%CHROM:%POS:%REF:%ALT' \
        "${hrpcvcf}" \
        -Oz \
        -o "${pathisec}/${name}.hrpc.idcomp.vcf.gz"
    
    hrpcvcfcomp="${pathisec}/${name}.hrpc.idcomp.vcf.gz"
    
    bcftools index "${hrpcvcfcomp}"


    
    # HLA-mapper
    
    bcftools \
        annotate --set-id '%CHROM:%POS:%REF:%ALT' \
        "${hlamvcf}" \
        -Oz \
        -o "${pathisec}/${name}.hlamapper.idcomp.vcf.gz"
    
    hlamvcfcomp="${pathisec}/${name}.hlamapper.idcomp.vcf.gz"
    
    bcftools index "${hlamvcfcomp}"



    # ISEC
    # equivalent expressions: -n=2, -n~11, -n+2
    bcftools \
        isec -n+2  \
        "${hrpcvcfcomp}" \
        "${hlamvcfcomp}" \
        -Oz \
        -p "${pathisec}"
    


    # Extract ID.comp, GT with the name of the sample
    bcftools \
        query -f "%ID[\t%SAMPLE=%GT]\n" \
        "${pathisec}/0000.vcf.gz" > \
        "${pathisec}/${name}.hrpc.tsv"

    bcftools \
        query -f "%ID[\t%SAMPLE=%GT]\n" \
        "${pathisec}/0001.vcf.gz" > \
        "${pathisec}/${name}.hlamapper.tsv"


    
    # merge by ID.comp
    awk '
        BEGIN {
            OFS="\t"
            print "idcomp", "hrpc.haplo", "hlam.haplo"
        }

        NR == FNR { hrpcidc[$1]=$2 ; next ; }

        {
            if ($1 in hrpcidc) {
                print $1, hrpcidc[$1], $2, (hrpcidc[$1]==$2 ? "MATCH" : "DIFF")
            }
        }' \
        "${pathisec}/${name}.hrpc.tsv" \
        "${pathisec}/${name}.hlamapper.tsv" > \
        "${pathisec}/${name}.hrpc.hlamapper.tsv"
    
    # WHATSHAP: common heterozygous variants
    grep -v "1|1" -v "." "${pathisec}/${name}.hrpc.hlamapper.tsv" > "${pathisec}/${name}.hrpc.hlamapper.heterozigous.tsv"

    inte_count=$(cat "${pathisec}/${name}.hrpc.hlamapper.heterozigous.tsv" | wc -l)
    diff_count=$(grep -c "DIFF" "${pathisec}/${name}.hrpc.hlamapper.heterozigous.tsv")
    echo -e "- sample: ${name} ; intersection variants: ${inte_count} ; diff count: ${diff_count} ; " >> "${log}"

done



# end