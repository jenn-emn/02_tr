#!/bin/bash



# OBSERVATION: there are missing genotypes in the original HRPC vcfs

# zgrep -w "31427785" /home/DATA/HRPC_PLUS/HG01109.f1_assembly_v2.dip.vcf.gz
# chr6	31427785	.	C	CTATATATATATTCTA	30	GAP1	.	GT:AD	.|1:0,1

# zgrep -w "31427785" /home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/HG01109.mapper.vcf.gz
# chr6	31427785	.	C	CTATATATATATTCTA	30	GAP1	.	GT:AD	.|1:0,1

# zgrep -w "31427785" /home/jennifer/02_datas/04_data_processing_trios/01_intermediate/graphic/isecHG01109/HG01109.hrpc.idcomp.vcf.gz | grep -e "\.|"
# chr6	31427785	chr6:31427785:C:CTATATATATATTCTA	C	CTATATATATATTCTA	30	GAP1	.	GT:AD	.|1:0,1




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

    hrpcvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hrpc.mhc/${name}.dip.reheaded.vcf.gz"
    hlamvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/${name}.mapper.vcf.gz"
    
    pathindiv="${pathgraph}/${name}"
    mkdir -p "${pathindiv}"


    
    # HRPC

    bcftools \
        query -f '%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT]' \
        "${hrpcvcf}" > \
        "${pathindiv}/${name}.hrpc.txt"
    
    hrpccomp="${pathindiv}/${name}.hrpc.txt"


    
    # HLA-mapper
    
    bcftools \
        query -f '%CHROM:%POS:%REF:%ALT\t%GT' \
        "${hlamvcf}" \
        -o "${pathindiv}/${name}.hlamapper.txt"
    
    hlamcomp="${pathindiv}/${name}.hlamapper.txt"


    
    # merge by ID.comp
    awk '
        BEGIN {
            OFS="\t"
            print "idcomp", "hrpc", "hlam", "status"
        }

        NR == FNR { hrpcidc[$1]=$2 ; next ; }

        {
            if ($1 in hrpcidc) {
                print $1, hrpcidc[$1], $2, (hrpcidc[$1]==$2 ? "MATCH" : "DIFF")
            }
        }' \
        "${hrpccomp}" \
        "${hlamcomp}" > \
        "${pathindiv}/${name}.hrpc.hlamapper.tsv"

    path_hrpc_hla="${pathindiv}/${name}.hrpc.hlamapper.tsv"
    
    
    # WHATSHAP just analyse common heterozygous variants

    echo -e "\nSAMPLE: ${name} ------------------------------------------------" >> "${log}"
    
    n_total_var=$(cat "${path_hrpc_hla}" | wc -l)
    echo -e "- Number of intersected variants (including missing genotypes): ${n_total_var} " >> "${log}"


    # OBS 1
    # the hrpc haplotypes have missing genotypes
    # we need to exclude them

    n_missing_var_truth=$(cut -f2 "${path_hrpc_hla}" | grep -e "\." | wc -l)
    echo -e "- Number of missing variants in HRPC: ${n_missing_var_truth} " >> "${log}"

    n_missing_var_phased=$(cut -f2 "${path_hrpc_hla}" | grep -e "\." | wc -l)
    echo -e "- Number of missing variants in HLA-mapper: ${n_missing_var_phased} " >> "${log}"

    # cleaning    
    grep -v -e "\." "${path_hrpc_hla}" > "${pathindiv}/${name}.hrpc.hlamapper.clean.tsv"
    path_hrpc_hla_clean="${pathindiv}/${name}.hrpc.hlamapper.clean.tsv"

    n_total_var=$(cat "${path_hrpc_hla_clean}" | wc -l)
    echo -e "- Number of intersected variants (without missing genotypes): ${n_total_var} " >> "${log}"


    # OBS 2
    # the unmatched ("DIFF") phased genotypes can not be excessive
    # if the proportion of "DIFF" > 5%
    # then we need to pair the first genotypes of HRPC with the second genotypes of the HLA-mapper
    n_diff=$(cut -f4 "${path_hrpc_hla_clean}" | grep "DIFF" | wc -l)


    # If the mismatch rate is greather than 0.05 (5%) the variable is equal 1
    is_high_diff=$(awk -v diff="$n_diff" -v total="$n_total_var" 'BEGIN { print ( (diff/total) > 0.05 ? 1 : 0 ) }')

    if [ "$is_high_diff" -eq 1 ]; then

        echo "- High mismatch rate: ($n_diff / $n_total_var)" >> "${log}"

        path_swapped="${pathindiv}/${name}.hrpc.hlamapper.clean.swapped.tsv"

        awk '
            BEGIN{
                OFS="\t"
                print "idcomp", "hrpc", "hlam", "status", "hlam_swapped", "status_swapped"
            }
            NR>1 {
                # Split HLAM ($3) by the pipe; "0|1": a[1]=0 e a[2]=1
                split($3, a, "|")
                swapped_gt = a[2]"|"a[1]
                
                # Comparing HRPC ($2) with HLAM_swapped
                if ($2 == swapped_gt) {
                    print $1, $2, $3, $4, swapped_gt, "MATCH"
                } else {
                    print $1, $2, $3, $4, swapped_gt, "DIFF"
                }
            }' "${path_hrpc_hla_clean}" > "${path_swapped}"
        
        n_diff_swapped=$(cut -f6 "${path_swapped}" | grep -c "DIFF")
        echo "- Number of unmatched genotypes with SWAPPED alleles: ${n_diff_swapped}" >> "${log}"

        n_homozigous_ref=$(cut -f5,6 "${path_swapped}" | grep "DIFF" | grep -F -c "0|0")
        echo "- Error due to reference homozigous '0\|0': ${n_homozigous_ref}" >> "${log}"

        n_heterozigous=$(cut -f5,6 "${path_swapped}" | grep "DIFF" | grep -F -v -e "0|0" -v -e "1|1" | wc -l)
        echo "- Error due to switch, i.e., '0\|1' instead '1\|0': ${n_heterozigous}" >> "${log}"

    else

        echo "- Low mismatch rate: ($n_diff / $n_total_var)" >> "${log}"
        echo "- Number of unmatched genotypes: ${n_diff}" >> "${log}"

        n_homozigous_ref=$(cut -f5 "${path_swapped}" | grep -c -e "0|0")
        echo "- Error due to reference homozigous '0\|0': ${n_homozigous_ref}" >> "${log}"

        n_heterozigous=$(cut -f5,6 "${path_swapped}" | grep "DIFF" | grep -F -v -e "0|0" -v -e "1|1" | wc -l)
        echo "- Error due to switch, i.e., '0\|1' instead '1\|0': ${n_heterozigous}" >> "${log}"

    fi


    


    #grep -v "1|1" -v -e "\." "${pathindiv}/${name}.hrpc.hlamapper.tsv" > "${pathindiv}/${name}.hrpc.hlamapper.heterozigous.tsv"
    #inte_count=$(cat "${pathindiv}/${name}.hrpc.hlamapper.heterozigous.tsv" | wc -l)
    #diff_count=$(grep -c "DIFF" "${pathindiv}/${name}.hrpc.hlamapper.heterozigous.tsv")
    #echo -e "- sample: ${name} ; intersection variants: ${inte_count} ; diff count: ${diff_count} ; " >> "${log}"

done



# end