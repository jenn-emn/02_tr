#!/bin/bash

# ===============================================
# Script to extract the MHC region from phased VCFs.
# For each sample in the trios list,
# extract the MHC region (chr6:29700000-33149972)
# from phased VCFs created by
# the hla-mapper (populational)
# and the HPRC (individual).
# ===============================================



# ===============================================
# Verifying the input call
# ===============================================

# Initialize variables for the options
path_truth=""
path_estimated_populational=""
path_out=""
name_job=""

Usage() {
    echo "Usage: $(basename "$0")"
    echo "path_truth <path to the folder containing the reference VCFs (one VCF per sample)>"
    echo "path_estimated_populational <path to the estimated population VCF containing phased haplotypes>"
    echo "path_out <path to the output folder where the extracted MHC VCFs will be saved>"
    echo "name <string to identify the job, e.g., 'mvn' or 'hla_mapper'>"
    echo ""
    echo "Example:"
    echo "$(basename "$0") --true /path/to/truth_vcfs --est /path/to/estimated_population.vcf --out /path/to/output --name trios_analysis"
    echo "02.processing.hprc.sh --true /dados/home/DATA/HPRC_PLUS --est /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz --out /home/jennifer/02_datas/04_data_processing_trios/01_intermediate --name test_hla_mapper"
    echo ""
}

# Get the absolute path of the directory where the script is located.
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# Check all arguments (flags and their values)
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --true)
            path_truth="$2"
            shift ; shift
            ;;
        --est)
            path_estimated_populational="$2"
            shift ; shift
            ;;
        --out)
            path_out="$2"
            shift ; shift
            ;;
        --name)
            name_job="$2"
            shift ; shift
            ;;
        *)
            echo "Error: Unknown option or invalid syntax: $1"
            Usage
            exit 1
            ;;
    esac
done

# Validate that all required options are provided
if [[ -z "$path_truth" || -z "$path_estimated_populational" || -z "$path_out" || -z "$name_job" ]]; then
    Usage
    exit 1
fi

# Validate directories and files
if [ ! -d "$path_truth" ]; then
    echo "The directory path ('$path_truth') is not an existing directory."
    exit 1
fi

if [ ! -d "$path_out" ]; then
    echo "The directory path ('$path_out') is not an existing directory."
    exit 1
fi

if [ ! -f "$path_estimated_populational" ]; then
    echo "The file path ('$path_estimated_populational') is not an existing file."
    exit 1
fi



# ===============================================
# Starting the analysis
# ===============================================

# MHC region from hla-mapper
# zgrep -v "^#" /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz | cut -f1-10 | head -n 1
# chr6	29700000	rs76597150;rs1773811957	TTGG	CTGG,T	.	.	AC=72;AF=0.00672646,0.000280269;CM=0,0	GT	0|0
# zgrep -v "^#" /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz | cut -f1-10 | tail -n 1
# chr6	33149972	rs546881640	T	G	.	.	AC=1;AF=9.3423e-05;CM=3.44997	GT	0|0
#    HLA coordinates = chr6:29700000-33149972

# MHC region from Kulski 2022
#    Human leukocyte antigen super-locus: nexus of genomic supergenes, SNPs, indels, transcripts, and haplotypes.
#    HLA genomic region that corresponds to the genomic coordinates of 29602228 (GABBR1) to 33410226 (KIFC1).



# intermediate paths
#path_int="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate"
path_int="${path_out}/${name_job}"
mkdir -p "${path_int}"

path_hprc="${path_int}/hprc.mhc"
path_estimated="${path_int}/${name_job}"

# to check names
samples_list="${path_int}/${name_job}.samples.in.hprc.txt"
# to .log
log="${path_int}/01.log"



# SAMPLES

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

bcftools \
    query -l \
    "${path_estimated_populational}" > \
    "${samples_list}"

for name in "${names[@]}"; do
    if ! grep -Fxq "${name}" "${samples_list}"; then
        echo "'${name}' is missing!" >> "${log}"
    fi
done
# cat /home/jennifer/02_datas/04_data_processing_trios/01_intermediate/01.log
    # 'HG002' is missing!
    # 'HG005' is missing!
    # 'HG02109' is missing!
    # 'NA21309' is missing!



# HPRC vcfs

# create folder to VCFs with the MHC region, one by sample (HPRC)
mkdir -p "${path_hprc}"

for name in "${names[@]}"; do
    
    # selecting chr6 from the phased VCFs
    bcftools view \
        -r chr6:29700000-33149972 \
        "${path_truth}/${name}.f1_assembly_v2.dip.vcf.gz" \
        -Oz \
        -o "${path_hprc}/${name}.dip.vcf.gz"
        #"/home/DATA/HPRC_PLUS/${name}.f1_assembly_v2.dip.vcf.gz" \
    
    # name.txt
    echo "${name}" > "${path_hprc}/${name}.to.reheader.txt"
    
    bcftools reheader \
        "${path_hprc}/${name}.dip.vcf.gz" \
        -s "${path_hprc}/${name}.to.reheader.txt" \
        -o "${path_hprc}/${name}.dip.reheaded.vcf.gz"

    bcftools index "${path_hprc}/${name}.dip.reheaded.vcf.gz"
    if [ -s "${path_hprc}/${name}.dip.reheaded.vcf.gz.csi" ]; then
        rm "${path_hprc}/${name}.dip.vcf.gz"
    fi

done



# ESTIMATED vcfs

# create folder to VCFs with the MHC region, one by sample (HLA-MAPPER)
mkdir -p "${path_estimated}"

# selecting MHC from the populational (SABE_1KGEN_HGDP) phased VCFs (hla-mapper), one VCF by sample.

for name in "${names[@]}"; do

    # selecting chr6 from the phased VCFs
    bcftools view \
        -s "${name}" \
        -r chr6:29700000-33149972 \
        "${path_estimated_populational}" \
        -Oz \
        -o "${path_estimated}/${name}.${name_job}.vcf.gz"
        #"/dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz" \

    bcftools index "${path_estimated}/${name}.${name_job}.vcf.gz"

done



# end