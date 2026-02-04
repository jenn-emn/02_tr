#!/bin/bash

# ===============================================
# Script to compare phased haplotypes in the MHC region between
# the hla-mapper (populational)
# and the HPRC (individual)
# using Whatshap compare
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
    echo "path_out <path to the output folder where the extracted MHC VCFs will be saved>"
    echo "name <string to identify the job, e.g., 'trios' or 'trios_hla-mapper'>"
    echo ""
    echo "Example:"
    echo "$(basename "$0") --out /path/to/output --name trios_analysis"
    echo "03.whatshap.sh --out /home/jennifer/02_datas/04_data_processing_trios/01_intermediate --name test_hla-mapper"
    echo ""
}

# Get the absolute path of the directory where the script is located.
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# Check all arguments (flags and their values)
while [[ "$#" -gt 0 ]]; do
    case "$1" in
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
if [[ -z "$path_out" || -z "$name_job" ]]; then
    Usage
    exit 1
fi

# Validate directories and files
if [ ! -d "$path_out" ]; then
    echo "The directory path ('$path_out') is not an existing directory."
    exit 1
fi



# ===============================================
# Starting the analysis
# ===============================================

# intermediate paths
path_int="${path_out}/${name_job}"
mkdir -p "${path_int}"

path_hprc="${path_int}/hprc.mhc"
path_hlamapper="${path_int}/hlamapper.mhc"

# SAMPLES
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

log="${path_int}/whatshap.all.samples.log"
> "${log}"

for name in "${names[@]}"; do

    echo -e "\nSAMPLE: ${name} ------------------------------------------------" >> "${log}"

    truthvcf="${path_hprc}/${name}.dip.reheaded.vcf.gz"
    phasedvcf="${path_hlamapper}/${name}.mapper.vcf.gz"
    tsvpairwise="${path_int}/${name}.tsv.pairwse.txt"
    echo -e "${name}\tchr6\ttruth\twhatshap\t${name}.dip.vcf.gz\t${name}.mapper.vcf.gz" > "${tsvpairwise}"

    # comparison hla.mapper vs trios
    whatshap compare \
        --names truth,phased \
        --tsv-pairwise "${tsvpairwise}" \
        "${truthvcf}" \
        "${phasedvcf}" &>> "${log}"

done



# end