#!/bin/bash

# ===============================================
# Script to:
# 1. extract the MHC region (chr6:29700000-33149972) from phased VCFs: hla-mapper (populational) and the HPRC (individual).
# 2. compare the extracted MHC regions using whatsHap compare.
# 3. compare the extracted MHC regions using genotype state comparison.
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
    echo "Usage of $(basename "$0")"
    echo "  --true <path to the folder containing the reference VCFs (one VCF per sample)>"
    echo "  --est  <path to the estimated population VCF containing phased haplotypes>"
    echo "  --out  <path to the output folder where the extracted MHC VCFs will be saved>"
    echo "  --name <name to identify the job, e.g., 'trios' or 'trios_hla-mapper'>"
    echo ""
    echo "Example:"
    echo "bash $(basename "$0") --true /path/to/truth_vcfs --est /path/to/estimated_population.vcf --out /path/to/output --name trios_analysis"
    echo "bash 00.main.sh --true /dados/home/DATA/HPRC_PLUS --est /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz --out /home/jennifer/02_datas/04_data_processing_trios/01_intermediate --name test_hla-mapper"
    echo ""
}

# Get the absolute path of the directory where the script is located.
path_script=$(dirname "$(readlink -f "$0")")

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

# Validate if the length of any required option is zero
if [[ -z "$path_truth" || -z "$path_estimated_populational" || -z "$path_out" || -z "$name_job" ]]; then
    Usage
    exit 1
fi

# Validate directories and files do not exist
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

# Extract the MHC region from phased VCFs
bash "${path_script}/02.processing.hprc.sh" \
    --true "${path_truth}" \
    --est "${path_estimated_populational}" \
    --out "${path_out}" \
    --name "${name_job}"

# Whatshap compare
bash "${path_script}/03.whatshap.sh" \
    --out "${path_out}" \
    --name "${name_job}"
exit
# switches
bash "${path_script}/04.switches.sh" \
    --out "${path_out}" \
    --name "${name_job}"

# report
cash "${path_script}/05.report.plot.sh" \
    --out "${path_out}" \
    --name "${name_job}"

# bash 00.main.sh --true /dados/home/DATA/HPRC_PLUS --est /dados/home/DATA/HLAcalls_1kgenHGDP_2024/SABE_1KGEN_HGDP/vcf_nay/whatshap/whatshap_bialelico_shapeit_multialelico_EDITADO7.vcf.gz --out /home/jennifer/02_datas/04_data_processing_trios/01_intermediate --name hla_mapper
# end
