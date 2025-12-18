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



# observations
echo -e "-------------------------- OBSERVATIONS --------------------------" >> "${log}"
echo -e "- Truth genotypes: HRPC, column 2." >> "${log}"
echo -e "- Phased genotypes: HLA-mapper, column 3 (and column 5 in swapped)." >> "${log}"
echo -e "- Unmatched tags in relation with the swapped HLA-mapper are in the column 4 (and column 6 in swapped)." >> "${log}"
echo -e "- Phasing Error: When the first allele of HRPC is equals to the second allele of the swapped HLA-mapped allele, and 2º hrpc == 1º hlam." >> "${log}"
echo -e "- Genotyping Error: When the phase is not inverted between HRPC and HLA-mapper, it is an error because is different by state." >> "${log}"
echo -e "                 - When the phase (HLA-mapper) is homozigous for one of the truth alleles or for a third extra allele." >> "${log}"
echo -e "                 - when the phase (HLA-mapper) is heterozygous and contains a different allele (by state) from the two truth alleles." >> "${log}"
echo -e "- Obs.: The allele 2 (2) represented by the asterisk in the ALT column of the HRPC vcf indicates a spanning deletion." >> "${log}"
echo -e "                 - It is when the physical position does not exist in this haplotype because it was removed by a larger structural deletion in this region.\n" >> "${log}"



# iterating between samples

for name in "${names[@]}"; do

    echo -e "SAMPLE: ${name} ------------------------------------------------"

    hrpcvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hrpc.mhc/${name}.dip.reheaded.vcf.gz"
    hlamvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/${name}.mapper.vcf.gz"
    
    pathindiv="${pathgraph}/${name}"
    mkdir -p "${pathindiv}"


    
    # HRPC

    bcftools \
        query -f '%CHROM:%POS:%REF:%ALT[\t%GT]' \
        "${hrpcvcf}" > \
        "${pathindiv}/${name}.hrpc.txt"
    
    hrpccomp="${pathindiv}/${name}.hrpc.txt"


    
    # HLA-mapper
    
    bcftools \
        query -f '%CHROM:%POS:%REF:%ALT[\t%GT]' \
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
    
    

    # Checking numbers
    echo -e "\nSAMPLE: ${name} ------------------------------------------------" >> "${log}"
    
    # Intersection by compoused ID variant in the "chr:pos:ref:alt" format
    # counting intersected variants
    n_total_var=$(cat "${path_hrpc_hla}" | wc -l)
    echo -e "- Number of intersected variants (including missing genotypes): ${n_total_var} " >> "${log}"


    
    # Before cleaning
    # The hrpc haplotypes have missing genotypes, we need to exclude them

    # counting hrpc
    n_missing_var_truth=$(cut -f2 "${path_hrpc_hla}" | grep -e "\." | wc -l)
    echo -e "- Number of missing variants in HRPC: ${n_missing_var_truth} " >> "${log}"
    # counting  hla-mapper
    n_missing_var_phased=$(cut -f3 "${path_hrpc_hla}" | grep -e "\." | wc -l)
    echo -e "- Number of missing variants in HLA-mapper: ${n_missing_var_phased} " >> "${log}"

    # cleaning
    path_hrpc_hla_clean="${pathindiv}/${name}.hrpc.hlamapper.clean.tsv"
    grep -v -e "\." "${path_hrpc_hla}" > "${path_hrpc_hla_clean}"
    
    
    
    # After cleaning

    # counting intersected variants
    n_total_var=$(cat "${path_hrpc_hla_clean}" | wc -l)
    echo -e "- Number of intersected variants (without missing genotypes): ${n_total_var} " >> "${log}"


    # Checking unmatched genotypes (tag = "DIFF")
    # If the proportion of "DIFF" > 3%, we need to invert the genotypes of the HLA-mapper
    
    # counting unmatched genotypes
    n_diff=$(cut -f4 "${path_hrpc_hla_clean}" | grep "DIFF" | wc -l)

    # calculating mismatch rate
    is_high_diff=$(awk -v diff="$n_diff" -v total="$n_total_var" 'BEGIN { print ( (diff/total) > 0.3 ? 1 : 0 ) }')
    echo -e "- Mismatch rate: ${is_high_diff}" >> "${log}"



    # Run taging

    if [ "$is_high_diff" -eq 1 ]; then

        # switch error rate before swap
        mismatch_rate=$((100 * "${n_diff}" / "${n_total_var}"))
        echo "- High mismatch rate (before swap): ${mismatch_rate}" >> "${log}"

        # creating a unswapped column
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
        
        # calculate the new Switch Eror Rate (and the subtype error)
        awk '
            $6 == "DIFF" {

                # split hrpc ($2) and swapped hla-mapper (%5)
                split($2, hrpc, "|")
                split($5, hlam, "|")
                
                # -- Phasing Error
                # The allele 1 of hrpc have to be equals to allele 2 of hlam, and the 2 of hrpc == 1 of hlam
                # HRPC is 0|1, 1|0, 0|2, 2|0
                # and HLa-mapper is 1|0, 0|1, 2|0, 0|2, respectly.

                if (hrpc[1] == hlam[2] && hrpc[2] == hlam[1]) {
                    phase_err++
                }
                
                # -- Genotyping Error
                # The alleles are not inverted, they are differents by status
                else {
                    geno_err++

                    # -- Genotyping Error: discrepant homozigous in the swapped HLA-mapper
                    #    when the swapped phase is homozigous for one of the truth alleles or for a third extra allele.
                    if (hlam[1] == hlam[2]) {
                        geno_hom_qry++
                    }

                    # -- Genotyping Error: discrepant heterozigous in the swapped HLA-mapper
                    #    when the swapped phase is heterozygous and contains a different allele (by state) from the two truth alleles
                    else {
                        geno_het_qry++
                    }
                }
            }

            END {
                
                total_err = geno_err + phase_err

                print "---------- AFTER SWAP ----------"
                print "- Total variants: " NR
                print "- Total unmatched variants (DIFF): " total_err
                print "- Switch (Phasing) Errors: " phase_err
                print "- Genotyping Errors: " geno_err
                print "- Genotyping Error by homozigous in phased: " geno_hom_qry
                print "- Genotyping Error by heterozigous in phased: " geno_het_qry
                print "- Total Error Rate: " (total_err * 100 / NR) "%"
                print "- Genotyping Error Rate: " (geno_err * 100 / NR) "%"
                print "- Switch (Phasing) Error Rate: " (phase_err * 100 / NR) "%"
            }
            ' "${path_swapped}" &>> "${log}"

    else

        # calculate the new Switch Eror Rate (and the subtype error)
        awk '
            $4 == "DIFF" {

                # split hrpc ($2) and swapped hla-mapper ($3)
                split($2, hrpc, "|")
                split($3, hlam, "|")
                
                # -- Phasing Error
                # The alleles in HLA-mapper are inverted
                if (hrpc[1] == hlam[2] && hrpc[2] == hlam[1]) {
                    phase_err++
                }
                
                # -- Genotyping Error
                # The alleles are not inverted, they are differents by status
                else {
                    geno_err++

                    # -- Genotyping Error: discrepant homozigous in the swapped HLA-mapper
                    #    when the swapped phase is homozigous for one of the truth alleles or for a third extra allele.
                    if (hlam[1] == hlam[2]) {
                        geno_hom_qry++
                    }

                    # -- Genotyping Error: discrepant heterozigous in the swapped HLA-mapper
                    #    when the swapped phase is heterozygous and contains a different allele (by state) from the two truth alleles
                    else {
                        geno_het_qry++
                    }
                }
            }

            END {
                
                total_err = geno_err + phase_err

                print "---------- WHITOUT SWAP ----------"
                print "- Total variants: " NR
                print "- Total unmatched variants (DIFF): " total_err
                print "- Switch (Phasing) Errors: " phase_err
                print "- Genotyping Errors: " geno_err
                print "- Genotyping Error by homozigous in phased: " geno_hom_qry
                print "- Genotyping Error by heterozigous in phased: " geno_het_qry
                print "- Total Error Rate: " (total_err * 100 / NR) "%"
                print "- Genotyping Error Rate: " (geno_err * 100 / NR) "%"
                print "- Switch (Phasing) Error Rate: " (phase_err * 100 / NR) "%"
            }
            ' "${path_hrpc_hla_clean}" &>> "${log}"

    fi

done



# end