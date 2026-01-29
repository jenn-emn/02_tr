#!/bin/bash



# OBSERVATION: there are missing genotypes in the original HPRC vcfs

# zgrep -w "31427785" /home/DATA/HPRC_PLUS/HG01109.f1_assembly_v2.dip.vcf.gz
# chr6	31427785	.	C	CTATATATATATTCTA	30	GAP1	.	GT:AD	.|1:0,1

# zgrep -w "31427785" /home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/HG01109.mapper.vcf.gz
# chr6	31427785	.	C	CTATATATATATTCTA	30	GAP1	.	GT:AD	.|1:0,1

# zgrep -w "31427785" /home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch/isecHG01109/HG01109.hprc.idcomp.vcf.gz | grep -e "\.|"
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

mkdir -p "/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch"
pathswitch="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch"

log="${pathswitch}/switches.inhouse.log"
> "${log}"



# observations
echo -e "-------------------------- OBSERVATIONS --------------------------" >> "${log}"
echo -e "- Truth genotypes: HPRC, column 2." >> "${log}"
echo -e "                 - Cleaning: exclusion of missingness genotypes and homozigous from HPRC." >> "${log}"
echo -e "- Phased genotypes: HLA-mapper, column 3 (and column 5 in swapped)." >> "${log}"
echo -e "- Tags: 'DIFF' to unmatched phases between HPRC and HLA-mapper phase." >> "${log}"
echo -e "- Phasing Error: When the first allele of HPRC is equals to the second allele of the swapped HLA-mapped allele, and 2º hprc == 1º hlam." >> "${log}"
echo -e "- Genotyping Error: When the phase is not inverted between HPRC and HLA-mapper, it is an error because is different by state." >> "${log}"
echo -e "                 - When the phase (HLA-mapper) is homozigous for one of the truth alleles or for a third extra allele." >> "${log}"
echo -e "                 - when the phase (HLA-mapper) is heterozygous and contains a different allele (by state) from the two truth alleles." >> "${log}"
echo -e "- Obs.: The allele 2 (2) represented by the asterisk in the ALT column of the HPRC vcf indicates a spanning deletion." >> "${log}"
echo -e "                 - It is when the physical position does not exist in this haplotype because it was removed by a larger structural deletion in this region.\n" >> "${log}"



# iterating between samples

for name in "${names[@]}"; do

    echo -e "SAMPLE: ${name} ------------------------------------------------"

    hprcvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hprc.mhc/${name}.dip.reheaded.vcf.gz"
    hlamvcf="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/hlamapper.mhc/${name}.mapper.vcf.gz"
    
    pathindiv="${pathswitch}/${name}"
    mkdir -p "${pathindiv}"


    
    # HPRC

    bcftools \
        query -f '%CHROM:%POS:%REF:%ALT[\t%GT]' \
        "${hprcvcf}" > \
        "${pathindiv}/${name}.hprc.txt"
    
    hprccomp="${pathindiv}/${name}.hprc.txt"


    
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
            print "idcomp", "hprc", "hlam", "status"
        }

        NR == FNR { hprcidc[$1]=$2 ; next ; }

        {
            if ($1 in hprcidc) {
                print $1, hprcidc[$1], $2, (hprcidc[$1]==$2 ? "MATCH" : "DIFF")
            }
        }' \
        "${hprccomp}" \
        "${hlamcomp}" > \
        "${pathindiv}/${name}.hprc.hlamapper.tsv"

    path_hprc_hla="${pathindiv}/${name}.hprc.hlamapper.tsv"
    
    

    # Checking numbers
    echo -e "\nSAMPLE: ${name} ------------------------------------------------" >> "${log}"
    
    # Intersection by compoused ID variant in the "chr:pos:ref:alt" format
    # counting intersected variants
    n_total_var=$(cat "${path_hprc_hla}" | wc -l)
    echo -e "- Number of intersected variants (pre-cleaning): ${n_total_var} " >> "${log}"


    
    # Before cleaning
    # The hprc haplotypes have missing genotypes, we need to exclude them

    # counting missingness in hprc
    n_missing_var_truth=$(cut -f2 "${path_hprc_hla}" | grep -e "\." | wc -l)
    echo -e "    - Number of missing variants in HPRC: ${n_missing_var_truth} " >> "${log}"
    
    # counting homozygous in hprc
    n_homozygous_var_truth=$(cut -f2 "${path_hprc_hla}" | grep -F "1|1" | wc -l)
    echo -e "    - Number of homozygous variants in HPRC: ${n_homozygous_var_truth} " >> "${log}"
    
    # counting  hla-mapper
    n_missing_var_phased=$(cut -f3 "${path_hprc_hla}" | grep -e "\." | wc -l)
    #echo -e "    - Number of missing variants in HLA-mapper: ${n_missing_var_phased} " >> "${log}"



    # Cleaning missing genotypes and heterozygous
    path_hprc_hla_clean="${pathindiv}/${name}.hprc.hlamapper.clean.tsv"
    grep -v -e "\." "${path_hprc_hla}" | \
        grep -v -F "1|1" > \
        "${path_hprc_hla_clean}"
    
    
    
    # After cleaning

    # counting intersected variants
    n_total_var=$(( $(wc -l < "${path_hrpc_hla_clean}") - 1 ))
    echo -e "- Number of common heterozygous variants (pos-cleaning): ${n_total_var} " >> "${log}"


    # Checking unmatched genotypes (tag = "DIFF")
    # If the proportion of "DIFF" > 3%, we need to invert the genotypes of the HLA-mapper
    
    # counting unmatched genotypes
    n_diff=$(cut -f4 "${path_hrpc_hla_clean}" | grep "DIFF" | wc -l)

    # calculating mismatch rate
    # Perform a floating-point calculation to determine if a mismatch rate exceeds a threshold of 30% (0.3).
    is_high_diff=$(awk -v diff="$n_diff" -v total="$n_total_var" 'BEGIN { print ( (diff/total) > 0.3 ? 1 : 0 ) }')
    echo -e "    - Mismatch rate > 0.3 (yes=1, no=0): ${is_high_diff}" >> "${log}"



    # Run

    # If needs swapping

    if [ "$is_high_diff" -eq 1 ]; then

        # switch error rate before swap
        mismatch_rate=$((100 * "${n_diff}" / "${n_total_var}"))
        echo "    - Mismatch rate before swap (> 30%)): ${mismatch_rate}%" >> "${log}"
        
        echo "- unswapped ${name}"

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
        
        path_diff="${path_swapped}"

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

                print "    ---------- after swap ----------"
                print "    - Genotyping Errors: " geno_err " (" (geno_err * 100 / ( NR - 1)) "%)"
                print "    - Final common heterozygous variants: " (NR - 1 - geno_err)
                print "    - Switch (Phasing) Errors: " phase_err " (" (phase_err * 100 / ( NR - 1 - geno_err)) "%)"

                #print "- Genotyping Error by homozigous in phased: " geno_hom_qry
                #print "- Genotyping Error by heterozigous in phased: " geno_het_qry
            }
            ' "${path_swapped}" &>> "${log}"

        
        # write genotyping errors
        path_geno_errors="${pathindiv}/${name}.geno.errors.tsv"
        awk '
            BEGIN{
                OFS="\t"
                print "idcomp", "hrpc", "hlam", "status", "hlam_swapped", "status_swapped"
            }
            NR>1 && $6 == "DIFF" {

                # split hrpc ($2) and swapped hla-mapper (%5)
                split($2, hrpc, "|")
                split($5, hlam, "|")

                # -- Phasing Error
                if (hrpc[1] == hlam[2] && hrpc[2] == hlam[1]) { next }
                
                # -- Genotyping Error
                else { print $0 }
            }
            ' "${path_swapped}" > "${path_geno_errors}"

        
        # tag genotyping errors
        path_diff_err="${pathindiv}/${name}.hrpc.hlamapper.switch.errors.tsv"
        awk '
            BEGIN{
                OFS="\t"
                print "idcomp", "hrpc", "hlam", "status", "hlam_swapped", "status_swapped", "status_error"
            }
            NR>1 {

                # split hrpc ($2) and swapped hla-mapper ($5)
                split($2, hrpc, "|")
                split($5, hlam, "|")

                # -- Genotyping Error
                if ($6 == "DIFF" && (hrpc[1] != hlam[2] || hrpc[2] != hlam[1])) { print $1, $2, $3, $4, $5, $6, "ERROR" }
                
                # -- Others
                else { print $1, $2, $3, $4, $5, $6, $6 }
            }
            ' "${path_diff}" > "${path_diff_err}"
        


    # If it's doesn't need

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

                #print "------------ NO SWAP ------------"
                print "    - Genotyping Errors: " geno_err " (" (geno_err * 100 / ( NR - 1)) "%)"
                print "    - Final common heterozygous variants: " (NR - 1 - geno_err)
                print "    - Switch (Phasing) Errors: " phase_err " (" (phase_err * 100 / ( NR - 1 - geno_err)) "%)"
                
                #print "    - Genotyping Error by homozigous in phased: " geno_hom_qry
                #print "    - Genotyping Error by heterozigous in phased: " geno_het_qry
            }
            ' "${path_hrpc_hla_clean}" &>> "${log}"

        
        # write genotyping errors
        path_geno_errors="${pathindiv}/${name}.geno.errors.tsv"
        awk '
            BEGIN{
                OFS="\t"
                print "idcomp", "hrpc", "hlam", "status"
            }
            NR>1 && $4 == "DIFF" {

                # split hrpc ($2) and hla-mapper ($3)
                split($2, hrpc, "|")
                split($3, hlam, "|")

                # -- Phasing Error
                if (hrpc[1] == hlam[2] && hrpc[2] == hlam[1]) { next }
                
                # -- Genotyping Error
                else { print $0 }
            }
            ' "${path_hrpc_hla_clean}" > "${path_geno_errors}"
        
        path_diff="${path_hrpc_hla_clean}"

        
        # tag genotyping errors
        path_diff_err="${pathindiv}/${name}.hrpc.hlamapper.switch.errors.tsv"
        awk '
            BEGIN{
                OFS="\t"
                print "idcomp", "hrpc", "hlam", "status", "status_error"
            }
            NR>1 {

                # split hrpc ($2) and swapped hla-mapper ($3)
                split($2, hrpc, "|")
                split($3, hlam, "|")

                # -- Genotyping Error
                if ($4 == "DIFF" && (hrpc[1] != hlam[2] || hrpc[2] != hlam[1])) { print $1, $2, $3, $4, "ERROR" }
                
                # -- Others
                else { print $1, $2, $3, $4, $4 }
            }
            ' "${path_diff}" > "${path_diff_err}"

    fi

    # calculate switch sizes
    path_switch_sizes="${pathindiv}/${name}.switch.sizes.tsv"
    grep -v "ERROR" "${path_diff_err}" | \
    cut -f1,4 | \
    awk '
        BEGIN { 
            OFS="\t"
            print "start", "end", "block_id", "block_size" 
        }
        $2 == "DIFF" {
            if (!in_block) {
                start = $1;
                in_block = 1;
                block_id++;
                count = 0;
            }
            end = $1;
            count++;
            next;
        }
        {
            if (in_block) {
                print start, end, block_id, count;
                in_block = 0;
            }
        }
        END {
            if (in_block) print start, end, block_id, count;
        }' > "${path_switch_sizes}"

echo "- File with tags ${path_diff_err}" >> "${log}"

echo "- File with switch blocks ${path_switch_sizes}" >> "${log}"

done



# end