
# Phasing comparison of the gold-standard with the pipeline optimized for HLA phasing using HPRC data

#
## Data

### HPRC

The Human Pangenome Reference Consortium was launched in May 2025 ([HPRC Release 2](https://humanpangenome.org/hprc-data-release-2/), [timeline](https://humanpangenome.org/release-timeline/), [ensemble](https://projects.ensembl.org/hprc/)).

- Several sequencing technologies were used: PacBio HiFi, ONT Ultralong, Dovetail/Illumina Hi-C, PacBio Kinnex RNA, and Illumina WGS. High coverage: 60X PacBio HiFi and 30X Oxford Nanopore in 100 kb.

There are [234 samples](https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables/sample) in total [metadata](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv):
- 216 from HPRC: 126 phased with parental information seq by Illumina and 90 unrelated phased with hifiasm-hic,
- 14 from HPP (phased with hifiasm)
- 4 are the gold-standard samples from the community (they are reference assemblies), sequenced and assembled by other projects: CHM13 is the first (haploid) telomere-to-telomere (T2T) assembly of the human genome, GRCh38 is a current "classic" reference. HG002 is a sample (trio) from Genome in a Bottle (GIAB) used for benchmarking.

- The [year1_freeze_assembly_v2](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG01109/assemblies/year1_freeze_assembly_v2/) is the curated and polished version. The consortium took the [raw data](https://github.com/human-pangenomics/HPP_Year1_Assemblies), cleaned it, checked the quality, and standardized it. This is the "Gold Standard" published in the Year 1 paper. They used dipcall ([dipcall github](https://github.com/lh3/dipcall/tree/master)) to convert the fasta files to vcf, writing each line perfectly phased (1|0 or 0|1).

- The [S3 root of HPRC_PLUS](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/) have 18 samples from the offspring, originating from the trios.

- Samples (14) into HLA-MAPPER vcf: HG00733, HG01109, HG01243, HG02055, HG02080, HG02145, HG02723, HG02818, HG03098, HG03486, HG03492, NA18906, NA19240 and NA20129.

- Samples (4) out HLA-MAPPER vcf: HG002, HG005, HG02109 and NA21309.

#
## Report

- [Report of mapping phasing errors](https://htmlpreview.github.io/?https://github.com/jenn-emn/02_tr/blob/main/report/phasing_report.html)
- [dasd](https://htmlpreview.github.io/?https://github.com/jenn-emn/02_tr/blob/main/reports/phasing_report.html)