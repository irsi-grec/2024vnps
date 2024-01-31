#!/bin/bash

## reference
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

-------------------------------------------------------------------------
# SNVs/short indels Calling


## combinegvcfs
GATK/4.1.8.0

output per chromosome: ${chr}.combined.g.vcf.gz

## genotypegvcfs

  	gatk GenotypeGVCFs \
	--reference hsapiens.GRCh38.hl.fasta \
	--variant ${chr}.combined.g.vcf.gz \
	-L chr${chr} \
	-D dbsnp_146.hg38.vcf.gz \
	-G StandardAnnotation \
	--output ${chr}.genotyped.vcf.gz

## normalize

	bcftools norm -m -both \
	-f hsapiens.GRCh38.hl.fasta \
	-Oz -o ${chr}.genotyped.norm.vcf.gz \
	${chr}.genotyped.vcf.gz

## annotate short variants with snpEff v5.0
for dbNSFP4.1a, gnomad v3.0, CADD v1.6 and ClinVar. For each preformatted database do (for example, dbnsfp):

 	java -jar SnpSift.jar dbnsfp \
  	-v -noCollapse -db dbNSFP4.1a_variant.vcf.gz \
	${chr}.genotyped.norm.vcf.gz | bgzip -c > ${chr}.genotyped.norm.annotated.vcf.gz

## gatk filters

 	gatk VariantFiltration \
	-R hsapiens.GRCh38.hl.fasta \
	-V ${chr}.norm.annot.snpEff.gnomAD.CADD.Clinvar.vcf.gz \
	--filter-expression "BaseQRankSum > 4.0 || BaseQRankSum < -4.0" --filter-name "BQB" \
	--filter-expression "FS > 60.000" --filter-name "SB60SNV" \
	--filter-expression "FS > 200.000" --filter-name "SB200INDEL" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPBSNV" \
	--filter-expression "ReadPosRankSum > 20.0" --filter-name "RPBINDEL" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQBSNV" \
	--filter-expression "QD < 2.0" --filter-name "QD2" \
	--filter-expression "MQ < 40.0" --filter-name "MQ40" \
	-O ${chr}.norm.annot.snpEff.gnomAD.CADD.Clinvar.GATKfilters.vcf.gz

## get supported variants
filter by at least 1 sample containing the variant has DP>9

## intervar annotation
1. convert variants from vcf into ANNOVAR's input file
2. predict pathogenicity
3. add intervar to vcf

## concat vcf's
per chromosome with bcftools concat

-------------------------------------------------------------------------
# HLA
-------------------------------------------------------------------------
singularity run \
xHLA_0.0.0.sif run.py \
--sample_id ${SAMPLE} \
--input_bam_path ${SAMPLE}.bam \
--output_path xHLA/ \
--full
