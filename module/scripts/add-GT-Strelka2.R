#!/usr/bin/env Rscript
# add-GT-Strelka2.R
####################################################################################################
#
# Adds GT to FORMAT fields
# 0/1 for called variants and ./. for missing data
#
####################################################################################################

suppressPackageStartupMessages({
    library(vcfR);
    library(data.table);
    library(argparse);
    });

###################################################################################################
# Input
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--input-vcf', type = 'character', help = 'Strelka2 VCF to add FORMAT/GT field to');
parser$add_argument('--output-vcf', type = 'character', help = 'Formatted output VCF path');
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

###################################################################################################
# Main
###################################################################################################
vcf <- read.vcfR(input.vcf);
format.fields <- unlist(strsplit(vcf@gt[1, 'FORMAT'], ':'));

if ('GT' %in% format.fields) {
    message('FORMAT/GT field already exists');
    write.vcf(vcf, output.vcf);
} else {
    vcf@meta <- append(vcf@meta, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', after = grep('##FORMAT', vcf@meta)[1] - 1);
    vcf@gt[, 'FORMAT'] <- paste0('GT:', vcf@gt[, 'FORMAT']);
    gt <- vcf@gt[, colnames(vcf@gt) != 'FORMAT'];
    vcf@gt[, colnames(vcf@gt) != 'FORMAT'] <- ifelse(startsWith(gt, '.'), paste0('./.:', gt), paste0('0/1:', gt));
    write.vcf(vcf, output.vcf);
    }
