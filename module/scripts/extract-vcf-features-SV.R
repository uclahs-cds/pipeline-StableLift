#!/usr/bin/env Rscript
# extract-vcf-features-SV.R
####################################################################################################
#
# Extract features from vcf
# Intersect and annotate with gnomAD-SV vcf
#
####################################################################################################

suppressPackageStartupMessages({
    library(vcfR);
    library(data.table);
    library(argparse);
    library(GenomicRanges);
    });

###################################################################################################
# Input
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--variant-caller', type = 'character', help = '');
parser$add_argument('--input-vcf', type = 'character', help = 'Delly2 vcf');
parser$add_argument('--output-rds', type = 'character', help = 'Rds output for use in RF model');
parser$add_argument('--gnomad-rds', type = 'character', help = 'gnomAD Rds file');
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

# Set parameters for interactive runs
if (interactive()) {
    variant.caller <- 'Delly2';
    input.vcf <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSV/stableLift/train_CPCG-40QC_Delly2/CPCG-40QC_Delly2_LiftOver-GRCh38.vcf.gz';
    gnomad.rds <- '/hot/code/nkwang/GitHub/uclahs-cds/project-method-AlgorithmEvaluation-BNCH-000142-GRCh37v38/report/manuscript/publish/data/gnomad.v4.0.sv.Rds';
    }

###################################################################################################
# Functions
###################################################################################################
vcf.info.to.dt <- function(vcf.info) {
    # Split each string by semicolon and convert to a list of key-value pairs
    vcf.info <- strsplit(vcf.info, ';');
    vcf.info <- lapply(vcf.info, function(x) {
        x <- strsplit(x, '=');
        as.list(stats::setNames(sapply(x, `[`, 2), sapply(x, `[`, 1)));
        })

    # Combine the list of key-value pairs into a data table
    rbindlist(vcf.info, fill = TRUE);
    }

calculate.VAF <- function(GT.row) {
    total <- sum(GT.row %in% c('0/0', '0/1', '1/1'), na.rm = TRUE) * 2;
    alt <- sum(GT.row == '0/1', na.rm = TRUE) + sum(GT.row == '1/1', na.rm = TRUE) * 2;
    return(alt / total);
    }

get.overlap <- function(start1, end1, start2, end2) {
    max.length <- pmax((end1 - start1), (end2 - start2));
    overlap.length <- pmin(end1, end2) - pmax(start1, start2);
    return(overlap.length / max.length);
    }

find.SV.match <- function(this.ID, input, reference, overlap, offset) {
    # Match SV type and CHR
    this.variant <- input[ID == this.ID];
    reference <- reference[SVTYPE == this.variant$SVTYPE & CHROM == this.variant$CHROM];

    if (this.variant$SVTYPE == 'BND') {
        # reference[, OFFSET := abs(POS - this.variant$POS) + abs(POS2 - this.variant$POS2)];
        reference[, OFFSET := abs(POS - this.variant$POS)];
        reference <- reference[OFFSET < offset & CHR2 == this.variant$CHR2][order(OFFSET)];
    } else {
        reference[, OVERLAP := get.overlap(POS, END, this.variant$POS, this.variant$END)];
        reference <- reference[OVERLAP > overlap][order(OVERLAP, decreasing = TRUE)];
        }

    return(list(gnomad.match.ID = reference[1, ID], gnomad.matches = nrow(reference)));
    }

###################################################################################################
# Load files
###################################################################################################
input.vcf <- read.vcfR(input.vcf);
features.dt.gnomad <- readRDS(gnomad.rds);

###################################################################################################
# Data preprocessing
###################################################################################################
# Convert variant information into dt
input.info <- vcf.info.to.dt(input.vcf@fix[, 'INFO']);
input.fix <- as.data.table(input.vcf@fix);
features.dt <- cbind(input.fix[, -c('INFO')], input.info);

# Format columns
features.dt[, CONSENSUS := NULL];
numeric.columns <- c('POS', 'QUAL', 'END', 'PE', 'MAPQ', 'SRMAPQ', 'INSLEN', 'HOMLEN', 'SR', 'SRQ', 'CE', 'RDRATIO', 'SVLEN', 'POS2');
features.dt[, (numeric.columns) := lapply(.SD, as.numeric), .SDcols = numeric.columns];

# Extract and aggregate per sample GT fields
gt.fields <- c('GQ', 'RC', 'RDCN', 'DR', 'DV', 'RR', 'RV');
for (field in gt.fields) {
    features.dt[, (field) := apply(extract.gt(input.vcf, element = ..field, as.numeric = TRUE), 1, mean, na.rm = TRUE)];
    }
features.dt[, COHORT_AF := apply(extract.gt(input.vcf, element = 'GT'), 1, calculate.VAF)];
features.dt[!SVTYPE %in% c('BND', 'INS'), SVLEN := END - POS + 1];
features.dt[, CIPOS := as.numeric(sapply(CIPOS, function(x) unlist(strsplit(x, ','))[2]))];

###################################################################################################
# Intersect variants with gnomAD SVs
###################################################################################################
start.time <- Sys.time();

# features.dt <- features.dt[1:100];
features.dt[, c('gnomad.match.ID', 'gnomad.matches') := rbindlist(lapply(ID, find.SV.match, input = features.dt, reference = features.dt.gnomad, overlap = 0.8, offset = 500))];

gnomad.features <- c('ID', 'AF', 'POPMAX_AF', 'NCR');
features.dt <- merge(features.dt, features.dt.gnomad[, ..gnomad.features], all.x = TRUE, by.x = 'gnomad.match.ID', by.y = 'ID');

cat(format(Sys.time() - start.time, nsmall = 2), '\n');

###################################################################################################
# Format features for RF
###################################################################################################
continuous.features <- c(
    'POS',
    'QUAL',
    'END',
    'PE',
    'MAPQ',
    'CIPOS',
    'SRMAPQ',
    'HOMLEN',
    'SR',
    'SRQ',
    'CE',
    'RDRATIO',
    'SVLEN',
    'GQ',
    'RC',
    'RDCN',
    'DR',
    'DV',
    'RR',
    'RV',
    'gnomad.matches',
    'AF',
    'POPMAX_AF',
    'NCR'
    );

categorical.features <- c(
    'CHROM',
    'SVTYPE',
    'CT'
    );

# Extract features and format
continuous.features <- continuous.features[continuous.features %in% names(features.dt)];
categorical.features <- categorical.features[categorical.features %in% names(features.dt)];
all.features <- c(continuous.features, categorical.features, 'ID');

features.dt <- features.dt[, ..all.features];
features.dt[, (continuous.features) := lapply(.SD, as.numeric), .SDcols = continuous.features];
features.dt[, (continuous.features) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = continuous.features];
features.dt[, (categorical.features) := lapply(.SD, function(x) ifelse(is.na(x), '', x)), .SDcols = categorical.features];
features.dt[, (categorical.features) := lapply(.SD, as.factor), .SDcols = categorical.features];
names(features.dt) <- make.names(names(features.dt));

# Remove rows with NA
features.dt.rows <- nrow(features.dt);
features.dt <- features.dt[apply(features.dt, 1, function(x) !any(is.na(x))), ];
cat('Removed', features.dt.rows - nrow(features.dt), 'rows with missing data\n');

###################################################################################################
# Save features.dt for input to RF
###################################################################################################
saveRDS(features.dt, output.rds);
