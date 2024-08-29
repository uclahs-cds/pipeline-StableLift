#!/usr/bin/env Rscript
# extract-VCF-features-SV.R
####################################################################################################
#
# LiftOver Delly2 structural variants and annotate with gnomAD population allele frequency
# Extract VCF features and save as Rds for input to predict-variant-stability.R
#
####################################################################################################

suppressPackageStartupMessages({
    library(vcfR);
    library(data.table);
    library(argparse);
    library(rtracklayer);
    library(GenomicRanges);
    });

###################################################################################################
# Input
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--input-vcf', type = 'character', help = 'Input Delly2 VCF');
parser$add_argument('--source-build', type = 'character', help = 'One of {GRCh37, GRCh38}');
parser$add_argument('--chain-file', type = 'character', help = 'Chain file for coordinate conversion');
parser$add_argument('--header-contigs', type = 'character', help = 'Resource file with VCF header for target build');
parser$add_argument('--gnomad-rds', type = 'character', help = 'gnomAD-SV v4 resource file');
parser$add_argument('--output-vcf', type = 'character', help = 'VCF output');
parser$add_argument('--output-rds', type = 'character', help = 'Rds output for input to RF model');
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

###################################################################################################
# Functions
###################################################################################################
# vcfR::getINFO() to data.table
vcf.info.to.dt <- function(vcf.info) {
    vcf.info <- lapply(vcf.info, function(x) vcf.info.string.to.list(x));
    feature.names <- unique(unlist(lapply(vcf.info, names)));
    vcf.info <- do.call(mapply, c(FUN = list, lapply(vcf.info, `[`, feature.names)));
    setNames(as.data.table(vcf.info), feature.names);
    }

# Split VCF info field to list
vcf.info.string.to.list <- function(vcf.info, keep.columns = NULL) {
    list.out <- strsplit(vcf.info, split = ';');
    list.out <- lapply(list.out, function(x) strsplit(x, split = '='));
    labels <- sapply(list.out[[1]], function(x) x[[1]]);
    values <- sapply(list.out[[1]], function(x) if (length(x) == 2) x[[2]] else x[[1]]);
    names(values) <- labels;
    if (is.null(keep.columns)) return(values);
    values <- values[labels %in% keep.columns];
    return(values);
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
        reference <- reference[CHR2 == this.variant$CHR2];
        reference[, OFFSET := abs(POS - this.variant$POS)];
        reference <- reference[OFFSET < offset][order(OFFSET)];
    } else {
        reference[, OVERLAP := get.overlap(POS, END, this.variant$POS, this.variant$END)];
        reference <- reference[OVERLAP > overlap][order(OVERLAP, decreasing = TRUE)];
        }

    return(list(gnomad.match.ID = reference[1, ID], gnomad.matches = nrow(reference)));
    }

annotate.gnomad.features <- function(features.dt, features.dt.gnomad) {
    gnomad.features <- c('ID', 'AF', 'POPMAX_AF');
    features.dt[, c('gnomad.match.ID', 'gnomad.matches') := rbindlist(lapply(ID, find.SV.match, input = features.dt, reference = features.dt.gnomad, overlap = 0.8, offset = 500))];
    features.dt <- merge(features.dt, features.dt.gnomad[, ..gnomad.features], all.x = TRUE, by.x = 'gnomad.match.ID', by.y = 'ID');
    }

###################################################################################################
# Load files
###################################################################################################
input.vcf <- read.vcfR(input.vcf);
header.contigs <- scan(header.contigs, character());
liftover.chain <- import.chain(chain.file);
features.dt.gnomad <- readRDS(gnomad.rds);

###################################################################################################
# Data preprocessing
###################################################################################################
# Convert variant information into dt
if (any(duplicated(input.vcf@fix[, 'ID']))) input.vcf@fix[, 'ID'] <- paste0(substr(input.vcf@fix[, 'ID'], 1, 3), sprintf('%08d', seq_len(nrow(input.vcf@fix))));
input.info <- vcf.info.to.dt(input.vcf@fix[, 'INFO']);
input.fix <- as.data.table(input.vcf@fix);
features.dt <- cbind(input.fix[, -c('INFO')], input.info);

# Format columns
features.dt[, CONSENSUS := NULL];
numeric.columns <- c('POS', 'QUAL', 'END', 'PE', 'MAPQ', 'SRMAPQ', 'INSLEN', 'HOMLEN', 'SR', 'SRQ', 'CE', 'RDRATIO', 'SVLEN', 'POS2');
character.columns <- names(features.dt)[!names(features.dt) %in% numeric.columns];
features.dt[, (numeric.columns) := lapply(.SD, as.numeric), .SDcols = numeric.columns];
features.dt[, (character.columns) := lapply(.SD, as.character), .SDcols = character.columns];

# Extract and aggregate per sample GT fields
gt.fields <- c('GQ', 'RC', 'RDCN', 'DR', 'DV', 'RR', 'RV');
for (field in gt.fields) {
    features.dt[, (field) := apply(extract.gt(input.vcf, element = ..field, as.numeric = TRUE), 1, mean, na.rm = TRUE)];
    }
features.dt[, COHORT_AF := apply(extract.gt(input.vcf, element = 'GT'), 1, calculate.VAF)];
features.dt[, CIPOS := as.numeric(sapply(CIPOS, function(x) unlist(strsplit(x, ','))[2]))];

if (source.build == 'GRCh37') {
    features.dt[, CHROM := paste0('chr', CHROM)];
    features.dt[, CHR2 := ifelse(is.na(CHR2), CHR2, paste0('chr', CHR2))];
} else if (source.build == 'GRCh38') {
    # Annotate with gnomAD before LiftOver if source build == GRCh38
    features.dt <- annotate.gnomad.features(features.dt, features.dt.gnomad);
    }

###################################################################################################
# LiftOver variants by breakpoint
###################################################################################################
# Create GRanges object
grange.source <- makeGRangesFromDataFrame(
    df = features.dt,
    seqnames.field = 'CHROM',
    start.field = 'POS',
    end.field = 'END',
    keep.extra.columns = TRUE
    );

# LiftOver using chain file
grange.target <- unlist(liftOver(grange.source, liftover.chain));
grange.target.dt <- as.data.table(grange.target);

# Create GRanges object using CHROM, CHR2, and POS2 from features.dt
grange.source.BND <- makeGRangesFromDataFrame(
    df = features.dt[SVTYPE == 'BND', ],
    seqnames.field = 'CHR2',
    start.field = 'POS2',
    end.field = 'POS2',
    keep.extra.columns = TRUE
    );
grange.target.BND <- as.data.table(unlist(liftOver(grange.source.BND, liftover.chain)));

# Remove multiple mappings
grange.target.dt <- grange.target.dt[!duplicated(ID)];
grange.target.BND <- grange.target.BND[!duplicated(ID)];
common <- intersect(grange.target.dt$ID, grange.target.BND$ID);

grange.target.dt[ID %in% common, c('CHR2', 'POS2') := grange.target.BND[ID %in% common, .(seqnames, start)]];

if (source.build == 'GRCh38') {
    grange.target.dt[, seqnames := sub('chr', '', seqnames)];
    grange.target.dt[, CHR2 := ifelse(is.na(CHR2), CHR2, sub('chr', '', CHR2))];
    }

###################################################################################################
# Write output VCF
###################################################################################################
pass.liftover <- as.data.table(input.vcf@fix)$ID %in% grange.target.dt$ID;
input.fix <- as.data.table(input.vcf@fix)[pass.liftover];
input.gt <- as.data.table(input.vcf@gt)[pass.liftover];
grange.target.dt <- grange.target.dt[match(input.fix$ID, grange.target.dt$ID)];

for (i in seq_len(nrow(input.fix))) {
    this.ID <- input.fix[i, ID];
    this.INFO <- vcf.info.string.to.list(input.fix[i, INFO]);
    this.INFO[['END']] <- grange.target.dt[i, end];
    if (this.INFO[['SVTYPE']] == 'BND') {
        this.INFO[['CHR2']] <- grange.target.dt[i, CHR2];
        this.INFO[['POS2']] <- grange.target.dt[i, POS2];
        }
    if (this.INFO[['SVTYPE']] %in% )
    this.INFO <- lapply(names(this.INFO), function(x) paste(x, this.INFO[[x]], sep = '='));
    this.INFO <- paste(this.INFO, collapse = ';');
    this.INFO <- gsub('IMPRECISE=IMPRECISE', 'IMPRECISE', this.INFO);
    this.INFO <- gsub('PRECISE=PRECISE', 'PRECISE', this.INFO);
    this.INFO <- gsub('SOMATIC=SOMATIC', 'SOMATIC', this.INFO);
    input.fix[i, c('CHROM', 'POS', 'INFO') := grange.target.dt[i, .(seqnames, start, ..this.INFO)]];
    }

lifted.vcf <- input.vcf;
lifted.vcf@fix <- as.matrix(input.fix);
lifted.vcf@gt <- as.matrix(input.gt);
lifted.vcf@meta <- lifted.vcf@meta[!grepl("^##(contig|reference)", lifted.vcf@meta)];
lifted.vcf@meta <- c(lifted.vcf@meta, header.contigs);

write.vcf(lifted.vcf, output.vcf);

###################################################################################################
# Format features for RF
###################################################################################################
features.dt <- features.dt[ID %in% grange.target.dt$ID];
features.dt <- features.dt[match(input.fix$ID, features.dt$ID)];
features.dt[, c('CHROM', 'POS', 'END', 'CHR2', 'POS2') := grange.target.dt[, .(seqnames, start, end, CHR2, POS2)]];
features.dt[!SVTYPE %in% c('BND', 'INS'), SVLEN := END - POS + 1];

if (source.build == 'GRCh37') {
    features.dt <- annotate.gnomad.features(features.dt, features.dt.gnomad);
    }

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
    'AF',
    'gnomad.matches',
    'POPMAX_AF'
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
