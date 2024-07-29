#!/usr/bin/env Rscript
# liftover-Delly2-vcf.R
###################################################################################################
#
#
#
###################################################################################################

suppressPackageStartupMessages({
    library(vcfR);
    library(data.table);
    library(argparse);
    library(rtracklayer);
    });

###################################################################################################
# Input
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--input-vcf', type = 'character', help = 'GRCh37 Delly2 vcf');
parser$add_argument('--header-contigs', type = 'character', help = 'Directory with vcf subsets');
parser$add_argument('--chain-file', type = 'character', help = 'hg19ToHg38.over.chain file');
parser$add_argument('--output', type = 'character', help = 'Where to write lifted vcf');
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

# Set parameters for interactive runs
if (interactive()) {
    # input.vcf <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/gSV/bcftools-merge/CPCG-40QC_GRCh37/CPCG-40QC_GRCh37_regenotype-gSV_delly_bcftools-merge_delly-filter-germline.vcf.gz';
    input.vcf <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSV/bcftools-merge/CPCG-40QC_GRCh37/CPCG-40QC_GRCh37_call-sSV_delly_bcftools-merge_somatic-only.vcf.gz';
    header.contigs <- '/hot/code/nkwang/GitHub/uclahs-cds/project-method-AlgorithmEvaluation-BNCH-000142-GRCh37v38/report/manuscript/publish/GRCh38-vcf-header-contigs.txt';
    chain.file <- '/hot/resource/genomics/liftover_chain_files/hg19ToHg38.over.chain';
    output <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSV/stableLift/train_CPCG-40QC_Delly2/CPCG-40QC_Delly2_LiftOver-GRCh38.vcf.gz';
    }

###################################################################################################
# Functions
###################################################################################################
vcf.fix.to.dt <- function(vcf.fix) {
    vcf.fix <- as.data.table(vcf.fix);
    vcf.info <- vcf.info.to.dt(vcf.fix$INFO);
    cbind(vcf.fix[, -'INFO'], vcf.info);
    }

# vcfR::getINFO() to data.table
vcf.info.to.dt <- function(vcf.info) {
    vcf.info <- lapply(vcf.info, function(x) vcf.info.string.to.list(x));
    feature.names <- unique(unlist(lapply(vcf.info, names)));
    vcf.info <- do.call(mapply, c(FUN = list, lapply(vcf.info, `[`, feature.names)));
    setNames(as.data.table(vcf.info), feature.names);
    }

# Split vcf info field to list
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

###################################################################################################
# Load files
###################################################################################################
input.vcf.path <- input.vcf;
input.vcf <- read.vcfR(input.vcf);
header.contigs <- scan(header.contigs, character());
liftover.chain <- import.chain(chain.file);

###################################################################################################
# Data preprocessing
###################################################################################################
if (any(duplicated(input.vcf@fix[, 'ID']))) input.vcf@fix[, 'ID'] <- paste0(substr(input.vcf@fix[, 'ID'], 1, 3), sprintf('%08d', seq_len(nrow(input.vcf@fix))));
# if (any(duplicated(input.vcf@fix[, 'ID']))) fix.dt[, ID := paste0(substr(ID, 1, 3), sprintf('%08d', seq_len(nrow(fix.dt))))];
fix.dt <- as.data.table(input.vcf@fix);
gt.dt <- as.data.table(input.vcf@gt);

fix.dt <- vcf.fix.to.dt(fix.dt);
fix.dt[, CHROM := paste0('chr', CHROM)];
fix.dt[, CHR2 := ifelse(is.na(CHR2), CHR2, paste0('chr', CHR2))];

fix.dt <- fix.dt[, -c('CONSENSUS')];
numeric.columns <- c('POS', 'QUAL', 'END', 'PE', 'MAPQ', 'SRMAPQ', 'INSLEN', 'HOMLEN', 'SR', 'SRQ', 'CE', 'SVLEN', 'POS2');
character.columns <- names(fix.dt)[!names(fix.dt) %in% numeric.columns];
fix.dt[, (numeric.columns) := lapply(.SD, as.numeric), .SDcols = numeric.columns];
fix.dt[, (character.columns) := lapply(.SD, as.character), .SDcols = character.columns];

###################################################################################################
# Liftover
###################################################################################################
# Create GRanges object
granges.37 <- makeGRangesFromDataFrame(
    df = fix.dt,
    seqnames.field = 'CHROM',
    start.field = 'POS',
    end.field = 'END',
    keep.extra.columns = TRUE
    );

# Liftover using chain file
granges.38 <- unlist(liftOver(granges.37, liftover.chain));
granges.38.dt <- as.data.table(granges.38);

# Create GRanges object using CHROM, CHR2, and POS2 from fix.dt
granges.37.BND <- makeGRangesFromDataFrame(
    df = fix.dt[SVTYPE == 'BND', ],
    seqnames.field = 'CHR2',
    start.field = 'POS2',
    end.field = 'POS2',
    keep.extra.columns = TRUE
    );
granges.38.BND <- as.data.table(unlist(liftOver(granges.37.BND, liftover.chain)));

# Remove multiple mappings
granges.38.dt <- granges.38.dt[!duplicated(ID)];
granges.38.BND <- granges.38.BND[!duplicated(ID)];
common <- intersect(granges.38.dt$ID, granges.38.BND$ID);

granges.38.dt[ID %in% common, c('CHR2', 'POS2') := granges.38.BND[ID %in% common, .(seqnames, start)]];

pass.liftover <- as.data.table(input.vcf@fix)$ID %in% granges.38.dt$ID;
fix.lifted <- as.data.table(input.vcf@fix)[pass.liftover];
gt.dt <- gt.dt[pass.liftover];
for (i in seq_len(nrow(fix.lifted))) {
    this.ID <- fix.lifted[i, ID];
    this.INFO <- vcf.info.string.to.list(fix.lifted[i, INFO]);
    this.INFO[['END']] <- granges.38.dt[i, end];
    if (this.INFO[['SVTYPE']] == 'BND') {
        this.INFO[['CHR2']] <- granges.38.dt[i, CHR2];
        this.INFO[['POS2']] <- granges.38.dt[i, POS2];
        }
    this.INFO <- lapply(names(this.INFO), function(x) paste(x, this.INFO[[x]], sep = '='));
    this.INFO <- paste(this.INFO, collapse = ';');
    this.INFO <- gsub('IMPRECISE=IMPRECISE', 'IMPRECISE', this.INFO);
    this.INFO <- gsub('PRECISE=PRECISE', 'PRECISE', this.INFO);
    this.INFO <- gsub('SOMATIC=SOMATIC', 'SOMATIC', this.INFO);
    fix.lifted[i, c('CHROM', 'POS', 'INFO') := granges.38.dt[ID == ..this.ID, .(seqnames, start, ..this.INFO)]];
    }

###################################################################################################
# Write output vcf
###################################################################################################
output.vcf <- input.vcf;
output.vcf@fix <- as.matrix(fix.lifted);
output.vcf@gt <- as.matrix(gt.dt);
output.vcf@meta <- output.vcf@meta[!grepl('^##(contig|reference)', output.vcf@meta)];
output.vcf@meta <- c(output.vcf@meta, header.contigs);

write.vcf(output.vcf, output);
