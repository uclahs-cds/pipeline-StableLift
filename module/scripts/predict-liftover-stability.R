#!/usr/bin/env Rscript
# predict-liftover-stability.R
####################################################################################################
#
# Apply random forest model to predict variant LiftOver stability
# Validate results and plot model performance with discordance file
#
####################################################################################################

suppressPackageStartupMessages({
    library(caret);
    library(ranger);
    library(argparse);
    library(data.table);
    });

###################################################################################################
# Set parameters
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--features-dt', type = 'character');
parser$add_argument('--rf-model', type = 'character');
parser$add_argument('--variant-caller', type = 'character');
parser$add_argument('--specificity', type = 'numeric', help = 'Target specificity, overrides `--threshold`');
parser$add_argument('--threshold', type = 'numeric', help = 'Stability score threshold', default = 0.5);
parser$add_argument('--output-tsv', type = 'character', help = 'TSV output file');
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

###################################################################################################
# Load data
###################################################################################################
features.dt.path <- features.dt;
features.dt <- readRDS(features.dt.path);
rf.model.path <- rf.model;
rf.model <- readRDS(rf.model);

####################################################################################################
# Feature engineering
####################################################################################################
if (variant.caller == 'HaplotypeCaller') {
    normalize.features <- c('FS', 'VQSLOD', 'QD', 'SOR');
    features.dt[, (normalize.features) := lapply(.SD, scale), .SDcols = normalize.features];
    features.dt[, c('QUAL', 'GQ', 'DP') := NULL];
} else if (variant.caller == 'Mutect2') {
    features.dt[, c('TRINUCLEOTIDE') := NULL];
    features.dt <- features.dt[Gencode_34_variantType != 'SNP'];
    features.dt[, c('Gencode_34_variantType') := NULL];
} else if (variant.caller == 'Muse2') {
    features.dt[, c('TRINUCLEOTIDE', 'Gencode_34_variantType') := NULL];
} else if (variant.caller == 'Strelka2') {
    features.dt[, c('TRINUCLEOTIDE_SEQ', 'DP', 'Gencode_34_variantType') := NULL];
} else if (variant.caller == 'Delly2') {
    normalize.features <- c('SR', 'SRQ', 'DV');
    features.dt[, (normalize.features) := lapply(.SD, scale), .SDcols = normalize.features];
    }

dim(features.dt);

###################################################################################################
# Apply random forest model
###################################################################################################
cat('\nPredicting liftover stability with', basename(rf.model.path), '\n');
stability <- predict(rf.model, data = features.dt);

if (!is.null(specificity) && is.numeric(specificity)) {
    cat('Target specificity =', specificity, '\n');
    operating.index <- max(which(unlist(rf.model$performance@x.values) < 1 - specificity));
    sensitivity <- unlist(rf.model$performance@y.values)[operating.index];
    cat('Projected sensitivity =', round(sensitivity, 3), '\n');
    threshold <- 1 - unlist(rf.model$performance@alpha.values)[operating.index];
    cat('Stability score threshold =', round(threshold, 3), '\n');
} else {
    cat('Target threshold =', threshold, '\n');
    operating.index <- min(which(unlist(rf.model$performance@alpha.values) <= 1 - threshold));
    specificity <- 1 - unlist(rf.model$performance@x.values)[operating.index];
    sensitivity <- unlist(rf.model$performance@y.values)[operating.index];
    cat('Projected specificity =', round(specificity, 3), '\n');
    cat('Projected sensitivity =', round(sensitivity, 3), '\n');
    }

stability.classification <- ifelse(stability$predictions[, 1] < threshold, 1, 0);
stability.classification <- as.factor(stability.classification);
cat('Proportion predicted unstable =', round(mean(as.numeric(as.character(stability.classification))), 3), '\n');

###################################################################################################
# Output stability scores
###################################################################################################
annotation.dt <- data.table(
    CHROM = features.dt$CHROM,
    POS = features.dt$POS,
    STABILITY_SCORE = format(round(stability$predictions[, 1], 4), nsmall = 4),
    STABILITY = ifelse(stability.classification == '1', 'UNSTABLE', 'STABLE')
    );
fwrite(annotation.dt, file = output.tsv, sep = '\t', col.names = TRUE);
