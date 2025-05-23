#!/usr/bin/env Rscript
# predict-liftover-stability.R
####################################################################################################
#
# Apply random forest model to predict variant stability
#
####################################################################################################

suppressPackageStartupMessages({
    library(ranger);
    library(argparse);
    library(ROCR);
    library(data.table);
    });

###################################################################################################
# Set parameters
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--variant-caller', type = 'character', help = 'One of {HaplotypeCaller, Mutect2, Strelka2, SomaticSniper, Muse2, Delly2-gSV, Delly2-sSV}');
parser$add_argument('--features-dt', type = 'character', help = 'Processed Rds file with variant info and annotations');
parser$add_argument('--rf-model', type = 'character', help = 'Pre-trained random forest model Rds file');
parser$add_argument('--specificity', type = 'numeric', help = 'Target specificity based on whole genome validation set, overrides `--threshold`');
parser$add_argument('--threshold', type = 'numeric', help = 'Stability score threshold, default based on maximizing F1-score in whole genome validation set');
parser$add_argument('--output-tsv', type = 'character', help = 'Output TSV with predicted Stability Scores');
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

###################################################################################################
# Apply random forest model
###################################################################################################
cat('\nPredicting variant stability with', basename(rf.model.path), '\n');
stability <- predict(rf.model, data = features.dt);
performance <- performance(rf.model$prediction, 'sens', 'spec');

if (!is.null(specificity) && is.numeric(specificity)) {
    cat('Target specificity =', specificity, '\n');
    operating.index <- max(which(unlist(performance@x.values) > specificity));
    threshold <- 1 - unlist(performance@alpha.values)[operating.index];
    sensitivity <- unlist(performance@y.values)[operating.index];

    cat(sprintf('Threshold = %.3f\n', threshold));
    cat(sprintf('Projected sensitivity = %.3f\n', sensitivity));
} else if (!is.null(threshold) && is.numeric(threshold)) {
    cat('Target threshold =', threshold, '\n');
    operating.index <- min(which(unlist(performance@alpha.values) <= 1 - threshold));
    specificity <- unlist(performance@x.values)[operating.index];
    sensitivity <- unlist(performance@y.values)[operating.index];

    cat(sprintf('Projected sensitivity = %.3f\n', sensitivity));
    cat(sprintf('Projected specificity = %.3f\n', specificity));
} else {
    performance.f <- performance(rf.model$prediction, measure = 'f');
    operating.index <- which.max(unlist(performance.f@y.values));
    threshold <- 1 - unlist(performance.f@x.values)[operating.index];
    sensitivity <- unlist(performance@y.values)[operating.index];
    specificity <- unlist(performance@x.values)[operating.index];

    cat(sprintf('Default threshold = %.3f\n', threshold));
    cat(sprintf('Projected sensitivity = %.3f\n', sensitivity));
    cat(sprintf('Projected specificity = %.3f\n', specificity));
    }

stability.classification <- ifelse(stability$predictions[, 1] < threshold, 1, 0);
cat(sprintf('Proportion predicted unstable = %.3f\n\n', mean(stability.classification)));
stability.classification <- as.factor(stability.classification);

###################################################################################################
# Output stability scores
###################################################################################################
annotation.dt <- data.table(
    CHROM = features.dt$CHROM,
    POS = features.dt$POS,
    STABILITY_SCORE = format(round(stability$predictions[, 1], 4), nsmall = 4),
    STABILITY = ifelse(stability.classification == '1', 'UNSTABLE', 'STABLE')
    );
setorder(annotation.dt, CHROM, POS);

fwrite(annotation.dt, file = output.tsv, sep = '\t', col.names = FALSE);
