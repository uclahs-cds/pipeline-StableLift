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
    library(ROCR);
    library(data.table);
    library(BoutrosLab.plotting.general);
    });

###################################################################################################
# Set parameters
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--features-dt', type = 'character');
parser$add_argument('--rf-model', type = 'character');
parser$add_argument('--variant-caller', type = 'character');
parser$add_argument('--discordance-file', type = 'character');
parser$add_argument('--specificity', type = 'numeric', help = 'Target specificity, overrides `--threshold`');
parser$add_argument('--threshold', type = 'numeric', help = 'Stability score threshold', default = 0.5);
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

# Set parameters for interactive runs
if (interactive()) {
    variant.caller <- 'Mutect2';
    # features.dt <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSNV/stableLift/validate_TCGA-SARC_WGS_Mutect2/TCGA-SARC_WGS_Mutect2_LiftOver-GRCh38_annotated.Rds';
    # discordance.file <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSNV/stableLift/validate_TCGA-SARC_WGS_Mutect2/TCGA-SARC_WGS_Mutect2_bcftools-stats-verbose_PSD.txt';
    features.dt <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSNV/stableLift/validate_TCGA-SARC_WXS_Mutect2/TCGA-SARC_WXS_Mutect2_LiftOver-GRCh38_annotated.Rds';
    discordance.file <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSNV/stableLift/validate_TCGA-SARC_WXS_Mutect2/TCGA-SARC_WXS_Mutect2_bcftools-stats-verbose_PSD.txt';
    rf.model <- '/hot/project/method/AlgorithmEvaluation/BNCH-000142-GRCh37v38/sSNV/stableLift/train_CPCG-40QC_Mutect2/RF-train_Mutect2_ntree1000_nodesize5_classratio0.Rds';
    specificity = 0.99;
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
# Combine with NRD for validation
###################################################################################################
if (!is.null(discordance.file)) {
    if (variant.caller != 'Delly2') {
        site.discordance <- fread(
            discordance.file,
            select = c(2, 3, 4, 5, 6)
            );

        # Format data
        names(site.discordance) <- c('CHROM', 'POS', 'n.match', 'n.mismatch', 'NRD');
        site.discordance[, NRD := 1 - (n.match / (n.match + n.mismatch))];
        site.discordance <- unique(site.discordance, by = c('CHROM', 'POS'));

        # Merge with all funcotator variants since bcftools stats only reports PSD for shared variants
        features.dt <- merge(
            x = features.dt,
            y = site.discordance[, c('CHROM', 'POS', 'NRD')],
            by = c('CHROM', 'POS'),
            all.x = TRUE
            );

        features.dt[is.na(NRD), NRD := 1];
    } else {
        site.discordance <- read.vcfR(discordance.file);
        site.discordance.info <- vcf.info.to.dt(site.discordance@fix[, 'INFO']);
        site.discordance.fix <- as.data.table(site.discordance@fix);
        site.discordance <- cbind(site.discordance.fix[, -c('INFO')], site.discordance.info);

        site.discordance[is.na(NON_REF_GENOTYPE_CONCORDANCE), NON_REF_GENOTYPE_CONCORDANCE := 0];
        site.discordance[, NRD := 1 - as.numeric(NON_REF_GENOTYPE_CONCORDANCE)];

        features.dt <- merge(
            x = features.dt,
            y = site.discordance[, c('ID', 'NRD')],
            by = c('ID'),
            all.x = TRUE
            );
        features.dt[, ID := NULL];
        }

    # Dichotomize NRD based on threshold
    NRD.threshold <- ifelse(variant.caller == 'Delly2', 0.2, 0.1);
    features.dt[, NRD := as.factor(as.numeric(NRD > NRD.threshold))];
    }

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
# Calculate and plot performance
###################################################################################################
if (!is.null(discordance.file)) {
    # Print confusion matrix
    print(confusionMatrix(stability.classification, features.dt$NRD, positive = '1'));

    # Calculate prediction performance
    stability.prediction <- prediction(stability$predictions[, 2], features.dt$NRD);
    stability.performance <- performance(stability.prediction, 'tpr', 'fpr');

    # Calculate AUC
    stability.auc <- performance(stability.prediction, 'auc')@y.values[[1]];
    cat('AUC =', round(stability.auc, 4), '\n');

    # Plot ROC curve
    this.data <- data.table(TPR = unlist(stability.performance@y.values), FPR = unlist(stability.performance@x.values));
    this.filename <- sub('_LiftOver-GRCh38_annotated.Rds', '_stableLift-GRCh38_ROC.png', features.dt.path);
    create.scatterplot(
        filename = this.filename,
        formula = TPR ~ FPR,
        data = this.data,
        type = 's',
        lwd = 3,
        xlimits = c(-0.03, 1.03),
        ylimits = c(-0.03, 1.03),
        add.xyline = TRUE,
        xyline.col = 'gray',
        xyline.lty = 2,
        add.text = TRUE,
        text.labels = paste0('AUC = ', format(round(stability.auc, 3), nsmall = 3)),
        text.x = 0.63,
        text.y = 0.07,
        text.cex = 1.7,
        height = 4,
        width = 4
        );

    # Plot stability score distribution
    this.data <- data.table(
        StabilityScore = stability$predictions[, 1],
        DiscordanceStatus = factor(ifelse(features.dt$NRD == 1, 'Discordant', 'Concordant'), levels = c('Concordant', 'Discordant'))
        );
    this.filename <- sub('_LiftOver-GRCh38_annotated.Rds', '_stableLift-GRCh38_stability-score-histogram.png', features.dt.path);
    bins <- seq(0, 1, 0.05);
    colors <- c(Concordant = 'darkgreen', Discordant = 'darkred');
    labels <- c(
        paste0('bold("Concordant (', round(mean(this.data$DiscordanceStatus == 'Concordant') * 100, 1), '%)")'),
        paste0('bold("Discordant (', round(mean(this.data$DiscordanceStatus == 'Discordant') * 100, 1), '%)")')
        );

    this.plot <- list();
    for (group in levels(this.data$DiscordanceStatus)) {
        this.axis <- cut(this.data[DiscordanceStatus == group, StabilityScore], breaks = bins, include.lowest = TRUE) |> table();
        this.axis <- auto.axis(this.axis, log.scaled = FALSE);
        if (group == 'Concordant') {
            this.legend <- list(
                inside = list(
                    fun = draw.key,
                    args = list(
                        key = list(
                            points = list(
                                pch = 22,
                                col = 'black',
                                fill = colors,
                                cex = 2.2
                                ),
                            text = list(
                                lab = parse(text = labels),
                                cex = 1.2
                                ),
                            padding.text = 2.2
                            )
                        ),
                    x = 0.02,
                    y = 0.95
                    )
                );
        } else {
            this.legend <- NULL;
            }

        this.plot[[group]] <- create.histogram(
            x = this.data[DiscordanceStatus == group, StabilityScore],
            xaxis.lab = if (group == 'Discordant') TRUE else NULL,
            xlab.lab = if (group == 'Discordant') 'Stability Score' else NULL,
            xlab.cex = 1.5,
            ylab.cex = 1.5,
            yaxis.cex = 1.5,
            xaxis.tck = c(1, 0),
            breaks = bins,
            yat = this.axis$at,
            yaxis.lab = this.axis$axis.lab,
            col = colors[group],
            type = 'count',
            legend = this.legend,
            height = 3,
            width = 5
            );
        }

    create.multipanelplot(
        filename = this.filename,
        plot.objects = this.plot,
        ylab.label = expression(bold('    Variant Count')),
        ylab.cex = 1.7,
        ylab.axis.padding = -5.5,
        xlab.axis.padding = 0.5,
        left.padding = 0.3,
        top.padding = -2.1,
        bottom.padding = -1.7,
        right.padding = -1.5,
        y.spacing = -1.5,
        height = 5.5,
        width = 5
        );
    }

###################################################################################################
# Output stability scores
###################################################################################################
this.filename <- sub('_LiftOver-GRCh38_annotated.Rds', '_stableLift-GRCh38_stability-scores.tsv', features.dt.path);
annotation.dt <- data.table(
    CHROM = features.dt$CHROM,
    POS = features.dt$POS,
    STABILITY_SCORE = format(round(stability$predictions[, 1], 4), nsmall = 4),
    STABILITY = ifelse(stability.classification == '1', 'UNSTABLE', 'STABLE')
    );
fwrite(annotation.dt, file = this.filename, sep = '\t', col.names = TRUE);
