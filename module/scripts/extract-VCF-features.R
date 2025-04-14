#!/usr/bin/env Rscript
# extract-VCF-features.R
####################################################################################################
#
# Extract VCF features and save as Rds for input to predict-variant-stability.R
#
####################################################################################################

suppressPackageStartupMessages({
    library(vcfR);
    library(data.table);
    library(argparse);
    library(doParallel);
    library(foreach);
    });

###################################################################################################
# Input
###################################################################################################
# Define command line arguments
parser <- ArgumentParser();
parser$add_argument('--input-vcf', type = 'character', help = 'Input VCF for feature extraction, mutually exclusive with --input-dir');
parser$add_argument('--input-dir', type = 'character', help = 'Directory with VCF subsets for parallelization, mutually exclusive with --input-vcf');
parser$add_argument('--output-rds', type = 'character', help = 'Rds output for input to RF model');
parser$add_argument('--variant-caller', type = 'character', help = 'One of {HaplotypeCaller, Mutect2, Strelka2, SomaticSniper, Muse2}');
parser$add_argument('--ncore', type = 'integer', help = 'Number of cores to use for processing VCF subsets in --input-dir', default = 1);
args <- parser$parse_args();

# Save command line arguments
for (arg in names(args)) {
    assign(gsub('_', '.', arg), args[[arg]]);
    }

if (!is.null(input.dir)) {
    vcf.subsets <- list.files(input.dir, full.names = TRUE, pattern = '(\\.vcf.gz|\\.vcf)$');
    output.path <- output.rds;
} else if (!is.null(input.vcf)) {
    vcf.subsets <- input.vcf;
    output.path <- ifelse(is.null(output.rds), sub('\\..*$', '.Rds', input.vcf), output.rds);
} else {
    stop('Either `--input-vcf` or `--input-dir` must be provided');
    }

###################################################################################################
# Functions
###################################################################################################
vcf.info.to.dt <- function(vcf.info) {
    # Split each string by semicolon and convert to a list of key-value pairs
    vcf.info <- strsplit(vcf.info, ';');
    vcf.info <- lapply(vcf.info, function(x) {
        x <- strsplit(x, '=');
        as.list(stats::setNames(
            sapply(x, function(pair) if (length(pair) == 1) pair[1] else pair[2]),
            sapply(x, `[`, 1)
            ));
        })

    # Combine the list of key-value pairs into a data table
    rbindlist(vcf.info, fill = TRUE);
    }

sum.MAF <- function(x) {
    x <- unlist(strsplit(x, ','));
    x <- 1 - as.numeric(x[1]);
    if (is.na(x)) 0 else x;
    }

###################################################################################################
# Extract RF features from Funcotator vcf
###################################################################################################
start.extract <- Sys.time();

registerDoParallel(cores = ncore);
features.dt.subsets <- foreach(vcf.subset = vcf.subsets) %dopar% {
    start.subset <- Sys.time();
    cat('Processing:', basename(vcf.subset), '\n');

    input.vcf <- read.vcfR(
        file = vcf.subset
        );

    # Parse info column
    info <- vcf.info.to.dt(input.vcf@fix[, 'INFO']);

    # Aggregate per sample metrics, average depth per GT called
    info$DP <- apply(extract.gt(input.vcf, element = 'DP', as.numeric = TRUE), 1, mean, na.rm = TRUE);
    if (variant.caller == 'HaplotypeCaller') {
        # Take mean of genotype quality and remove cohort allele frequency
        info$GQ <- apply(extract.gt(input.vcf, element = 'GQ', as.numeric = TRUE), 1, mean, na.rm = TRUE);
        info[, AF := NULL];
    } else if (variant.caller == 'Mutect2') {
        # Get metric for ALT allele
        info$MBQ <- sapply(extract.info(input.vcf, element = 'MBQ'), function(x) as.numeric(unlist(strsplit(x, ','))[2]));
        info$MFRL <- sapply(extract.info(input.vcf, element = 'MFRL'), function(x) as.numeric(unlist(strsplit(x, ','))[2]));
        info$MMQ <- sapply(extract.info(input.vcf, element = 'MMQ'), function(x) as.numeric(unlist(strsplit(x, ','))[2]));
        # Take max of somatic variant allele frequencies
        info$AF <- apply(extract.gt(input.vcf, element = 'AF', as.numeric = TRUE), 1, max, na.rm = TRUE);
    } else if (variant.caller == 'Muse2') {
        # Calculate VAF from allelic depths
        info$AF <- apply(extract.gt(input.vcf, element = 'AD'), 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[2]) / (as.numeric(y[1]) + as.numeric(y[2]))), na.rm = TRUE));
        info$BQ <- apply(extract.gt(input.vcf, element = 'BQ'), 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[2])), na.rm = TRUE));
    } else if (variant.caller == 'Strelka2') {
        # Calculate VAF from allelic depths
        info[input.vcf@fix[, 'REF'] == 'A', REFCOUNTS := apply(extract.gt(input.vcf, element = 'AU')[input.vcf@fix[, 'REF'] == 'A', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'ALT'] == 'A', ALTCOUNTS := apply(extract.gt(input.vcf, element = 'AU')[input.vcf@fix[, 'ALT'] == 'A', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'REF'] == 'T', REFCOUNTS := apply(extract.gt(input.vcf, element = 'TU')[input.vcf@fix[, 'REF'] == 'T', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'ALT'] == 'T', ALTCOUNTS := apply(extract.gt(input.vcf, element = 'TU')[input.vcf@fix[, 'ALT'] == 'T', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'REF'] == 'C', REFCOUNTS := apply(extract.gt(input.vcf, element = 'CU')[input.vcf@fix[, 'REF'] == 'C', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'ALT'] == 'C', ALTCOUNTS := apply(extract.gt(input.vcf, element = 'CU')[input.vcf@fix[, 'ALT'] == 'C', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'REF'] == 'G', REFCOUNTS := apply(extract.gt(input.vcf, element = 'GU')[input.vcf@fix[, 'REF'] == 'G', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[input.vcf@fix[, 'ALT'] == 'G', ALTCOUNTS := apply(extract.gt(input.vcf, element = 'GU')[input.vcf@fix[, 'ALT'] == 'G', ], 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[1])), na.rm = TRUE))];
        info[, AF := ALTCOUNTS / (REFCOUNTS + ALTCOUNTS)];
    } else if (variant.caller == 'SomaticSniper') {
        # Calculate VAF from allelic depths
        info$AF <- apply(extract.gt(input.vcf, element = 'DP4'), 1, function(x) mean(
            sapply(strsplit(x, ','), function(y) {
                y <- as.numeric(y);
                return((y[3] + y[4]) / sum(y));
                }),
            na.rm = TRUE
            ));
        info$AMQ <- apply(extract.gt(input.vcf, element = 'AMQ'), 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[2])), na.rm = TRUE));
        info$BQ <- apply(extract.gt(input.vcf, element = 'BQ'), 1, function(x) mean(sapply(strsplit(x, ','), function(y) as.numeric(y[2])), na.rm = TRUE));
        info$GQ <- apply(extract.gt(input.vcf, element = 'GQ', as.numeric = TRUE), 1, mean, na.rm = TRUE);
        info$MQ <- apply(extract.gt(input.vcf, element = 'GQ', as.numeric = TRUE), 1, mean, na.rm = TRUE);
        info$SSC <- apply(extract.gt(input.vcf, element = 'SSC', as.numeric = TRUE), 1, mean, na.rm = TRUE);
        info$VAQ <- apply(extract.gt(input.vcf, element = 'VAQ', as.numeric = TRUE), 1, mean, na.rm = TRUE);
        }

    # Get funcotation fields
    funcotator.fields <- grep('ID=FUNCOTATION', input.vcf@meta, value = TRUE);
    funcotator.fields <- sub('.*Funcotation fields are: (.*)\\".*', '\\1', funcotator.fields);
    funcotator.fields <- strsplit(funcotator.fields, split = '\\|')[[1]];

    # Parse funcotation field from info column
    funcotation <- sapply(info$FUNCOTATION, function(x) unlist(strsplit(x, split = ','))[1]);
    names(funcotation) <- NULL;
    funcotation <- gsub('\\[|\\]', '', funcotation);
    funcotation <- lapply(funcotation, function(x) {
        x <- unlist(strsplit(x, split = '\\|'));
        x <- as.list(x[seq_along(funcotator.fields)]);
        names(x) <- funcotator.fields;
        return(x);
        });
    features.dt <- rbindlist(funcotation, use.names = TRUE, fill = TRUE);
    rm(funcotation);

    # names(features.dt) <- funcotator.fields;
    funcotator.fields.keep <- c(
        'GC Content' = 'Gencode_34_gcContent',
        'GENCODE Variant Classification' = 'Gencode_34_variantClassification',
        'GENCODE Variant Type' = 'Gencode_34_variantType',
        'HGNC Locus Group' = 'HGNC_Locus_Group',
        '1KG AF' = 'dbSNP_CAF',
        'TOPMED AF' = 'dbSNP_TOPMED'
        );
    features.dt <- features.dt[, ..funcotator.fields.keep];

    # Replace ascii character codes
    features.dt <- as.data.table(lapply(features.dt, function(x) gsub('_%7C_', ';', x, fixed = TRUE)));
    features.dt <- as.data.table(lapply(features.dt, function(x) gsub('_%2C_', ',', x, fixed = TRUE)));
    features.dt <- as.data.table(lapply(features.dt, function(x) gsub('_%20_', ' ', x, fixed = TRUE)));
    features.dt <- as.data.table(lapply(features.dt, function(x) gsub('_%3D_', '=', x, fixed = TRUE)));

    # Aggregate minor allele frequencies
    features.dt[, dbSNP_CAF := sapply(dbSNP_CAF, sum.MAF)];
    features.dt[, dbSNP_TOPMED := sapply(dbSNP_TOPMED, sum.MAF)];

    # Aggregate flattened fields
    features.dt <- cbind(input.vcf@fix[, 1:7], info[, -c('FUNCOTATION')], features.dt);

    # Collapse trinucleotide reverse complements
    features.dt[Gencode_34_variantType != 'SNP', TRINUCLEOTIDE := ''];
    features.dt[Gencode_34_variantType == 'SNP', TRINUCLEOTIDE := sapply(TRINUCLEOTIDE, function(x) {
        x <- as.character(x);
        if (substr(x, 2, 2) %in% c('A', 'G')) {
            x <- chartr('ATGC', 'TACG', x);
            paste(rev(unlist(strsplit(x, ''))), collapse = '');
        } else {
            x;
            }
        })];

    cat(format(Sys.time() - start.subset, nsmall = 2), '\n');
    return(features.dt);
    }

features.dt <- rbindlist(features.dt.subsets, use.names = TRUE, fill = TRUE);

# Remove variants with identical CHROM and POS
features.dt <- features.dt[!duplicated(features.dt[, .(CHROM, POS)])];

rm(features.dt.subsets);
gc();
cat(format(Sys.time() - start.extract, nsmall = 2), '\n');

###################################################################################################
# Format features for RF
###################################################################################################
continuous.features <- c(
    'Chromosome Position (POS)' = 'POS',
    'Local GC Content' = 'Gencode_34_gcContent',
    '1k Genomes Pop AF' = 'dbSNP_CAF',
    'TOPMED Pop AF' = 'dbSNP_TOPMED',
    'Sequencing Depth (DP)' = 'DP',
    'Variant Frequency (AF)' = 'AF'
    );
if (variant.caller == 'HaplotypeCaller') continuous.features <- c(continuous.features,
    'Variant Quality (QUAL)' = 'QUAL',
    'ExcessHet' = 'ExcessHet',
    "Fisher's Strand Bias (FS)" = 'FS',
    'Mapping Quality (MQ)' = 'MQ',
    'Quality by Depth (QD)' = 'QD',
    'Strand Bias Odds Ratio (SOR)' = 'SOR',
    'VQSR Log Odds (VQSLOD)' = 'VQSLOD',
    'Genotype Quality (GQ)' = 'GQ'
    );
if (variant.caller == 'Mutect2') continuous.features <- c(continuous.features,
    'Tumor Variant Log Odds (TLOD)' = 'TLOD',
    'Events in Haplotype (ECNT)' = 'ECNT',
    'Germline Quality (GERMQ)' = 'GERMQ',
    'Median Distance to Read End (MPOS)' = 'MPOS',
    'Normal Artifact Log Odds (NALOD)' = 'NALOD',
    'Normal Variant Log Odds (NLOD)' = 'NLOD',
    'Read Orientation Artifact (ROQ)' = 'ROQ',
    'Median Base Quality (MBQ)' = 'MBQ',
    'Median Fragment Length (MFRL)' = 'MFRL',
    'Median Mapping Quality (MMQ)' = 'MMQ'
    );
if (variant.caller == 'Strelka2') continuous.features <- c(continuous.features,
    'Mapping Quality (MQ)' = 'MQ',
    'MQ0 Reads (MQ0)' = 'MQ0',
    'Empirical Variant Score (SomaticEVS)' = 'SomaticEVS',
    'Quality Score Somatic (QSS)' = 'QSS',
    'Quality Score Joint (QSS_NT)' = 'QSS_NT',
    'ReadPosRankSum'
    );
if (variant.caller == 'Muse2') continuous.features <- c(continuous.features,
    'Variant Base Quality (BQ)' = 'BQ'
    );
if (variant.caller == 'SomaticSniper') continuous.features <- c(continuous.features,
    'Variant Mapping Quality (AMQ)' = 'AMQ',
    'Base Quality (BQ)' = 'BQ',
    'Genotype Quality (GQ)' = 'GQ',
    'Mapping Quality (MQ)' = 'MQ',
    'Somatic Score (SSC)' = 'SSC',
    'Variant Allele Quality (VAQ)' = 'VAQ'
    );

categorical.features <- c(
    'Chromosome (CHR)' = 'CHROM',
    'GENCODE Variant Type' = 'Gencode_34_variantType',
    'GENCODE Variant Type' = 'Gencode_34_variantClassification',
    'HGNC Locus Group' = 'HGNC_Locus_Group',
    'Repeat Class' = 'REPEAT',
    'Trinucleotide Context' = 'TRINUCLEOTIDE',
    'VQSR Culprit' = 'culprit',
    'Somatic Genotype (SGT)' = 'SGT',
    'LiftOver Allele Flip' = 'FLIP',
    'LiftOver Allele Swap' = 'SWAP'
    );

# Extract features and format
continuous.features <- continuous.features[continuous.features %in% names(features.dt)];
categorical.features <- categorical.features[categorical.features %in% names(features.dt)];
all.features <- c(continuous.features, categorical.features);

features.dt <- features.dt[, ..all.features];
features.dt[, (continuous.features) := lapply(.SD, as.numeric), .SDcols = continuous.features];
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
saveRDS(features.dt, output.path);
