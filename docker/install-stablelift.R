install.packages('renv', lib = .Library)

options(
    renv.settings.bioconductor.version = Sys.getenv("BIOC_VERSION")
)

renv::init(
    project = "/tmp/stablelift",
    bare = TRUE,
    bioconductor = Sys.getenv("BIOC_VERSION")
)

renv::install(c(
    'ROCR@1.0-11',
    'argparse@2.2.2',
    'caret@6.0-94',
    'data.table@1.14.8',
    'doParallel@1.0.17',
    'foreach@1.5.2',
    'ranger@0.15.1',
    'vcfR@1.14.0',
    'bioc::rtracklayer@1.62.0'
))
