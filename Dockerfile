ARG R_VERSION=4.3.1

# FROM rocker/tidyverse:${R_VERSION} AS build
FROM rocker/r-ver:${R_VERSION} AS build

RUN R --quiet --no-save <<EOF
# Install the remotes package to the library
install.packages('remotes', lib = .Library)

# Make a temporary directory to hold all of the installed packages
localdir <- '/tmp/userlib'
dir.create(localdir)

dependencies <- c(
    'argparse' = '2.2.2',
    'caret' = '6.0-94',
    'data.table' = '1.14.8',
    'doParallel' = '1.0.17',
    'foreach' = '1.5.2',
    'ranger' = '0.15.1',
    'vcfR' = '1.14.0'
)

# Unfortunately, this will install the dependencies multiple times
for (name in names(dependencies)) {
    remotes::install_version(
        name,
        unname(dependencies[name]),
        lib = localdir
    )
}
EOF

FROM rocker/r-ver:${R_VERSION}
COPY --from=build /tmp/userlib /usr/local/lib/R/site-library

LABEL maintainer="Nicholas Wiltsie <nwiltsie@mednet.ucla.edu" \
      org.opencontainers.image.source=https://github.com/uclahs-cds/pipeline-StableLift
