ARG R_VERSION=4.3.1

FROM rocker/r-ver:${R_VERSION} AS build

COPY docker/install-stablelift.R /tmp
RUN Rscript /tmp/install-stablelift.R

FROM rocker/r-ver:${R_VERSION}

# Overwrite the site library with just the desired packages. By default rocker
# only bundles docopt and littler in that directory.
COPY --from=build /tmp/userlib /usr/local/lib/R/site-library

# Install python (required for argparse). The version is not important, but
# let's pin it for stability.
ARG PYTHON_VERSION=3.10.6-1~22.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3=${PYTHON_VERSION} \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Add a new user/group called bldocker
RUN groupadd -g 500001 bldocker && \
    useradd -l -r -u 500001 -g bldocker bldocker

# Change the default user to bldocker from root
USER bldocker

LABEL maintainer="Nicholas Wiltsie <nwiltsie@mednet.ucla.edu" \
      org.opencontainers.image.source=https://github.com/uclahs-cds/pipeline-StableLift
