ARG R_VERSION="4.3.1"

ARG LIBBZ2_VERSION="1.0.8-5build1"
ARG LIBCURL_VERSION="7.81.0-1ubuntu1.16"
ARG LIBLZMA_VERSION="5.2.5-2ubuntu1"
ARG LIBXML2_VERSION="2.9.13+dfsg-1ubuntu0.4"
ARG PYTHON_VERSION="3.10.6-1~22.04"
ARG ZLIB_VERSION="1:1.2.11.dfsg-2ubuntu9.2"

FROM rocker/r-ver:${R_VERSION} AS build

ARG LIBBZ2_VERSION
ARG LIBCURL_VERSION
ARG LIBLZMA_VERSION
ARG LIBXML2_VERSION
ARG ZLIB_VERSION

# Install build-time dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libbz2-dev=${LIBBZ2_VERSION} \
        libcurl4-openssl-dev=${LIBCURL_VERSION} \
        liblzma-dev=${LIBLZMA_VERSION} \
        libxml2-dev=${LIBXML2_VERSION} \
        zlib1g-dev=${ZLIB_VERSION} \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ARG BIOC_VERSION="3.18"
ENV BIOC_VERSION=${BIOC_VERSION}

COPY docker/install-stablelift.R /tmp
RUN Rscript /tmp/install-stablelift.R

FROM rocker/r-ver:${R_VERSION}

# Overwrite the site library with just the desired packages. By default rocker
# only bundles docopt and littler in that directory.
COPY --from=build \
    /tmp/stablelift/renv/library/R-4.3/aarch64-unknown-linux-gnu \
    /usr/local/lib/R/site-library

# Install python (required for argparse). The version is not important, but
# let's pin it for stability.
ARG PYTHON_VERSION

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
