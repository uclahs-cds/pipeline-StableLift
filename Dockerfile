ARG R_VERSION="4.3.1"

ARG LIBBZ2_VERSION="1.0.8-5build1"
ARG LIBCURL_VERSION="7.81.0-1ubuntu1.16"
ARG LIBLZMA_VERSION="5.2.5-2ubuntu1"
ARG LIBXML2_VERSION="2.9.13+dfsg-1ubuntu0.4"
ARG PYTHON_VERSION="3.10.6-1~22.04"
ARG ZLIB_VERSION="1:1.2.11.dfsg-2ubuntu9.2"
ARG RLIBDIR="/usr/local/stablelift-R"

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

ARG RLIBDIR
ENV RENV_PATHS_CACHE ${RLIBDIR}/.cache

RUN mkdir -p ${RENV_PATHS_CACHE}

WORKDIR ${RLIBDIR}

COPY docker/install-stablelift.R /tmp
RUN Rscript /tmp/install-stablelift.R

# renv prints to stdout, so we need to change directories
WORKDIR /
RUN echo ".libPaths( c( .libPaths(), \"/usr/local/stablelift-R/renv/library/R-4.3/$(Rscript -e "cat(unname(unlist(R.version['platform'])))")\" ) )" >> /usr/local/lib/R/etc/Rprofile.site

FROM rocker/r-ver:${R_VERSION}

ARG RLIBDIR
COPY --from=build ${RLIBDIR} ${RLIBDIR}
COPY --from=build \
    /usr/local/lib/R/etc/Rprofile.site \
    /usr/local/lib/R/etc/Rprofile.site

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
