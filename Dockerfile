# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:2

LABEL org.opencontainers.image.title="phyling" \
      org.opencontainers.image.description="A lightweight phylogenetic tree builder from annotated genomes." \
      org.opencontainers.image.source="https://github.com/stajichlab/Phyling" \
      org.opencontainers.image.licenses="MIT"

USER root
RUN mkdir -p /src && chown $MAMBA_USER:$MAMBA_USER /src
WORKDIR /src
COPY --chown=$MAMBA_USER:$MAMBA_USER . .
USER $MAMBA_USER

# Runtime environment: python, pip/git (needed to install and version phyling
# from this checkout), and the external phylogenetics binaries phyling shells
# out to for alignment and tree inference.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.13 \
        pip \
        git \
        "aster>=1.19" \
        "fasttree>=2.1.1" \
        "iqtree>=3.1.2" \
        muscle \
        raxml \
        "raxml-ng>=2.0.2" \
    && micromamba clean --all --yes

# versioningit derives the package version from git metadata, so keep the .git
# directory available at build time. git refuses to read a repo owned by
# another user by default, so mark it safe for $MAMBA_USER first.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN git config --global --add safe.directory /src \
    && pip install --no-cache-dir .

# Fail the build early if any binary phyling actually looks for is missing or
# was installed under an unexpected name (e.g. bioconda iqtree3 vs iqtree/iqtree2).
RUN set -eux; \
    command -v astral; \
    command -v FastTree; \
    command -v muscle; \
    command -v raxml-ng; \
    command -v iqtree || command -v iqtree2; \
    phyling --help

USER root
RUN rm -rf /src && mkdir -p /data && chown $MAMBA_USER:$MAMBA_USER /data
USER $MAMBA_USER
WORKDIR /data

# _entrypoint.sh (from the base image) activates the micromamba env before
# exec'ing the given command, so `phyling` resolves from PATH at run time.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "phyling"]
CMD ["--help"]
