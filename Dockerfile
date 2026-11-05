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

# Release tag expected from the build (e.g. v2.6.0), checked against `phyling --version` below; empty skips the
# check (manual builds of a branch). versioningit itself derives the version from the .git directory in the context.
ARG PHYLING_VERSION=""

# Runtime environment: python, pip/git (needed to install and version phyling
# from this checkout), the external phylogenetics binaries phyling shells
# out to for alignment and tree inference, and procps-ng (provides `ps`, which
# Nextflow runs inside the container to collect task metrics).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.13 \
        pip \
        git \
        "aster>=1.19" \
        "fasttree>=2.1.1" \
        "iqtree>=3.1.2" \
        muscle \
        "procps-ng" \
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
    command -v ps; \
    command -v raxml-ng; \
    command -v iqtree || command -v iqtree2; \
    phyling --help; \
    if [ -n "${PHYLING_VERSION}" ]; then \
        python -c 'import sys; from importlib.metadata import version; from pip._vendor.packaging.version import Version; v = version("phyling"); sys.exit(0 if Version(v) == Version(sys.argv[1]) else f"phyling reports {v}, expected {sys.argv[1]}")' "${PHYLING_VERSION}"; \
    fi

USER root
# /etc/profile resets PATH in login shells, so re-add the environment there too
RUN rm -rf /src && mkdir -p /data && chown $MAMBA_USER:$MAMBA_USER /data \
    && echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/phyling.sh
USER $MAMBA_USER
WORKDIR /data

# Put the environment on PATH without the entrypoint, for runtimes that bypass it: `apptainer exec`, `docker run
# --entrypoint` and Nextflow/SLURM task scripts (`bash -l`).
ENV PATH=/opt/conda/bin:$PATH

# _entrypoint.sh (from the base image) activates the micromamba env before
# exec'ing the given command, so `phyling` resolves from PATH at run time.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "phyling"]
CMD ["--help"]
