# TEMPORARY hand-off image A2-19c-fastapair2 (second run: more DNA sets): MAFFT 7.505-1 (Debian 12 "bookworm" binary package, the version the port is verified against) plus the program fasta34 that
# `mafft --fastapair` / `--fastswpair` call (dndfast7; the audit machine has only fasta36, which is a different program and cannot serve as the oracle).
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12-slim
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates gcc libc6-dev make python3 zip file procps mafft \
 && rm -rf /var/lib/apt/lists/* \
 && mafft --version 2>&1 | head -1

# fasta34: the FASTA 3.4 package of W. R. Pearson (release 34.26.5, inputs/fasta-34.26.5.shar.Z, parsed not executed by fetch_fasta34.py) and builds `fasta34` with the Makefile of the
# package (Makefile.linux).
COPY fetch_fasta34.py /opt/fetch_fasta34.py
COPY inputs /opt/inputs
RUN mkdir -p /opt/fasta34src && cd /opt/fasta34src && python3 -I /opt/fetch_fasta34.py /opt/inputs \
    > /opt/fasta34.build.log 2>&1 || { echo FASTA34_BUILD_FAILED > /opt/FASTA34_BUILD_FAILED; }; tail -3 /opt/fasta34.build.log
WORKDIR /work
