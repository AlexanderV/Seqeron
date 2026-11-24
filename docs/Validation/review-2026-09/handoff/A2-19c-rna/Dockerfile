# TEMPORARY hand-off image A2-19c-rna: MAFFT 7.505 built WITH its extensions (MXSCARNA modified for MAFFT: `mxscarnamod`, and the other programs of the `extensions` directory) plus ViennaRNA's RNAalifold, so that
# `mafft --mccaskill`, `--scarnapair`, `--rnaalifold`, `--dafs`, `--contrafold` (the RNA options of the script mafft that the Debian binary package cannot run: "Install MXSCARNA", "RNAalifold: not found") can be run.
# Every step that can fail is logged and never fails the build. Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12-slim
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl g++ make zlib1g-dev python3 zip file procps mafft vienna-rna ruby \
 && rm -rf /var/lib/apt/lists/* \
 && mafft --version 2>&1 | head -1 ; RNAalifold --version 2>&1 | head -1

COPY build_ext.sh /opt/build_ext.sh
RUN chmod +x /opt/build_ext.sh && /opt/build_ext.sh > /opt/build_ext.log 2>&1 ; tail -30 /opt/build_ext.log; ls -l /opt/mafft-ext/libexec/mafft 2>&1 | head -40
WORKDIR /work
