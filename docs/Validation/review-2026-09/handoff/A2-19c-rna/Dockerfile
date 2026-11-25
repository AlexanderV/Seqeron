# TEMPORARY hand-off image A2-19c-rna: MAFFT 7.505 built WITH its extensions (MXSCARNA modified for MAFFT: `mxscarnamod`, and the other programs of the `extensions` directory) plus ViennaRNA's RNAalifold, so that
# `mafft --mccaskill`, `--scarnapair`, `--rnaalifold`, `--dafs`, `--contrafold` (the RNA options of the script mafft that the Debian binary package cannot run: "Install MXSCARNA", "RNAalifold: not found") can be run.
# RUN 1 FAILED because `vienna-rna` has no installation candidate on debian:12 and that aborted the single `apt-get install`: the base packages are now installed one by one (a missing package is logged, not fatal), ViennaRNA is built from
# the upstream release (build_vienna.sh), the extensions by build_ext.sh. Every step that can fail is logged and never fails the build. Build and run with `--platform linux/amd64` (Apple Silicon: emulated).
# No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12-slim
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update; \
    for p in ca-certificates curl g++ gcc make zlib1g-dev python3 zip file procps mafft ruby; do \
      DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends $p > /tmp/apt-$p.log 2>&1 && echo "installed $p" || { echo "NOT INSTALLED: $p"; tail -3 /tmp/apt-$p.log; }; \
    done; \
    rm -rf /var/lib/apt/lists/*; \
    mafft --version 2>&1 | head -1; true

COPY build_vienna.sh /opt/build_vienna.sh
RUN chmod +x /opt/build_vienna.sh && (apt-get update > /dev/null 2>&1; /opt/build_vienna.sh) > /opt/build_vienna.log 2>&1 ; tail -30 /opt/build_vienna.log; rm -rf /var/lib/apt/lists/*; true

COPY build_ext.sh /opt/build_ext.sh
RUN chmod +x /opt/build_ext.sh && /opt/build_ext.sh > /opt/build_ext.log 2>&1 ; tail -30 /opt/build_ext.log; ls -l /opt/mafft-ext/libexec/mafft 2>&1 | head -40; true
WORKDIR /work
