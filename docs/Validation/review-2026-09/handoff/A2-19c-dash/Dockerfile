# TEMPORARY hand-off image A2-19c-dash: the MAFFT 7.505 source package (upstream, "without extensions") built WITH the Go program dash_client (core/Makefile: DASH_CLIENT = dash_client), which the Debian
# binary package does not contain. Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; every download uses the neutral User-Agent "Seqeron-review/1.0".
# The tarball is the official one (https://mafft.cbrc.jp/alignment/software/); its SHA-256 is recorded in out/versions.txt (the audit machine's copy: 9b9ee0dea1d493f0fad6f9c227bcc291ddedd045b5558b7fea30dc5932db33b8,
# a different value is reported, not fatal).
ARG BASE=debian:12-slim
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl gcc libc6-dev make golang-go python3 zip file procps gawk \
 && rm -rf /var/lib/apt/lists/* \
 && go version

ARG MAFFT_URL=https://mafft.cbrc.jp/alignment/software/mafft-7.505-without-extensions-src.tgz
RUN mkdir -p /opt/src && cd /opt/src \
 && curl -fsSL -A "Seqeron-review/1.0" -o mafft-src.tgz "$MAFFT_URL" \
 && sha256sum mafft-src.tgz | tee /opt/mafft-src.sha256 \
 && tar xzf mafft-src.tgz \
 && cd mafft-7.505-without-extensions/core \
 && sed -i 's/^#DASH_CLIENT = dash_client/DASH_CLIENT = dash_client/' Makefile \
 && grep -n "^DASH_CLIENT" Makefile \
 && make ENABLE_MULTITHREAD= PREFIX=/usr/local BINDIR=/usr/local/bin LIBDIR=/usr/local/libexec/mafft > /opt/mafft.build.log 2>&1 \
 && make ENABLE_MULTITHREAD= PREFIX=/usr/local BINDIR=/usr/local/bin LIBDIR=/usr/local/libexec/mafft install >> /opt/mafft.build.log 2>&1 \
 ; tail -5 /opt/mafft.build.log; ls -l /usr/local/libexec/mafft/dash_client || { echo DASH_CLIENT_BUILD_FAILED > /opt/DASH_CLIENT_BUILD_FAILED; }
WORKDIR /work
