# TEMPORARY hand-off image: NCBI ORFfinder 0.4.3 (official Linux x86-64 binary) on a pinned Debian 12.7 slim base.
# The binary is x86-64 only: build and run with `--platform linux/amd64` (Apple Silicon: Docker Desktop emulates it).
# No account data, no secrets; every download uses the neutral User-Agent "Seqeron-review/1.0".
#
# Base: debian:12.7-slim (= bookworm-20241016-slim), linux/amd64 manifest pinned by digest.
ARG BASE=debian:12.7-slim@sha256:d83056144b2dd301730d2739635c8cbdeaaae20d6887146434184f8c060f03ce
FROM --platform=linux/amd64 ${BASE}

# Official NCBI FTP location of the standalone program (USAGE.txt documents the options, CHANGELOG.txt the versions).
ARG ORFFINDER_URL=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
ARG USAGE_URL=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/USAGE.txt
ARG CHANGELOG_URL=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/CHANGELOG.txt
# sha256 of the files as downloaded on 2026-11-18 (ORFfinder.gz Last-Modified 2021-05-04, 8386303 bytes).
# If NCBI ever replaces the file the build stops here on purpose: re-run with --build-arg ORFFINDER_GZ_SHA256=<new value>
# and say so when sending the results back.
ARG ORFFINDER_GZ_SHA256=f30a0806b984931a967d8024066a53dbb811035ca14564f36f7f7045a6a509f2
ARG USAGE_SHA256=368d0fe3c98b3383f04747b17ce2cd57031fa8e0ece4e56ef139ee0df8f7b409
ARG CHANGELOG_SHA256=a36fb2a240f9510b7a42ca71ed5df8256c7feb03b8cf7938f79be0f481724316

# Runtime libraries of the binary (readelf -d: libz, libbz2, libuv, libnghttp2, libdw, libstdc++, libgcc_s, libc),
# python3 for the committed input / job generators, zip for the result archive.
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ca-certificates curl python3-minimal zip file procps \
        zlib1g libbz2-1.0 libuv1 libnghttp2-14 libdw1 libstdc++6 libgcc-s1 \
 && rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    mkdir -p /opt/orffinder; cd /opt/orffinder; \
    curl -fsSL -A "Seqeron-review/1.0" -o ORFfinder.gz "$ORFFINDER_URL"; \
    curl -fsSL -A "Seqeron-review/1.0" -o USAGE.txt "$USAGE_URL"; \
    curl -fsSL -A "Seqeron-review/1.0" -o CHANGELOG.txt "$CHANGELOG_URL"; \
    echo "$ORFFINDER_GZ_SHA256  ORFfinder.gz" | sha256sum -c -; \
    echo "$USAGE_SHA256  USAGE.txt" | sha256sum -c -; \
    echo "$CHANGELOG_SHA256  CHANGELOG.txt" | sha256sum -c -; \
    gzip -dc ORFfinder.gz > ORFfinder; chmod 0755 ORFfinder; \
    sha256sum ORFfinder.gz ORFfinder USAGE.txt CHANGELOG.txt > SHA256SUMS; \
    ldd ./ORFfinder | tee ldd.txt; \
    if grep -q "not found" ldd.txt; then echo "missing shared library for ORFfinder" >&2; exit 1; fi

ENV PATH=/opt/orffinder:$PATH
WORKDIR /work
CMD ["/work/in_container.sh"]
