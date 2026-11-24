# TEMPORARY hand-off image A2-19c-last: MAFFT 7.505-1 (Debian 12 "bookworm" binary package, the version the port is verified against) plus the LAST programs lastdb / lastal that `mafft --lastpair` / `--lastmultipair` call (pairlocalalign -r / -R, "Calling last").
# Several LAST versions are tried: the Debian package (last-align) and release sources from gitlab.com/mcfrith/last (built in the image); a smoke test runs the exact mafft call on a small DNA set for each and
# /opt/last/bin points to the first that answers within 60 seconds (the log of every candidate is kept in /opt/last-candidates.txt and copied to out/). The step never fails the build.
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12-slim
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl g++ make zlib1g-dev python3 zip file procps mafft \
 && rm -rf /var/lib/apt/lists/* \
 && mafft --version 2>&1 | head -1

COPY last-log /opt/last-log
COPY build_last.sh /opt/build_last.sh
COPY smoke.fa /opt/smoke.fa
ARG LAST_VERSIONS="1256 1250 1282 1411 1542"
RUN chmod +x /opt/last-log /opt/build_last.sh && LAST_VERSIONS="$LAST_VERSIONS" /opt/build_last.sh > /opt/build_last.log 2>&1 ; tail -25 /opt/build_last.log; ls -l /opt/last 2>&1 | head
WORKDIR /work
