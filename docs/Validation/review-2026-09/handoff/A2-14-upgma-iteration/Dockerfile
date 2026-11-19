# TEMPORARY hand-off image A2-14-upgma-iteration: ClustalW 2.1 (Debian source package clustalw 2.1+lgpl-7, src/ of the unmodified upstream sources) built at -O0 with g++.
# WHY: the installed Debian binary /usr/bin/clustalw (2.1+lgpl-7) and an -O2 build of the same source SEGMENTATION-FAULT on every `-CLUSTERING=UPGMA` run; the -O0 build runs
# (handoff A2-14-clustalw-upgma: 600 / 600 sets equal to the port). This package runs the UPGMA code paths that package did not reach: `-CLUSTERING=UPGMA` together with
# `-ITERATION=alignment|tree` (TreeInterface::getWeightsForProfileAlignUPGMA) and with the profile commands `-PROFILE` / `-SEQUENCES`. No patch of any kind.
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12.7-slim@sha256:d83056144b2dd301730d2739635c8cbdeaaae20d6887146434184f8c060f03ce
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates g++ make python3 zip file procps \
 && rm -rf /var/lib/apt/lists/*

COPY fetch_src.py /opt/fetch_src.py
RUN mkdir -p /opt/cwsrc && cd /opt/cwsrc && python3 -I /opt/fetch_src.py

# The build never aborts the image: its log goes to /opt/cw/build_O0.log, a failed build leaves /opt/cw/BUILD_FAILED_O0 (in_container.sh then returns the log).
RUN mkdir -p /opt/cw; cd /opt/cwsrc/src; \
    INC="-I. -Igeneral -Ialignment -Itree -Itree/UPGMA -IsubstitutionMatrix -Ipairwise -ImultipleAlign -IfileInput -Iinterface -DHAVE_CONFIG_H -DOS_UNIX=1"; \
    ( g++ -std=gnu++98 -O0 -g -w $INC -o /opt/cw/clustalw_O0 $(find . -name '*.cpp') ) > /opt/cw/build_O0.log 2>&1 \
      || echo BUILD_FAILED > /opt/cw/BUILD_FAILED_O0; \
    sha256sum /opt/cw/clustalw_* > /opt/cw/binaries.sha256 2>/dev/null; cp /opt/cwsrc/fetch_src.sha256 /opt/cw/; ls -l /opt/cw
WORKDIR /work
