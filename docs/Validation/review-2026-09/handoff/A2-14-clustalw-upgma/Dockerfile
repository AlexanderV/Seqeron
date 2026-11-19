# TEMPORARY hand-off image A2-14-clustalw-upgma: ClustalW 2.1 (Debian source package clustalw 2.1+lgpl-7, src/ of the unmodified upstream sources) built TWICE from source with g++:
#   clustalw_O0   : -O0 -g          (the build that worked in the audit sandbox with -fsanitize=address)
#   clustalw_O2   : -O2             (what a distribution build looks like)
# WHY: the installed Debian binary /usr/bin/clustalw (2.1+lgpl-7) SEGMENTATION-FAULTS on every `-CLUSTERING=UPGMA -ALIGN` / `-TREE` run (2, 3, 4 sequences alike), so it cannot be the oracle of the
# UPGMA port; the same source compiled at -O0 with AddressSanitizer ran the UPGMA alignment to the end. No patch of any kind. Build and run with `--platform linux/amd64` (Apple Silicon: emulated).
# No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12.7-slim@sha256:d83056144b2dd301730d2739635c8cbdeaaae20d6887146434184f8c060f03ce
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates g++ make python3 zip file procps \
 && rm -rf /var/lib/apt/lists/*

COPY fetch_src.py /opt/fetch_src.py
RUN mkdir -p /opt/cwsrc && cd /opt/cwsrc && python3 -I /opt/fetch_src.py

# The build never aborts the image: its log goes to /opt/cw/build.log, a failed build leaves /opt/cw/BUILD_FAILED_<name> (in_container.sh then returns the log).
RUN mkdir -p /opt/cw; cd /opt/cwsrc/src; \
    INC="-I. -Igeneral -Ialignment -Itree -Itree/UPGMA -IsubstitutionMatrix -Ipairwise -ImultipleAlign -IfileInput -Iinterface -DHAVE_CONFIG_H -DOS_UNIX=1"; \
    for v in O0 O2; do \
      FL="-O2"; [ "$v" = O0 ] && FL="-O0 -g"; \
      ( g++ -std=gnu++98 $FL -w $INC -o /opt/cw/clustalw_$v $(find . -name '*.cpp') ) > /opt/cw/build_$v.log 2>&1 \
        || echo BUILD_FAILED > /opt/cw/BUILD_FAILED_$v; \
    done; \
    sha256sum /opt/cw/clustalw_* > /opt/cw/binaries.sha256 2>/dev/null; cp /opt/cwsrc/fetch_src.sha256 /opt/cw/; ls -l /opt/cw
WORKDIR /work
