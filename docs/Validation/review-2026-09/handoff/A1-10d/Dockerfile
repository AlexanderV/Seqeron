# TEMPORARY hand-off image A1-10d: ClustalW 2.1+lgpl (Debian source 2.1+lgpl-7) built from source with a 3-line patch that makes the command-line
# `-ENDGAPS` handler read the environment variables CW_USEENDGAPS / CW_ENDGAPPEN (so that UserParameters::useEndGaps / endGapPenalties, which
# the real command line never sets to true, can be switched on). Nothing else is changed; the control run (no variables) must reproduce the
# rows of the installed clustalw 2.1 (inputs/sets.tsv holds them). Build and run with `--platform linux/amd64` (Apple Silicon: emulated).
# No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12.7-slim@sha256:d83056144b2dd301730d2739635c8cbdeaaae20d6887146434184f8c060f03ce
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl g++ make python3-minimal zip file procps \
 && rm -rf /var/lib/apt/lists/*

# Source: sources.debian.org clustalw 2.1+lgpl-7, directory src/ (the files are fetched one by one through the sources.debian.org JSON API).
COPY fetch_src.py /opt/fetch_src.py
RUN mkdir -p /opt/cwsrc && cd /opt/cwsrc && python3 -I /opt/fetch_src.py

# The patch: -ENDGAPS (command-line: setUseEndGaps(false)) now reads CW_USEENDGAPS / CW_ENDGAPPEN.
RUN set -eux; cd /opt/cwsrc/src; \
    grep -n "userParameters->setUseEndGaps(false);" interface/CommandLineParser.cpp; \
    sed -i 's|userParameters->setUseEndGaps(false);|userParameters->setUseEndGaps(getenv("CW_USEENDGAPS") != NULL); userParameters->setEndGapPenalties(getenv("CW_ENDGAPPEN") != NULL);|' interface/CommandLineParser.cpp; \
    grep -n "CW_USEENDGAPS" interface/CommandLineParser.cpp; \
    sed -e 's/@CLUSTALW_NAME@/clustalw/;s/@CLUSTALW_VERSION@/2.1/' clustalw_version.h.in > clustalw_version.h; : > config.h; \
    mkdir -p /opt/clustalw /tmp/obj; \
    for f in $(find . -name '*.cpp'); do g++ -O2 -w -DHAVE_CONFIG_H -DUNIX=1 -DOS_UNIX=1 -I. -c $f -o /tmp/obj/$(echo $f | tr '/' '_' | sed 's/^\._//').o; done; \
    g++ -O2 /tmp/obj/*.o -o /opt/clustalw/clustalw; \
    sha256sum /opt/clustalw/clustalw > /opt/clustalw/clustalw.sha256; \
    cp -r /opt/cwsrc/src/interface/CommandLineParser.cpp /opt/clustalw/CommandLineParser.patched.cpp
WORKDIR /work
