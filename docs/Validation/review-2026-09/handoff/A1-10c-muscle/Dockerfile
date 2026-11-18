# TEMPORARY hand-off image A1-10c-muscle: MUSCLE 3.8.1551 (the final v3.8 release; Debian source package muscle3 3.8.1551-2) built from source with the upstream
# Makefile flags (-O3 -funroll-loops -DNDEBUG=1). No patch of any kind. It is the oracle of SequenceAligner.MultipleAlignMuscle (the port of the same source);
# the installed muscle 5.x of the audit machine is a different program (probcons-like PPP pipeline) and cannot serve.
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=debian:12.7-slim@sha256:d83056144b2dd301730d2739635c8cbdeaaae20d6887146434184f8c060f03ce
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates g++ make python3 zip file procps \
 && rm -rf /var/lib/apt/lists/*

COPY fetch_src.py /opt/fetch_src.py
RUN mkdir -p /opt/muscle3src && cd /opt/muscle3src && python3 -I /opt/fetch_src.py && ls | wc -l

# Build exactly as the upstream Makefile does (CFLAGS += -O3 -funroll-loops -Winline -DNDEBUG=1; LDLIBS += -lm); fall back to the permissive mode of
# newer compilers only if the plain build fails (the fallback does not change floating-point semantics).
RUN set -eux; cd /opt/muscle3src; mkdir -p /opt/muscle3 /tmp/obj; \
    (for f in *.cpp; do g++ -c -O3 -funroll-loops -w -DNDEBUG=1 $f -o /tmp/obj/${f%.cpp}.o; done) \
 || (for f in *.cpp; do g++ -c -O3 -funroll-loops -w -fpermissive -std=gnu++98 -DNDEBUG=1 $f -o /tmp/obj/${f%.cpp}.o; done); \
    g++ -O3 -funroll-loops -o /opt/muscle3/muscle3 /tmp/obj/*.o -lm; \
    sha256sum /opt/muscle3/muscle3 > /opt/muscle3/muscle3.sha256; \
    cp /opt/muscle3src/fetch_src.sha256 /opt/muscle3/
WORKDIR /work
