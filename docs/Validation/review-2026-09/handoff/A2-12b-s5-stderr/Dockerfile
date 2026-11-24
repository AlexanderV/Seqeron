# TEMPORARY hand-off image A2-12b-s5-stderr (follow-up of A2-12b run 10, integrated by B09 WP71): re-runs ONLY the 15 dali_pair sets that ended with exit status 1 in run 10 (set positions 5 8 10 15 16 17 18 20 21 23 24 25 29 32 38 of the
# 60-set plan) with t_coffee's standard error recorded (TC_STDERR=1), on the SAME real DaliLite 3.3 (full-backbone structures, same seeds).
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
FROM --platform=linux/amd64 ubuntu:24.04

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl git build-essential gfortran perl python3 tcsh zip file strace binutils \
        t-coffee tm-align mustang ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

# ----- DaliLite 3.3 (the version whose command line `DaliLite -pairwise` t_coffee uses): from the live server of its authors (http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite_3.3.tar.gz, sha256 85f6dd8d...),
#       else (as in runs 4 - 10) from the Internet Archive copy of its download page; fetch_dali3.py builds it, repairs its precompiled g77 programs (libg2c shim) and installs the logging wrapper. A failure never aborts the image.
COPY fetch_dali3.py /opt/fetch_dali3.py
RUN mkdir -p /opt/dalilite3 && python3 -I /opt/fetch_dali3.py /opt/dalilite3 > /opt/dalilite3/fetch.log 2>&1 || echo FETCH_FAILED >> /opt/dalilite3/fetch.log; tail -5 /opt/dalilite3/fetch.log

WORKDIR /work
