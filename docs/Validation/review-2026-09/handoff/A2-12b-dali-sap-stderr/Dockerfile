# TEMPORARY hand-off image A2-12b (follow-up of A2-12-sap-dali, integrated 2026-11-23): (1) WHY t_coffee 13.41 stops (exit status 1) with the real SAP when a pair of structures has no SAP result (stderr of the 5 stopped sets);
# (2) t_coffee's dali_pair with DaliLite v5 on the PATH (the official tarball has no `DaliLite` command: the glue tc_generic_method.pl calls `DaliLite -pairwise pdb1 pdb2`, the command line of DaliLite v3) and, if the archived
# DaliLite v3 can be fetched, with the real DaliLite v3 behind a logging wrapper; (3) -mode 3dcoffee with the real SAP (the user's machine has a network: the name look-ups of the profile templates may work).
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
FROM --platform=linux/amd64 ubuntu:24.04

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl git build-essential gfortran autoconf automake libtool m4 perl python3 tcsh zip file strace binutils \
        t-coffee tm-align mustang ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

# ----- SAP 1.1.3 exactly as in A2-12 (sha256 of the archive pinned)
RUN mkdir -p /opt/sap /opt/src && cd /opt/src \
 && ( curl -sSL -A "Seqeron-review/1.0" -o sap.tar.gz https://github.com/mathbio-nimr-mrc-ac-uk/SAP/archive/v.1.1.3.tar.gz \
      && echo "1ee5025f8a900cd9d9c490f7038b98d80a619e3015f2dc97b869ea3033c459e0  sap.tar.gz" | sha256sum -c - \
      && tar xzf sap.tar.gz && cd SAP-v.1.1.3 && ./bootstrap && ./configure --prefix=/opt/sap && make && make install \
      && sha256sum /opt/sap/bin/sap > /opt/sap/sap.sha256 ) > /opt/sap/build.log 2>&1 \
 || echo BUILD_FAILED > /opt/sap/BUILD_FAILED; tail -3 /opt/sap/build.log

# ----- DaliLite v5 from the OFFICIAL tarball (http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz); its bin/ is put on the PATH of the dali_pair run to show that no `DaliLite` command exists there
RUN mkdir -p /opt/dali5 && cd /opt/dali5 \
 && ( curl -sSL -A "Seqeron-review/1.0" -o DaliLite.v5.tar.gz http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz && sha256sum DaliLite.v5.tar.gz > tarball.sha256 \
      && tar xzf DaliLite.v5.tar.gz && ls DaliLite.v5/bin > bin.list && cd DaliLite.v5/bin && make ) > build.log 2>&1 \
 || echo BUILD_FAILED > BUILD_FAILED; tail -3 build.log; ls /opt/dali5/DaliLite.v5/bin | head -80 > /opt/dali5/bin_after_build.list || true

# ----- DaliLite v3 (the version whose command line `DaliLite -pairwise` t_coffee uses): the download page of the paper (Hasegawa & Holm 2008) is gone from the live server (404); the Internet Archive keeps the page.
#       fetch_dali3.py lists the snapshot, downloads every .tar.gz / .zip link into /opt/dalilite3/src and records everything; a failure never aborts the image.
COPY fetch_dali3.py /opt/fetch_dali3.py
RUN mkdir -p /opt/dalilite3 && python3 -I /opt/fetch_dali3.py /opt/dalilite3 > /opt/dalilite3/fetch.log 2>&1 || echo FETCH_FAILED >> /opt/dalilite3/fetch.log; tail -5 /opt/dalilite3/fetch.log

WORKDIR /work
