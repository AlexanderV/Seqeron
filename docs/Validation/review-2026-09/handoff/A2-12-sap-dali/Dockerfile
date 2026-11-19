# TEMPORARY hand-off image A2-12: the REAL structure aligners behind t_coffee's `sap_pair` / `-mode 3dcoffee` (SAP 1.1.3) and `dali_pair` (DaliLite v5) next to the same
# t_coffee 13.41 package the audit machine uses (Ubuntu 24.04 `t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1`).
# Why: on the audit machine SAP and DaliLite are not installed; t_coffee then runs TMalign instead of SAP ("SAP is not installed, TMalign will be used instead") and stops on dali_pair
# (ERROR: UNSPECIFIED). Both situations are reproduced and tested; THIS image provides the real programs so that the real SAP / Dali glue of t_coffee can be ported and verified.
# Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
FROM --platform=linux/amd64 ubuntu:24.04

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl git build-essential gfortran autoconf automake libtool m4 perl python3 tcsh zip file \
        t-coffee tm-align mustang ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

# ----- SAP 1.1.3 (Taylor & Orengo, Protein Sci 8:654, 1999), the source of the bioconda package `sap`: https://github.com/mathbio-nimr-mrc-ac-uk/SAP/archive/v.1.1.3.tar.gz
# (the recipe's build is `./bootstrap && ./configure --prefix=... && make && make install`). The build never aborts the image: log in /opt/sap/build.log, failure marker /opt/sap/BUILD_FAILED.
RUN mkdir -p /opt/sap /opt/src && cd /opt/src \
 && ( curl -sSL -A "Seqeron-review/1.0" -o sap.tar.gz https://github.com/mathbio-nimr-mrc-ac-uk/SAP/archive/v.1.1.3.tar.gz \
      && echo "1ee5025f8a900cd9d9c490f7038b98d80a619e3015f2dc97b869ea3033c459e0  sap.tar.gz" | sha256sum -c - \
      && tar xzf sap.tar.gz && cd SAP-v.1.1.3 && ./bootstrap && ./configure --prefix=/opt/sap && make && make install \
      && sha256sum /opt/sap/bin/sap > /opt/sap/sap.sha256 ) > /opt/sap/build.log 2>&1 \
 || echo BUILD_FAILED > /opt/sap/BUILD_FAILED; tail -3 /opt/sap/build.log

# ----- DaliLite v5 (Holm, Bioinformatics 35:5326, 2019): https://github.com/liisaholm/DaliLite.v5 (`cd bin && make`, "ignore warnings"). The commit that was built is recorded in /opt/dali/COMMIT.
RUN mkdir -p /opt/dali \
 && ( git clone --depth 1 https://github.com/liisaholm/DaliLite.v5 /opt/dali/DaliLite.v5 && cd /opt/dali/DaliLite.v5 && git rev-parse HEAD > /opt/dali/COMMIT && cd bin && make ) > /opt/dali/build.log 2>&1 \
 || echo BUILD_FAILED > /opt/dali/BUILD_FAILED; tail -3 /opt/dali/build.log

WORKDIR /work
