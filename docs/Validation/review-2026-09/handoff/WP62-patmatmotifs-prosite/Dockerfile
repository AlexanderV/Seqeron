# TEMPORARY hand-off image WP62-patmatmotifs-prosite: EMBOSS 6.6.0 (Ubuntu 22.04 package emboss 6.6.0+dfsg-12 = the version the port is built against) with the real PROSITE database
# (prosite.dat, prosite.doc of the SIB FTP site) turned into the data files of patmatmotifs by prosextract. It is the oracle of the PROSITE input of `patmatmotifs` (the sandbox has no PROSITE data and no
# network access to it). Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
ARG BASE=ubuntu:22.04
FROM --platform=linux/amd64 ${BASE}

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates emboss emboss-data python3 curl zip file procps valgrind \
 && rm -rf /var/lib/apt/lists/*

# the PROSITE release (the SIB site; the files are text); the checksums are written to the output so that the release is known
RUN mkdir -p /opt/prosite && cd /opt/prosite \
 && curl -fsSL -A "Seqeron-review/1.0" -O https://ftp.expasy.org/databases/prosite/prosite.dat \
 && curl -fsSL -A "Seqeron-review/1.0" -O https://ftp.expasy.org/databases/prosite/prosite.doc \
 && sha256sum prosite.dat prosite.doc > SHA256SUMS && head -3 prosite.dat && wc -l prosite.dat prosite.doc

# prosextract writes prosite.lines (the index patmatmotifs reads) and the documentation files into the PROSITE directory of the EMBOSS data
RUN mkdir -p /opt/emboss-data/PROSITE && cd /opt/emboss-data/PROSITE && cp /opt/prosite/prosite.dat /opt/prosite/prosite.doc . \
 && ( printf '%s\n' /opt/emboss-data/PROSITE/ | prosextract -auto -prositedir /opt/emboss-data/PROSITE/ ) > /opt/emboss-data/prosextract.log 2>&1; ls -l /opt/emboss-data/PROSITE | head -20
WORKDIR /work
