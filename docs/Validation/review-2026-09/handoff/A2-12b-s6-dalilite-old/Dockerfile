# TEMPORARY hand-off image A2-12b-s6-dalilite-old (B09 WP71): is an OLDER DaliLite (v2.x / v1) obtainable whose `DaliLite -pairwise` prints the alignment as `Query` / `Sbjct` lines to standard output (the glue of t_coffee 13.41, tc_generic_method.pl
# seq2pdb_pair, reads only that)? DaliLite 3.3 (sha256 85f6dd8d...) does not (272 calls of run 10: 0 lines). This image only FETCHES and INSPECTS (file lists, CHANGES / README / the `DaliLite` script as text); it compiles and runs nothing that it downloads.
# No account data, no secrets; downloads use the neutral User-Agent "Seqeron-review/1.0".
FROM --platform=linux/amd64 ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl python3 file zip && rm -rf /var/lib/apt/lists/*
COPY fetch_old_dalilite.py /opt/fetch_old_dalilite.py
WORKDIR /work
