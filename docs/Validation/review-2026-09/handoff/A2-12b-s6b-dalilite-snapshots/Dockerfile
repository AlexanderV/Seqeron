# TEMPORARY hand-off image A2-12b-s6b-dalilite-snapshots (B09 WP76): direct Internet Archive snapshot probes for OLDER DaliLite releases (1.0 ... 3.2, DaliLite.tar.gz / .zip) under every known distribution place; keeps what is a real archive (magic bytes) and
# reads CHANGES / README / INSTALL / Makefile / the `DaliLite` script as TEXT (python tarfile / zipfile in memory). Nothing that is downloaded is extracted to disk, compiled or run. No account data, no secrets; User-Agent "Seqeron-review/1.0".
FROM --platform=linux/amd64 ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl python3 file zip gzip && rm -rf /var/lib/apt/lists/*
COPY fetch_dalilite_snapshots.py /opt/fetch_dalilite_snapshots.py
WORKDIR /work
