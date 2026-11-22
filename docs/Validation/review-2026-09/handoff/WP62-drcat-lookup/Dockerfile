# TEMPORARY hand-off image WP62-drcat-lookup: EMBOSS 6.6.0 (Ubuntu 22.04 package emboss 6.6.0+dfsg-12) with network access, the oracle of the DRCAT catalogue lookup (`queryDbfind` of ajquery.c: a USA `db:id` of a database that
# the configuration does not define is looked up in DRCAT.dat and fetched from the URL / method the entry names). The sandbox of the audit has no network access to the databases. Build and run with `--platform linux/amd64`.
# No account data, no secrets.
ARG BASE=ubuntu:22.04
FROM --platform=linux/amd64 ${BASE}
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates emboss emboss-data python3 zip file procps \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /work
