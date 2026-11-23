# Seqeron review package A2-29b: runs the published ksw2 / minimap2 extension functions (ksw_extz, ksw_extz2_sse, ksw_extd, ksw_extd2_sse) on the review cases.
# x86_64 only (SSE intrinsics). On Apple Silicon build and run with --platform linux/amd64 (see run.sh).
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq && apt-get install -y -qq gcc libc6-dev curl ca-certificates python3 coreutils valgrind && rm -rf /var/lib/apt/lists/*
COPY harness.c cases.txt inside.sh /in/
RUN chmod +x /in/inside.sh
CMD ["/in/inside.sh"]
