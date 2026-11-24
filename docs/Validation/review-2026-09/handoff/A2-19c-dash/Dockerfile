# TEMPORARY hand-off image A2-19c-dash: the MAFFT 7.505 source package (upstream, "without extensions") built together with the Go program dash_client (core/dash_client.go),
# which the Debian binary package does not contain.  Run 1 (2026-11-24) failed in the build ("dash_client: No such file", no mafft): Go compiled INSIDE the emulated amd64 container.
# Now: stage 1 compiles dash_client on the HOST architecture (cross-compiled to linux/amd64, static, CGO off: a pure HTTP client, no emulation); stage 2 builds the C programs of
# MAFFT as before (linux/amd64, emulated on Apple Silicon) WITHOUT the Go step and copies dash_client in.  Every step logs to /opt/*.log and a failure never aborts the image build.
# No account data, no secrets; every download uses the neutral User-Agent "Seqeron-review/1.0".  The tarball is the official one (https://mafft.cbrc.jp/alignment/software/); its SHA-256
# is recorded in out/versions.txt (the audit machine's copy: 2cdc797068b4b52dfd2a6a46502c832911cacc1cdd8a20c8e26b3b96a6208896).
ARG MAFFT_URL=https://mafft.cbrc.jp/alignment/software/mafft-7.505-without-extensions-src.tgz

FROM --platform=$BUILDPLATFORM golang:1.19-bullseye AS gobuild
ARG MAFFT_URL
RUN mkdir -p /opt/src /out && cd /opt/src \
 && curl -fsSL -A "Seqeron-review/1.0" -o mafft-src.tgz "$MAFFT_URL" \
 && sha256sum mafft-src.tgz | tee /out/mafft-src.sha256 \
 && tar xzf mafft-src.tgz \
 && cd mafft-7.505-without-extensions/core \
 && ( GO111MODULE=off GOOS=linux GOARCH=amd64 CGO_ENABLED=0 go build -o /out/dash_client dash_client.go > /out/go-build.log 2>&1 ; echo "go build exit $?" >> /out/go-build.log ) \
 ; cat /out/go-build.log; ls -l /out

FROM debian:12-slim
ARG MAFFT_URL
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates curl gcc libc6-dev make python3 zip file procps gawk \
 && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /opt/src && cd /opt/src \
 && curl -fsSL -A "Seqeron-review/1.0" -o mafft-src.tgz "$MAFFT_URL" \
 && sha256sum mafft-src.tgz | tee /opt/mafft-src.sha256 \
 && tar xzf mafft-src.tgz
RUN cd /opt/src/mafft-7.505-without-extensions/core \
 && ( make ENABLE_MULTITHREAD= PREFIX=/usr/local BINDIR=/usr/local/bin LIBDIR=/usr/local/libexec/mafft > /opt/mafft.build.log 2>&1 ; echo "make exit $?" >> /opt/mafft.build.log \
   ; make ENABLE_MULTITHREAD= PREFIX=/usr/local BINDIR=/usr/local/bin LIBDIR=/usr/local/libexec/mafft install >> /opt/mafft.build.log 2>&1 ; echo "make install exit $?" >> /opt/mafft.build.log ) \
 ; tail -8 /opt/mafft.build.log
COPY --from=gobuild /out/ /opt/gobuild/
RUN mkdir -p /usr/local/libexec/mafft \
 && { [ -x /opt/gobuild/dash_client ] && cp /opt/gobuild/dash_client /usr/local/libexec/mafft/dash_client && chmod 755 /usr/local/libexec/mafft/dash_client || true; } \
 ; cp /opt/gobuild/go-build.log /opt/go-build.log 2> /dev/null; ls -l /usr/local/libexec/mafft/dash_client /usr/local/bin/mafft 2>&1 \
 ; { [ -x /usr/local/libexec/mafft/dash_client ] && [ -x /usr/local/bin/mafft ]; } || echo DASH_CLIENT_BUILD_FAILED > /opt/DASH_CLIENT_BUILD_FAILED
WORKDIR /work
