# syntax=docker/dockerfile:1
# E8 hand-off image: the B-DCS pipeline (ADR-004, command e6, PH0-PH5) plus the generator, the verifier, the runner and the sampler of the E8 protocol, for Linux on x86-64 and arm64.
#
# STATUS: UNBUILT AND UNTESTED. The session that wrote this file had no Docker daemon. The build context closure was checked by building the program from exactly the files that Dockerfile.dockerignore lets in,
# offline, with an empty NuGet cache (see E8-RUNBOOK.md); the image itself, the apt step, the scripts inside it, macOS and arm64 behaviour are NOT verified.
#
# Build from the repository root (the build context), natively on the architecture of the machine that will run it:
#   docker build -f docs/refactoring/suffix-tree-upgrade/analysis/e8/Dockerfile --build-arg COMMIT=<git hash> -t seqeron-e8 .
# both architectures:
#   docker buildx build --platform linux/amd64,linux/arm64 -f docs/refactoring/suffix-tree-upgrade/analysis/e8/Dockerfile .
# e8-docker.sh does the first of these with the right --platform.
#
# What the build fetches: the two public .NET base images from mcr.microsoft.com and the Debian packages listed in the apt line (from the distribution's own mirrors). No NuGet package (the substitute Directory.Build.props
# leaves the analyzer package out), no source from the network, no credential, no build argument that could carry one.

ARG DOTNET_TAG=10.0

# ---- build stage: runs on the BUILD platform; a framework-dependent publish is architecture neutral (IL only; no apphost, no RID-specific pack)
FROM --platform=$BUILDPLATFORM mcr.microsoft.com/dotnet/sdk:${DOTNET_TAG} AS build
WORKDIR /src
COPY . .
# the repository's Directory.Build.props (not sent, see Dockerfile.dockerignore) is replaced by a file without the analyzer package
COPY docs/refactoring/suffix-tree-upgrade/analysis/e8/docker/Directory.Build.props /src/Directory.Build.props
# the projects ask for the trimming / AOT analyzers (an implicit ILLink package): switched off from the command line so that the restore needs no package at all; a newer SDK than the one the code was
# developed with must not fail the build over a new warning (warnings are errors in the project files, so that is switched off here, too)
RUN dotnet publish apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj -c Release -o /out --nologo \
    -p:UseAppHost=false -p:IsTrimmable=false -p:IsAotCompatible=false -p:EnableTrimAnalyzer=false -p:EnableAotAnalyzer=false \
    -p:GenerateDocumentationFile=false -p:TreatWarningsAsErrors=false -p:EnforceCodeStyleInBuild=false

# ---- runtime stage: per TARGET architecture
FROM mcr.microsoft.com/dotnet/runtime:${DOTNET_TAG} AS runtime
ARG COMMIT=unknown
ARG TARGETARCH
LABEL org.opencontainers.image.title="E8 hand-off image: B-DCS pipeline at 1-3 Gbp" \
      org.opencontainers.image.description="ADR-004 E8: generator, e6 pipeline, verifier, runner and sampler. Unbuilt and untested when this file was written." \
      org.opencontainers.image.revision="${COMMIT}"
# only what the scripts need: bash, python3 (sampler, tables, preflight I/O test), procps, coreutils, bc, util-linux (mountpoint), findutils for the cleanups
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash python3 procps coreutils bc util-linux findutils \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /out /opt/e8/bin
COPY docs/refactoring/suffix-tree-upgrade/analysis/e8/*.sh docs/refactoring/suffix-tree-upgrade/analysis/e8/*.py /opt/e8/scripts/
RUN chmod +x /opt/e8/scripts/*.sh
ENV E8_DLL=/opt/e8/bin/SuffixTree.EmProto.dll \
    E8_DATA=/data \
    E8_IMAGE_ARCH=${TARGETARCH} \
    E8_IMAGE_COMMIT=${COMMIT} \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    DOTNET_NOLOGO=1
VOLUME ["/data"]
WORKDIR /data
ENTRYPOINT ["/opt/e8/scripts/e8.sh"]
CMD ["help"]
