# syntax=docker/dockerfile:1
# Hand-off image of the package s1b-handoff (ADR-004 section 5a, 1b.9f): the production chain of the SESA builder (the Bench application: sesa resume build, sesa handoff), the prototype's verifier and
# generator (e8-verify, e8-gen), the scripts, and (target "full") the .NET SDK with the repository's sources and warm builds of the test projects, for Linux on x86-64 and arm64.
#
# STATUS: built and run on linux/amd64 where it was written (a Docker daemon, started by hand, native architecture, builds through a proxy); NOT built, NOT run on linux/arm64 or on Docker Desktop: the arm64 build is the
# owner's first run. See S1B-RUNBOOK.md section 9.
#
# Build from the repository root (the build context), natively on the architecture of the machine that will run it:
#   docker build --target full  -f docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff/Dockerfile --build-arg COMMIT=<git hash> -t seqeron-s1b .
#   docker build --target chain -f ...                                                                                             -t seqeron-s1b .     (no SDK, no test sets: 0.4 GB)
# both architectures:  docker buildx build --platform linux/amd64,linux/arm64 --target full -f ... .
# s1b-docker.sh image does the first of these with the right --platform.
#
# What the build fetches: the public .NET base images from mcr.microsoft.com, the Ubuntu packages of the apt line (from the distribution's mirrors), and, in the "full" target only, the NuGet packages of the test
# projects and of the repository's analyzer (from nuget.org: the test frameworks NUnit, CsCheck, coverlet, the SonarAnalyzer). The chain programs need no NuGet package at all. Everything is DOWNLOADED; nothing is
# uploaded; no credential, no build argument that could carry one. The containers run with no network needed (the test sets work with --network none).

ARG DOTNET_TAG=10.0

# ---- the chain programs: a framework-dependent publish is architecture neutral (IL only; no apphost, no RID-specific pack); the repository's root props (the analyzer package) are replaced by the E8 substitute
FROM --platform=$BUILDPLATFORM mcr.microsoft.com/dotnet/sdk:${DOTNET_TAG} AS chain-build
WORKDIR /chain
COPY . .
COPY docs/refactoring/suffix-tree-upgrade/analysis/e8/docker/Directory.Build.props /chain/Directory.Build.props
RUN rm -f /chain/Directory.Packages.props \
 && for proj in apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj; do \
      dotnet publish "$proj" -c Release -o "/out/$(basename "$(dirname "$proj")")" --nologo -p:UseAppHost=false -p:IsTrimmable=false -p:IsAotCompatible=false -p:EnableTrimAnalyzer=false \
        -p:EnableAotAnalyzer=false -p:GenerateDocumentationFile=false -p:TreatWarningsAsErrors=false -p:EnforceCodeStyleInBuild=false || exit 1; \
    done

# ---- the chain image: the runtime, the programs, the scripts
FROM mcr.microsoft.com/dotnet/runtime:${DOTNET_TAG} AS chain
ARG COMMIT=unknown
ARG TARGETARCH
LABEL org.opencontainers.image.title="S1B hand-off image: the production chain of the SESA builder" \
      org.opencontainers.image.description="ADR-004 1b.9f: sesa resume build, sesa handoff, e8-verify, e8-gen, the scripts, the sampler, the form and the checker." \
      org.opencontainers.image.revision="${COMMIT}"
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash python3 procps coreutils util-linux findutils \
 && rm -rf /var/lib/apt/lists/*
COPY --from=chain-build /out/SuffixTree.Esa.Bench /opt/s1b/bin/bench
COPY --from=chain-build /out/SuffixTree.EmProto /opt/s1b/bin/emproto
COPY docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff /opt/s1b/analysis/s1b-handoff
COPY docs/refactoring/suffix-tree-upgrade/analysis/e8/e8lib.sh /opt/s1b/analysis/e8/e8lib.sh
COPY docs/refactoring/suffix-tree-upgrade/analysis/bench/gen/gen_inputs.py /opt/s1b/analysis/bench/gen/gen_inputs.py
RUN chmod +x /opt/s1b/analysis/s1b-handoff/*.sh
ENV S1B_DATA=/data \
    S1B_BENCH=/opt/s1b/bin/bench/SuffixTree.Esa.Bench.dll \
    S1B_EMPROTO=/opt/s1b/bin/emproto/SuffixTree.EmProto.dll \
    S1B_IMAGE_ARCH=${TARGETARCH} \
    S1B_COMMIT=${COMMIT} \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    DOTNET_NOLOGO=1
VOLUME ["/data"]
WORKDIR /data
ENTRYPOINT ["/opt/s1b/analysis/s1b-handoff/s1b.sh"]
CMD ["help"]

# ---- the test sets: the SDK, the repository's sources as they are (its own props), the NuGet restore and warm Debug and Release builds of the three test projects, so that `docker run --network none` works
FROM mcr.microsoft.com/dotnet/sdk:${DOTNET_TAG} AS full
ARG COMMIT=unknown
ARG TARGETARCH
# the flags are the ones that s1b.sh uses at run time (S1B_MSBUILD_FLAGS): the same flags keep the warm builds valid. The audit of known vulnerabilities needs the network; a newer SDK than the one the code was
# developed with must not fail the build over a new analyzer warning (warnings are errors in the project files, so that is switched off here, too; the compiled code is the same)
ARG MSBUILD_FLAGS="-p:NuGetAudit=false -p:TreatWarningsAsErrors=false -p:EnforceCodeStyleInBuild=false"
LABEL org.opencontainers.image.title="S1B hand-off image with the test sets" \
      org.opencontainers.image.description="ADR-004 1b.9f: the chain image plus the .NET SDK, the sources and the warm builds of the Esa, EmProto and Persistent test projects." \
      org.opencontainers.image.revision="${COMMIT}"
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash python3 procps coreutils util-linux findutils \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN set -e; for p in Esa EmProto Persistent; do for c in Debug Release; do \
      dotnet build "tests/SuffixTree/SuffixTree.$p.Tests" -c "$c" --nologo -v:q $MSBUILD_FLAGS; \
    done; done
COPY --from=chain-build /out/SuffixTree.Esa.Bench /opt/s1b/bin/bench
COPY --from=chain-build /out/SuffixTree.EmProto /opt/s1b/bin/emproto
RUN mkdir -p /opt/s1b/analysis/bench/gen /opt/s1b/analysis/e8 \
 && cp -r /src/docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff /opt/s1b/analysis/s1b-handoff \
 && cp /src/docs/refactoring/suffix-tree-upgrade/analysis/e8/e8lib.sh /opt/s1b/analysis/e8/e8lib.sh \
 && cp /src/docs/refactoring/suffix-tree-upgrade/analysis/bench/gen/gen_inputs.py /opt/s1b/analysis/bench/gen/gen_inputs.py \
 && chmod +x /opt/s1b/analysis/s1b-handoff/*.sh
ENV S1B_DATA=/data \
    S1B_BENCH=/opt/s1b/bin/bench/SuffixTree.Esa.Bench.dll \
    S1B_EMPROTO=/opt/s1b/bin/emproto/SuffixTree.EmProto.dll \
    S1B_SRC=/src \
    S1B_MSBUILD_FLAGS="${MSBUILD_FLAGS}" \
    S1B_IMAGE_ARCH=${TARGETARCH} \
    S1B_COMMIT=${COMMIT} \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    DOTNET_NOLOGO=1
VOLUME ["/data"]
WORKDIR /data
ENTRYPOINT ["/opt/s1b/analysis/s1b-handoff/s1b.sh"]
CMD ["help"]
