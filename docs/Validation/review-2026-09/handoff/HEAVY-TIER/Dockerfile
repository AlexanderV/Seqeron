# TEMPORARY hand-off image HEAVY-TIER: runs the .NET test tiers of the Seqeron repository on the user's machine (usually faster than the cloud sandbox).
# Official Microsoft .NET 10 SDK image (multi-arch: native arm64 on Apple Silicon, no emulation needed). No secrets, no account data.
FROM mcr.microsoft.com/dotnet/sdk:10.0
ENV DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_gcServer=0 DOTNET_TieredPGO=0
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3 rsync zip ca-certificates procps \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /work
COPY in_container.sh /usr/local/bin/in_container.sh
RUN chmod +x /usr/local/bin/in_container.sh
