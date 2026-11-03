#!/usr/bin/env bash
# Checked-pointer configuration of the ESA tests (ST-S1 1o.0, review MAJ-1): compiles the Esa assembly with ESA_BOUNDS, so that every Buf<T> access is
# bounds-checked, and runs the test project in that mode. ESA_EXPECT_BOUNDS=true makes a run that silently fell back to the Release accessor FAIL.
#   usage: run_esa_checked.sh [extra dotnet-test args, e.g. --filter Category=Exhaustive]
# The Release (raw pointer) configuration is the plain `dotnet test ... -c Release`; both must be green. Afterwards the default build is restored.
set -eu
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../.." && pwd)"
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH DOTNET_CLI_TELEMETRY_OPTOUT=1
P=$ROOT/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj
ESA_EXPECT_BOUNDS=true dotnet test "$P" -c Release --nologo -p:EsaBounds=true "$@"
status=$?
dotnet build "$P" -c Release --nologo -p:EsaBounds=false >/dev/null   # restore the raw-pointer build in bin/
exit $status
