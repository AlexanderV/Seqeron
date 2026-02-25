param(
    [switch]$IncludeExtra
)

$ErrorActionPreference = "Stop"

function Run-Test {
    param([string]$Command)
    Write-Host ">> $Command"
    Invoke-Expression $Command
}

Run-Test "dotnet test tests/SuffixTree/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --nologo"
Run-Test "dotnet test tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj -c Release --nologo"
Run-Test "dotnet test tests/SuffixTree/SuffixTree.Mcp.Core.Tests/SuffixTree.Mcp.Core.Tests.csproj -c Release --nologo"

if ($IncludeExtra) {
    Run-Test "dotnet test tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj -c Release --nologo --filter Category=Parity"
    Run-Test "dotnet test tests/SuffixTree/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --nologo --filter Category=Performance"
}

Write-Host "All suffix-tree refactor gates passed."
