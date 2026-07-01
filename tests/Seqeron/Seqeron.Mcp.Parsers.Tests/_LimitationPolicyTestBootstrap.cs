using System.Runtime.CompilerServices;
using Seqeron.Genomics.Core;

namespace Seqeron.Mcp.Parsers.Tests;

/// <summary>
/// The library default is <see cref="LimitationMode.Strict"/> (limitation branches throw). These MCP
/// parser tests exercise FASTQ quality-encoding auto-detection on overlap-ambiguous input — a guarded
/// branch (PARSE-FASTQ-001) — so this assembly runs under <see cref="LimitationMode.Permissive"/>.
/// The Strict default is covered by <c>LimitationPolicy_Strict_Tests</c> in Seqeron.Genomics.Tests.
/// </summary>
internal static class LimitationPolicyTestBootstrap
{
    [ModuleInitializer]
    internal static void Init() => LimitationPolicy.DefaultMode = LimitationMode.Permissive;
}
