using System.Runtime.CompilerServices;
using Seqeron.Genomics.Core;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// The Seqeron.Genomics library default is <see cref="LimitationMode.Strict"/>, under which the
/// guarded disorder units (DISORDER-REGION-001: <c>PredictDisorder</c>, <c>PredictMoRFs</c>) throw
/// <see cref="SeqeronLimitationException"/> instead of returning their best-effort result. The
/// Analysis MCP tools <c>predict_disorder</c> and <c>predict_morfs</c> wrap those units, so this
/// test assembly runs under <see cref="LimitationMode.Permissive"/> to exercise the real algorithm
/// output (mirrors <c>Seqeron.Genomics.Tests</c>' bootstrap).
/// </summary>
internal static class LimitationPolicyTestBootstrap
{
    [ModuleInitializer]
    internal static void Init() => LimitationPolicy.DefaultMode = LimitationMode.Permissive;
}
