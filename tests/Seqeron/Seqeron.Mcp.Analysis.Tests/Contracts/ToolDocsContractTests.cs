using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;
using Seqeron.Mcp.Tests.Contracts;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Audit A5-5 (B09 F170): the documents of the Analysis tools that expose B09's classes (<c>methodId</c> = SequenceAligner / GenomicAnalyzer / MultipleAlignment* / ...: the same list as
/// <c>scripts/docs/sync_mcp_tool_docs.py</c>) equal the tool metadata of the code. The other batches' tools of this server are theirs.
/// </summary>
[TestFixture]
public sealed class ToolDocsContractTests : ToolDocsContractTestsBase
{
    protected override Type ToolsType => typeof(AnalysisTools);

    protected override string DocsServer => "analysis";

    protected override string[] OwnedMethodIdPrefixes { get; } =
        ["SequenceAligner.", "GenomicAnalyzer.", "MultipleAlignment", "AnchorBasedAligner.", "DotMatcher.", "Minimap2Chaining.", "NucmerClusters.", "Ksw2Aligner.", "AlignmentFormats.", "MashCommands.", "MotifFinder.GenerateEmbossAmbiguousConsensus", "MummerDeltaFilter.", "MummerShow", "EmbossEst2genomePrograms."];
}
