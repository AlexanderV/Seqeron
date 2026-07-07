using System.Text.Json.Serialization;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;

namespace Seqeron.Mcp.MolTools;

[JsonSourceGenerationOptions(WriteIndented = false)]
// Primitive / shared collection types
[JsonSerializable(typeof(string))]
[JsonSerializable(typeof(string[]))]
[JsonSerializable(typeof(List<string>))]
[JsonSerializable(typeof(IReadOnlyList<string>))]
[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(int[]))]
[JsonSerializable(typeof(List<int>))]
[JsonSerializable(typeof(IReadOnlyList<int>))]
[JsonSerializable(typeof(double))]
[JsonSerializable(typeof(bool))]

// PrimerDesigner DTOs (from Seqeron.Genomics.MolTools)
[JsonSerializable(typeof(PrimerParameters))]
[JsonSerializable(typeof(PrimerParameters?))]
[JsonSerializable(typeof(PrimerCandidate))]
[JsonSerializable(typeof(PrimerCandidate[]))]
[JsonSerializable(typeof(List<PrimerCandidate>))]
[JsonSerializable(typeof(IReadOnlyList<PrimerCandidate>))]
[JsonSerializable(typeof(PrimerPairResult))]

// PrimerDesigner result wrappers
[JsonSerializable(typeof(TmResult))]
[JsonSerializable(typeof(HomopolymerLengthResult))]
[JsonSerializable(typeof(DinucleotideRepeatResult))]
[JsonSerializable(typeof(HairpinPotentialResult))]
[JsonSerializable(typeof(PrimerDimerResult))]
[JsonSerializable(typeof(ThreePrimeStabilityResult))]
[JsonSerializable(typeof(PrimerCandidateListResult))]

// RestrictionAnalyzer DTOs (from Seqeron.Genomics.MolTools) — fully qualified to avoid
// any potential collision with future server-side types of the same simple name.
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.RestrictionEnzyme))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.RestrictionEnzyme[]))]
[JsonSerializable(typeof(List<global::Seqeron.Genomics.MolTools.RestrictionEnzyme>))]
[JsonSerializable(typeof(IReadOnlyList<global::Seqeron.Genomics.MolTools.RestrictionEnzyme>))]
[JsonSerializable(typeof(OverhangType))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.RestrictionSite))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.RestrictionSite[]))]
[JsonSerializable(typeof(List<global::Seqeron.Genomics.MolTools.RestrictionSite>))]
[JsonSerializable(typeof(IReadOnlyList<global::Seqeron.Genomics.MolTools.RestrictionSite>))]
[JsonSerializable(typeof(DigestFragment))]
[JsonSerializable(typeof(DigestFragment[]))]
[JsonSerializable(typeof(List<DigestFragment>))]
[JsonSerializable(typeof(IReadOnlyList<DigestFragment>))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.DigestSummary))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.RestrictionMap))]
[JsonSerializable(typeof(Dictionary<string, List<int>>))]
[JsonSerializable(typeof(IReadOnlyDictionary<string, List<int>>))]

// RestrictionAnalyzer result wrappers
[JsonSerializable(typeof(EnzymeLookupResult))]
[JsonSerializable(typeof(EnzymeListResult))]
[JsonSerializable(typeof(RestrictionSiteListResult))]
[JsonSerializable(typeof(DigestResult))]
[JsonSerializable(typeof(EnzymeCompatibilityPair))]
[JsonSerializable(typeof(EnzymeCompatibilityPair[]))]
[JsonSerializable(typeof(List<EnzymeCompatibilityPair>))]
[JsonSerializable(typeof(IReadOnlyList<EnzymeCompatibilityPair>))]
[JsonSerializable(typeof(CompatibleEnzymesResult))]
[JsonSerializable(typeof(EnzymeCompatibilityResult))]

// Shared dictionary types for codon usage / RSCU / frequency tables
[JsonSerializable(typeof(Dictionary<string, int>))]
[JsonSerializable(typeof(IReadOnlyDictionary<string, int>))]
[JsonSerializable(typeof(Dictionary<string, double>))]
[JsonSerializable(typeof(IReadOnlyDictionary<string, double>))]
[JsonSerializable(typeof(Dictionary<string, string>))]
[JsonSerializable(typeof(IReadOnlyDictionary<string, string>))]

// CodonUsageAnalyzer DTOs
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CodonUsageStatistics))]
[JsonSerializable(typeof(CodonCountsResult))]
[JsonSerializable(typeof(RscuResult))]
[JsonSerializable(typeof(CaiResult))]
[JsonSerializable(typeof(EncResult))]

// CodonOptimizer DTOs (source types are nested in CodonOptimizer; fully qualified to avoid ambiguity)
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CodonOptimizer.CodonUsageTable))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CodonOptimizer.CodonUsageTable?))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CodonOptimizer.OptimizationResult))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CodonOptimizer.OptimizationStrategy))]
[JsonSerializable(typeof(CodonUsageTableInput))]
[JsonSerializable(typeof(CodonUsageTableDto))]
[JsonSerializable(typeof(CodonChange))]
[JsonSerializable(typeof(CodonChange[]))]
[JsonSerializable(typeof(List<CodonChange>))]
[JsonSerializable(typeof(IReadOnlyList<CodonChange>))]
[JsonSerializable(typeof(OptimizationResultDto))]
[JsonSerializable(typeof(OptimizedSequenceResult))]
[JsonSerializable(typeof(RareCodon))]
[JsonSerializable(typeof(RareCodon[]))]
[JsonSerializable(typeof(List<RareCodon>))]
[JsonSerializable(typeof(IReadOnlyList<RareCodon>))]
[JsonSerializable(typeof(RareCodonsResult))]
[JsonSerializable(typeof(SimilarityResult))]

// CrisprDesigner DTOs (top-level records in Seqeron.Genomics.MolTools)
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CrisprSystemType))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.CrisprSystem))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.PamSite))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.PamSite[]))]
[JsonSerializable(typeof(List<global::Seqeron.Genomics.MolTools.PamSite>))]
[JsonSerializable(typeof(IReadOnlyList<global::Seqeron.Genomics.MolTools.PamSite>))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.GuideRnaParameters))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.GuideRnaParameters?))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.GuideRnaCandidate))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.GuideRnaCandidate[]))]
[JsonSerializable(typeof(List<global::Seqeron.Genomics.MolTools.GuideRnaCandidate>))]
[JsonSerializable(typeof(IReadOnlyList<global::Seqeron.Genomics.MolTools.GuideRnaCandidate>))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.OffTargetSite))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.OffTargetSite[]))]
[JsonSerializable(typeof(List<global::Seqeron.Genomics.MolTools.OffTargetSite>))]
[JsonSerializable(typeof(IReadOnlyList<global::Seqeron.Genomics.MolTools.OffTargetSite>))]
[JsonSerializable(typeof(PamSitesResult))]
[JsonSerializable(typeof(GuideRnasResult))]
[JsonSerializable(typeof(OffTargetsResult))]
[JsonSerializable(typeof(SpecificityResult))]

// ProbeDesigner DTOs (records nested in ProbeDesigner static class)
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.ProbeParameters))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.ProbeParameters?))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.ProbeType))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe?))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe[]))]
[JsonSerializable(typeof(List<global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe>))]
[JsonSerializable(typeof(IReadOnlyList<global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe>))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.TilingProbeSet))]
[JsonSerializable(typeof(global::Seqeron.Genomics.MolTools.ProbeDesigner.ProbeValidation))]
[JsonSerializable(typeof(ProbesResult))]
[JsonSerializable(typeof(MolecularBeaconResult))]
[JsonSerializable(typeof(OligoAnalysisResult))]
[JsonSerializable(typeof(ExtinctionCoefficientResult))]
[JsonSerializable(typeof(ConcentrationResult))]
public partial class AppJsonContext : JsonSerializerContext
{
}
