using System.Collections.Generic;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Mcp.MolTools.Models;

// PrimerDesigner result wrappers
public sealed record TmResult(double Tm);
public sealed record HomopolymerLengthResult(int Length);
public sealed record DinucleotideRepeatResult(int Repeats);
public sealed record HairpinPotentialResult(bool HasHairpin);
public sealed record PrimerDimerResult(bool HasDimer, int ComplementaryBases);
public sealed record ThreePrimeStabilityResult(double DeltaG);
public sealed record PrimerCandidateListResult(IReadOnlyList<PrimerCandidate> Candidates);

// RestrictionAnalyzer result wrappers
public sealed record EnzymeLookupResult(RestrictionEnzyme? Enzyme);
public sealed record EnzymeListResult(IReadOnlyList<RestrictionEnzyme> Enzymes);
public sealed record RestrictionSiteListResult(IReadOnlyList<RestrictionSite> Sites);
public sealed record DigestResult(IReadOnlyList<DigestFragment> Fragments);
public sealed record EnzymeCompatibilityPair(string Enzyme1, string Enzyme2, string CompatibleEnd);
public sealed record CompatibleEnzymesResult(IReadOnlyList<EnzymeCompatibilityPair> Pairs);
public sealed record EnzymeCompatibilityResult(bool Compatible);

// CodonUsageAnalyzer result wrappers
public sealed record CodonCountsResult(Dictionary<string, int> Counts);
public sealed record RscuResult(Dictionary<string, double> Rscu);
public sealed record CaiResult(double Cai);
public sealed record EncResult(double Enc);

// CodonOptimizer input/output wrappers
/// <summary>
/// Caller-supplied codon usage table for codon-optimizer tools. Either supply a
/// preset id (<c>EColiK12</c>, <c>Yeast</c>, <c>Human</c>) or a custom inline table
/// with RNA-alphabet codons (U, not T). When both are present, <c>Preset</c> wins.
/// </summary>
public sealed record CodonUsageTableInput(
    string? Preset = null,
    string? OrganismName = null,
    Dictionary<string, double>? CodonFrequencies = null);

public sealed record CodonChange(int Position, string Original, string Optimized);

public sealed record OptimizationResultDto(
    string OriginalSequence,
    string OptimizedSequence,
    string ProteinSequence,
    double OriginalCAI,
    double OptimizedCAI,
    double GcContentOriginal,
    double GcContentOptimized,
    int ChangedCodons,
    IReadOnlyList<CodonChange> Changes);

public sealed record OptimizedSequenceResult(string OptimizedSequence);

public sealed record RareCodon(int Position, string Codon, string AminoAcid, double Frequency);
public sealed record RareCodonsResult(IReadOnlyList<RareCodon> RareCodons);

public sealed record SimilarityResult(double Similarity);

public sealed record CodonUsageTableDto(
    string OrganismName,
    Dictionary<string, double> CodonFrequencies,
    Dictionary<string, string> CodonToAminoAcid);

// CrisprDesigner result wrappers
public sealed record PamSitesResult(IReadOnlyList<global::Seqeron.Genomics.MolTools.PamSite> Sites);
public sealed record GuideRnasResult(IReadOnlyList<global::Seqeron.Genomics.MolTools.GuideRnaCandidate> Guides);
public sealed record OffTargetsResult(IReadOnlyList<global::Seqeron.Genomics.MolTools.OffTargetSite> OffTargets);
public sealed record SpecificityResult(double Specificity);

// ProbeDesigner result wrappers
public sealed record ProbesResult(IReadOnlyList<global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe> Probes);
public sealed record MolecularBeaconResult(global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe? Probe);
public sealed record OligoAnalysisResult(double Tm, double GcContent, double MolecularWeight, double ExtinctionCoefficient);
public sealed record ExtinctionCoefficientResult(double ExtinctionCoefficient);
public sealed record ConcentrationResult(double ConcentrationMicromolar);
