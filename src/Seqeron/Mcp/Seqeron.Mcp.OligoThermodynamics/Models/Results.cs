// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.OligoThermodynamics, licensed under the terms stated in LICENSE (this directory).
namespace Seqeron.Mcp.OligoThermodynamics.Tools;

public sealed record PrimerCandidateListResult(IReadOnlyList<PrimerCandidate> Candidates);

/// <summary>
/// design_primers output: the best pair (same fields as <see cref="PrimerPairResult"/>; Primer3 PRIMER_PAIR_0_*)
/// plus <see cref="Pairs"/>, the ranked valid pairs (PRIMER_NUM_RETURN; empty when no pair qualifies). With a mispriming
/// library, <c>LibraryMispriming</c>/<c>LibraryMisprimingName</c> are PRIMER_PAIR_0_LIBRARY_MISPRIMING (score, entry).
/// </summary>
public sealed record DesignPrimersResult(
    PrimerCandidate? Forward,
    PrimerCandidate? Reverse,
    bool IsValid,
    string Message,
    int ProductSize,
    double? PairPenalty,
    double? ProductTm,
    double? ComplAnyTh,
    double? ComplEndTh,
    ProbeDesigner.Primer3Probe? InternalOligo,
    IReadOnlyList<PrimerPairResult> Pairs,
    double? ComplAny = null,
    double? ComplEnd = null,
    double? LibraryMispriming = null,
    string? LibraryMisprimingName = null,
    double? TemplateMispriming = null);

public sealed record ProbesResult(IReadOnlyList<global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe> Probes);

public sealed record Primer3ProbesResult(IReadOnlyList<global::Seqeron.Genomics.MolTools.ProbeDesigner.Primer3Probe> Probes);

public sealed record MolecularBeaconResult(global::Seqeron.Genomics.MolTools.ProbeDesigner.Probe? Probe);

public sealed record OligoAnalysisResult(double? Tm, double GcContent, double MolecularWeight, double ExtinctionCoefficient);

/// <summary>primer_dimer output: Primer3 alignment-mode PRIMER_PAIR_COMPL_END (flag, floor, exact) and PRIMER_PAIR_COMPL_ANY.</summary>
public sealed record PrimerDimerResult(bool HasDimer, int ComplementaryBases, double ComplEndScore, double ComplAnyScore = 0.0);
