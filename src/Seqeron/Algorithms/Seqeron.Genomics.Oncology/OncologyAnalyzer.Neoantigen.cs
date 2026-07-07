using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region GenerateNeoantigenPeptides

    /// <summary>
    /// Shortest MHC class I peptide length (8-mer) over which candidate neoantigen windows are enumerated.
    /// Source: Jurtz et al. (2017) NetMHCpan-4.0, <i>J. Immunol.</i> 199(9):3360–3368 — class I predictions
    /// are made for peptides of length 8–14; the NetMHCpan-4.1 web service offers 8/9/10/11/12/13/14-mer
    /// peptide options (https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/). pVACtools restricts the
    /// canonical class I neoantigen search to lengths 8–11 (Hundal et al. 2020, <i>Cancer Immunol. Res.</i>
    /// 8(3):409–420 — "8–11-mer for Class I MHC").
    /// </summary>
    public const int MhcClassIMinPeptideLength = 8;

    /// <summary>
    /// Longest MHC class I peptide length (14-mer) over which candidate neoantigen windows are enumerated by
    /// default. Source: the NetMHCpan-4.1 class I service offers 8/9/10/11/12/13/14-mer peptide options
    /// (Reynisson et al. 2020, <i>Nucleic Acids Res.</i> 48(W1):W449–W454, doi:10.1093/nar/gkaa379;
    /// https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), and the MHCflurry pan-allele peptide
    /// encoding admits up to 15-mers (<see cref="MhcflurryAffinityPredictor.PeptideMaxLength"/> = 15; O'Donnell
    /// et al. 2020, <i>Cell Systems</i> 11(1):42–48). The default upper bound follows NetMHCpan-4.1's
    /// class I window (8–14); callers may pass <c>maxLength</c> up to 15 to reach the MHCflurry encoding ceiling.
    /// (The earlier 8–11 default mirrored pVACtools' canonical search — Hundal et al. 2020, <i>Cancer Immunol.
    /// Res.</i> 8(3):409–420 — and remains reachable by passing <c>maxLength: 11</c>.)
    /// </summary>
    public const int MhcClassIMaxPeptideLength = 14;

    /// <summary>
    /// A candidate neoantigen peptide: a fixed-length window of the mutant protein that spans the mutated
    /// residue, paired with the wild-type peptide occupying the same coordinates (the agretope). Mutant and
    /// wild-type peptides differ only at the substituted residue (Hundal et al. 2020, <i>Cancer Immunol.
    /// Res.</i> 8(3):409–420; agretopicity is the differential binding of these two — Wells et al. 2020,
    /// <i>Cell</i> 183(3):818–834).
    /// </summary>
    /// <param name="Length">Peptide length k (number of residues), in [<see cref="MhcClassIMinPeptideLength"/>,
    /// <see cref="MhcClassIMaxPeptideLength"/>].</param>
    /// <param name="StartPosition">1-based position in the protein of the first residue of the window.</param>
    /// <param name="MutantPeptide">The k-mer taken from the mutant protein (carries the substituted residue).</param>
    /// <param name="WildTypePeptide">The k-mer at the same coordinates in the wild-type protein (the agretope).</param>
    /// <param name="MutationOffset">0-based offset of the mutated residue within the peptide window.</param>
    public readonly record struct NeoantigenPeptide(
        int Length,
        int StartPosition,
        string MutantPeptide,
        string WildTypePeptide,
        int MutationOffset);

    /// <summary>
    /// Generates the candidate MHC class I neoantigen peptide windows arising from a single somatic missense
    /// (amino-acid substitution) mutation. For each peptide length k the method enumerates every length-k
    /// window of the mutant protein that <b>spans</b> the mutated residue, and pairs it with the wild-type
    /// window at the same coordinates (the agretope). This is the windowing step of neoantigen prediction as
    /// defined by pVACtools (Hundal et al. 2020) and the 21-mer ± 10-flank construction of ProGeo-neo (Li et
    /// al. 2020, <i>BMC Med. Genomics</i> 13:52): a 21-mer with 10 residues flanking the substitution on each
    /// side contains exactly the 8–11-mer windows that overlap the mutation.
    /// <para>
    /// Binding affinity / IC50 is NOT computed here — that requires a trained MHC-binding model. The bundled
    /// <see cref="MhcflurryAffinityPredictor"/> port scores any window of length 8–15 (with caller-supplied
    /// weights); the default upper bound here is 14, matching the NetMHCpan-4.1 class I peptide window.
    /// </para>
    /// </summary>
    /// <param name="wildTypeProtein">The wild-type (reference) protein sequence (one-letter amino-acid codes).</param>
    /// <param name="mutantResidue">The substituted (mutant) amino acid (one-letter code).</param>
    /// <param name="mutationPosition">1-based position of the substituted residue within the protein.</param>
    /// <param name="minLength">Minimum peptide length (default <see cref="MhcClassIMinPeptideLength"/>).</param>
    /// <param name="maxLength">Maximum peptide length (default <see cref="MhcClassIMaxPeptideLength"/>).</param>
    /// <returns>
    /// All candidate peptides, ordered by length ascending then by start position ascending. Empty only if no
    /// window of any requested length fits within the protein bounds while spanning the mutation.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="wildTypeProtein"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// The protein is empty; <paramref name="mutantResidue"/> is not a single character; the wild-type residue
    /// at the mutation position already equals the mutant residue (not a substitution); or the length range is
    /// invalid (min &gt; max, or min &lt; 1).
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="mutationPosition"/> is outside [1, protein length].
    /// </exception>
    public static IReadOnlyList<NeoantigenPeptide> GenerateNeoantigenPeptides(
        string wildTypeProtein,
        char mutantResidue,
        int mutationPosition,
        int minLength = MhcClassIMinPeptideLength,
        int maxLength = MhcClassIMaxPeptideLength)
    {
        ArgumentNullException.ThrowIfNull(wildTypeProtein);

        if (wildTypeProtein.Length == 0)
        {
            throw new ArgumentException("Protein sequence must be non-empty.", nameof(wildTypeProtein));
        }

        if (minLength < 1)
        {
            throw new ArgumentException($"Minimum peptide length must be at least 1; got {minLength}.", nameof(minLength));
        }

        if (maxLength < minLength)
        {
            throw new ArgumentException(
                $"Maximum peptide length ({maxLength}) must be ≥ minimum ({minLength}).", nameof(maxLength));
        }

        if (mutationPosition < 1 || mutationPosition > wildTypeProtein.Length)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutationPosition),
                mutationPosition,
                $"Mutation position must be in [1, {wildTypeProtein.Length}].");
        }

        int mutationIndex = mutationPosition - 1; // 0-based index of the substituted residue.
        char wildTypeResidue = wildTypeProtein[mutationIndex];
        if (wildTypeResidue == mutantResidue)
        {
            throw new ArgumentException(
                $"Mutant residue '{mutantResidue}' equals the wild-type residue at position {mutationPosition}; " +
                "a missense mutation requires a different amino acid.", nameof(mutantResidue));
        }

        // The mutant protein differs from the wild type only at the substituted residue.
        char[] mutantChars = wildTypeProtein.ToCharArray();
        mutantChars[mutationIndex] = mutantResidue;
        string mutantProtein = new(mutantChars);

        var peptides = new List<NeoantigenPeptide>();
        int proteinLength = wildTypeProtein.Length;

        // For each requested length k, a length-k window starts at 0-based index s and covers [s, s+k-1].
        // It spans the mutation iff s ≤ mutationIndex ≤ s+k-1, i.e. s ∈ [mutationIndex-k+1, mutationIndex],
        // additionally clamped to the protein bounds s ∈ [0, proteinLength-k]. This is exactly the set of
        // 8–11-mers contained in the 21-mer ±10-flank window of pVACtools / ProGeo-neo.
        for (int k = minLength; k <= maxLength; k++)
        {
            if (k > proteinLength)
            {
                continue; // No window of this length fits in the protein.
            }

            int firstStart = Math.Max(0, mutationIndex - k + 1);
            int lastStart = Math.Min(mutationIndex, proteinLength - k);
            for (int start = firstStart; start <= lastStart; start++)
            {
                string mutantPeptide = mutantProtein.Substring(start, k);
                string wildTypePeptide = wildTypeProtein.Substring(start, k);
                int offset = mutationIndex - start; // 0-based offset of the mutation within the window.
                peptides.Add(new NeoantigenPeptide(k, start + 1, mutantPeptide, wildTypePeptide, offset));
            }
        }

        return peptides;
    }

    #endregion


    #region ClassifyMhcBinding

    /// <summary>
    /// MHC molecule class for peptide-binding classification. Class I and class II have different accepted
    /// peptide-length ranges and different %Rank cutoffs (Reynisson et al. 2020, <i>Nucleic Acids Res.</i>
    /// 48(W1):W449–W454).
    /// </summary>
    public enum MhcClass
    {
        /// <summary>MHC class I (HLA-A/B/C). Presented peptide length 8–14 (NetMHCpan-4.1 class I window).</summary>
        ClassI,

        /// <summary>MHC class II (HLA-DR/DQ/DP). Peptide length 13–25.</summary>
        ClassII
    }

    /// <summary>
    /// Binding-strength category assigned to a peptide–MHC pair from a caller-supplied predicted affinity
    /// (IC50) or %Rank. Categories follow the IEDB / NetMHCpan strong-/weak-binder convention.
    /// </summary>
    public enum BindingStrength
    {
        /// <summary>Strong binder (IC50 &lt; 50 nM, or class I %Rank &lt; 0.5% / class II %Rank &lt; 2%).</summary>
        Strong,

        /// <summary>Weak (intermediate) binder (IC50 &lt; 500 nM, or class I %Rank &lt; 2% / class II %Rank &lt; 10%).</summary>
        Weak,

        /// <summary>Not a binder (above the weak-binder cutoff).</summary>
        NonBinder
    }

    /// <summary>
    /// IC50 (nM) below which a peptide–MHC pair is a strong (high-affinity) binder. Source: Sette et al.
    /// (1994), <i>J. Immunol.</i> 153(12):5586–5592 — "an affinity threshold of approximately 500 nM
    /// (preferably 50 nM or less) apparently determines the capacity" to elicit a CTL response; the IEDB
    /// states "Peptides with IC50 values &lt;50 nM are considered high affinity". Strict inequality.
    /// </summary>
    public const double StrongBinderIc50Nm = 50.0;

    /// <summary>
    /// IC50 (nM) below which a peptide–MHC pair is at least a weak (intermediate-affinity) binder. Source:
    /// IEDB — "&lt;500 nM intermediate affinity"; Sette et al. (1994) (≈500 nM threshold); corroborated by
    /// Roomp, Antes &amp; Lengauer (2010), <i>BMC Bioinformatics</i> 11:90 (500 nM binder demarcation). Strict
    /// inequality.
    /// </summary>
    public const double WeakBinderIc50Nm = 500.0;

    /// <summary>
    /// Class I %Rank below which a peptide is a strong binder. Source: Reynisson et al. (2020) — "by default,
    /// %Rank &lt; 0.5% and %Rank &lt; 2% thresholds are considered for detecting SBs and WBs for class I".
    /// Strict inequality.
    /// </summary>
    public const double ClassIStrongBinderRankPercent = 0.5;

    /// <summary>
    /// Class I %Rank below which a peptide is at least a weak binder. Source: Reynisson et al. (2020) (class I
    /// WB &lt; 2%). Strict inequality.
    /// </summary>
    public const double ClassIWeakBinderRankPercent = 2.0;

    /// <summary>
    /// Class II %Rank below which a peptide is a strong binder. Source: Reynisson et al. (2020) — "%Rank &lt; 2%
    /// and %Rank &lt; 10%, for SBs and WBs for class II". Strict inequality.
    /// </summary>
    public const double ClassIIStrongBinderRankPercent = 2.0;

    /// <summary>
    /// Class II %Rank below which a peptide is at least a weak binder. Source: Reynisson et al. (2020) (class
    /// II WB &lt; 10%). Strict inequality.
    /// </summary>
    public const double ClassIIWeakBinderRankPercent = 10.0;

    /// <summary>
    /// Minimum accepted peptide length for MHC class II. Source: IEDB MHC class II tool description —
    /// "Peptides binding to MHC class II molecules ... typically range between 13 and 25 amino acids long".
    /// </summary>
    public const int MhcClassIIMinPeptideLength = 13;

    /// <summary>
    /// Maximum accepted peptide length for MHC class II. Source: IEDB MHC class II tool description (13–25).
    /// </summary>
    public const int MhcClassIIMaxPeptideLength = 25;

    /// <summary>
    /// Largest finite %Rank value (a %Rank is a percentile and must lie in [0, 100]). Source: Reynisson et
    /// al. (2020) — %Rank is "the top X% scores from random natural peptides".
    /// </summary>
    private const double MaxRankPercent = 100.0;

    /// <summary>
    /// Classifies a caller-supplied predicted peptide–MHC binding affinity (IC50 in nanomolar) into
    /// <see cref="BindingStrength.Strong"/> (IC50 &lt; 50 nM), <see cref="BindingStrength.Weak"/>
    /// (IC50 &lt; 500 nM), or <see cref="BindingStrength.NonBinder"/> (IC50 ≥ 500 nM). The cutoffs are the
    /// standard IEDB / NetMHCpan convention (Sette et al. 1994; IEDB) and the boundaries are strict
    /// inequalities, so 50 nM is weak (not strong) and 500 nM is a non-binder (not weak).
    /// <para>
    /// This method does NOT predict the IC50 — that requires a trained MHC-binding model (e.g. NetMHCpan) and
    /// is caller-supplied / out of scope (ONCO-MHC-001). It only classifies a supplied value.
    /// </para>
    /// </summary>
    /// <param name="ic50Nm">Predicted half-maximal inhibitory concentration in nM (must be &gt; 0).</param>
    /// <returns>The binding-strength category.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="ic50Nm"/> is not a finite value greater than 0 (IC50 is a positive concentration).
    /// </exception>
    public static BindingStrength ClassifyBindingAffinity(double ic50Nm)
    {
        if (double.IsNaN(ic50Nm) || double.IsInfinity(ic50Nm) || ic50Nm <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(ic50Nm), ic50Nm, "IC50 must be a finite concentration greater than 0 nM.");
        }

        if (ic50Nm < StrongBinderIc50Nm)
        {
            return BindingStrength.Strong;
        }

        return ic50Nm < WeakBinderIc50Nm ? BindingStrength.Weak : BindingStrength.NonBinder;
    }

    /// <summary>
    /// Classifies a caller-supplied predicted %Rank into <see cref="BindingStrength"/> using the NetMHCpan-4.1
    /// default cutoffs (Reynisson et al. 2020): class I — strong &lt; 0.5%, weak &lt; 2%; class II — strong
    /// &lt; 2%, weak &lt; 10%. The boundaries are strict inequalities (a value exactly at a cutoff falls into
    /// the weaker category).
    /// <para>
    /// This method does NOT predict the %Rank — the trained model is caller-supplied / out of scope.
    /// </para>
    /// </summary>
    /// <param name="percentRank">Predicted %Rank as a percentile in [0, 100].</param>
    /// <param name="mhcClass">MHC class selecting the cutoff set.</param>
    /// <returns>The binding-strength category.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="percentRank"/> is NaN or outside [0, 100] (a %Rank is a percentile).
    /// </exception>
    public static BindingStrength ClassifyBindingRank(double percentRank, MhcClass mhcClass)
    {
        if (double.IsNaN(percentRank) || percentRank < 0.0 || percentRank > MaxRankPercent)
        {
            throw new ArgumentOutOfRangeException(
                nameof(percentRank), percentRank, "%Rank must be a percentile in [0, 100].");
        }

        double strongCutoff = mhcClass == MhcClass.ClassI
            ? ClassIStrongBinderRankPercent
            : ClassIIStrongBinderRankPercent;
        double weakCutoff = mhcClass == MhcClass.ClassI
            ? ClassIWeakBinderRankPercent
            : ClassIIWeakBinderRankPercent;

        if (percentRank < strongCutoff)
        {
            return BindingStrength.Strong;
        }

        return percentRank < weakCutoff ? BindingStrength.Weak : BindingStrength.NonBinder;
    }

    /// <summary>
    /// Determines whether <paramref name="length"/> is a valid presented-peptide length for the given MHC
    /// class: class I 8–14 (the NetMHCpan-4.1 class I peptide window — Reynisson et al. 2020, <i>Nucleic Acids
    /// Res.</i> 48(W1):W449–W454; matching <see cref="MhcClassIMinPeptideLength"/>/<see
    /// cref="MhcClassIMaxPeptideLength"/>), class II 13–25 (IEDB MHC class II tool description). Both bounds
    /// are inclusive.
    /// </summary>
    /// <param name="length">Peptide length (residue count).</param>
    /// <param name="mhcClass">MHC class selecting the accepted length range.</param>
    /// <returns><see langword="true"/> iff <paramref name="length"/> is within the class's accepted range.</returns>
    public static bool IsValidPeptideLength(int length, MhcClass mhcClass)
    {
        return mhcClass == MhcClass.ClassI
            ? length >= MhcClassIMinPeptideLength && length <= MhcClassIMaxPeptideLength
            : length >= MhcClassIIMinPeptideLength && length <= MhcClassIIMaxPeptideLength;
    }

    /// <summary>
    /// Classifies a candidate peptide–MHC pair end-to-end: a peptide whose length is not valid for the MHC
    /// class is not a presentable candidate and is classified <see cref="BindingStrength.NonBinder"/>
    /// regardless of affinity; otherwise the supplied IC50 is classified by
    /// <see cref="ClassifyBindingAffinity(double)"/>. Thin convenience wrapper over the length gate
    /// (<see cref="IsValidPeptideLength(int, MhcClass)"/>) and the affinity classifier.
    /// </summary>
    /// <param name="peptideLength">Peptide length (residue count).</param>
    /// <param name="ic50Nm">Caller-supplied predicted IC50 in nM (must be &gt; 0).</param>
    /// <param name="mhcClass">MHC class.</param>
    /// <returns>The binding-strength category, or <see cref="BindingStrength.NonBinder"/> for invalid length.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ic50Nm"/> is not finite and &gt; 0.</exception>
    public static BindingStrength ClassifyMhcBinding(int peptideLength, double ic50Nm, MhcClass mhcClass)
    {
        if (!IsValidPeptideLength(peptideLength, mhcClass))
        {
            return BindingStrength.NonBinder;
        }

        return ClassifyBindingAffinity(ic50Nm);
    }

    #endregion


    #region Matrix-based pMHC binding prediction (ONCO-MHC-001)

    /// <summary>
    /// Scoring convention for a position-specific peptide–MHC binding matrix. Selects how a per-residue
    /// matrix is combined into a prediction and what the prediction means.
    /// </summary>
    public enum PmhcScoringMethod
    {
        /// <summary>
        /// BIMAS / Parker et al. (1994) "stabilized matrix" convention: the prediction is the
        /// <b>product</b> of the position-specific coefficients times a final constant, and (for HLA-A2) it
        /// estimates the <b>half-time of dissociation</b> (T½, arbitrary BIMAS units) of the HLA–peptide
        /// complex at 37 °C, pH 6.5. Higher T½ = stronger binder. Source: BIMAS HLA peptide-motif-search
        /// scoring documentation; Parker, Bednarek &amp; Coligan (1994), <i>J. Immunol.</i> 152(1):163–175.
        /// </summary>
        BimasHalfLife,

        /// <summary>
        /// SMM (Peters &amp; Sette 2005) / IEDB convention: the prediction is the <b>sum</b> of the
        /// position-specific values plus an intercept, giving a <c>log50k</c> score, which is converted to
        /// an IC50 (nM) by <c>IC50 = 50000^(1 − score)</c>. Lower IC50 = stronger binder. Source:
        /// Peters &amp; Sette (2005), <i>BMC Bioinformatics</i> 6:132; IEDB <c>log50k = 1 − log(IC50)/log(50000)</c>.
        /// </summary>
        SmmIc50
    }

    /// <summary>
    /// A position-specific scoring matrix (PSSM) for peptide–MHC binding: one per-position row mapping each
    /// amino-acid residue (single-letter code) to a numeric value, plus a single scalar
    /// <see cref="FinalConstant"/>. The interpretation of the values and the constant depends on the
    /// <see cref="PmhcScoringMethod"/> used to score with it:
    /// <list type="bullet">
    /// <item><description><see cref="PmhcScoringMethod.BimasHalfLife"/>: each value is a multiplicative
    /// coefficient (default 1.0 ≡ neutral), and <see cref="FinalConstant"/> is the BIMAS per-allele final
    /// multiplier; the score is their product (Parker 1994 / BIMAS).</description></item>
    /// <item><description><see cref="PmhcScoringMethod.SmmIc50"/>: each value is an additive log50k
    /// contribution (default 0.0 ≡ neutral), and <see cref="FinalConstant"/> is the SMM intercept; the score
    /// is their sum (Peters &amp; Sette 2005).</description></item>
    /// </list>
    /// <para>
    /// <b>This matrix is caller-supplied, not bundled.</b> No redistributable, cross-verifiable trained
    /// HLA coefficient matrix was obtainable for embedding: the public BIMAS coefficient files are served
    /// only by a now-defunct dynamic CGI (not archived) and the Parker (1994) table is behind a paywall;
    /// the IEDB SMM matrices carry a non-commercial / no-redistribution licence (like CIBERSORT LM22).
    /// The library therefore implements the published <i>scoring rules</i> and a <see cref="LoadScoringMatrix"/>
    /// loader, and the caller provides the coefficient values under their own licence (ONCO-MHC-001).
    /// </para>
    /// </summary>
    /// <param name="Rows">Per-position residue→value maps; <c>Rows[i]</c> is position <c>i</c> (0-based).</param>
    /// <param name="FinalConstant">The BIMAS final multiplier (product convention) or SMM intercept (sum convention).</param>
    public readonly record struct PmhcScoringMatrix(
        IReadOnlyList<IReadOnlyDictionary<char, double>> Rows,
        double FinalConstant);

    /// <summary>
    /// Largest possible IC50 (nM) produced by the SMM transform, reached at score 0:
    /// <c>50000^(1 − 0) = 50000</c>. Source: IEDB <c>log50k = 1 − log(IC50)/log(50000)</c> (Peters &amp;
    /// Sette 2005, <i>BMC Bioinformatics</i> 6:132), so IC50 = 50000^(1 − score).
    /// </summary>
    public const double SmmIc50Base = 50000.0;

    /// <summary>
    /// Default neutral BIMAS coefficient: an amino acid "not known to make either a favorable or unfavorable
    /// contribution" has a coefficient of exactly 1.0 (multiplicative identity). Source: BIMAS scoring
    /// documentation — ambiguous/unlisted residues "given a coefficient of 1.00 … leaves the score unchanged".
    /// </summary>
    private const double BimasNeutralCoefficient = 1.0;

    /// <summary>
    /// Default neutral SMM contribution: 0.0 (additive identity), so an unlisted residue does not move the
    /// log50k score. Source: SMM is an additive position-specific model (Peters &amp; Sette 2005).
    /// </summary>
    private const double SmmNeutralContribution = 0.0;

    /// <summary>
    /// Loads a caller-supplied position-specific peptide–MHC scoring matrix from text rows. Each non-blank,
    /// non-comment line is one matrix position and lists residue values as whitespace- or comma-separated
    /// <c>RESIDUE=VALUE</c> tokens (e.g. <c>L=2.5 M=1.8 V=1.1</c>); the very first token of the form
    /// <c>CONST=VALUE</c> (or a line <c>CONST VALUE</c>) sets <see cref="PmhcScoringMatrix.FinalConstant"/>.
    /// Lines beginning with <c>#</c> are comments. Residues are upper-cased.
    /// <para>
    /// <b>Licence / provenance.</b> The matrix VALUES are NOT supplied by this library — see
    /// <see cref="PmhcScoringMatrix"/>. Obtain coefficients from a source you are licensed to use (the
    /// public-domain BIMAS/Parker 1994 HLA-A2 table, or an IEDB SMM matrix under its non-commercial licence)
    /// and pass them here. The library bundles only the published scoring <i>rules</i>.
    /// </para>
    /// </summary>
    /// <param name="lines">Matrix text lines: optional <c>CONST=…</c> plus one line per position.</param>
    /// <returns>The parsed scoring matrix.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="lines"/> is null.</exception>
    /// <exception cref="FormatException">A token is malformed or a value is non-numeric.</exception>
    public static PmhcScoringMatrix LoadScoringMatrix(IEnumerable<string> lines)
    {
        ArgumentNullException.ThrowIfNull(lines);

        var rows = new List<IReadOnlyDictionary<char, double>>();
        double finalConstant = BimasNeutralCoefficient; // 1.0: neutral for product; harmless if overridden.
        bool constantSet = false;

        foreach (string raw in lines)
        {
            if (string.IsNullOrWhiteSpace(raw))
            {
                continue;
            }

            string line = raw.Trim();
            if (line.StartsWith('#'))
            {
                continue;
            }

            string[] tokens = line.Split(new[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
            var row = new Dictionary<char, double>(tokens.Length);
            foreach (string token in tokens)
            {
                int eq = token.IndexOf('=');
                string keyPart;
                string valuePart;
                if (eq >= 0)
                {
                    keyPart = token[..eq];
                    valuePart = token[(eq + 1)..];
                }
                else
                {
                    // Allow a bare "CONST 1234" pair only when it is the line's sole content handled below.
                    throw new FormatException($"Malformed matrix token '{token}': expected RESIDUE=VALUE or CONST=VALUE.");
                }

                if (!double.TryParse(valuePart, System.Globalization.NumberStyles.Float,
                        System.Globalization.CultureInfo.InvariantCulture, out double value))
                {
                    throw new FormatException($"Non-numeric value in token '{token}'.");
                }

                if (string.Equals(keyPart, "CONST", StringComparison.OrdinalIgnoreCase))
                {
                    finalConstant = value;
                    constantSet = true;
                    continue;
                }

                if (keyPart.Length != 1)
                {
                    throw new FormatException($"Residue key '{keyPart}' must be a single amino-acid letter.");
                }

                row[char.ToUpperInvariant(keyPart[0])] = value;
            }

            if (row.Count > 0)
            {
                rows.Add(row);
            }
        }

        // If no CONST token was given, leave the additive/multiplicative identity (1.0); callers using SMM
        // can pass CONST=<intercept>. constantSet is retained only to make the intent explicit.
        _ = constantSet;
        return new PmhcScoringMatrix(rows, finalConstant);
    }

    /// <summary>
    /// Predicts the BIMAS half-time of dissociation (T½) of a peptide–MHC complex as the
    /// <b>product</b> of the position-specific coefficients times the matrix's final constant:
    /// <c>T½ = FinalConstant · ∏_i matrix.Rows[i][peptide[i]]</c>. A residue absent from a position's row
    /// contributes the neutral coefficient 1.0 (no effect), exactly as BIMAS treats ambiguous residues.
    /// Source: BIMAS scoring documentation — "The initial (running) score is set to 1.0 … the running score
    /// is then multiplied by the coefficient for that amino acid type, at that position … The resulting
    /// running score is multiplied by a final constant to yield an estimate of the half time of
    /// disassociation"; Parker, Bednarek &amp; Coligan (1994), <i>J. Immunol.</i> 152(1):163–175 ("calculated
    /// by multiplying together the corresponding coefficients").
    /// </summary>
    /// <param name="peptide">The peptide; its length must equal <c>matrix.Rows.Count</c>.</param>
    /// <param name="matrix">A caller-supplied BIMAS-convention coefficient matrix (see <see cref="PmhcScoringMatrix"/>).</param>
    /// <returns>The predicted half-time of dissociation (BIMAS units; higher = stronger binding).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="peptide"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix has no rows, or the peptide length ≠ row count.</exception>
    public static double PredictBindingHalfLifeBimas(string peptide, PmhcScoringMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(peptide);
        ValidateMatrixAgainstPeptide(peptide, matrix);

        // Matrix score from a caller-supplied matrix — not a bundled validated matrix nor the trained
        // NetMHCpan/MHCflurry predictor. Strict mode throws; Moderate/Permissive allow it.
        Seqeron.Genomics.Core.LimitationPolicy.Enforce("ONCO-MHC-001");

        double score = 1.0; // BIMAS: "The initial (running) score is set to 1.0."
        for (int i = 0; i < peptide.Length; i++)
        {
            char residue = char.ToUpperInvariant(peptide[i]);
            double coefficient = matrix.Rows[i].TryGetValue(residue, out double c) ? c : BimasNeutralCoefficient;
            score *= coefficient;
        }

        return score * matrix.FinalConstant;
    }

    /// <summary>
    /// Predicts the IC50 (nM) of a peptide–MHC pair under the SMM / IEDB convention: the log50k score is the
    /// <b>sum</b> of the position-specific values plus the matrix intercept, and the IC50 is recovered by
    /// <c>IC50 = 50000^(1 − score)</c>. A residue absent from a position's row contributes 0 (no effect).
    /// Source: IEDB linearisation <c>log50k = 1 − log(IC50)/log(50000)</c> (inverted gives
    /// <c>IC50 = 50000^(1 − log50k)</c>); Peters &amp; Sette (2005), <i>BMC Bioinformatics</i> 6:132 (SMM is
    /// an additive position-specific model fitted on log-transformed IC50).
    /// </summary>
    /// <param name="peptide">The peptide; its length must equal <c>matrix.Rows.Count</c>.</param>
    /// <param name="matrix">A caller-supplied SMM-convention matrix (values = additive log50k contributions; <see cref="PmhcScoringMatrix.FinalConstant"/> = intercept).</param>
    /// <returns>The predicted IC50 in nM (lower = stronger binding).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="peptide"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix has no rows, or the peptide length ≠ row count.</exception>
    public static double PredictIc50Smm(string peptide, PmhcScoringMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(peptide);
        ValidateMatrixAgainstPeptide(peptide, matrix);

        // Matrix score from a caller-supplied matrix — not a bundled validated matrix nor the trained
        // NetMHCpan/MHCflurry predictor. Strict mode throws; Moderate/Permissive allow it.
        Seqeron.Genomics.Core.LimitationPolicy.Enforce("ONCO-MHC-001");

        double score = matrix.FinalConstant; // SMM intercept.
        for (int i = 0; i < peptide.Length; i++)
        {
            char residue = char.ToUpperInvariant(peptide[i]);
            score += matrix.Rows[i].TryGetValue(residue, out double v) ? v : SmmNeutralContribution;
        }

        // IC50 = 50000^(1 - score)  ⇔  log50k = 1 - log(IC50)/log(50000).
        return Math.Pow(SmmIc50Base, 1.0 - score);
    }

    /// <summary>
    /// End-to-end SMM prediction → classification: predicts the IC50 with <see cref="PredictIc50Smm"/> and
    /// classifies it into <see cref="BindingStrength"/> via the existing <see cref="ClassifyBindingAffinity"/>
    /// (strong &lt; 50 nM, weak &lt; 500 nM, else non-binder). This chains the matrix-based predictor into the
    /// established threshold classifier so prediction and classification compose end-to-end (ONCO-MHC-001).
    /// </summary>
    /// <param name="peptide">The peptide; length must equal <c>matrix.Rows.Count</c>.</param>
    /// <param name="matrix">A caller-supplied SMM-convention matrix.</param>
    /// <returns>A tuple of the predicted IC50 (nM) and its binding-strength category.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="peptide"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix has no rows, or the peptide length ≠ row count.</exception>
    public static (double Ic50Nm, BindingStrength Strength) PredictAndClassifySmm(
        string peptide, PmhcScoringMatrix matrix)
    {
        double ic50 = PredictIc50Smm(peptide, matrix);
        return (ic50, ClassifyBindingAffinity(ic50));
    }

    private static void ValidateMatrixAgainstPeptide(string peptide, PmhcScoringMatrix matrix)
    {
        if (matrix.Rows is null || matrix.Rows.Count == 0)
        {
            throw new ArgumentException("Scoring matrix has no position rows.", nameof(matrix));
        }

        if (peptide.Length != matrix.Rows.Count)
        {
            throw new ArgumentException(
                $"Peptide length {peptide.Length} does not match the matrix position count {matrix.Rows.Count}.",
                nameof(peptide));
        }
    }

    #endregion

}
