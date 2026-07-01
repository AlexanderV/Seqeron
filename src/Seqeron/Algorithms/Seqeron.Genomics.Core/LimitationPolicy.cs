using System;
using System.Collections.Generic;
using System.Text;
using System.Threading;

namespace Seqeron.Genomics.Core;

/// <summary>
/// Library-wide policy controlling what happens when a caller reaches a documented limitation —
/// a code path that would otherwise hand back a value that is <b>not</b> the fully-resolved,
/// calibrated, validated "ideal" result (a defaulted FASTQ encoding, an uncalibrated disorder
/// confidence, a partial context++ score, an approximate miRNA* span, an unresolved SF1/SF2 call).
///
/// <para>In <see cref="LimitationMode.Strict"/> (the default) such a branch throws a
/// <see cref="SeqeronLimitationException"/> that names the limitation, what it is related to, and
/// how to obtain the result another way. In <see cref="LimitationMode.Permissive"/> the library
/// returns the honest best-effort value instead (the historical behaviour).</para>
///
/// <para>The flag is a process-wide singleton (<see cref="DefaultMode"/>) with an optional
/// async-local scoped override (<see cref="Use"/>/<see cref="UsePermissive"/>/<see cref="UseStrict"/>)
/// for a single region of code.</para>
/// </summary>
public static class LimitationPolicy
{
    private static readonly AsyncLocal<LimitationMode?> ScopedMode = new();

    /// <summary>
    /// The process-wide policy. Defaults to <see cref="LimitationMode.Moderate"/> — non-ideal-output
    /// branches throw, while correct-but-incomplete / narrower-contract results are allowed. Set to
    /// <see cref="LimitationMode.Strict"/> for only ideal-and-complete results, or
    /// <see cref="LimitationMode.Permissive"/> for the historical best-effort behaviour everywhere.
    /// </summary>
    public static LimitationMode DefaultMode { get; set; } = LimitationMode.Moderate;

    /// <summary>The effective mode for the current async flow: the scoped override if one is active,
    /// otherwise <see cref="DefaultMode"/>.</summary>
    public static LimitationMode CurrentMode => ScopedMode.Value ?? DefaultMode;

    /// <summary>True when the effective mode is <see cref="LimitationMode.Strict"/>.</summary>
    public static bool IsStrict => CurrentMode == LimitationMode.Strict;

    /// <summary>
    /// Whether a call guarded by <paramref name="limitationId"/> is allowed under the effective mode —
    /// i.e. the effective mode is at least as permissive as the limitation's
    /// <see cref="LimitationInfo.MinimumMode"/>.
    /// </summary>
    public static bool IsAllowed(string limitationId)
        => CurrentMode >= LimitationCatalog.Get(limitationId).MinimumMode;

    /// <summary>
    /// Pushes <paramref name="mode"/> as the effective mode for the current async flow until the
    /// returned token is disposed (restoring the previous effective mode). Nestable.
    /// </summary>
    public static IDisposable Use(LimitationMode mode) => new Scope(mode);

    /// <summary>Scopes <see cref="LimitationMode.Permissive"/> for a region of code.</summary>
    public static IDisposable UsePermissive() => Use(LimitationMode.Permissive);

    /// <summary>Scopes <see cref="LimitationMode.Strict"/> for a region of code.</summary>
    public static IDisposable UseStrict() => Use(LimitationMode.Strict);

    /// <summary>
    /// Enforces the policy for a guarded branch. In <see cref="LimitationMode.Strict"/> this throws
    /// the <see cref="SeqeronLimitationException"/> built from the catalog entry for
    /// <paramref name="limitationId"/>; in <see cref="LimitationMode.Permissive"/> it is a no-op and
    /// the caller proceeds to return the best-effort value.
    /// </summary>
    /// <param name="limitationId">The limitation unit id, e.g. <c>"PARSE-FASTQ-001"</c>.</param>
    /// <exception cref="SeqeronLimitationException">Thrown when the effective mode is more restrictive
    /// than the limitation's <see cref="LimitationInfo.MinimumMode"/>.</exception>
    /// <exception cref="KeyNotFoundException">Thrown when <paramref name="limitationId"/> is not in the catalog.</exception>
    public static void Enforce(string limitationId)
    {
        var info = LimitationCatalog.Get(limitationId);
        if (CurrentMode < info.MinimumMode)
            throw new SeqeronLimitationException(info);
    }

    private sealed class Scope : IDisposable
    {
        private readonly LimitationMode? _previous;
        private bool _disposed;

        public Scope(LimitationMode mode)
        {
            _previous = ScopedMode.Value;
            ScopedMode.Value = mode;
        }

        public void Dispose()
        {
            if (_disposed)
                return;
            _disposed = true;
            ScopedMode.Value = _previous;
        }
    }
}

/// <summary>
/// How permissive the library is about documented-limitation branches. Ordered from least to most
/// permissive: <see cref="Strict"/> &lt; <see cref="Moderate"/> &lt; <see cref="Permissive"/>. Each
/// limitation declares the minimum mode at which it is allowed (<see cref="LimitationInfo.MinimumMode"/>);
/// a call throws when the effective mode is more restrictive than that minimum.
/// </summary>
public enum LimitationMode
{
    /// <summary>
    /// Only the ideal <b>and complete</b> result. Throws on every documented-limitation branch —
    /// both the non-ideal-output branches (defaulted / uncalibrated / partial / approximate /
    /// unresolved) and the correct-but-incomplete / narrower-contract ones.
    /// </summary>
    Strict,

    /// <summary>
    /// The default. Throws on the non-ideal-output branches (a value that is not correct as returned),
    /// but <b>allows</b> the correct-but-incomplete / narrower-contract results (e.g. domain-level
    /// CheckM, a matrix score from a caller-supplied matrix, qualitative MGB rules).
    /// </summary>
    Moderate,

    /// <summary>Allows everything — the honest best-effort value (the historical behaviour).</summary>
    Permissive
}

/// <summary>Kind of limitation, mirroring the sections of <c>docs/Validation/LIMITATIONS.md</c>.</summary>
public enum LimitationCategory
{
    /// <summary>No algorithm or model can close it (physics / information theory).</summary>
    Irreducible,

    /// <summary>Needs a trained model / matrix / database that is gated, non-redistributable, or never measured.</summary>
    DataBlocked,

    /// <summary>A deliberate out-of-scope boundary; use the named reference tool, or supply the input.</summary>
    Scope
}

/// <summary>
/// One catalog entry describing a guarded limitation branch: the unit id, its category, the specific
/// branch, a one-line summary, what it is related to (the blocker), exhaustive workaround guidance,
/// and the path to the per-unit validation report.
/// </summary>
public sealed record LimitationInfo(
    string Id,
    LimitationCategory Category,
    LimitationMode MinimumMode,
    string Branch,
    string Summary,
    string RelatedTo,
    string Workaround,
    string ReportPath);

/// <summary>
/// Single source of truth for the runtime-guarded limitations. Text mirrors
/// <c>docs/Validation/LIMITATIONS.md</c> and the per-unit reports. Only branches that would return a
/// non-ideal value are listed here; irreducible / exact-narrower-contract limitations are
/// documentation-only and are not enforced at runtime.
/// </summary>
public static class LimitationCatalog
{
    private static readonly IReadOnlyDictionary<string, LimitationInfo> ByIdInternal = BuildCatalog();

    /// <summary>All guarded limitation entries, keyed by unit id.</summary>
    public static IReadOnlyDictionary<string, LimitationInfo> Entries => ByIdInternal;

    /// <summary>Gets the catalog entry for <paramref name="limitationId"/>.</summary>
    /// <exception cref="KeyNotFoundException">Thrown when the id is not a guarded limitation.</exception>
    public static LimitationInfo Get(string limitationId)
    {
        if (ByIdInternal.TryGetValue(limitationId, out var info))
            return info;
        throw new KeyNotFoundException(
            $"'{limitationId}' is not a runtime-guarded limitation. Known ids: {string.Join(", ", ByIdInternal.Keys)}.");
    }

    private static IReadOnlyDictionary<string, LimitationInfo> BuildCatalog()
    {
        var entries = new[]
        {
            new LimitationInfo(
                Id: "PARSE-FASTQ-001",
                Category: LimitationCategory.Irreducible,
                MinimumMode: LimitationMode.Permissive,
                Branch: "Phred+33 vs Phred+64 auto-detection on overlap-confined input",
                Summary: "The quality encoding cannot be auto-detected: every quality character lies in the " +
                         "Phred+33/Phred+64 overlap range (ASCII 64-74), so the input is information-theoretically " +
                         "ambiguous; the detector would otherwise default to Phred+33.",
                RelatedTo: "The overlapping ASCII ranges of the Sanger (Phred+33) and Illumina 1.3+ (Phred+64) " +
                           "encodings (Cock et al., 2010).",
                Workaround: "Decode with an explicit offset instead of Auto: pass QualityEncoding.Phred33 or " +
                            ".Phred64 to ParseQualityString / FastqParser.Parse. If the producing instrument is " +
                            "known the offset is fixed (modern Illumina 1.8+ is Phred+33). A single non-overlap read " +
                            "elsewhere in the file also resolves it via DetectEncoding(IEnumerable<string>).",
                ReportPath: "docs/Validation/reports/PARSE-FASTQ-001.md"),

            new LimitationInfo(
                Id: "CHROM-CENT-001",
                Category: LimitationCategory.DataBlocked,
                MinimumMode: LimitationMode.Permissive,
                Branch: "SF1-vs-SF2 separation (dimeric alpha-satellite)",
                Summary: "The array is dimeric (A->B) but SF1 (J1.J2) and SF2 (D1.D2) cannot be separated from the " +
                         "bundled CC0 reference — the result is Sf1OrSf2Dimeric, not a single SF.",
                RelatedTo: "SF1 and SF2 share the identical A->B dimeric box pattern; resolving them needs an " +
                           "SF-resolved consensus-monomer library that is not CC0/redistributable " +
                           "(HumAS-HMMER / logsdon-lab).",
                Workaround: "Pass an SF-resolved consensus-monomer reference to " +
                            "AssignSuprachromosomalFamily(sequence, reference) to obtain SF1 vs SF2; or classify with " +
                            "HumAS-HMMER / the T2T-CHM13 cenSat AS-HOR-SF annotation. The dimeric call itself, the HOR " +
                            "period, and the box pattern are exact and available under Permissive.",
                ReportPath: "docs/Validation/reports/CHROM-CENT-001.md"),

            new LimitationInfo(
                Id: "DISORDER-REGION-001",
                Category: LimitationCategory.DataBlocked,
                MinimumMode: LimitationMode.Permissive,
                Branch: "calibrated per-residue / per-region disorder confidence",
                Summary: "A calibrated disorder CONFIDENCE value is not available; PredictDisorder otherwise returns " +
                         "a heuristic per-region Confidence that is not a published calibrated standard.",
                RelatedTo: "No disorder predictor publishes a calibrated confidence standard; only the region " +
                           "boundaries follow the validated TOP-IDP threshold.",
                Workaround: "Use PredictDisorderRegions(...) for the validated TOP-IDP region boundaries and " +
                            "per-residue calls (no confidence) — the ideal result. For a calibrated confidence run an " +
                            "external predictor (IUPred2A / MobiDB-lite / AlphaFold pLDDT-derived). Enable Permissive " +
                            "to receive the heuristic Confidence as-is.",
                ReportPath: "docs/Validation/reports/DISORDER-REGION-001.md"),

            new LimitationInfo(
                Id: "MIRNA-TARGET-001",
                Category: LimitationCategory.Scope,
                MinimumMode: LimitationMode.Permissive,
                Branch: "full context++ score (one or more optional inputs not supplied)",
                Summary: "The context++ score is PARTIAL: one or more features (TA_3UTR / SPS / PCT / Len_ORF / " +
                         "ORF-8mer) were not supplied, so ContextScorePartial is not the full Agarwal (2015) " +
                         "context++ score (see OmittedFeatures).",
                RelatedTo: "TargetScan's compiled SPS table, the per-family PCT sigmoid parameters, and a default " +
                           "transcriptome are citation-required / not bundled.",
                Workaround: "Supply the missing inputs via ContextPlusPlusInputs: derive TA_3UTR with " +
                            "ComputeTa3Utr(miRna, 3'UTR set); provide SPS, ORF features, and a PctConservation " +
                            "(alignment + tree + sigmoid parameters) for PCT. With every feature supplied " +
                            "OmittedFeatures is empty and the full score is returned. Or use TargetScan directly. " +
                            "Permissive returns the partial score with OmittedFeatures listed.",
                ReportPath: "docs/Validation/reports/MIRNA-TARGET-001.md"),

            new LimitationInfo(
                Id: "MIRNA-CLEAVAGE-001",
                Category: LimitationCategory.Scope,
                MinimumMode: LimitationMode.Permissive,
                Branch: "3p-arm (miRNA*/star) cleavage boundary",
                Summary: "The opposite-arm (3p / star) boundary is a linear 2-nt-3'-overhang approximation, not the " +
                         "exact miRBase-annotated cut; only the 5p Drosha/Dicer cut reproduces miRBase exactly.",
                RelatedTo: "miRBase mature boundaries encode the dominant sequencing-read cut sites, not a " +
                           "deterministic fold + fixed-overhang rule.",
                Workaround: "Supply the miRBase mature-3p coordinates (MIMAT) or small-RNA-seq read pileups to define " +
                            "the star arm. The 5p mature / Drosha cut is exact and available. Permissive returns the " +
                            "approximate 3p span.",
                ReportPath: "docs/Validation/reports/MIRNA-CLEAVAGE-001.md"),

            // ── correct-but-incomplete / narrower-contract group (allowed in Moderate, blocked in Strict) ──

            new LimitationInfo(
                Id: "ONCO-MHC-001",
                Category: LimitationCategory.DataBlocked,
                MinimumMode: LimitationMode.Moderate,
                Branch: "matrix (SMM / BIMAS) pMHC binding prediction",
                Summary: "Matrix-based pMHC binding (BIMAS half-life / SMM IC50) is computed from a caller-supplied " +
                         "scoring matrix; the library bundles no validated, cross-verifiable matrix, and this is not " +
                         "the trained NetMHCpan-4.1 / MHCflurry predictor.",
                RelatedTo: "No redistributable, cross-verifiable SMM/BIMAS matrix is obtainable (the BIMAS server is a " +
                           "defunct CGI, the Parker 1994 table is paywalled, IEDB SMM matrices are non-commercial); " +
                           "NetMHCpan-4.1 / the MHCflurry presentation models are the vendors' trained models.",
                Workaround: "Supply a validated scoring matrix to PredictBindingHalfLifeBimas / PredictIc50Smm and call " +
                            "under Moderate/Permissive; or use the bundled MHCflurry NN predictor " +
                            "(MhcflurryAffinityPredictor); or run NetMHCpan-4.1 directly.",
                ReportPath: "docs/Validation/reports/ONCO-MHC-001.md"),

            new LimitationInfo(
                Id: "ONCO-IMMUNE-001",
                Category: LimitationCategory.DataBlocked,
                MinimumMode: LimitationMode.Moderate,
                Branch: "immune-cell deconvolution / ESTIMATE purity (no CIBERSORT-LM22 parity)",
                Summary: "Deconvolution runs against the bundled ABIS signature (or a caller-supplied matrix), NOT the " +
                         "CIBERSORT-LM22-identical matrix; the ESTIMATE->ABSOLUTE purity transform is calibrated on " +
                         "Affymetrix / ABSOLUTE data.",
                RelatedTo: "LM22 is distributed by Stanford under a no-redistribution licence; the ESTIMATE " +
                           "absolute-purity transform (Yoshihara 2013) is platform-calibrated.",
                Workaround: "Use the bundled ABIS matrix (LoadBundledAbisSignatureMatrix) or supply your own signature " +
                            "matrix to DeconvoluteImmuneCells / DeconvoluteImmuneCellsNuSvr under Moderate/Permissive; " +
                            "for exact CIBERSORT parity run CIBERSORT with LM22; apply EstimateTumorPurity only to " +
                            "Affymetrix-scale input.",
                ReportPath: "docs/Validation/reports/ONCO-IMMUNE-001.md"),

            new LimitationInfo(
                Id: "META-BIN-001",
                Category: LimitationCategory.DataBlocked,
                MinimumMode: LimitationMode.Moderate,
                Branch: "domain-level CheckM completeness/contamination (no lineage-specific refinement)",
                Summary: "Bin quality is estimated from the bundled DOMAIN-level marker sets only; per-lineage-specific " +
                         "CheckM marker refinement and the reference genome tree are not bundled.",
                RelatedTo: "The gated checkm_data DB and the TIGRFAM-defined markers (CC BY-SA 4.0, not redistributable " +
                           "here) are caller-supplied.",
                Workaround: "Pass lineage-specific marker sets / HMMs to EstimateBinQualityFromMarkers(Counts) under " +
                            "Moderate/Permissive; or run CheckM with its gated lineage DB for lineage-resolved " +
                            "completeness/contamination.",
                ReportPath: "docs/Validation/reports/META-BIN-001.md"),

            new LimitationInfo(
                Id: "PROBE-DESIGN-001",
                Category: LimitationCategory.Scope,
                MinimumMode: LimitationMode.Moderate,
                Branch: "MGB probe design (qualitative rules; no quantitative ΔTm)",
                Summary: "EvaluateMgbProbeDesign returns the citable qualitative MGB design rules (3'-MGB placement, " +
                         "12-20mer window) only; the quantitative MGB ΔTm is not computed.",
                RelatedTo: "The MGB ΔTm model (Kutyavin 2000 / MGB-Eclipse) is empirical/proprietary with no published " +
                           "closed-form expression.",
                Workaround: "Use the qualitative MGB rules under Moderate/Permissive; for a quantitative MGB ΔTm use a " +
                            "chemistry-specific tool (e.g. the manufacturer's MGB-Eclipse model).",
                ReportPath: "docs/Validation/reports/PROBE-DESIGN-001.md"),
        };

        var dict = new Dictionary<string, LimitationInfo>(StringComparer.Ordinal);
        foreach (var e in entries)
            dict[e.Id] = e;
        return dict;
    }
}

/// <summary>
/// Thrown (in <see cref="LimitationMode.Strict"/>) when a caller reaches a documented limitation — a
/// branch that would otherwise return a value that is not the ideal, confirmed result. Carries the
/// limitation id, category, branch, what it is related to, and exhaustive workaround guidance.
/// </summary>
public sealed class SeqeronLimitationException : InvalidOperationException
{
    /// <summary>The limitation unit id (e.g. <c>"PARSE-FASTQ-001"</c>).</summary>
    public string LimitationId { get; }

    /// <summary>The limitation category.</summary>
    public LimitationCategory Category { get; }

    /// <summary>The minimum <see cref="LimitationMode"/> at which the guarded call is allowed.</summary>
    public LimitationMode MinimumMode { get; }

    /// <summary>The specific guarded branch.</summary>
    public string Branch { get; }

    /// <summary>What blocks the ideal result.</summary>
    public string RelatedTo { get; }

    /// <summary>Exhaustive guidance on how to obtain the result another way.</summary>
    public string Workaround { get; }

    /// <summary>Path to the per-unit validation report.</summary>
    public string ReportPath { get; }

    /// <summary>Builds the exception from a catalog entry.</summary>
    public SeqeronLimitationException(LimitationInfo info)
        : base(Compose(info ?? throw new ArgumentNullException(nameof(info))))
    {
        LimitationId = info.Id;
        Category = info.Category;
        MinimumMode = info.MinimumMode;
        Branch = info.Branch;
        RelatedTo = info.RelatedTo;
        Workaround = info.Workaround;
        ReportPath = info.ReportPath;
    }

    private static string Compose(LimitationInfo info)
    {
        var sb = new StringBuilder();
        sb.Append('[').Append(info.Id).Append(" / ").Append(info.Category).Append("] ")
          .Append(info.Branch)
          .AppendLine(" is gated in Strict mode — Seqeron.Genomics returns only the ideal, confirmed result.");
        sb.AppendLine(info.Summary);
        sb.Append("Related to: ").AppendLine(info.RelatedTo);
        sb.Append("How to obtain the result: ").AppendLine(info.Workaround);
        sb.Append("See: ").AppendLine(info.ReportPath);
        sb.Append("Minimum mode that allows this call: ").Append(info.MinimumMode)
          .Append(". Set LimitationPolicy.DefaultMode to at least that, or wrap the call in ")
          .Append("`using (LimitationPolicy.Use(LimitationMode.").Append(info.MinimumMode).Append(")) { ... }`.");
        return sb.ToString();
    }
}
