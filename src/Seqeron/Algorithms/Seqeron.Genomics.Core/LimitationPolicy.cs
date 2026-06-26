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
    /// The process-wide policy. Defaults to <see cref="LimitationMode.Strict"/> — the library
    /// returns only the ideal, confirmed result and throws on any non-ideal branch.
    /// </summary>
    public static LimitationMode DefaultMode { get; set; } = LimitationMode.Strict;

    /// <summary>The effective mode for the current async flow: the scoped override if one is active,
    /// otherwise <see cref="DefaultMode"/>.</summary>
    public static LimitationMode CurrentMode => ScopedMode.Value ?? DefaultMode;

    /// <summary>True when the effective mode is <see cref="LimitationMode.Strict"/>.</summary>
    public static bool IsStrict => CurrentMode == LimitationMode.Strict;

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
    /// <exception cref="SeqeronLimitationException">Thrown when the effective mode is Strict.</exception>
    /// <exception cref="KeyNotFoundException">Thrown when <paramref name="limitationId"/> is not in the catalog.</exception>
    public static void Enforce(string limitationId)
    {
        if (!IsStrict)
            return;
        throw new SeqeronLimitationException(LimitationCatalog.Get(limitationId));
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

/// <summary>How the library handles a documented-limitation branch.</summary>
public enum LimitationMode
{
    /// <summary>Throw <see cref="SeqeronLimitationException"/> — return only the ideal result. Default.</summary>
    Strict,

    /// <summary>Return the honest best-effort value (the historical behaviour).</summary>
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
                Branch: "3p-arm (miRNA*/star) cleavage boundary",
                Summary: "The opposite-arm (3p / star) boundary is a linear 2-nt-3'-overhang approximation, not the " +
                         "exact miRBase-annotated cut; only the 5p Drosha/Dicer cut reproduces miRBase exactly.",
                RelatedTo: "miRBase mature boundaries encode the dominant sequencing-read cut sites, not a " +
                           "deterministic fold + fixed-overhang rule.",
                Workaround: "Supply the miRBase mature-3p coordinates (MIMAT) or small-RNA-seq read pileups to define " +
                            "the star arm. The 5p mature / Drosha cut is exact and available. Permissive returns the " +
                            "approximate 3p span.",
                ReportPath: "docs/Validation/reports/MIRNA-CLEAVAGE-001.md"),
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
        sb.Append("To receive the best-effort value instead, set ")
          .Append("LimitationPolicy.DefaultMode = LimitationMode.Permissive, or wrap the call in ")
          .Append("`using (LimitationPolicy.UsePermissive()) { ... }`.");
        return sb.ToString();
    }
}
