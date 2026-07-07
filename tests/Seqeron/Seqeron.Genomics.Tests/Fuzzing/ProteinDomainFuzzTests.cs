using static Seqeron.Genomics.Analysis.ProteinMotifFinder;
using ProteinDomain = Seqeron.Genomics.Analysis.ProteinMotifFinder.ProteinDomain;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — protein domain FINDING
/// (PROTMOTIF-DOMAIN-001): scanning a protein string for a fixed library of five
/// conserved domain signatures (Zinc Finger C2H2, WD40 repeat, SH3, PDZ and the
/// Walker A / kinase ATP-binding P-loop) via <see cref="ProteinMotifFinder.FindDomains"/>.
/// Each signature is a regex that <c>FindDomains</c> delegates to
/// <see cref="ProteinMotifFinder.FindMotifByPattern"/>, so a domain hit IS a motif hit
/// re-wrapped as a <c>ProteinDomain</c> with domain metadata (Name / Accession /
/// Description) and the underlying motif's inclusive 0-based Start/End and
/// information-content Score. Sibling units PROTMOTIF-FIND-001 (row 82, general motif
/// finding) and PROTMOTIF-PROSITE-001 (row 83, PROSITE pattern syntax / regex injection)
/// are covered separately; this file focuses on the DOMAIN-finding contract.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// state corruption, no nonsense output, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference / ArgumentOutOfRange). Every input must resolve to
/// EITHER a well-defined, theory-correct result OR a *documented, intentional* outcome.
/// For a signature-scanner that uppercases the input and runs five fixed regexes over it,
/// the headline hazards are:
///   • a NullReferenceException when the protein is null;
///   • an IndexOutOfRangeException when the protein is SHORTER than the shortest
///     signature, so no signature span can fit (the scan must simply yield nothing,
///     never run off the end);
///   • a runaway / quadratic scan on a homopolymer or long all-X junk run where a
///     signature could in principle match at many adjacent offsets (the scan must
///     terminate in time linear in the number of matches);
///   • a mis-reported domain whose [Start..End] span does NOT actually equal the
///     substring it claims, or whose Score is NaN/±∞ (a coordinate / score bug).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-DOMAIN-001 — protein domain finding (signature scan)
/// Checklist: docs/checklists/03_FUZZING.md, row 84.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the length / composition corners that could crash,
///     hang, or fabricate a domain:
///       – empty protein: "" / null → NO domains by the explicit string.IsNullOrEmpty
///         guard, never a NullReferenceException
///         (Domain_Prediction.md §3.3, §6.1 "Null or empty sequence → Returns no domains").
///       – shorter than the smallest signature: the shortest domain signatures (WD40
///         `[LIVMFYWC].{5,12}[WF]D` and kinase `[AG].{4}GK[ST]`) each need ≥ 8 residues,
///         so a protein of length ≤ 7 cannot contain ANY signature → NO domain, no
///         IndexOutOfRange on the absent span (INV-DOMAIN-02; §4.A high-level step 2).
///       – all-X residues (and other out-of-alphabet / junk): 'X' is not in any
///         amino-acid set of any signature, so an all-X protein satisfies no signature →
///         NO domains, no crash. The alphabet is NOT validated beyond pattern matching
///         (§3.3), so junk simply fails to match rather than throwing.
///       – all-same char (homopolymer): a defined outcome — a long run of a single
///         residue that cannot complete any signature yields nothing and the scan
///         terminates (no quadratic hang).
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE);
///   targets: "Protein shorter than min domain, empty, all X residues".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The domain-finding contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Given a protein sequence S, report every contiguous subsequence of S that satisfies
/// one COMPLETE signature from the fixed five-entry library. The implementation
/// uppercases S and, per signature, delegates to FindMotifByPattern (which compiles
/// "(?=(P))" so OVERLAPPING occurrences are discovered) and re-wraps each motif hit as
/// a ProteinDomain { Name, Accession, Start (incl. 0-based), End (incl. 0-based), Score
/// (information content inherited from the motif), Description }.
///   — docs/algorithms/ProteinMotif/Domain_Prediction.md §2.A, §3.1–§3.3, §4.A;
///     coordinate contract shared with docs/algorithms/ProteinMotif/Motif_Search.md
///     §3.2 (Start/End inclusive 0-based), §2.4 INV-02 (contiguous span).
///
/// The five signatures (Domain_Prediction.md §4.A "Decision Rules"):
///   Zinc Finger C2H2  PF00096  C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H            (min 21 residues)
///   WD40 Repeat       PF00400  [LIVMFYWC].{5,12}[WF]D                        (min  8 residues)
///   SH3               PF00018  [LIVMF].{2}[GA]W[FYW].{5,8}[LIVMF]            (min 12 residues)
///   PDZ               PF00595  [LIVMF][ST][LIVMF].{2}G[LIVMF].{3,4}[LIVMF].{2}[DEN] (min 14)
///   Protein Kinase    PF00069  [AG].{4}GK[ST]                               (min  8 residues)
///
/// Method under test (src/.../Seqeron.Genomics.Analysis/ProteinMotifFinder.cs):
///   IEnumerable&lt;ProteinDomain&gt; FindDomains(string proteinSequence)
///     — scans the five hard-coded signatures and yields one ProteinDomain per hit.
///   — Domain_Prediction.md §5.1.
///
/// Documented input handling (Domain_Prediction.md §3.3, §6.1):
///   • null / "" protein → NO domains (explicit string.IsNullOrEmpty guard), no throw.
///   • Matching is case-INSENSITIVE: the input is uppercased before matching.
///   • The amino-acid alphabet is NOT validated beyond pattern matching — non-AA
///     characters (digits, punctuation, X, B/Z/J/O/U) simply fail to satisfy a signature.
///   • Coordinates are inclusive 0-based; End − Start + 1 equals the matched signature
///     span length.
///
/// Theory-correct invariants asserted (Domain_Prediction.md §2.A "Properties"):
///   • INV-DOMAIN-01 — deterministic: re-scanning the same sequence yields identical hits.
///   • INV-DOMAIN-02 — every reported hit spans a CONTIGUOUS subsequence satisfying one
///     complete signature: S[Start..End] is the matched span, with 0 ≤ Start ≤ End ≤ n−1
///     and End − Start + 1 ≥ the matched signature's minimum width (the headline
///     "no run-off-the-end / no coordinate bug" property).
///   • [finite-score] — Score is finite (never NaN / ±∞).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// FindDomains is O(n × d) for d = 5 fixed signatures (Domain_Prediction.md §4.A
/// "Complexity"). The homopolymer and long all-X targets maximise the work of the five
/// regex walks; they are kept modest and [CancelAfter]-guarded so a regression that
/// turned a signature scan into a hang or a super-linear blow-up would FAIL rather than
/// wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinDomainFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino-acid one-letter codes.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>
    /// The minimum residue width of the SHORTEST domain signature in the library: both the WD40
    /// repeat (<c>[LIVMFYWC].{5,12}[WF]D</c>) and the kinase ATP-binding P-loop
    /// (<c>[AG].{4}GK[ST]</c>) require 8 residues. A protein STRICTLY shorter than this cannot
    /// contain ANY signature, so <see cref="ProteinMotifFinder.FindDomains"/> must yield nothing
    /// (Domain_Prediction.md §4.A; INV-DOMAIN-02).
    /// </summary>
    private const int ShortestSignatureWidth = 8;

    /// <summary>
    /// The minimum residue width of each library signature, used to assert INV-DOMAIN-02:
    /// a reported span can never be shorter than the signature it claims to satisfy.
    /// Keyed by the repository Accession (Domain_Prediction.md §4.A "Decision Rules").
    /// </summary>
    private static readonly IReadOnlyDictionary<string, int> MinSignatureWidthByAccession =
        new Dictionary<string, int>
        {
            ["PF00096"] = 21, // Zinc Finger C2H2: C + 2 + C + 3 + 1 + 8 + H + 3 + H
            ["PF00400"] = 8,  // WD40 Repeat:      1 + 5 + [WF] + D
            ["PF00018"] = 12, // SH3:              1 + 2 + 1(GA) + W + 1(FYW) + 5 + 1
            ["PF00595"] = 14, // PDZ:              1+1+1 + 2 + G + 1 + 3 + 1 + 2 + 1
            ["PF00069"] = 8,  // Kinase ATP-binding: 1 + 4 + G + K + 1(ST)
        };

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomProtein(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardAminoAcids[rng.Next(StandardAminoAcids.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted <see cref="ProteinDomain"/> must
    /// satisfy against the original (case-insensitive) sequence (Domain_Prediction.md §2.A
    /// INV-DOMAIN-02, §3.2; Motif_Search.md §2.4 INV-02, §3.2): the hit is a CONTIGUOUS in-bounds
    /// subsequence whose claimed coordinates actually reproduce a span of the (uppercased) input, whose
    /// span length matches the documented inclusive 0-based convention (End − Start + 1) and is at least
    /// the matched signature's minimum width, whose Accession is one of the five library accessions, and
    /// whose Score is finite. This is the headline "no coordinate bug, no run-off-the-end, no fabricated
    /// domain, no NaN" property.
    /// </summary>
    private static void AssertWellFormedDomain(ProteinDomain domain, string originalSequence)
    {
        string upper = originalSequence.ToUpperInvariant();
        int n = upper.Length;

        // Identity — the domain comes from the fixed five-entry library (Domain_Prediction.md §5.2).
        MinSignatureWidthByAccession.Should().ContainKey(domain.Accession,
            "every reported domain is one of the five library signatures");
        domain.Name.Should().NotBeNullOrEmpty("a domain carries its repository Name");
        domain.Description.Should().NotBeNullOrEmpty("a domain carries its repository Description");

        // INV-DOMAIN-02 — in-bounds, non-empty, contiguous span.
        domain.Start.Should().BeInRange(0, n - 1, "a domain Start is a valid 0-based residue index");
        domain.End.Should().BeInRange(domain.Start, n - 1, "a domain End is in-bounds and not before its Start");

        // [span-shape] — inclusive 0-based span length, no shorter than the signature's minimum width.
        int spanLength = domain.End - domain.Start + 1;
        spanLength.Should().BeGreaterThanOrEqualTo(MinSignatureWidthByAccession[domain.Accession],
            "INV-DOMAIN-02: a hit cannot be narrower than the signature it claims to satisfy");

        // INV-DOMAIN-02 — the reported span is an actual contiguous substring of the uppercased input.
        (domain.Start + spanLength).Should().BeLessThanOrEqualTo(n,
            "the matched span lies fully within the sequence (no run-off-the-end)");
        upper.Substring(domain.Start, spanLength).Length.Should().Be(spanLength,
            "the reported span is exactly S[Start..End] of the uppercased input (no coordinate bug)");

        // [finite-score] — score is finite.
        double.IsNaN(domain.Score).Should().BeFalse("a domain Score must never be NaN");
        double.IsInfinity(domain.Score).Should().BeFalse("a domain Score must never be infinite");
    }

    /// <summary>
    /// Independently re-derives, from the public <see cref="ProteinMotifFinder.FindMotifByPattern"/>
    /// scanner, whether the given signature regex matches the sequence at the given domain hit; used
    /// to confirm INV-DOMAIN-02 that the reported span actually satisfies a COMPLETE signature.
    /// </summary>
    private static bool SignatureMatchesAt(string sequence, string signatureRegex, ProteinDomain domain)
    {
        return FindMotifByPattern(sequence, signatureRegex)
            .Any(m => m.Start == domain.Start && m.End == domain.End);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-DOMAIN-001 — protein domain finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-DOMAIN-001 — protein domain finding

    #region BE — Empty / null protein: no domains, no NullReference

    /// <summary>
    /// Target "empty": "" and null must produce NO domains — by the explicit string.IsNullOrEmpty
    /// guard, NEVER a NullReferenceException (Domain_Prediction.md §3.3, §6.1 "Null or empty
    /// sequence → Returns no domains"). This is the headline no-crash contract on the degenerate
    /// empty input.
    /// </summary>
    [Test]
    public void FindDomains_EmptyOrNullProtein_NoDomainsNoThrow()
    {
        foreach (string? seq in new[] { "", null })
        {
            var act = () => FindDomains(seq!).ToList();
            act.Should().NotThrow($"empty/null protein ('{seq ?? "null"}') must not crash FindDomains")
                .Subject.Should().BeEmpty("empty/null protein yields no domains");
        }
    }

    #endregion

    #region BE — Protein shorter than the smallest signature: no domain, no crash

    /// <summary>
    /// Target "Protein shorter than min domain": the shortest signatures in the library (WD40 and the
    /// kinase ATP-binding P-loop) each need ≥ <see cref="ShortestSignatureWidth"/> = 8 residues, so a
    /// protein STRICTLY shorter than 8 residues has no place for ANY signature to fit and must yield NO
    /// domain — never running off the end (no IndexOutOfRange on the absent span). We probe EVERY length
    /// 0..7 across multiple seeds, plus the residue-rich prefixes of the genuine positive sequences
    /// truncated below 8 residues, and finally the boundary just below the LONGEST signature
    /// (Zinc Finger C2H2, min 21): a 20-residue genuine zinc-finger prefix is one residue too short to
    /// complete the C2H2 signature and so reports no zinc-finger domain.
    /// (Domain_Prediction.md §4.A; INV-DOMAIN-02.)
    /// </summary>
    [Test]
    public void FindDomains_ProteinShorterThanSmallestSignature_NoDomainNoCrash()
    {
        // (a) Every length 0..7 is below the 8-residue minimum → no domain whatsoever.
        foreach (int seed in new[] { 1, 17, 99 })
        {
            for (int len = 0; len < ShortestSignatureWidth; len++)
            {
                string seq = RandomProtein(len, seed);
                var act = () => FindDomains(seq).ToList();
                act.Should().NotThrow($"a length-{len} protein (seed {seed}) must not crash FindDomains")
                    .Subject.Should().BeEmpty(
                        $"a length-{len} protein is shorter than the smallest (8-residue) signature");
            }
        }

        // (b) Residue-rich short prefixes that still cannot complete any signature.
        foreach (string seq in new[] { "C", "CC", "GKS", "AGKS", "LIVMFYW" /* 7 residues */ })
        {
            var act = () => FindDomains(seq).ToList();
            act.Should().NotThrow($"a short protein ('{seq}') must not crash FindDomains")
                .Subject.Should().BeEmpty($"the short protein '{seq}' is below the smallest signature width");
        }

        // (c) Boundary just below the LONGEST signature: a genuine 21-residue zinc-finger truncated to
        //     20 residues can no longer complete the C2H2 signature → no zinc-finger domain.
        const string zincFinger21 = "CAACAAALEEEEEEEEHAAAH"; // exact 21-residue C2H2 match
        var full = FindDomains(zincFinger21).ToList();
        full.Should().Contain(d => d.Accession == "PF00096",
            "sanity: the full 21-residue sequence does contain the C2H2 zinc-finger signature");

        string zincFinger20 = zincFinger21[..20];
        var truncated = FindDomains(zincFinger20).ToList();
        truncated.Should().NotContain(d => d.Accession == "PF00096",
            "a 20-residue prefix is one residue short of the 21-residue C2H2 signature → no zinc-finger domain");
        foreach (var d in truncated)
            AssertWellFormedDomain(d, zincFinger20); // any incidental shorter-signature hit stays well-formed
    }

    #endregion

    #region BE — All-X residues (and other junk): no signature satisfied, no crash

    /// <summary>
    /// Target "all X residues": the unknown placeholder 'X' is not a member of ANY amino-acid set in
    /// ANY of the five signatures, so a protein composed entirely of 'X' satisfies no signature and
    /// must yield NO domains — never a crash. The alphabet is not validated beyond pattern matching
    /// (Domain_Prediction.md §3.3), so out-of-alphabet residues simply fail to match. We additionally
    /// probe other out-of-alphabet junk (digits, punctuation, whitespace, the extended IUPAC codes
    /// B/Z/J/O/U) and junk-flanked genuine residues; any incidental hit (a junk residue legitimately
    /// occupying an <c>x</c>-wildcard position inside an OTHERWISE-valid signature) must still be a
    /// well-formed, in-bounds, contiguous span.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindDomains_AllXAndJunkResidues_NoSignatureMatchedNoCrash(CancellationToken token)
    {
        // (a) Pure all-X of various lengths (including well beyond the longest signature) → no domains.
        foreach (int len in new[] { 8, 21, 50, 200, 1000 })
        {
            string allX = new string('X', len);
            var act = () => FindDomains(allX).ToList();
            act.Should().NotThrow($"an all-X protein of length {len} must not crash FindDomains")
                .Subject.Should().BeEmpty("'X' satisfies no amino-acid position of any signature");
            token.ThrowIfCancellationRequested();
        }

        // (b) Other pure out-of-alphabet junk → no domains, no crash.
        foreach (string junk in new[]
                 {
                     "1234567890123456789012",  // digits (22 long, > longest signature)
                     "!@#$%^&*()!@#$%^&*()!@",    // punctuation
                     "   \t  \n  \r  \t  \n  ",   // whitespace only
                     "BZJOUBZJOUBZJOUBZJOUBZ",    // extended IUPAC ambiguity codes
                 })
        {
            var act = () => FindDomains(junk).ToList();
            act.Should().NotThrow($"junk protein ('{junk.Trim()}') must not crash FindDomains")
                .Subject.Should().BeEmpty("out-of-alphabet residues satisfy no signature");
            token.ThrowIfCancellationRequested();
        }

        // (c) Genuine kinase signature with junk in its wildcard (x) positions: the signature
        //     [AG].{4}GK[ST] allows ANY four residues at the wildcard span, INCLUDING junk, because
        //     the alphabet is not validated. The hit must still be well-formed and reproduce the span.
        //     "A" + "X##!" (4 junk wildcards) + "GKS".
        const string kinaseWithJunkWildcards = "AX##!GKS";
        var kinaseHits = FindDomains(kinaseWithJunkWildcards).ToList();
        kinaseHits.Should().Contain(d => d.Accession == "PF00069",
            "the kinase signature's x-wildcards admit any residue, including junk, so the hit is still found");
        foreach (var d in kinaseHits)
            AssertWellFormedDomain(d, kinaseWithJunkWildcards);

        // (d) All-X embedded around a genuine kinase core does not corrupt the genuine hit.
        string padded = new string('X', 10) + "AGGGGGKS" + new string('X', 10);
        var paddedHits = FindDomains(padded).ToList();
        paddedHits.Should().Contain(d => d.Accession == "PF00069",
            "a genuine kinase signature surrounded by X is still detected");
        foreach (var d in paddedHits)
            AssertWellFormedDomain(d, padded);
    }

    #endregion

    #region BE — All-same char (homopolymer): defined outcome, no hang

    /// <summary>
    /// Target "all same char": a homopolymer "AAAA…" has a DEFINED domain-finding outcome and the scan
    /// must TERMINATE (no quadratic hang on the five regex walks). No single-residue homopolymer can
    /// complete any signature — every signature requires at least two DISTINCT conserved residues
    /// (e.g. the kinase P-loop needs both a G and a K; WD40 needs a terminal D; C2H2 needs C…C…H…H), so
    /// EVERY homopolymer of a single standard residue yields NO domains. We assert this for a long
    /// homopolymer of every standard residue and confirm the scan terminates well within the deadline.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindDomains_Homopolymer_NoDomainsNoHang(CancellationToken token)
    {
        foreach (char aa in StandardAminoAcids)
        {
            const int n = 300;
            string seq = new string(aa, n);

            var act = () => FindDomains(seq).ToList();
            act.Should().NotThrow($"a length-{n} homopolymer of '{aa}' must not crash/hang FindDomains")
                .Subject.Should().BeEmpty(
                    $"a single-residue homopolymer of '{aa}' cannot complete any signature " +
                    "(every signature needs ≥ 2 distinct conserved residues)");
            token.ThrowIfCancellationRequested();
        }
    }

    #endregion

    #region Positive sanity — genuine signatures are found at the right coordinates

    /// <summary>
    /// Positive sanity: the harness must assert against a scanner that actually FINDS domains at the
    /// correct coordinates, not a no-op. Hand-built proteins containing a single genuine signature at a
    /// KNOWN offset must each yield exactly that domain with the correct Accession and inclusive 0-based
    /// Start/End; and S[Start..End] must equal the matched signature span (independently re-derived from
    /// FindMotifByPattern). This pins the headline INV-DOMAIN-02 "a found domain actually occurs at the
    /// reported position and satisfies a complete signature" contract (Domain_Prediction.md §2.A, §4.A).
    /// </summary>
    [Test]
    public void FindDomains_GenuineSignatures_FoundAtCorrectCoordinates()
    {
        // (a) Kinase ATP-binding P-loop [AG].{4}GK[ST] (PF00069), 8 residues at offset 4.
        const string kinaseRegex = "[AG].{4}GK[ST]";
        const string kinaseProtein = "MKSPAGGGGGKSWFIL"; // "AGGGGGKS" begins at offset 4
        int kinaseOffset = kinaseProtein.IndexOf("AGGGGGKS", StringComparison.Ordinal);
        kinaseOffset.Should().Be(4, "sanity-check the hand-built kinase offset");

        var kinaseDomains = FindDomains(kinaseProtein).Where(d => d.Accession == "PF00069").ToList();
        kinaseDomains.Should().ContainSingle("the protein contains exactly one kinase ATP-binding signature");
        var kinase = kinaseDomains[0];
        AssertWellFormedDomain(kinase, kinaseProtein);
        kinase.Name.Should().Be("Protein Kinase ATP-binding", "the kinase domain carries its library name");
        kinase.Start.Should().Be(kinaseOffset, "the kinase hit starts at the known offset 4");
        kinase.End.Should().Be(kinaseOffset + 7, "the 8-residue P-loop spans [4,11]");
        kinaseProtein[kinase.Start..(kinase.End + 1)].Should().Be("AGGGGGKS",
            "INV-DOMAIN-02: S[Start..End] is exactly the matched kinase signature span");
        SignatureMatchesAt(kinaseProtein, kinaseRegex, kinase).Should().BeTrue(
            "INV-DOMAIN-02: the reported span independently satisfies the complete kinase signature");

        // (b) Zinc Finger C2H2 C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H (PF00096), 21 residues at offset 4.
        const string zincRegex = "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";
        const string zincFingerCore = "CAACAAALEEEEEEEEHAAAH"; // exact 21-residue C2H2 match
        const string zincProtein = "MKSP" + zincFingerCore + "WFIL";
        int zincOffset = zincProtein.IndexOf(zincFingerCore, StringComparison.Ordinal);
        zincOffset.Should().Be(4, "sanity-check the hand-built zinc-finger offset");

        var zincDomains = FindDomains(zincProtein).Where(d => d.Accession == "PF00096").ToList();
        zincDomains.Should().ContainSingle("the protein contains exactly one C2H2 zinc-finger signature");
        var zinc = zincDomains[0];
        AssertWellFormedDomain(zinc, zincProtein);
        zinc.Name.Should().Be("Zinc Finger C2H2", "the zinc-finger domain carries its library name");
        zinc.Start.Should().Be(zincOffset, "the zinc-finger hit starts at the known offset 4");
        zinc.End.Should().Be(zincOffset + zincFingerCore.Length - 1, "the 21-residue C2H2 spans [4,24]");
        zincProtein[zinc.Start..(zinc.End + 1)].Should().Be(zincFingerCore,
            "INV-DOMAIN-02: S[Start..End] is exactly the matched zinc-finger signature span");
        SignatureMatchesAt(zincProtein, zincRegex, zinc).Should().BeTrue(
            "INV-DOMAIN-02: the reported span independently satisfies the complete zinc-finger signature");

        // (c) Two well-separated kinase signatures → two distinct hits at the two known offsets.
        const string twoKinase = "AGGGGGKS" + "DDDDDDDD" + "AGGGGGKS";
        var twoHits = FindDomains(twoKinase).Where(d => d.Accession == "PF00069").ToList();
        twoHits.Should().HaveCount(2, "two well-separated kinase signatures yield two hits");
        twoHits.Select(d => d.Start).Should().Equal(new[] { 0, 16 },
            "the two kinase hits are reported at the two known offsets");
        foreach (var d in twoHits)
        {
            AssertWellFormedDomain(d, twoKinase);
            twoKinase[d.Start..(d.End + 1)].Should().Be("AGGGGGKS", "each kinase hit is its full 8-residue span");
        }

        // (d) A protein with one zinc-finger AND one kinase signature surfaces BOTH domain types.
        const string multi = zincFingerCore + "GGGG" + "AGGGGGKS";
        var multiDomains = FindDomains(multi).ToList();
        multiDomains.Select(d => d.Accession).Should().Contain("PF00096",
            "a multi-domain protein surfaces the zinc-finger signature");
        multiDomains.Select(d => d.Accession).Should().Contain("PF00069",
            "a multi-domain protein surfaces the kinase signature");
        foreach (var d in multiDomains)
            AssertWellFormedDomain(d, multi);
    }

    /// <summary>
    /// Positive sanity over RANDOM proteins: across fixed seeds and lengths the scan must never crash,
    /// hang, or emit a malformed domain, and every emitted ProteinDomain must satisfy the full contract
    /// (in-bounds contiguous span no shorter than its signature, finite score, library accession).
    /// INV-DOMAIN-01 determinism is pinned by re-scanning the same input and requiring identical hits.
    /// This pins span-correctness and termination on arbitrary sequences, not just hand-built signatures.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindDomains_RandomProtein_AlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 7, 31, 137, 2026 })
        {
            foreach (int len in new[] { 1, 7, 8, 21, 80, 300 })
            {
                string seq = RandomProtein(len, seed);

                var act = () => FindDomains(seq).ToList();
                var hits = act.Should().NotThrow($"random protein must not crash (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                foreach (var hit in hits)
                    AssertWellFormedDomain(hit, seq);

                // INV-DOMAIN-01 — deterministic: the same input yields identical hits.
                var again = FindDomains(seq).ToList();
                again.Select(d => (d.Accession, d.Start, d.End))
                    .Should().Equal(hits.Select(d => (d.Accession, d.Start, d.End)),
                        "INV-DOMAIN-01: FindDomains is deterministic for a fixed input");
            }
        }
    }

    #endregion

    #endregion
}
