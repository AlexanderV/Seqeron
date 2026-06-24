using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PROTMOTIF-DOMAIN-001: Protein Domain Identification (FindDomains).
/// Evidence: docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md
/// Spec:     tests/TestSpecs/PROTMOTIF-DOMAIN-001.md
/// Source:   PROSITE PS00028 (zinc finger C2H2), PS00678 (WD-repeats), PS00017 (Walker A);
///           SH3 (PS50002) and PDZ (PS50106) are PROFILE-only — no deterministic pattern exists.
/// </summary>
[TestFixture]
[Category("PROTMOTIF-DOMAIN-001")]
public class ProteinMotifFinder_DomainPrediction_Tests
{
    #region Evidence-sourced constants

    // --- Zinc Finger C2H2 (PROSITE PATTERN PS00028, Pfam PF00096) ---
    // Pattern verbatim: C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
    // Source: PROSITE PS00028; Krishna SS et al. (2003) NAR 31:532–550
    private const string ZincFingerSequence = "AAAACAACAAALEEEEEEEEHAAAHAAAA";
    private const int ZincFingerExpectedStart = 4;
    private const int ZincFingerExpectedEnd = 24;
    private const string ZincFingerExpectedName = "Zinc Finger C2H2";
    private const string ZincFingerExpectedAccession = "PF00096";

    // --- Walker A / P-loop (PROSITE PATTERN PS00017, Pfam PF00069) ---
    // Pattern verbatim: [AG]-x(4)-G-K-[ST]
    // Source: PROSITE PS00017; Walker JE et al. (1982) EMBO J 1:945–951
    private const string KinaseSequence = "AAAAGAAEAGKSAAAA";
    private const int KinaseExpectedStart = 4;
    private const int KinaseExpectedEnd = 11;
    private const string KinaseExpectedName = "Protein Kinase ATP-binding";
    private const string KinaseExpectedAccession = "PF00069";

    // --- WD40 Repeat (EXACT PROSITE PATTERN PS00678 = WD_REPEATS_1, Pfam PF00400) ---
    // Pattern verbatim:
    //   [LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]
    // Translated to regex (PROSITE syntax: x→'.', x(2)→'.{2}', {P}→'[^P]', {DP}→'[^DP]', '-' dropped):
    //   [LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]
    // Source: PROSITE PS00678 (https://prosite.expasy.org/PS00678)
    //
    // Real positive: a WD repeat of GBB1_HUMAN (P62873, β-transducin / Gβ1) embedded in Ala padding.
    // PS00678 is a 14-element signature spanning 15 residues. Segment "LVSASQDGKLIIWDS" is a WD repeat
    // of P62873 and matches PS00678 (verified by ScanProsite-style regex; hand-trace below). Padded so
    // the 15-residue match starts at 0-based index 4 and ends at index 18.
    private const string Wd40Sequence = "AAAALVSASQDGKLIIWDSAAAA";
    private const int Wd40ExpectedStart = 4;
    private const int Wd40ExpectedEnd = 18;
    private const string Wd40ExpectedName = "WD40 Repeat";
    private const string Wd40ExpectedAccession = "PF00400";
    private const string Wd40RealSegment = "LVSASQDGKLIIWDS"; // exact 15-mer matching PS00678

    // Near-miss: the conserved Trp (literal 'W', the 12th element of PS00678) is the most diagnostic
    // residue of the WD signature; replacing W→A destroys the match (PS00678 mandates literal W).
    private const string Wd40NearMissNoTrp = "AAAALVSASQDGKLIIADSAAAA";

    // Full GBB1_HUMAN (UniProt P62873) sequence — a canonical WD40 β-propeller (7 blades).
    // Source: https://rest.uniprot.org/uniprotkb/P62873.fasta
    private const string Gbb1HumanSequence =
        "MSELDQLRQEAEQLKNQIRDARKACADATLSQITNNIDPVGRIQMRTRRTLRGHLAKIYAMHWGTDSRLLVSASQDGKLIIWDSY" +
        "TTNKVHAIPLRSSWVMTCAYAPSGNYVACGGLDNICSIYNLKTREGNVRVSRELAGHTGYLSCCRFLDDNQIVTSSGDTTCALWDI" +
        "ETGQQTTTFTGHTGDVMSLSLAPDTRLFVSGACDASAKLWDVREGMCRQTFTGHESDINAICFFPNGNAFATGSDDATCRLFDLRA" +
        "DQELMTYSHDNIICGITSVSFSKSGRLLLAGYDDFNCNVWDALKADRAQGVLAGHDNRVSCLGVTDDGMAVATGSWDSFLKIWN";

    // --- Multi-domain sequence (zinc finger + kinase) ---
    private const string MultiDomainSequence = "AAAACAACAAALEEEEEEEEHAAAHAAAGAAEAGKSAAAA";

    // SH3 (PROSITE PS50002) and PDZ (PROSITE PS50106) are weight-matrix PROFILES, not patterns:
    // they have NO deterministic signature and are intentionally NOT detected by FindDomains.
    // Signal-peptide prediction is covered by its own Test Unit PROTMOTIF-SP-001.

    #endregion

    #region MUST: FindDomains Tests

    /// <summary>
    /// M1: Zinc finger C2H2 pattern detection.
    /// Evidence: PROSITE PS00028 — C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
    /// </summary>
    [Test]
    public void FindDomains_ZincFinger_MatchesPrositeConsensus()
    {
        var domains = FindDomains(ZincFingerSequence).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Exactly one Zinc Finger C2H2 domain expected");
        Assert.Multiple(() =>
        {
            Assert.That(domains[0].Name, Is.EqualTo(ZincFingerExpectedName),
                "Domain name must match PROSITE PS00028 classification");
            Assert.That(domains[0].Accession, Is.EqualTo(ZincFingerExpectedAccession),
                "Accession must be Pfam PF00096");
            Assert.That(domains[0].Start, Is.EqualTo(ZincFingerExpectedStart),
                "Match start at first Cys of C2H2 motif");
            Assert.That(domains[0].End, Is.EqualTo(ZincFingerExpectedEnd),
                "Match end at second His of C2H2 motif");
        });
    }

    /// <summary>
    /// M2: Walker A / P-loop (kinase ATP-binding) pattern detection.
    /// Evidence: PROSITE PS00017 — [AG]-x(4)-G-K-[ST]; Walker JE et al. (1982) EMBO J
    /// </summary>
    [Test]
    public void FindDomains_WalkerA_MatchesKinaseDomain()
    {
        var domains = FindDomains(KinaseSequence).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Exactly one Kinase ATP-binding domain expected");
        Assert.Multiple(() =>
        {
            Assert.That(domains[0].Name, Is.EqualTo(KinaseExpectedName),
                "Domain name must match P-loop/Walker A classification");
            Assert.That(domains[0].Accession, Is.EqualTo(KinaseExpectedAccession),
                "Accession must be Pfam PF00069");
            Assert.That(domains[0].Start, Is.EqualTo(KinaseExpectedStart),
                "Match start at [AG] position of P-loop");
            Assert.That(domains[0].End, Is.EqualTo(KinaseExpectedEnd),
                "Match end at [ST] position of P-loop");
        });
    }

    /// <summary>
    /// M3: Empty sequence produces no domains.
    /// </summary>
    [Test]
    public void FindDomains_EmptySequence_ReturnsEmpty()
    {
        var domains = FindDomains("").ToList();

        Assert.That(domains, Is.Empty, "Empty input must yield no domains");
    }

    /// <summary>
    /// M4: Null sequence produces no domains.
    /// </summary>
    [Test]
    public void FindDomains_NullSequence_ReturnsEmpty()
    {
        var domains = FindDomains(null!).ToList();

        Assert.That(domains, Is.Empty, "Null input must yield no domains");
    }

    /// <summary>
    /// M5: Every returned domain has non-empty Name, Accession, and Description across all
    /// exact-PROSITE-pattern domain types (zinc finger PS00028, WD40 PS00678, Walker A PS00017).
    /// </summary>
    [Test]
    public void FindDomains_DomainMetadata_HasCorrectFields()
    {
        var testSequences = new[] { ZincFingerSequence, KinaseSequence, Wd40Sequence };

        foreach (var sequence in testSequences)
        {
            var domains = FindDomains(sequence).ToList();
            Assert.That(domains, Is.Not.Empty,
                $"At least one domain must be found in '{sequence[..Math.Min(20, sequence.Length)]}...'");

            foreach (var domain in domains)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(domain.Name, Is.Not.Null.And.Not.Empty,
                        $"Domain at {domain.Start}-{domain.End}: Name must not be null or empty");
                    Assert.That(domain.Accession, Is.Not.Null.And.Not.Empty,
                        $"Domain at {domain.Start}-{domain.End}: Accession must not be null or empty");
                    Assert.That(domain.Description, Is.Not.Null.And.Not.Empty,
                        $"Domain at {domain.Start}-{domain.End}: Description must not be null or empty");
                    Assert.That(domain.Score, Is.GreaterThan(0),
                        $"Domain at {domain.Start}-{domain.End}: Score must be positive");
                });
            }
        }
    }

    /// <summary>
    /// M6: INV-3 — Every domain has Start ≤ End.
    /// </summary>
    [Test]
    public void FindDomains_StartLessOrEqualEnd()
    {
        var testSequences = new[] { ZincFingerSequence, KinaseSequence, Wd40Sequence, Gbb1HumanSequence };

        foreach (var sequence in testSequences)
        {
            var domains = FindDomains(sequence).ToList();
            foreach (var domain in domains)
            {
                Assert.That(domain.Start, Is.LessThanOrEqualTo(domain.End),
                    $"Domain '{domain.Name}': Start ({domain.Start}) must be ≤ End ({domain.End})");
            }
        }
    }

    #endregion

    #region MUST: WD40 exact PROSITE pattern PS00678

    /// <summary>
    /// M7: WD40 repeat detection by EXACT PROSITE PATTERN PS00678 on a real β-transducin segment.
    /// Evidence: PROSITE PS00678; positive = GBB1_HUMAN (P62873) WD repeat "LVSASQDGKLIIWDSY",
    /// padded so the 16-residue match spans 0-based 4..19.
    /// </summary>
    [Test]
    public void FindDomains_WD40Repeat_MatchesPrositePS00678()
    {
        var domains = FindDomains(Wd40Sequence).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Exactly one WD40 Repeat domain expected");
        Assert.Multiple(() =>
        {
            Assert.That(domains[0].Name, Is.EqualTo(Wd40ExpectedName),
                "Domain name must be WD40 Repeat");
            Assert.That(domains[0].Accession, Is.EqualTo(Wd40ExpectedAccession),
                "Accession must be Pfam PF00400");
            Assert.That(domains[0].Start, Is.EqualTo(Wd40ExpectedStart),
                "Match start at first residue of the 15-mer PS00678 signature");
            Assert.That(domains[0].End, Is.EqualTo(Wd40ExpectedEnd),
                "Match end at the last residue ([LIVMFSTAGCN]) of the PS00678 signature");
            Assert.That(domains[0].End - domains[0].Start + 1, Is.EqualTo(Wd40RealSegment.Length),
                "PS00678 is a fixed-length 15-residue signature");
        });
    }

    /// <summary>
    /// M8: WD40 near-miss — replacing the conserved Trp (literal W in PS00678) with Ala
    /// abolishes the match. Confirms the exact pattern is enforced, not a loose heuristic.
    /// Evidence: PS00678 position 12 is the invariant 'W' (the W of "WD"/Trp-Asp).
    /// </summary>
    [Test]
    public void FindDomains_WD40NearMiss_NoConservedTrp_ReturnsNoWD40()
    {
        var domains = FindDomains(Wd40NearMissNoTrp).ToList();

        Assert.That(domains.Any(d => d.Accession == Wd40ExpectedAccession), Is.False,
            "Removing the invariant Trp of WD must abolish the PS00678 match");
        Assert.That(domains, Is.Empty,
            "No other exact-pattern domain should match this near-miss sequence");
    }

    /// <summary>
    /// M9: Real WD40 β-propeller protein GBB1_HUMAN (P62873) yields multiple PS00678 hits.
    /// Evidence: PS00678 is a per-repeat signature; the regex matches the β-transducin repeats
    /// at 0-based starts 69, 156, 284 (verified with the ScanProsite-syntax regex).
    /// </summary>
    [Test]
    public void FindDomains_GBB1Human_DetectsMultipleWD40Repeats()
    {
        var wd40 = FindDomains(Gbb1HumanSequence)
            .Where(d => d.Accession == Wd40ExpectedAccession)
            .OrderBy(d => d.Start)
            .ToList();

        Assert.That(wd40.Select(d => d.Start), Is.EqualTo(new[] { 69, 156, 284 }),
            "GBB1_HUMAN WD repeats must be detected at the PS00678 match positions");
        Assert.That(wd40, Has.All.Matches<ProteinDomain>(d => d.End - d.Start + 1 == 15),
            "Each PS00678 hit is a fixed 15-residue window");
    }

    /// <summary>
    /// M10: Translation verification — the PS00678→regex mapping is reproduced exactly on a
    /// hand-traced 15-mer. Each element of "LVSASQDGKLIIWDS" satisfies its PROSITE element:
    ///  L∈[LIVMSTAC], V∈[LIVMFYWSTAGC], S∈[LIMSTAG], A∈[LIVMSTAGC], (SQ)=x(2), D∈[DN], (G)=x,
    ///  K∈{P}≡[^P], L∈[LIVMWSTAC], I∈{DP}≡[^DP], I∈[LIVMFSTAG], W=W, D∈[DEN], S∈[LIVMFSTAGCN].
    /// </summary>
    [Test]
    public void FindMotifByProsite_PS00678_Translation_MatchesHandTracedSegment()
    {
        const string ps00678 =
            "[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]";

        var matches = FindMotifByProsite(Wd40RealSegment, ps00678, "WD40 Repeat").ToList();

        Assert.That(matches, Has.Count.EqualTo(1),
            "PS00678 translated from PROSITE syntax must match the hand-traced WD repeat exactly once");
        Assert.Multiple(() =>
        {
            Assert.That(matches[0].Start, Is.EqualTo(0), "Match starts at residue 0 of the 15-mer");
            Assert.That(matches[0].End, Is.EqualTo(14), "Match ends at residue 14 (fixed 15-residue signature)");
            Assert.That(matches[0].Sequence, Is.EqualTo(Wd40RealSegment),
                "Captured substring must equal the full PS00678 signature window");
        });
    }

    #endregion

    #region SHOULD: FindDomains Tests

    /// <summary>
    /// S4: Random short sequence with no domain signatures returns empty.
    /// </summary>
    [Test]
    public void FindDomains_NoMatchingDomains_ReturnsEmpty()
    {
        var domains = FindDomains("AAAEEE").ToList();

        Assert.That(domains, Is.Empty, "Short random peptide must yield no domain matches");
    }

    /// <summary>
    /// S5: SH3 / PDZ are PROFILE-only (PS50002 / PS50106) — FindDomains must NOT report them.
    /// Evidence: PROSITE has no deterministic PATTERN for SH3 or PDZ; only weight-matrix profiles.
    /// A canonical SH3 core ("ALYDYEARTEDDLSFKKGERLQI", chicken Src SH3) must not yield an SH3 domain.
    /// </summary>
    [Test]
    public void FindDomains_NoFabricatedSH3OrPDZ_ProfileOnlyDomains()
    {
        // Src SH3 core region; there is no PROSITE pattern for it, so it must not be reported as SH3/PDZ.
        const string sh3Core = "ALYDYEARTEDDLSFKKGERLQIVNNTE";

        var domains = FindDomains(sh3Core).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(domains.Any(d => d.Name == "SH3"), Is.False,
                "SH3 (PROFILE-only PS50002) must not be reported as a fabricated pattern match");
            Assert.That(domains.Any(d => d.Name == "PDZ"), Is.False,
                "PDZ (PROFILE-only PS50106) must not be reported as a fabricated pattern match");
        });
    }

    /// <summary>
    /// S7: INV-10 — FindDomains produces identical results for upper and lower case input.
    /// Evidence: Implementation uses ToUpperInvariant + RegexOptions.IgnoreCase.
    /// </summary>
    [Test]
    public void FindDomains_CaseInsensitive()
    {
        string lower = ZincFingerSequence.ToLowerInvariant();

        var domainsUpper = FindDomains(ZincFingerSequence).ToList();
        var domainsLower = FindDomains(lower).ToList();

        Assert.That(domainsUpper, Has.Count.EqualTo(1), "Upper case must detect zinc finger");
        Assert.That(domainsLower, Has.Count.EqualTo(1), "Lower case must detect zinc finger");
        Assert.Multiple(() =>
        {
            Assert.That(domainsLower[0].Name, Is.EqualTo(domainsUpper[0].Name),
                "Domain name must be identical regardless of case");
            Assert.That(domainsLower[0].Start, Is.EqualTo(domainsUpper[0].Start),
                "Start position must be identical regardless of case");
            Assert.That(domainsLower[0].End, Is.EqualTo(domainsUpper[0].End),
                "End position must be identical regardless of case");
            Assert.That(domainsLower[0].Score, Is.EqualTo(domainsUpper[0].Score).Within(1e-10),
                "Score must be identical regardless of case");
        });
    }

    #endregion

    #region COULD Tests

    /// <summary>
    /// C1: Sequence containing both zinc finger and kinase motifs detects both.
    /// </summary>
    [Test]
    public void FindDomains_MultipleDomainTypes()
    {
        var domains = FindDomains(MultiDomainSequence).ToList();

        Assert.That(domains, Has.Count.EqualTo(2),
            "Multi-domain sequence must yield both Zinc Finger C2H2 and Protein Kinase ATP-binding");
        Assert.Multiple(() =>
        {
            Assert.That(domains.Any(d => d.Name == ZincFingerExpectedName), Is.True,
                "Zinc Finger C2H2 must be detected in multi-domain sequence");
            Assert.That(domains.Any(d => d.Name == KinaseExpectedName), Is.True,
                "Protein Kinase ATP-binding must be detected in multi-domain sequence");
        });
    }

    #endregion
}
