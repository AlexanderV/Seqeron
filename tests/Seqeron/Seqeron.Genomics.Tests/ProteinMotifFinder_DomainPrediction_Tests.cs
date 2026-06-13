using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PROTMOTIF-DOMAIN-001: Domain Prediction &amp; Signal Peptide Prediction.
/// Evidence: docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md
/// Spec:     tests/TestSpecs/PROTMOTIF-DOMAIN-001.md
/// </summary>
[TestFixture]
[Category("PROTMOTIF-DOMAIN-001")]
public class ProteinMotifFinder_DomainPrediction_Tests
{
    #region Evidence-sourced constants

    // --- Zinc Finger C2H2 (PROSITE PS00028, Pfam PF00096) ---
    // Consensus: C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
    // Source: Krishna SS et al. (2003) NAR 31:532–550; PROSITE PS00028
    private const string ZincFingerSequence = "AAAACAACAAALEEEEEEEEHAAAHAAAA";
    private const int ZincFingerExpectedStart = 4;
    private const int ZincFingerExpectedEnd = 24;
    private const string ZincFingerExpectedName = "Zinc Finger C2H2";
    private const string ZincFingerExpectedAccession = "PF00096";

    // --- Walker A / P-loop (PROSITE PS00017, Pfam PF00069) ---
    // Consensus: [AG]-x(4)-G-K-[ST]
    // Source: Walker JE et al. (1982) EMBO J 1:945–951
    private const string KinaseSequence = "AAAAGAAEAGKSAAAA";
    private const int KinaseExpectedStart = 4;
    private const int KinaseExpectedEnd = 11;
    private const string KinaseExpectedName = "Protein Kinase ATP-binding";
    private const string KinaseExpectedAccession = "PF00069";

    // --- WD40 Repeat (Pfam PF00400) ---
    // Simplified pattern: [LIVMFYWC]-x(5,12)-[WF]-D
    private const string Wd40Sequence = "AAAALAAAAAAAAWDAAAA";
    private const int Wd40ExpectedStart = 4;
    private const int Wd40ExpectedEnd = 14;
    private const string Wd40ExpectedName = "WD40 Repeat";
    private const string Wd40ExpectedAccession = "PF00400";

    // --- SH3 Domain (Pfam PF00018) ---
    // Simplified pattern: [LIVMF]-x(2)-[GA]-W-[FYW]-x(5,8)-[LIVMF]
    // Sequence: L(4) AA(5-6) A(7) W(8) F(9) AAAAAA(10-15) L(16) → .{6} between FYW and LIVMF
    private const string Sh3Sequence = "AAAALAAAWFAAAAAALAAAA";
    private const int Sh3ExpectedStart = 4;
    private const int Sh3ExpectedEnd = 16;
    private const string Sh3ExpectedName = "SH3";
    private const string Sh3ExpectedAccession = "PF00018";

    // --- PDZ Domain (Pfam PF00595) ---
    // Simplified pattern: [LIVMF]-[ST]-[LIVMF]-x(2)-G-[LIVMF]-x(3,4)-[LIVMF]-x(2)-[DEN]
    private const string PdzSequence = "AAAALSLAAGLAAAALAADAAAA";
    private const int PdzExpectedStart = 4;
    private const int PdzExpectedEnd = 18;
    private const string PdzExpectedName = "PDZ";
    private const string PdzExpectedAccession = "PF00595";

    // --- Multi-domain sequence (zinc finger + kinase) ---
    private const string MultiDomainSequence = "AAAACAACAAALEEEEEEEEHAAAHAAAGAAEAGKSAAAA";

    // Signal-peptide prediction (PredictSignalPeptide) is covered by its own Test Unit
    // PROTMOTIF-SP-001 — see ProteinMotifFinder_PredictSignalPeptide_Tests.cs.

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
    /// M5: Every returned domain has non-empty Name, Accession, and Description.
    /// Evidence: Pfam domain families require metadata — PF00096, PF00069, PF00400, PF00018, PF00595
    /// Verifies all 5 domain pattern types, not just one.
    /// </summary>
    [Test]
    public void FindDomains_DomainMetadata_HasCorrectFields()
    {
        // INV-2: Verify metadata across all 5 domain pattern types
        var testSequences = new[] { ZincFingerSequence, KinaseSequence, Wd40Sequence, Sh3Sequence, PdzSequence };

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
        // Test across multiple pattern types
        var testSequences = new[] { ZincFingerSequence, KinaseSequence, Wd40Sequence, Sh3Sequence, PdzSequence };

        foreach (var sequence in testSequences)
        {
            var domains = FindDomains(sequence).ToList();
            foreach (var domain in domains)
            {
                Assert.That(domain.Start, Is.LessThanOrEqualTo(domain.End),
                    $"Domain '{domain.Name}' in sequence: Start ({domain.Start}) must be ≤ End ({domain.End})");
            }
        }
    }

    #endregion

    #region SHOULD: FindDomains Tests

    /// <summary>
    /// S1: WD40 repeat pattern detection.
    /// Evidence: Pfam PF00400 — simplified pattern [LIVMFYWC]-x(5,12)-[WF]-D
    /// </summary>
    [Test]
    public void FindDomains_WD40Repeat()
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
                "Match start at [LIVMFYWC] position");
            Assert.That(domains[0].End, Is.EqualTo(Wd40ExpectedEnd),
                "Match end at D position of WD motif");
        });
    }

    /// <summary>
    /// S2: SH3 domain pattern detection.
    /// Evidence: Pfam PF00018 — simplified pattern [LIVMF]-x(2)-[GA]-W-[FYW]-x(5,8)-[LIVMF]
    /// </summary>
    [Test]
    public void FindDomains_SH3Domain()
    {
        var domains = FindDomains(Sh3Sequence).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Exactly one SH3 domain expected");
        Assert.Multiple(() =>
        {
            Assert.That(domains[0].Name, Is.EqualTo(Sh3ExpectedName),
                "Domain name must be SH3");
            Assert.That(domains[0].Accession, Is.EqualTo(Sh3ExpectedAccession),
                "Accession must be Pfam PF00018");
            Assert.That(domains[0].Start, Is.EqualTo(Sh3ExpectedStart),
                "Match start at [LIVMF] anchor");
            Assert.That(domains[0].End, Is.EqualTo(Sh3ExpectedEnd),
                "Match end at terminal [LIVMF]");
        });
    }

    /// <summary>
    /// S3: PDZ domain pattern detection.
    /// Evidence: Pfam PF00595 — simplified pattern [LIVMF]-[ST]-[LIVMF]-x(2)-G-[LIVMF]-x(3,4)-[LIVMF]-x(2)-[DEN]
    /// </summary>
    [Test]
    public void FindDomains_PDZDomain()
    {
        var domains = FindDomains(PdzSequence).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Exactly one PDZ domain expected");
        Assert.Multiple(() =>
        {
            Assert.That(domains[0].Name, Is.EqualTo(PdzExpectedName),
                "Domain name must be PDZ");
            Assert.That(domains[0].Accession, Is.EqualTo(PdzExpectedAccession),
                "Accession must be Pfam PF00595");
            Assert.That(domains[0].Start, Is.EqualTo(PdzExpectedStart),
                "Match start at [LIVMF] anchor");
            Assert.That(domains[0].End, Is.EqualTo(PdzExpectedEnd),
                "Match end at [DEN] position");
        });
    }

    /// <summary>
    /// S4: Random short sequence with no domain signatures returns empty.
    /// </summary>
    [Test]
    public void FindDomains_NoMatchingDomains_ReturnsEmpty()
    {
        // Short peptide with no domain pattern matches
        var domains = FindDomains("AAAEEE").ToList();

        Assert.That(domains, Is.Empty, "Short random peptide must yield no domain matches");
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
