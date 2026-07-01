// MOTIF-REGULATORY-001 — Regulatory Element Scan (known regulatory consensus motifs)
// Evidence: docs/Evidence/MOTIF-REGULATORY-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-REGULATORY-001.md
// Source: Bucher (1990) J Mol Biol 212:563; Harley & Reynolds (1987) Nucleic Acids Res 15:2343;
//         Lundin/Nehlin/Ronne (1994) Mol Cell Biol 14:1979; Kozak (1987) Nucleic Acids Res 15:8125;
//         Proudfoot & Brownlee (1976) Nature 263:211; Massari & Murre (2000) Mol Cell Biol 20:429;
//         Lee/Mitchell/Tjian (1987) Cell 49:741; Sen & Baltimore (1986) Cell 46:705;
//         Montminy et al. (1986) PNAS 83:6682.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Canonical test class for MOTIF-REGULATORY-001. Verifies
/// <see cref="MotifFinder.FindRegulatoryElements(DnaSequence)"/> against the published consensus
/// sequences of 12 regulatory elements. Each probe embeds one consensus at a known 0-based offset;
/// expected positions/patterns are derived from the source consensus strings, not from the code.
/// </summary>
[TestFixture]
public class MotifFinder_FindRegulatoryElements_Tests
{
    private static RegulatoryElement[] Scan(string sequence) =>
        MotifFinder.FindRegulatoryElements(new DnaSequence(sequence)).ToArray();

    private static RegulatoryElement Single(string sequence, string name)
    {
        var hits = Scan(sequence).Where(e => e.Name == name).ToArray();
        Assert.That(hits.Length, Is.EqualTo(1),
            $"Expected exactly one '{name}' hit in '{sequence}' but found {hits.Length}.");
        return hits[0];
    }

    #region FindRegulatoryElements — MUST (each element detected at exact position)

    // M1 — TATA box consensus TATAAA (Bucher 1990). Probe "GGGTATAAAGGG" -> pos 3.
    [Test]
    public void FindRegulatoryElements_TataBoxProbe_DetectedAtConsensusPosition()
    {
        var e = Single("GGGTATAAAGGG", "TATA Box");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(3), "TATAAA begins at 0-based index 3 in GGGTATAAAGGG.");
            Assert.That(e.Pattern, Is.EqualTo("TATAAA"), "TATA box consensus per Bucher (1990).");
            Assert.That(e.Sequence, Is.EqualTo("TATAAA"), "Matched substring equals the consensus.");
        });
    }

    // M2 — -10 Pribnow box consensus TATAAT (Harley & Reynolds 1987). Probe "CCTATAATCC" -> pos 2.
    [Test]
    public void FindRegulatoryElements_MinusTenBoxProbe_DetectedAtConsensusPosition()
    {
        var e = Single("CCTATAATCC", "-10 Box");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "TATAAT begins at index 2 in CCTATAATCC.");
            Assert.That(e.Pattern, Is.EqualTo("TATAAT"), "-10 hexamer per Harley & Reynolds (1987).");
        });
    }

    // M3 — -35 box consensus TTGACA (Harley & Reynolds 1987). Probe "AATTGACAGG" -> pos 2.
    [Test]
    public void FindRegulatoryElements_MinusThirtyFiveBoxProbe_DetectedAtConsensusPosition()
    {
        var e = Single("AATTGACAGG", "-35 Box");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "TTGACA begins at index 2 in AATTGACAGG.");
            Assert.That(e.Pattern, Is.EqualTo("TTGACA"), "-35 hexamer per Harley & Reynolds (1987).");
        });
    }

    // M4 — CCAAT box pentanucleotide CCAAT (Bucher 1990). Probe "GGCCAATGG" -> pos 2.
    [Test]
    public void FindRegulatoryElements_CaatBoxProbe_DetectedAtConsensusPosition()
    {
        var e = Single("GGCCAATGG", "CAAT Box");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "CCAAT begins at index 2 in GGCCAATGG.");
            Assert.That(e.Pattern, Is.EqualTo("CCAAT"), "CCAAT pentanucleotide per Bucher (1990).");
        });
    }

    // M5 — GC box consensus GGGCGG (Lundin et al. 1994). Probe "AAGGGCGGTT" -> pos 2.
    [Test]
    public void FindRegulatoryElements_GcBoxProbe_DetectedAtConsensusPosition()
    {
        var e = Single("AAGGGCGGTT", "GC Box");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "GGGCGG begins at index 2 in AAGGGCGGTT.");
            Assert.That(e.Pattern, Is.EqualTo("GGGCGG"), "GC box consensus per Lundin et al. (1994).");
        });
    }

    // M6 — Kozak optimal context GCCGCCACCATGG (Kozak 1987). Probe "TTGCCGCCACCATGGAA" -> pos 2.
    [Test]
    public void FindRegulatoryElements_KozakProbe_DetectedAtConsensusPosition()
    {
        var e = Single("TTGCCGCCACCATGGAA", "Kozak");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "GCCGCCACCATGG begins at index 2.");
            Assert.That(e.Pattern, Is.EqualTo("GCCGCCACCATGG"), "Kozak optimal context per Kozak (1987).");
        });
    }

    // M7 — Shine-Dalgarno consensus AGGAGG. Probe "TTAGGAGGTTT" -> pos 2.
    [Test]
    public void FindRegulatoryElements_ShineDalgarnoProbe_DetectedAtConsensusPosition()
    {
        var e = Single("TTAGGAGGTTT", "Shine-Dalgarno");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "AGGAGG begins at index 2 in TTAGGAGGTTT.");
            Assert.That(e.Pattern, Is.EqualTo("AGGAGG"), "Shine-Dalgarno consensus (3' 16S rRNA complement).");
        });
    }

    // M8 — Poly(A) signal hexamer AATAAA (Proudfoot & Brownlee 1976). Probe "CCAATAAACC" -> pos 2.
    [Test]
    public void FindRegulatoryElements_PolyASignalProbe_DetectedAtConsensusPosition()
    {
        var e = Single("CCAATAAACC", "Poly(A) Signal");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "AATAAA begins at index 2 in CCAATAAACC.");
            Assert.That(e.Pattern, Is.EqualTo("AATAAA"), "Poly(A) signal per Proudfoot & Brownlee (1976).");
        });
    }

    // M9 — E-box CANNTG matches a degenerate occurrence CACGTG (Massari & Murre 2000). Probe "GGCACGTGGG" -> pos 2.
    [Test]
    public void FindRegulatoryElements_EboxDegenerateProbe_MatchesViaIupac()
    {
        var e = Single("GGCACGTGGG", "E-box");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "CACGTG begins at index 2 in GGCACGTGGG.");
            Assert.That(e.Pattern, Is.EqualTo("CANNTG"), "E-box IUPAC consensus per Massari & Murre (2000).");
            Assert.That(e.Sequence, Is.EqualTo("CACGTG"), "Matched substring is the canonical palindromic E-box.");
        });
    }

    // M10 — AP-1 corrected consensus TGACTCA (Lee/Mitchell/Tjian 1987). Probe "AATGACTCAGG" -> pos 2.
    [Test]
    public void FindRegulatoryElements_Ap1CorrectConsensus_DetectedAtConsensusPosition()
    {
        var e = Single("AATGACTCAGG", "AP-1");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "TGACTCA begins at index 2 in AATGACTCAGG.");
            Assert.That(e.Pattern, Is.EqualTo("TGACTCA"), "AP-1 recognition motif per Lee/Mitchell/Tjian (1987).");
        });
    }

    // M11 — Regression: the old WRONG pattern TGAGTCA must NOT be reported as AP-1.
    [Test]
    public void FindRegulatoryElements_Ap1OldWrongPattern_NotReported()
    {
        // "AATGAGTCAGG" contains TGAGTCA (the prior buggy constant), which is NOT the AP-1 consensus.
        var ap1 = Scan("AATGAGTCAGG").Where(e => e.Name == "AP-1").ToArray();
        Assert.That(ap1.Length, Is.EqualTo(0),
            "TGAGTCA is not the AP-1 consensus (TGACTCA per Lee/Mitchell/Tjian 1987); it must not match.");
    }

    // M12 — NF-κB reference κB site GGGACTTTCC (Sen & Baltimore 1986). Probe "AAGGGACTTTCCAA" -> pos 2.
    [Test]
    public void FindRegulatoryElements_NfKbProbe_DetectedAtConsensusPosition()
    {
        var e = Single("AAGGGACTTTCCAA", "NF-κB");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "GGGACTTTCC begins at index 2 in AAGGGACTTTCCAA.");
            Assert.That(e.Pattern, Is.EqualTo("GGGACTTTCC"), "NF-κB reference κB site per Sen & Baltimore (1986).");
        });
    }

    // M13 — CREB CRE palindrome TGACGTCA (Montminy et al. 1986). Probe "CCTGACGTCAGG" -> pos 2.
    [Test]
    public void FindRegulatoryElements_CrebProbe_DetectedAtConsensusPosition()
    {
        var e = Single("CCTGACGTCAGG", "CREB");
        Assert.Multiple(() =>
        {
            Assert.That(e.Position, Is.EqualTo(2), "TGACGTCA begins at index 2 in CCTGACGTCAGG.");
            Assert.That(e.Pattern, Is.EqualTo("TGACGTCA"), "CREB CRE palindrome per Montminy et al. (1986).");
        });
    }

    // M14 — Null sequence -> ArgumentNullException.
    [Test]
    public void FindRegulatoryElements_NullSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => MotifFinder.FindRegulatoryElements(null!).ToList(),
            "Null sequence must throw ArgumentNullException per the contract.");
    }

    // M15 — Empty sequence -> no elements.
    [Test]
    public void FindRegulatoryElements_EmptySequence_ReturnsNoElements()
    {
        var hits = Scan("");
        Assert.That(hits, Is.Empty, "An empty sequence has no scan window, so no elements are reported.");
    }

    // M16 — All consensus constants equal their cited source strings (locks INV-04, no fabricated values).
    [Test]
    public void KnownMotifs_Constants_EqualSourceConsensusStrings()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.KnownMotifs.TataBox, Is.EqualTo("TATAAA"), "Bucher (1990).");
            Assert.That(MotifFinder.KnownMotifs.CaatBox, Is.EqualTo("CCAAT"), "Bucher (1990).");
            Assert.That(MotifFinder.KnownMotifs.GcBox, Is.EqualTo("GGGCGG"), "Lundin et al. (1994).");
            Assert.That(MotifFinder.KnownMotifs.MinusTenBox, Is.EqualTo("TATAAT"), "Harley & Reynolds (1987).");
            Assert.That(MotifFinder.KnownMotifs.MinusThirtyFiveBox, Is.EqualTo("TTGACA"), "Harley & Reynolds (1987).");
            Assert.That(MotifFinder.KnownMotifs.Kozak, Is.EqualTo("GCCGCCACCATGG"), "Kozak (1987).");
            Assert.That(MotifFinder.KnownMotifs.ShineDalgarno, Is.EqualTo("AGGAGG"), "Shine-Dalgarno.");
            Assert.That(MotifFinder.KnownMotifs.PolyASignal, Is.EqualTo("AATAAA"), "Proudfoot & Brownlee (1976).");
            Assert.That(MotifFinder.KnownMotifs.EBox, Is.EqualTo("CANNTG"), "Massari & Murre (2000).");
            Assert.That(MotifFinder.KnownMotifs.Ap1, Is.EqualTo("TGACTCA"), "Lee/Mitchell/Tjian (1987).");
            Assert.That(MotifFinder.KnownMotifs.NfKb, Is.EqualTo("GGGACTTTCC"), "Sen & Baltimore (1986).");
            Assert.That(MotifFinder.KnownMotifs.Creb, Is.EqualTo("TGACGTCA"), "Montminy et al. (1986).");
        });
    }

    #endregion

    #region FindRegulatoryElements — SHOULD

    // S1 — Multiple occurrences of one element are all reported with their positions (INV-03).
    [Test]
    public void FindRegulatoryElements_RepeatedPolyASignal_ReportsAllPositions()
    {
        // "AATAAACGAATAAA": AATAAA at index 0 and index 8.
        var polyA = Scan("AATAAACGAATAAA")
            .Where(e => e.Name == "Poly(A) Signal")
            .Select(e => e.Position)
            .OrderBy(p => p)
            .ToArray();
        Assert.That(polyA, Is.EqualTo(new[] { 0, 8 }),
            "Both AATAAA occurrences (index 0 and 8) must be reported (exhaustive scan, INV-03).");
    }

    // S2 — Multiple distinct elements in one sequence are all reported.
    [Test]
    public void FindRegulatoryElements_MultipleDistinctElements_AllReported()
    {
        // "TATAAAAGGAGG": TATAAA at index 0 (TATA Box); AGGAGG at index 6 (Shine-Dalgarno).
        var hits = Scan("TATAAAAGGAGG");
        var tata = hits.Where(e => e.Name == "TATA Box").ToArray();
        var sd = hits.Where(e => e.Name == "Shine-Dalgarno").ToArray();
        Assert.Multiple(() =>
        {
            Assert.That(tata.Length, Is.EqualTo(1), "Exactly one TATA box (TATAAA at index 0).");
            Assert.That(tata[0].Position, Is.EqualTo(0), "TATAAA begins at index 0.");
            Assert.That(sd.Length, Is.EqualTo(1), "Exactly one Shine-Dalgarno (AGGAGG at index 6).");
            Assert.That(sd[0].Position, Is.EqualTo(6), "AGGAGG begins at index 6 in TATAAAAGGAGG.");
        });
    }

    // S3 — A sequence with no library consensus yields no elements at all.
    [Test]
    public void FindRegulatoryElements_NoConsensusPresent_ReturnsNoElements()
    {
        // GGGGCCCCGCGC contains NONE of the 12 consensus strings (no GGGCGG; no AT-rich
        // TATAAA/TATAAT/AATAAA/TTGACA; no CANNTG since there is no "CA..TG"; etc.), so the
        // whole scan must be empty, not merely the two probes checked before.
        var hits = Scan("GGGGCCCCGCGC");
        Assert.That(hits, Is.Empty,
            "A sequence containing none of the 12 library consensus strings yields an empty result.");
    }

    #endregion

    #region FindRegulatoryElements — COULD

    // C1 — Invariant INV-01: matched Sequence length equals Pattern length for every hit.
    [Test]
    public void FindRegulatoryElements_EveryHit_SequenceLengthEqualsPatternLength()
    {
        var hits = Scan("GGGTATAAAGGGCACGTGGGCCAATGG");
        Assert.That(hits.Length, Is.GreaterThan(0), "Sanity: the probe contains at least one element.");
        Assert.That(hits.All(e => e.Sequence.Length == e.Pattern.Length), Is.True,
            "Each matched substring must have the same length as its consensus pattern (INV-01).");
    }

    #endregion
}
