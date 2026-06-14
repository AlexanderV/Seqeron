// SEQ-CODON-FREQ-001 — Codon Frequencies (non-overlapping in-frame triplet usage)
// Evidence: docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-CODON-FREQ-001.md
// Source: Nakamura Y, Gojobori T, Ikemura T (2000). Nucleic Acids Res 28(1):292. DOI 10.1093/nar/28.1.292.
//         Kazusa Codon Usage Database (CUTG) README, https://www.kazusa.or.jp/codon/readme_codon.html.
//         EMBOSS cusp documentation, https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html.

using System.Linq;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_CalculateCodonFrequencies_Tests
{
    // Expected values are exact rationals derived from the Kazusa CUTG count/total
    // definition (frequency = count(codon) / total counted codons), computed
    // independently of the implementation.
    private const double Tolerance = 1e-10;

    #region CalculateCodonFrequencies

    // M1 — frame 0 of ATGATGAAA: non-overlapping triplets ATG, ATG, AAA -> ATG=2/3, AAA=1/3.
    // Evidence: Kazusa CUTG count/total over non-overlapping in-frame triplets.
    [Test]
    public void CalculateCodonFrequencies_Frame0_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(freq["ATG"], Is.EqualTo(2.0 / 3.0).Within(Tolerance),
                "ATG occurs 2 of 3 counted codons in frame 0 -> 2/3 (count/total)");
            Assert.That(freq["AAA"], Is.EqualTo(1.0 / 3.0).Within(Tolerance),
                "AAA occurs 1 of 3 counted codons in frame 0 -> 1/3");
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "ATG", "AAA" }),
                "only ATG and AAA are observed in frame 0 (INV-01: only observed codons are keys)");
        });
    }

    // M2 — frame 1 of ATGATGAAA: from index 1 the triplets are TGA, TGA (then 'AA' leftover) -> TGA=1.0.
    // Evidence: CUTG non-overlapping triplets read from the frame start.
    [Test]
    public void CalculateCodonFrequencies_Frame1_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 1);

        Assert.Multiple(() =>
        {
            Assert.That(freq["TGA"], Is.EqualTo(1.0).Within(Tolerance),
                "TGA is the only codon in frame 1 (both triplets) -> 1.0");
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "TGA" }),
                "trailing 'AA' is ignored; only TGA is counted (frame offset changes the multiset)");
        });
    }

    // M3 — INV-02: codon frequencies over all counted codons sum to 1.0.
    // Evidence: Kazusa CUTG count/total normalization.
    [Test]
    public void CalculateCodonFrequencies_Frame0_SumsToOne()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 0);

        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
            "Sum of codon frequencies = total/total = 1.0 (INV-02)");
    }

    // M4 — INV-03: non-ACGT triplets excluded. ATGNNNAAA frame 0 -> ATG, NNN(excluded), AAA -> each 1/2.
    // Evidence: Kazusa CUTG "codons containing ambiguous base were excluded from count".
    [Test]
    public void CalculateCodonFrequencies_NonAcgtTriplet_ExcludedFromCounts()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGNNNAAA", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(freq.ContainsKey("NNN"), Is.False,
                "NNN contains non-ACGT bases and is excluded from the table (INV-03)");
            Assert.That(freq["ATG"], Is.EqualTo(0.5).Within(Tolerance),
                "ATG is 1 of 2 valid codons (NNN not counted in total) -> 1/2");
            Assert.That(freq["AAA"], Is.EqualTo(0.5).Within(Tolerance),
                "AAA is 1 of 2 valid codons -> 1/2");
        });
    }

    // M5 — INV-05: count/total equals Kazusa per-thousand frequency / 1000. Reproduce the EMBOSS cusp
    // sample dataset exactly (Sum of Number column = 386 codons) and check two codons against the
    // published per-thousand values: CGC Number=22 -> 56.995, GGC Number=23 -> 59.585.
    // Evidence: EMBOSS cusp sample output (22/386*1000 = 56.995; 23/386*1000 = 59.585).
    [Test]
    public void CalculateCodonFrequencies_CuspDataset_MatchesPerThousandDividedBy1000()
    {
        // EMBOSS cusp sample output: each codon repeated by its "Number" column value.
        var cuspNumbers = new (string Codon, int Number)[]
        {
            ("GCA", 3), ("GCC", 18), ("GCG", 18), ("GCT", 0),
            ("TGC", 4), ("TGT", 0), ("GAC", 19), ("GAT", 3),
            ("GAA", 7), ("GAG", 19), ("TTC", 11), ("TTT", 0),
            ("GGA", 2), ("GGC", 23), ("GGG", 4), ("GGT", 3),
            ("CAC", 8), ("CAT", 3), ("ATA", 0), ("ATC", 16),
            ("ATT", 4), ("AAA", 0), ("AAG", 2), ("CTA", 0),
            ("CTC", 7), ("CTG", 15), ("CTT", 0), ("TTA", 0),
            ("TTG", 4), ("ATG", 6), ("AAC", 11), ("AAT", 0),
            ("CCA", 2), ("CCC", 6), ("CCG", 17), ("CCT", 2),
            ("CAA", 1), ("CAG", 15), ("AGA", 0), ("AGG", 1),
            ("CGA", 0), ("CGC", 22), ("CGG", 11), ("CGT", 1),
            ("AGC", 7), ("AGT", 2), ("TCA", 0), ("TCC", 6),
            ("TCG", 7), ("TCT", 1), ("ACA", 0), ("ACC", 11),
            ("ACG", 4), ("ACT", 0), ("GTA", 1), ("GTC", 13),
            ("GTG", 19), ("GTT", 0), ("TGG", 5), ("TAC", 13),
            ("TAT", 8), ("TAA", 0), ("TAG", 0), ("TGA", 1),
        };

        var sb = new System.Text.StringBuilder();
        foreach (var (codon, number) in cuspNumbers)
            for (int i = 0; i < number; i++)
                sb.Append(codon);
        string sequence = sb.ToString();

        int totalCodons = cuspNumbers.Sum(c => c.Number);
        var freq = SequenceStatistics.CalculateCodonFrequencies(sequence, readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(totalCodons, Is.EqualTo(386),
                "EMBOSS cusp sample: sum of the Number column is 386 codons");
            Assert.That(freq["CGC"], Is.EqualTo(22.0 / 386.0).Within(Tolerance),
                "CGC = 22/386 count/total = 0.056995 = cusp 56.995 per-thousand / 1000 (INV-05)");
            Assert.That(freq["CGC"] * 1000.0, Is.EqualTo(56.995).Within(1e-3),
                "CGC count/total x 1000 reproduces the cusp per-thousand value 56.995");
            Assert.That(freq["GGC"], Is.EqualTo(23.0 / 386.0).Within(Tolerance),
                "GGC = 23/386 count/total = cusp 59.585 / 1000 (INV-05)");
            Assert.That(freq["GGC"] * 1000.0, Is.EqualTo(59.585).Within(1e-3),
                "GGC count/total x 1000 reproduces the cusp per-thousand value 59.585");
        });
    }

    // S1 — trailing 1-2 bases are ignored: ATGAA frame 0 -> only ATG; ATG=1.0.
    // Evidence: CUTG non-overlapping triplets; remainder is no codon.
    [Test]
    public void CalculateCodonFrequencies_TrailingBases_Ignored()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGAA", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(freq["ATG"], Is.EqualTo(1.0).Within(Tolerance),
                "only the full ATG triplet is counted -> 1.0");
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "ATG" }),
                "trailing 'AA' adds no partial codon");
        });
    }

    // S2 — INV-04: case-insensitive. Lowercase 'atgaaa' equals uppercase 'ATGAAA'.
    // Evidence: codons are case-independent; implementation upper-cases the input.
    [Test]
    public void CalculateCodonFrequencies_LowercaseInput_EqualsUppercase()
    {
        var lower = SequenceStatistics.CalculateCodonFrequencies("atgaaa", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(lower["ATG"], Is.EqualTo(0.5).Within(Tolerance),
                "lowercase 'atg' is counted as ATG -> 1 of 2 codons (INV-04)");
            Assert.That(lower["AAA"], Is.EqualTo(0.5).Within(Tolerance),
                "lowercase 'aaa' is counted as AAA -> 1 of 2 codons");
            Assert.That(lower.Keys, Is.EquivalentTo(new[] { "ATG", "AAA" }),
                "keys are upper-cased; case does not change the result");
        });
    }

    // S3 — all-ambiguous sequence: total = 0 -> empty table, no division by zero.
    // Evidence: CUTG ambiguous codons excluded; count/total undefined for total=0 -> empty.
    [Test]
    public void CalculateCodonFrequencies_AllAmbiguous_ReturnsEmpty()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("NNNNNN", readingFrame: 0);

        Assert.That(freq, Is.Empty,
            "every triplet is ambiguous, so total = 0 and the table is empty (no division by zero)");
    }

    // S4 — INV-01: all frequencies are in (0, 1].
    // Evidence: count(x) >= 1 and count(x) <= total for every key.
    [Test]
    public void CalculateCodonFrequencies_AllValues_AreInUnitInterval()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 0);

        Assert.That(freq.Values, Is.All.InRange(double.Epsilon, 1.0),
            "every reported frequency is in (0, 1] (INV-01)");
    }

    // C1 — guards: null, empty, and length < 3 return an empty dictionary.
    // Evidence: no full codon to count -> empty (shared guard contract).
    [Test]
    public void CalculateCodonFrequencies_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateCodonFrequencies(null!), Is.Empty,
                "null input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateCodonFrequencies(string.Empty), Is.Empty,
                "empty input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateCodonFrequencies("AT"), Is.Empty,
                "length < 3 yields empty dictionary (no full codon)");
        });
    }

    #endregion
}
