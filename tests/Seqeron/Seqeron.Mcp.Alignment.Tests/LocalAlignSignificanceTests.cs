using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-7 (B09 F300): the <c>significance</c> option of <c>local_align_affine</c> / <c>local_align</c> / <c>water_align</c>: bit score, E-value, P-value and the Karlin-Altschul parameters of the score. Truth: the installed BLAST+
/// 2.12.0 on the same pairs (protein: <c>blastp -comp_based_stats 0 -seg no -subject</c>, 150 x 160, score 256: bit score 103.219, E 1.24008e-33; nucleotide: <c>blastn -task blastn -dust no</c>, 430 x 880, score 145:
/// bit score 132.031, E 1.21e-34), the same numbers the library test locks over 2 283 HSPs.
/// </summary>
[TestFixture]
public class LocalAlignSignificanceTests
{
    private const string ProteinQuery = "MDLNLSMGSSGCKANPAVQNPWARCGYHEIRNTNTKREWNYKNTARICGENSIPVEWIAHQKGPGDFYYRFFAAHHGGLMHVHGHPLANQGFKDMLYWAYLCQDHMTYNFMKVDLMLGDFLSGCDYVPCIYNKRQFCCSMHFWFQEGQNF";
    private const string ProteinSubject = "CQLFRYGTRSMSKLSPFEPVGSMGDSKTVTNDNWCLNVKSKLMGWASVKMKRLTNNKNQNGRNMTFTGHNYKNTAWICGENSIWVEGISHFGGPGDFHNPFFAAHHGGLMHVHQHPLYNQFFKDELYWAYWPQQISLSPPGYYKLSKQAMLSLFSAEYRI";
    private const string DnaQuery = "TATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGATTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTCCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTATTCCCCCCGCGGCCCACCCAGTATTCCTAA";
    private const string DnaSubject = "CGGAGCATAAATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAACGTTACTAAGGGGTATAATCTTCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAGTCTTGTGCGCCTAGGGGCGCCCCAAAGGTAAACGAACCGTTGCGGTCAATCTTGTCGCGGCTGATGAATTTGAAGCAGTGGCCGGGAGTGTGTGCTCAGGAGTTCGTCCCATGACACGATAGAGAGAGAACATCCTGTTGTTTCCTCATGCAATTCAAAACCATGTCCGTAGTGTAGGTGAAATAGTAAACCATTTTACGGAGGATACCAACTTCCTCCTGCGCGACCACCCAGGATTAGACTCATCATTCGGGTAGTAGACATTATATTCGATACCGTGGTAGCCTAGGGTGTTAACACCCCTATAACACATTAGTCCCTTGTATGCAGGCGGTATCGGACGGCGCCCACACCTTGGAGGTATCCAGCGCAAGGCGCCATATCCGTACCTTACTATCGCGCGAACTTATGTTGTTTTAAGTTAGAGTTGGACATCTATACGTCAGTCCTAAACATAGCGAGCATTTCGCAGATGGGTCTCCGACGGTACCCCAAGGGTCGTTACCGACGCCGGGACGCCGCATATAAAGGTACGCCCGACCATTATACAGGTAGCCATCTGCGTCTGACATCGCATTTGAAACCCAGTAGGTACTGCCTTAGTTGCACTCCTAACTCATGTTAACGGACTTACGGGCACTAGCTTCTTACTGCCCTCTCTGTTTCTCTTAAGGGACGTCGAGACGCCAAGTTATGGAGTCTACCCACGTTTCGGTTCCG";

    [Test]
    public void LocalAlignAffine_Protein_Blosum62_11_1_EqualsBlastp()
    {
        var r = AlignmentTools.LocalAlignAffine(ProteinQuery, ProteinSubject, gapOpen: -11, gapExtend: -1, matrix: "BLOSUM62", significance: true);
        Assert.That(r.Score, Is.EqualTo(256));
        Assert.That(r.Significance, Is.Not.Null);
        var s = r.Significance!;
        Assert.That(s.Method, Is.EqualTo("SpougeFiniteSizeCorrection"));
        Assert.That(s.BitScore, Is.EqualTo(103.219).Within(0.001));
        Assert.That(s.EValue, Is.EqualTo(1.24008e-33).Within(1e-38));
        Assert.That(s.PValue, Is.EqualTo(s.EValue).Within(1e-40));
        Assert.That((s.Lambda, s.K, s.H), Is.EqualTo((0.267, 0.041, 0.14)));
    }

    [Test]
    public void LocalAlignAffine_Dna_Blastn_2_m3_5_2_EqualsBlastn()
    {
        var r = AlignmentTools.LocalAlignAffine(DnaQuery, DnaSubject, match: 2, mismatch: -3, gapOpen: -5, gapExtend: -2, significance: true);
        Assert.That(r.Score, Is.EqualTo(145));
        var s = r.Significance!;
        Assert.That(s.Method, Is.EqualTo("LengthAdjustedKarlinAltschul"));
        Assert.That(s.BitScore, Is.EqualTo(132.031).Within(0.001));
        Assert.That(s.EValue, Is.EqualTo(1.21e-34).Within(0.01e-34), "blastn prints 1.21e-34 (odd score of an even-only table: rounded down to 144 for the E-value)");
        Assert.That(s.LengthAdjustment, Is.EqualTo(13));
    }

    [Test]
    public void DbLength_ScalesTheEValue_AndOffByDefault()
    {
        var plain = AlignmentTools.LocalAlignAffine(ProteinQuery, ProteinSubject, gapOpen: -11, gapExtend: -1, matrix: "BLOSUM62");
        Assert.That(plain.Significance, Is.Null);
        var one = AlignmentTools.LocalAlignAffine(ProteinQuery, ProteinSubject, gapOpen: -11, gapExtend: -1, matrix: "BLOSUM62", significance: true).Significance!;
        var db = AlignmentTools.LocalAlignAffine(ProteinQuery, ProteinSubject, gapOpen: -11, gapExtend: -1, matrix: "BLOSUM62", significance: true, dbLength: 1600).Significance!;
        Assert.That(db.EValue / one.EValue, Is.EqualTo(10.0).Within(1e-9), "a database 10 times as long: 10 times the E-value");
    }

    [Test]
    public void LocalAlign_LinearGaps_UseTheLinearRowOfTheNucleotideTable()
    {
        // blastn_values_1_2 row { 0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89 }: reward 1 / penalty -2 with linear gaps (existence 0, extension 2); the row { 0, 0, ... } is the non-affine (megablast) row
        var r = AlignmentTools.LocalAlign("ACGTACGTACGTACGTACGTACGTACGTACGT", "TTACGTACGTACGTACGTACGTACGTACGTACGTTT", match: 1, mismatch: -2, gapExtend: -2, significance: true);
        Assert.That(r.Significance, Is.Not.Null);
        Assert.That((r.Significance!.Lambda, r.Significance.K), Is.EqualTo((1.19, 0.34)));
    }

    [Test]
    public void WaterAlign_Protein_IntegerGapCosts_UseExistenceMinusExtension()
    {
        // EMBOSS gapopen 12 / gapextend 1 = BLAST existence 11 extension 1 (a gap of length k costs 12 + (k - 1)): the same row of the BLOSUM62 table
        var r = AlignmentTools.WaterAlign(ProteinQuery, ProteinSubject, matrix: "EBLOSUM62", gapOpen: 12, gapExtend: 1, significance: true);
        Assert.That(r.Significance, Is.Not.Null);
        Assert.That((r.Significance!.Lambda, r.Significance.K), Is.EqualTo((0.267, 0.041)));
        Assert.That(r.Significance.Method, Is.EqualTo("SpougeFiniteSizeCorrection"));
    }

    [Test]
    public void SignificanceOutsideTheTables_IsRefused()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignAffine(ProteinQuery, ProteinSubject, gapOpen: -7, gapExtend: -1, matrix: "BLOSUM62", significance: true));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignAffine(DnaQuery, DnaSubject, match: 7, mismatch: -3, gapOpen: -5, gapExtend: -2, significance: true));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignAffine(ProteinQuery, ProteinSubject, gapOpen: -10, gapExtend: -1, matrix: "BLOSUM80", significance: true), "the EMBOSS BLOSUM80 is not the NCBI matrix");
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.LocalAlignAffine("ACGT", "ACGT", significance: true, dbLength: 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignAffine("ACGT", "ACGT", match: 2, mismatch: -3, gapOpen: -2.5, gapExtend: -2, significance: true), "fractional gap costs");
    }
}

/// <summary>R14-7 (B09 F300): <c>local_align_zscore</c>, the shuffle Z-score of a local score (EMBOSS prss / FASTA idea; seeded, so reproducible).</summary>
[TestFixture]
public class LocalAlignZScoreTests
{
    private const string Core = "ACGTTGCAAGGCTTAACCGGTTAACGTACGGTACCATGCA";

    [Test]
    public void ZScore_OfASharedMotif_IsLarge_AndReproducible()
    {
        var a = AlignmentTools.LocalAlignZScore("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC", shuffles: 200, seed: 7);
        var b = AlignmentTools.LocalAlignZScore("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC", shuffles: 200, seed: 7);
        Assert.That(a, Is.EqualTo(b));
        Assert.That(a.Score, Is.EqualTo(AlignmentTools.LocalAlign("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC").Score), "the score of the pair is the local_align score");
        Assert.That(a.ZScore, Is.GreaterThan(5));
        Assert.That(a.StandardDeviation, Is.GreaterThan(0));
        Assert.That((a.Shuffles, a.Seed), Is.EqualTo((200, 7)));
        Assert.That(AlignmentTools.LocalAlignZScore("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC", shuffles: 200, seed: 8).Mean, Is.Not.EqualTo(a.Mean), "another seed, other shuffles");
    }

    [Test]
    public void ZScore_Validation()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignZScore(null!, "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.LocalAlignZScore("ACGT", "ACGT", shuffles: 1));
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => AlignmentTools.LocalAlignZScore("ACGTACGT", "ACGTACGT", cancellationToken: cts.Token));
    }
}
