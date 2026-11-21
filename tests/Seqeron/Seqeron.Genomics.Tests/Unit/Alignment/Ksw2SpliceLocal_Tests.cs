using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A3-9c2, B09 F187): the ksw2 kernels that minimap2 reaches through <c>align.c</c> (the two new ones and the two ported in F125), against the C functions of minimap2 v2.26 compiled with <c>-msse4.1</c> and run on random inputs
/// (TestData/Alignment/ksw2_c_oracle.tsv, scripts/alignment/ksw2/generate_ksw2_oracle.py + ksw2_oracle_driver.c): <see cref="Ksw2Aligner.LocalScore"/> = <c>ksw_ll_qinit</c> + <c>ksw_ll_i16</c> (score and end coordinates, including the zero-scored padding lanes
/// of the striped layout) on 500 pairs, <see cref="Ksw2Aligner.ExtendS2"/> = <c>ksw_exts2_sse</c> (donor / acceptor scores of the miniprot and the old splice model, both strands, reversed CIGAR, junction bonus, extension only, right-aligned gaps, approximate maximum,
/// Z-drop) on 500 spliced pairs, and <see cref="Ksw2Aligner.ExtendZ2"/> / <see cref="Ksw2Aligner.ExtendD2"/> = <c>ksw_extz2_sse</c> / <c>ksw_extd2_sse</c> (band, Z-drop, end bonus, extension only, left / right gaps, reversed CIGAR, approximate maximum / drop, score only,
/// dual affine cost with a second extension penalty 0 or 1) on 500 pairs each: every field of <c>ksw_extz_t</c> and the CIGAR with its <c>N</c> operations.
/// </summary>
[TestFixture]
public class Ksw2SpliceLocal_Tests
{
    private static readonly Lazy<List<string[]>> Rows = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.ksw2_c_oracle.tsv";
        using var stream = typeof(Ksw2SpliceLocal_Tests).Assembly.GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                rows.Add(line.Split('\t'));
        return rows;
    });

    private static int I(string s) => int.Parse(s, CultureInfo.InvariantCulture);

    private static readonly sbyte[] Matrix = Ksw2Aligner.SimpleMatrix(5, 2, -4, 1);

    [Test]
    public void LocalScore_EqualsKswLlI16_OnAllCases()
    {
        var rows = Rows.Value.Where(r => r[0] == "LL").ToList();
        Assert.That(rows, Has.Count.EqualTo(500));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            var got = Ksw2Aligner.LocalScore(Ksw2Aligner.EncodeNucleotides(r[1]), Ksw2Aligner.EncodeNucleotides(r[2]), 5, Matrix, I(r[3]), I(r[4]));
            if (got.Score != I(r[5]) || got.QueryEnd != I(r[6]) || got.TargetEnd != I(r[7]))
                failures.Add($"{r[1]} / {r[2]} gap {r[3]},{r[4]}: expected {r[5]} {r[6]} {r[7]}, actual {got}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(3)));
    }

    [Test]
    public void ExtendS2_EqualsKswExts2Sse_OnAllCases()
    {
        var rows = Rows.Value.Where(r => r[0] == "S2").ToList();
        Assert.That(rows, Has.Count.EqualTo(500));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            string target = r[2];
            byte[]? junc = r[10] == "-" ? null : [.. r[10].Select(c => (byte)(c - '0'))];
            var got = Ksw2Aligner.ExtendS2(Ksw2Aligner.EncodeNucleotides(r[1]), Ksw2Aligner.EncodeNucleotides(target), 5, Matrix, I(r[3]), I(r[4]), I(r[5]), I(r[6]), I(r[7]), I(r[8]), (Ksw2Option)I(r[9]), junc);
            string expected = string.Join(' ', r[11..22]) + (r.Length > 22 ? " " + string.Join(' ', r[22..]) : "");
            string actual = $"{got.Max} {(got.ZDropped ? 1 : 0)} {got.MaxQ} {got.MaxT} {got.Mqe} {got.MqeT} {got.Mte} {got.MteQ} {got.Score} {(got.ReachEnd ? 1 : 0)} {got.Cigar.Count}"
                + (got.Cigar.Count > 0 ? " " + string.Join(' ', got.Cigar.Select(c => c.Length.ToString(CultureInfo.InvariantCulture) + c.Operation)) : "");
            if (expected != actual)
                failures.Add($"flag {r[9]} {r[1]} / {r[2]}:\n   expected {expected}\n   actual   {actual}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(3)));
    }

    private static string Fields(Ksw2Result got) =>
        $"{got.Max} {(got.ZDropped ? 1 : 0)} {got.MaxQ} {got.MaxT} {got.Mqe} {got.MqeT} {got.Mte} {got.MteQ} {got.Score} {(got.ReachEnd ? 1 : 0)} {got.Cigar.Count}"
        + (got.Cigar.Count > 0 ? " " + string.Join(' ', got.Cigar.Select(c => c.Length.ToString(CultureInfo.InvariantCulture) + c.Operation)) : "");

    private static void CheckExtend(string kind, int columns, Func<string[], Ksw2Result> run)
    {
        var rows = Rows.Value.Where(r => r[0] == kind).ToList();
        Assert.That(rows, Has.Count.EqualTo(500));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            string expected = string.Join(' ', r[columns..]);
            string actual = Fields(run(r));
            if (expected != actual)
                failures.Add($"{string.Join(' ', r[3..columns])}  {r[1]} / {r[2]}:\n   expected {expected}\n   actual   {actual}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(3)));
    }

    [Test]
    public void ExtendZ2_EqualsKswExtz2Sse_OnAllCases() => CheckExtend("Z2", 9, r =>
        Ksw2Aligner.ExtendZ2(Ksw2Aligner.EncodeNucleotides(r[1]), Ksw2Aligner.EncodeNucleotides(r[2]), 5, Matrix, I(r[3]), I(r[4]), I(r[5]), I(r[6]), I(r[7]), (Ksw2Option)I(r[8])));

    [Test]
    public void ExtendD2_EqualsKswExtd2Sse_OnAllCases() => CheckExtend("D2", 11, r =>
        Ksw2Aligner.ExtendD2(Ksw2Aligner.EncodeNucleotides(r[1]), Ksw2Aligner.EncodeNucleotides(r[2]), 5, Matrix, I(r[3]), I(r[4]), I(r[5]), I(r[6]), I(r[7]), I(r[8]), I(r[9]), (Ksw2Option)I(r[10])));

    [Test]
    public void ExtendS2_SplicesAnIntronAsANSkip()
    {
        // exon 1 (30 bp) + GT...AG intron (100 bp) + exon 2 (30 bp) on the target, the two exons joined on the query
        string e1 = "ACGTTGCAGGCTAAGCTTCGAGGATCCTAG", e2 = "TTGACCGATAGGCTCATCGGATCCAGTCAA";
        string intron = "GT" + new string('A', 96) + "AG";
        var r = Ksw2Aligner.ExtendS2(Ksw2Aligner.EncodeNucleotides(e1 + e2), Ksw2Aligner.EncodeNucleotides(e1 + intron + e2), 5, Matrix, 4, 2, 24, 9, -1, 9, Ksw2Option.SpliceForward | Ksw2Option.SpliceComplex);
        Assert.That(r.CigarString, Is.EqualTo("30M100N30M"));
    }

    [Test]
    public void LocalScore_OfIdenticalSequences_IsTwiceTheLength()
    {
        var q = Ksw2Aligner.EncodeNucleotides("ACGTACGTTT");
        var r = Ksw2Aligner.LocalScore(q, q, 5, Matrix, 4, 2);
        Assert.That(r.Score, Is.EqualTo(20));
        Assert.That(r.TargetEnd, Is.EqualTo(9));
    }
}
