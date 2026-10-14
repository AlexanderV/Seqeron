// TRANS-SIXFRAME-001 — FindOrfs with the EMBOSS getorf -find modes 0–6 (B02 FIN WP11, audit A1-12)
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-SIXFRAME-001.md
// Source: EMBOSS getorf.c (getorf_FindORFs / getorf_WriteORF / getorf_AppORF) and acd/getorf.acd
//         (raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/); reference values are the
//         output of the EMBOSS 6.6.0 getorf binary (TestData/Getorf/getorf_find_modes.tsv).

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class Translator_GetOrfModes_Tests
{
    // Worked example (EMBOSS 6.6.0: getorf -find N -minsize 9 -flanking 2 -table 1), L = 18.
    private const string Worked = "GGGATGAAACCCTAAGGG";

    #region Worked examples (getorf binary output)

    // M16 — -find 0: STOP→STOP translations, leading regions open at the sequence start, trailing
    // regions open at the last complete codon; getorf headers [1 - 12] … [18 - 1] (REVERSE SENSE).
    [Test]
    public void FindOrfs_Find0_WorkedExample_MatchesGetorf()
    {
        var dna = new DnaSequence(Worked);

        var result = Translator.FindOrfs(dna, OrfFindMode.TranslationStopToStop, minSize: 9).ToList();

        Assert.That(result.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated, r.Sequence)), Is.EqualTo(new[]
        {
            (0, 14, 1, true, "GMKP"),     // getorf [1 - 12]  (stop TAA at 12..14 excluded)
            (7, 15, 2, false, "NPK"),     // getorf [8 - 16]
            (2, 16, 3, false, "DETLR"),   // getorf [3 - 17]
            (7, 15, -2, false, "GFI"),    // getorf [11 - 3] (REVERSE SENSE)
            (2, 16, -3, false, "LRVSS"),  // getorf [16 - 2] (REVERSE SENSE)
            (0, 17, -1, false, "PLGFHP"), // getorf [18 - 1] (REVERSE SENSE)
        }));
        Assert.That(result.Select(r => r.ToGetOrfRange(Worked.Length)), Is.EqualTo(new[]
        {
            (1, 12), (8, 16), (3, 17), (11, 3), (16, 2), (18, 1),
        }));
    }

    // M17 — -find 2: the same regions as nucleotides (STOP excluded).
    [Test]
    public void FindOrfs_Find2_WorkedExample_MatchesGetorf()
    {
        var result = Translator.FindOrfs(new DnaSequence(Worked), OrfFindMode.NucleotidesStopToStop, minSize: 9)
            .Select(r => r.Sequence).ToList();

        Assert.That(result, Is.EqualTo(new[]
        {
            "GGGATGAAACCC", "AACCCTAAG", "GATGAAACCCTAAGG", "GGTTTCATC", "CTTAGGGTTTCATCC", "CCCTTAGGGTTTCATCCC",
        }));
    }

    // M18 — -find 3: START→STOP nucleotides; getorf [4 - 12] ATGAAACCC.
    [Test]
    public void FindOrfs_Find3_WorkedExample_MatchesGetorf()
    {
        var result = Translator.FindOrfs(new DnaSequence(Worked), OrfFindMode.NucleotidesStartToStop, minSize: 9).ToList();

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That((result[0].StartPosition, result[0].EndPosition, result[0].Frame, result[0].IsTerminated, result[0].Sequence),
            Is.EqualTo((3, 14, 1, true, "ATGAAACCC")));
        Assert.That(result[0].ToGetOrfRange(Worked.Length), Is.EqualTo((4, 12)));
    }

    // M19 — -find 4/5/6 -flanking 2: windows [codon − 2, codon + 4]; getorf "Around codon at 4/5/13".
    [Test]
    public void FindOrfs_FlankingModes_WorkedExample_MatchGetorf()
    {
        var dna = new DnaSequence(Worked);

        var start = Translator.FindOrfs(dna, OrfFindMode.FlankingStart, minSize: 9, flanking: 2).ToList();
        var initial = Translator.FindOrfs(dna, OrfFindMode.FlankingInitialStop, minSize: 9, flanking: 2).ToList();
        var ending = Translator.FindOrfs(dna, OrfFindMode.FlankingEndingStop, minSize: 9, flanking: 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(start.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.CodonPosition, r.Sequence)),
                Is.EqualTo(new[] { (1, 7, 1, 3, "GGATGAA") }));                      // [2 - 8], codon at 4
            Assert.That(initial.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.CodonPosition, r.Sequence)),
                Is.EqualTo(new[] { (2, 8, 2, 4, "GATGAAA"), (2, 8, -2, 4, "CTTAGGG") })); // [3 - 9], [16 - 10] (REVERSE SENSE)
            Assert.That(initial.Select(r => r.ToGetOrfRange(Worked.Length)), Is.EqualTo(new[] { (3, 9), (16, 10) }));
            Assert.That(ending.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.CodonPosition, r.Sequence)),
                Is.EqualTo(new[] { (10, 16, 1, 12, "CCTAAGG") }));                   // [11 - 17], codon at 13
        });
    }

    #endregion

    #region getorf binary oracle (TSV)

    // M20 — every C line of the oracle: the records (getorf [s - e], strand, sequence) equal the
    // EMBOSS 6.6.0 getorf output, in order, except the documented one-codon open-ORF divergence.
    [TestCaseSource(nameof(OracleConfigIds))]
    public void FindOrfs_AllModes_MatchGetorfBinaryOracle(string configId)
    {
        var oracle = Oracle.Value;
        var c = oracle.Configs[configId];
        var code = GeneticCode.GetByTableNumber(c.Table);

        Assert.Multiple(() =>
        {
            foreach (var (name, seq) in oracle.Sequences)
            {
                var actual = Translator.FindOrfs(new DnaSequence(seq), c.Mode, code, c.MinSize, c.MaxSize,
                        c.Flanking, c.Reverse, c.Methionine)
                    .Where(r => !IsDocumentedDivergence(r, seq.Length))
                    .Select(r => { var (s, e) = r.ToGetOrfRange(seq.Length); return (s, e, r.IsReverseStrand, r.Sequence); })
                    .ToList();
                var expected = oracle.Records.TryGetValue((configId, name), out var list) ? list : [];
                Assert.That(actual, Is.EqualTo(expected), $"{configId} {name} {seq}");
            }
        });
    }

    // M21 — the START→STOP mode through the new overload is the legacy START→STOP overload:
    // same ORF set (minSize = 3 · minLength; legacy order is frame by frame).
    [Test]
    public void FindOrfs_Find1_EqualsLegacyStartToStopOverload()
    {
        Assert.Multiple(() =>
        {
            foreach (var table in new[] { 1, 2, 11 })
            {
                var code = GeneticCode.GetByTableNumber(table);
                foreach (var (name, seq) in Oracle.Value.Sequences)
                {
                    var dna = new DnaSequence(seq);
                    var legacy = Translator.FindOrfs(dna, code, minLength: 2)
                        .Select(o => (o.StartPosition, o.EndPosition, o.Frame, o.Protein.Sequence)).Order().ToList();
                    var getorf = Translator.FindOrfs(dna, OrfFindMode.TranslationStartToStop, code, minSize: 6)
                        .Select(o => (o.StartPosition, o.EndPosition, o.Frame, o.Sequence)).Order().ToList();
                    Assert.That(getorf, Is.EqualTo(legacy), $"table {table} {name}");
                }
            }
        });
    }

    #endregion

    #region Conventions and divergence

    // M22 — documented divergence: "ATG" alone is a one-codon open ORF; getorf -find 3 -minsize 0
    // -table 1 writes nothing (oracle config c28, sequence s4), this API reports it (as the legacy
    // overload does).
    [Test]
    public void FindOrfs_Find3_StartAsLastCodon_ReportedOneCodonOpenOrf()
    {
        var result = Translator.FindOrfs(new DnaSequence("ATG"), OrfFindMode.NucleotidesStartToStop, minSize: 0,
            searchBothStrands: false).ToList();

        Assert.That(result.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated, r.Sequence)),
            Is.EqualTo(new[] { (0, 2, 1, false, "ATG") }));
        Assert.That(Oracle.Value.Configs["c28"], Is.EqualTo(new OracleConfig(OrfFindMode.NucleotidesStartToStop, 1, 0,
            1_000_000, 0, true, true)));
        Assert.That(Oracle.Value.Sequences[4].Sequence, Is.EqualTo("ATG"));
        Assert.That(Oracle.Value.Records.ContainsKey(("c28", "s4")), Is.False,
            "getorf -find 3 -minsize 0 writes no ORF for 'ATG'");
    }

    // M23 — getorf -minsize 0 writes empty STOP→STOP regions (adjacent STOPs): "TAATAA" -find 0
    // -reverse N gives, in getorf order, the empty leading region [1 - 0], the open frame-2/3 regions
    // [2 - 4] N and [3 - 5] I, and the empty region between the STOPs [4 - 3].
    [Test]
    public void FindOrfs_Find0_MinSizeZero_ReportsEmptyRegions()
    {
        var result = Translator.FindOrfs(new DnaSequence("TAATAA"), OrfFindMode.TranslationStopToStop, minSize: 0,
            searchBothStrands: false).ToList();

        Assert.That(result.Select(r => (r.StartPosition, r.EndPosition, r.Sequence, r.ToGetOrfRange(6))), Is.EqualTo(new[]
        {
            (0, 2, "", (1, 0)),
            (1, 3, "N", (2, 4)),
            (2, 4, "I", (3, 5)),
            (3, 5, "", (4, 3)),
        }));
    }

    // M24 — getorf.acd -flanking: "If the region of flanking nucleotides crosses the start or end of
    // the sequence, no output is given for this codon"; the leading STOP-free region has no initial STOP.
    [Test]
    public void FindOrfs_FlankingWindowCrossingSequenceEnd_IsNotReported()
    {
        var dna = new DnaSequence(Worked);

        Assert.Multiple(() =>
        {
            Assert.That(Translator.FindOrfs(dna, OrfFindMode.FlankingStart, minSize: 9, flanking: 4), Is.Empty);
            Assert.That(Translator.FindOrfs(dna, OrfFindMode.FlankingEndingStop, minSize: 9, flanking: 2,
                searchBothStrands: false).Select(r => r.CodonPosition), Is.EqualTo(new[] { 12 }));
            Assert.That(Translator.FindOrfs(dna, OrfFindMode.FlankingEndingStop, minSize: 9, flanking: 6), Is.Empty);
        });
    }

    // M25 — getorf -methionine N: the initiator keeps its own residue (TTG = L in table 1).
    [Test]
    public void FindOrfs_Find1_MethionineOff_InitiatorKeepsItsResidue()
    {
        var dna = new DnaSequence("TTGAAAGGGTAA");

        var on = Translator.FindOrfs(dna, OrfFindMode.TranslationStartToStop, minSize: 0, searchBothStrands: false).Single();
        var off = Translator.FindOrfs(dna, OrfFindMode.TranslationStartToStop, minSize: 0, searchBothStrands: false,
            methionine: false).Single();

        Assert.That((on.Sequence, off.Sequence), Is.EqualTo(("MKG", "LKG")));
    }

    // M26 — argument validation; dual-coding-stop tables rejected as by the legacy overload.
    [Test]
    public void FindOrfs_GetOrfModes_InvalidArguments_Throw()
    {
        var dna = new DnaSequence("ATGTGATAAGGG");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => Translator.FindOrfs(null!, OrfFindMode.TranslationStopToStop));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.FindOrfs(dna, (OrfFindMode)7));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.FindOrfs(dna, OrfFindMode.TranslationStopToStop, minSize: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.FindOrfs(dna, OrfFindMode.TranslationStopToStop, maxSize: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.FindOrfs(dna, OrfFindMode.FlankingStart, flanking: -1));
            Assert.Throws<ArgumentException>(() => Translator.FindOrfs(dna, OrfFindMode.NucleotidesStopToStop,
                GeneticCode.GetByTableNumber(27)));
        });
    }

    #endregion

    #region Oracle loading

    // A one-codon open ORF opened at the last complete codon of its frame (documented divergence).
    private static bool IsDocumentedDivergence(OrfFeature r, int length) =>
        r.Mode is OrfFindMode.TranslationStartToStop or OrfFindMode.NucleotidesStartToStop or OrfFindMode.FlankingStart
        && !r.IsTerminated
        && (r.Mode == OrfFindMode.FlankingStart
            ? r.CodonPosition + 6 > length
            : r.EndPosition - r.StartPosition == 2 && r.StartPosition + 6 > length);

    private sealed record OracleConfig(OrfFindMode Mode, int Table, int MinSize, int MaxSize, int Flanking,
        bool Reverse, bool Methionine);

    private sealed record OracleData(
        List<(string Name, string Sequence)> Sequences,
        Dictionary<string, OracleConfig> Configs,
        Dictionary<(string, string), List<(int, int, bool, string)>> Records);

    private static readonly Lazy<OracleData> Oracle = new(LoadOracle);

    private static IEnumerable<string> OracleConfigIds() => Oracle.Value.Configs.Keys;

    private static OracleData LoadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Getorf.getorf_find_modes.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var data = new OracleData([], [], []);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line.StartsWith('#')) continue;
            var f = line.Split('\t');
            switch (f[0])
            {
                case "S":
                    data.Sequences.Add((f[1], f[2]));
                    break;
                case "C":
                    data.Configs[f[1]] = new OracleConfig((OrfFindMode)int.Parse(f[2]), int.Parse(f[3]), int.Parse(f[4]),
                        int.Parse(f[5]), int.Parse(f[6]), f[7] == "Y", f[8] == "Y");
                    break;
                case "R":
                    var key = (f[1], f[2]);
                    if (!data.Records.TryGetValue(key, out var list))
                        data.Records[key] = list = [];
                    list.Add((int.Parse(f[4]), int.Parse(f[5]), f[6] == "R", f.Length > 7 ? f[7] : ""));
                    break;
            }
        }
        return data;
    }

    #endregion
}
