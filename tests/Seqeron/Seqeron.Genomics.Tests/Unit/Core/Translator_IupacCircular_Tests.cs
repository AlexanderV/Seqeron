// TRANS-SIXFRAME-001 — IUPAC-aware string overloads of TranslateSixFrames / FindOrfs (B02 FIN WP12,
// audit A1-13) and getorf -circular (audit A1-20).
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-SIXFRAME-001.md
// Sources: Biopython 1.88 Seq.translate / reverse_complement (six frames); EMBOSS getorf.c
//          (getorf_FindORFs circular branch, getorf_WriteORF) and ajtranslate.c (ajTrnCodonstrTypeC,
//          trnComplete, ajTrnSeqDangleC) from raw.githubusercontent.com/kimrutherford/EMBOSS/master/;
//          reference values are the output of Biopython 1.88 and of the EMBOSS 6.6.0 getorf binary
//          (TestData/Getorf/getorf_iupac_circular.tsv).

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class Translator_IupacCircular_Tests
{
    #region Six frames on IUPAC strings (Biopython 1.88)

    // M27 — frame +k = Seq(s)[k-1:].translate(table), frame −k = Seq(s).reverse_complement()[k-1:]
    // .translate(table); values printed by Biopython 1.88 (reverse complement CARTTARTYNGCYTANNNCAT).
    [TestCase("ATGNNNTARGCNRAYTAAYTG", 1, "MX*AB*L", "XXXXXN", "XXXXLX", "QLXXLXH", "X*XAXX", "XXXXXX")]
    [TestCase("TTGTANSARGAYMGR", 11, "LXZDR", "CXXX", "VXXX", "XXXXQ", "XSXY", "XXXT")]
    public void TranslateSixFrames_IupacString_MatchesBiopython(string sequence, int table,
        string f1, string f2, string f3, string r1, string r2, string r3)
    {
        var frames = Translator.TranslateSixFrames(sequence, GeneticCode.GetByTableNumber(table));

        Assert.That(new[] { 1, 2, 3, -1, -2, -3 }.Select(k => frames[k].Sequence),
            Is.EqualTo(new[] { f1, f2, f3, r1, r2, r3 }));
    }

    // M28 — case-insensitive DNA/RNA input (U read as T): same frames as the upper-case DNA spelling,
    // and on A/C/G/T input the string overload equals the DnaSequence overload.
    [Test]
    public void TranslateSixFrames_StringOverload_CaseAndRnaInsensitive_EqualsDnaOverload()
    {
        static string[] Frames(IReadOnlyDictionary<int, ProteinSequence> d) =>
            new[] { 1, 2, 3, -1, -2, -3 }.Select(k => d[k].Sequence).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(Frames(Translator.TranslateSixFrames("augnnnuargcnrayuaayug")),
                Is.EqualTo(Frames(Translator.TranslateSixFrames("ATGNNNTARGCNRAYTAAYTG"))));
            const string acgt = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";
            foreach (var table in new[] { 1, 2, 11 })
            {
                var code = GeneticCode.GetByTableNumber(table);
                Assert.That(Frames(Translator.TranslateSixFrames(acgt, code)),
                    Is.EqualTo(Frames(Translator.TranslateSixFrames(new DnaSequence(acgt), code))), $"table {table}");
            }
        });
    }

    [Test]
    public void TranslateSixFrames_StringOverload_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => Translator.TranslateSixFrames((string)null!));
            Assert.Throws<ArgumentException>(() => Translator.TranslateSixFrames("ATG-CCC"));
            Assert.Throws<ArgumentException>(() => Translator.TranslateSixFrames("ATGXCC"));
        });
    }

    #endregion

    #region FindOrfs on IUPAC strings (getorf ambiguity rules)

    // M29 — getorf 6.6.0 (-find 1 -table 1 -minsize 3 -reverse N): an N codon reads through as X
    // ("ATGNNNCCCTAA" → [1 - 9] MXP); a possible stop (RTA = ATA/GTA) is X, not a STOP
    // ("atgaaacccrtaa" → open ORF [1 - 12] MKPX); U is read as T ("AUGAAACCCUAAGG" → [1 - 9] MKP).
    [TestCase("ATGNNNCCCTAA", 0, 11, true, "MXP", 1, 9)]
    [TestCase("atgaaacccrtaa", 0, 11, false, "MKPX", 1, 12)]
    [TestCase("AUGAAACCCUAAGG", 0, 11, true, "MKP", 1, 9)]
    public void FindOrfs_IupacString_Find1_MatchesGetorf(string sequence, int start, int end, bool terminated,
        string protein, int getorfStart, int getorfEnd)
    {
        var orf = Translator.FindOrfs(sequence, OrfFindMode.TranslationStartToStop, minSize: 3,
            searchBothStrands: false).Single();

        Assert.That((orf.StartPosition, orf.EndPosition, orf.Frame, orf.IsTerminated, orf.Sequence),
            Is.EqualTo((start, end, 1, terminated, protein)));
        Assert.That(orf.ToGetOrfRange(sequence.Length), Is.EqualTo((getorfStart, getorfEnd)));
    }

    // M30 — EMBOSS ajTrnCodonstrTypeC: the Starts table holds only A/C/G/T codons, so YTG (CTG/TTG, both
    // table-1 starts) does not open an ORF; getorf -find 1 writes nothing for "YTGAAACCCTAA".
    [Test]
    public void FindOrfs_IupacString_AmbiguousStartCodon_DoesNotOpenOrf()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Translator.FindOrfs("YTGAAACCCTAA", OrfFindMode.TranslationStartToStop, minSize: 3,
                searchBothStrands: false), Is.Empty);
            Assert.That(Translator.FindOrfs("YTGAAACCCTAA", minLength: 1, searchBothStrands: false), Is.Empty);
            Assert.That(Translator.FindOrfs("TTGAAACCCTAA", minLength: 1, searchBothStrands: false).Single().Protein.Sequence,
                Is.EqualTo("MKP"));
        });
    }

    // M31 — documented difference: TRA (TAA/TGA, every expansion a stop) ends the ORF here (Biopython:
    // TRA → '*'); EMBOSS's precomputed ambiguous table gives TRA → X, so getorf -find 1 prints
    // [1 - 12] MKPX for "ATGAAACCCTRA" — the ORF without the STOP is the same [1 - 9] here.
    [Test]
    public void FindOrfs_IupacString_AllStopAmbiguousCodon_TerminatesOrf()
    {
        var orf = Translator.FindOrfs("ATGAAACCCTRA", OrfFindMode.TranslationStartToStop, minSize: 3,
            searchBothStrands: false).Single();

        Assert.That((orf.StartPosition, orf.EndPosition, orf.IsTerminated, orf.Sequence), Is.EqualTo((0, 11, true, "MKP")));
        Assert.That(orf.ToGetOrfRange(12), Is.EqualTo((1, 9)));
    }

    // M32 — the legacy START→STOP overload on a string: equal to the DnaSequence overload on A/C/G/T
    // input; on IUPAC input it applies the same codon rules (N codon → X).
    [Test]
    public void FindOrfs_LegacyStringOverload_EqualsDnaOverload_AndReadsIupac()
    {
        const string acgt = "CCATGAAACCCGGGTTTTAGCCATGGCATTTAAACCCTTAGG";
        foreach (var table in new[] { 1, 2, 11 })
        {
            var code = GeneticCode.GetByTableNumber(table);
            var fromString = Translator.FindOrfs(acgt.ToLowerInvariant(), code, minLength: 1).ToList();
            var fromDna = Translator.FindOrfs(new DnaSequence(acgt), code, minLength: 1).ToList();
            Assert.That(fromString, Is.EqualTo(fromDna), $"table {table}");
        }

        var iupac = Translator.FindOrfs("ATGNNNCCCTAA", minLength: 1, searchBothStrands: false).Single();
        Assert.That((iupac.StartPosition, iupac.EndPosition, iupac.Frame, iupac.Protein.Sequence),
            Is.EqualTo((0, 11, 1, "MXP")));
    }

    [Test]
    public void FindOrfs_StringOverloads_InvalidInput_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => Translator.FindOrfs((string)null!).ToList());
            Assert.Throws<ArgumentNullException>(() => Translator.FindOrfs((string)null!, OrfFindMode.TranslationStopToStop));
            Assert.Throws<ArgumentException>(() => Translator.FindOrfs("ATG*TAA"));
            Assert.Throws<ArgumentException>(() => Translator.FindOrfs("ATGXTAA", OrfFindMode.TranslationStopToStop));
            Assert.Throws<ArgumentException>(() => Translator.FindOrfs("ATGTAA", GeneticCode.GetByTableNumber(27)));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.FindOrfs("ATGTAA", (OrfFindMode)9));
        });
    }

    #endregion

    #region getorf -circular (EMBOSS 6.6.0 binary)

    private const string Circle = "CCCTAAGGGATGAAA"; // L = 15

    // M33 — getorf -find 1 -minsize 3 -circular Y -reverse N: "[10 - 18] (ORF crosses the breakpoint) MKP";
    // linear: "[10 - 15] MK" (open at the end). EndPosition 20 ≥ L: the STOP TAA lies past the origin.
    [Test]
    public void FindOrfs_Circular_Find1_OrfCrossesOrigin()
    {
        var circular = Translator.FindOrfs(Circle, OrfFindMode.TranslationStartToStop, minSize: 3,
            searchBothStrands: false, circular: true).Single();
        var linear = Translator.FindOrfs(Circle, OrfFindMode.TranslationStartToStop, minSize: 3,
            searchBothStrands: false).Single();

        Assert.Multiple(() =>
        {
            Assert.That((circular.StartPosition, circular.EndPosition, circular.Frame, circular.IsTerminated, circular.Sequence),
                Is.EqualTo((9, 20, 1, true, "MKP")));
            Assert.That(circular.ToGetOrfRange(Circle.Length), Is.EqualTo((10, 18)));
            Assert.That((linear.Sequence, linear.ToGetOrfRange(Circle.Length)), Is.EqualTo(("MK", (10, 15))));
            Assert.That(Translator.FindOrfs(Circle, OrfFindMode.NucleotidesStartToStop, minSize: 3,
                    searchBothStrands: false, circular: true).Single().Sequence, Is.EqualTo("ATGAAACCC"));
        });
    }

    // M34 — getorf -find 0 / -find 6 -flanking 2 -circular Y: no leading region (every frame starts closed);
    // "[7 - 18] GMKP", "[14 - 25] NPKG" (both cross the breakpoint); -find 6: "[17 - 23] Around codon at 19.
    // CCTAAGG", "[24 - 30] Around codon at 26. GATGAAA". Frame 3 holds no STOP at all: its full-length record
    // [3 - 45] LRDETLRDETLRDEX (dangling base → X) is written here, while getorf 6.6.0 omits it (getorf.c
    // tests "ajSeqGetLen(pep) >= maxsize" for frame 3 only — documented divergence).
    [Test]
    public void FindOrfs_Circular_StopModes_MatchGetorf()
    {
        var find0 = Translator.FindOrfs(Circle, OrfFindMode.TranslationStopToStop, minSize: 3,
            searchBothStrands: false, circular: true).ToList();
        var find6 = Translator.FindOrfs(Circle, OrfFindMode.FlankingEndingStop, minSize: 3, flanking: 2,
            searchBothStrands: false, circular: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(find0.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated, r.Sequence)), Is.EqualTo(new[]
            {
                (6, 20, 1, true, "GMKP"),
                (13, 27, 2, true, "NPKG"),
                (2, 44, 3, false, "LRDETLRDETLRDEX"), // divergence: getorf omits the frame-3 full-length record
            }));
            Assert.That(find0.Select(r => r.ToGetOrfRange(Circle.Length)), Is.EqualTo(new[] { (7, 18), (14, 25), (3, 45) }));
            Assert.That(find6.Select(r => (r.ToGetOrfRange(Circle.Length), r.CodonPosition, r.Sequence)), Is.EqualTo(new[]
            {
                ((17, 23), 18, "CCTAAGG"),
                ((24, 30), 25, "GATGAAA"),
            }));
        });
    }

    // M35 — a STOP-free circular frame (getorf "full length of a circular genome" block, WriteORF(f, seqlen-1)):
    // getorf -find 0 -minsize 3 -circular Y -reverse N on GCGCGCCCGGG (L = 11) prints [1 - 33] ARPGRARGAPG and
    // [2 - 33] RARGAPGARPG; with -maxsize 30 only the frame-3 record [3 - 33] APGARPGRARX (two dangling bases
    // would be read as b1 b2 N, one as X). Mode 2 prints the same proteins in getorf; here it reports the
    // frame's nucleotides (documented divergence).
    [Test]
    public void FindOrfs_Circular_StopFreeFrame_FullLengthRecord()
    {
        const string s = "GCGCGCCCGGG";
        var find0 = Translator.FindOrfs(s, OrfFindMode.TranslationStopToStop, minSize: 3, searchBothStrands: false,
            circular: true).ToList();
        var find2 = Translator.FindOrfs(s, OrfFindMode.NucleotidesStopToStop, minSize: 3, searchBothStrands: false,
            circular: true).ToList();
        string tripled = s + s + s;

        Assert.Multiple(() =>
        {
            Assert.That(find0.Select(r => (r.ToGetOrfRange(s.Length), r.Frame, r.IsTerminated, r.Sequence)), Is.EqualTo(new[]
            {
                ((1, 33), 1, false, "ARPGRARGAPG"),
                ((2, 33), 2, false, "RARGAPGARPG"),
                ((3, 33), 3, false, "APGARPGRARX"),
            }));
            Assert.That(Translator.FindOrfs(s, OrfFindMode.TranslationStopToStop, minSize: 3, maxSize: 30,
                    searchBothStrands: false, circular: true).Select(r => (r.ToGetOrfRange(s.Length), r.Sequence)),
                Is.Empty, "-maxsize 30 → 10 residues: all three 11-residue records exceed it here");
            Assert.That(find2.Select(r => r.Sequence), Is.EqualTo(new[] { tripled, tripled[1..], tripled[2..] }));
        });
    }

    // M36 — circular: the DnaSequence and string overloads agree; default (circular: false) is unchanged.
    [Test]
    public void FindOrfs_Circular_DnaAndStringOverloadsAgree()
    {
        const string acgt = "ATGAAATTTGGGCCCTAGGCATGCCCAAATTTGA";
        foreach (OrfFindMode mode in Enum.GetValues<OrfFindMode>())
        {
            var fromDna = Translator.FindOrfs(new DnaSequence(acgt), mode, minSize: 3, flanking: 3, circular: true).ToList();
            var fromString = Translator.FindOrfs(acgt, mode, minSize: 3, flanking: 3, circular: true).ToList();
            Assert.That(fromString, Is.EqualTo(fromDna), mode.ToString());
            Assert.That(Translator.FindOrfs(acgt, mode, minSize: 3, flanking: 3).ToList(),
                Is.EqualTo(Translator.FindOrfs(new DnaSequence(acgt), mode, minSize: 3, flanking: 3).ToList()), mode.ToString());
        }
    }

    #endregion

    #region getorf binary oracle (TSV): IUPAC sequences, linear and -circular

    // M37 — every C line of the oracle: getorf [s - e], strand and nucleotide output (upper-cased; getorf keeps
    // the input case and writes U as T) equal, in order; protein residues equal getorf's except where getorf
    // prints X for an ambiguous codon that the genetic code still resolves (e.g. RAY → B). The oracle sequences
    // avoid the codons whose EMBOSS ambiguous translation is inconsistent (see the file header) and the
    // circular full-length record; the linear one-codon open-ORF divergence is excluded as in M20.
    [TestCaseSource(nameof(OracleConfigIds))]
    public void FindOrfs_IupacAndCircular_MatchGetorfBinaryOracle(string configId)
    {
        var oracle = Oracle.Value;
        var c = oracle.Configs[configId];
        var code = GeneticCode.GetByTableNumber(c.Table);
        bool protein = c.Mode is OrfFindMode.TranslationStopToStop or OrfFindMode.TranslationStartToStop;

        Assert.Multiple(() =>
        {
            foreach (var (name, seq) in oracle.Sequences)
            {
                var actual = Translator.FindOrfs(seq, c.Mode, code, c.MinSize, c.MaxSize, c.Flanking, c.Reverse,
                        c.Methionine, c.Circular)
                    .Where(r => c.Circular || !IsDocumentedDivergence(r, seq.Length))
                    .ToList();
                var expected = oracle.Records.TryGetValue((configId, name), out var list) ? list : [];
                Assert.That(actual.Select(r => r.ToGetOrfRange(seq.Length)).Zip(actual, (g, r) => (g.Start, g.End, r.IsReverseStrand)),
                    Is.EqualTo(expected.Select(e => (e.Start, e.End, e.Reverse))), $"{configId} {name} {seq}");
                if (actual.Count != expected.Count)
                    continue;
                for (int i = 0; i < actual.Count; i++)
                {
                    string want = expected[i].Sequence.ToUpperInvariant();
                    string got = actual[i].Sequence;
                    bool same = protein
                        ? got.Length == want.Length && got.Zip(want).All(p => p.First == p.Second || p.Second == 'X')
                        : got == want;
                    Assert.That(same, Is.True, $"{configId} {name} #{i + 1}: getorf {want}, here {got}");
                }
            }
        });
    }

    // M38 — the oracle exercises what it is meant to: ambiguity codes, lower case, U, and origin-crossing ORFs.
    [Test]
    public void IupacCircularOracle_CoversAmbiguityAndOriginCrossing()
    {
        var oracle = Oracle.Value;
        var lengths = oracle.Sequences.ToDictionary(s => s.Name, s => s.Sequence.Length);
        var circularIds = oracle.Configs.Where(kv => kv.Value.Circular).Select(kv => kv.Key).ToHashSet();

        Assert.Multiple(() =>
        {
            Assert.That(oracle.Sequences.Count(s => s.Sequence.Any(ch => "RYSWKMBDHVN".Contains(char.ToUpperInvariant(ch)))),
                Is.GreaterThanOrEqualTo(20));
            Assert.That(oracle.Sequences.Any(s => s.Sequence.Any(char.IsLower)), Is.True);
            Assert.That(oracle.Sequences.Any(s => s.Sequence.Contains('U')), Is.True);
            Assert.That(circularIds, Has.Count.EqualTo(43));
            Assert.That(oracle.Records.Where(kv => circularIds.Contains(kv.Key.Item1))
                .Sum(kv => kv.Value.Count(r => Math.Max(r.Start, r.End) > lengths[kv.Key.Item2])), Is.GreaterThan(100));
        });
    }

    #endregion

    #region Oracle loading

    // A one-codon open ORF opened at the last complete codon of its frame (linear; documented divergence).
    private static bool IsDocumentedDivergence(OrfFeature r, int length) =>
        r.Mode is OrfFindMode.TranslationStartToStop or OrfFindMode.NucleotidesStartToStop or OrfFindMode.FlankingStart
        && !r.IsTerminated
        && (r.Mode == OrfFindMode.FlankingStart
            ? r.CodonPosition + 6 > length
            : r.EndPosition - r.StartPosition == 2 && r.StartPosition + 6 > length);

    private sealed record OracleConfig(OrfFindMode Mode, int Table, int MinSize, int MaxSize, int Flanking,
        bool Reverse, bool Methionine, bool Circular);

    private sealed record OracleData(
        List<(string Name, string Sequence)> Sequences,
        Dictionary<string, OracleConfig> Configs,
        Dictionary<(string, string), List<(int Start, int End, bool Reverse, string Sequence)>> Records);

    private static readonly Lazy<OracleData> Oracle = new(LoadOracle);

    private static IEnumerable<string> OracleConfigIds() => Oracle.Value.Configs.Keys;

    private static OracleData LoadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Getorf.getorf_iupac_circular.tsv";
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
                        int.Parse(f[5]), int.Parse(f[6]), f[7] == "Y", f[8] == "Y", f[9] == "Y");
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
