// GENOMIC-ORF-001 (B09 A1-18 / A1-19) - FindOpenReadingFrames against the REAL NCBI ORFfinder 0.4.3 standalone program.
// Evidence: docs/Evidence/GENOMIC-ORF-001-Evidence.md  TestSpec: tests/TestSpecs/GENOMIC-ORF-001.md
// Oracle: scripts/analysis/generate_orffinder_oracle.py -> TestData/Analysis/orffinder_cli_oracle.tsv. The oracle is the unmodified NCBI
//         binary (ORFfinder.gz sha256 f30a0806...09f2, `ORFfinder -version` = 0.4.3, Debian 12 x86_64) run by the repository owner
//         in Docker over 629 deterministic inputs / 1750 invocations (docs/Validation/review-2026-09/handoff/A1-18-orffinder/,
//         returned archive returned/A1-18-orffinder-out.zip). Every row is one (option set, sequence) result: ORFs as
//         strand,lo,hi,flag (0-based inclusive forward coordinates, stop included, ORFfinder output order) and a hash of the proteins
//         exactly as ORFfinder printed them.

using System.Security.Cryptography;
using System.Text;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_FindOpenReadingFrames_OrfFinderCli_Tests
{
    private sealed record Config(string Id, int Table, int StartMode, int MinLength, bool IgnoreNested, string Strand, string Circular,
        int Begin, int End);

    private sealed record Row(Config Config, string SeqId, string Orfs, string ProteinHash);

    private sealed class Oracle
    {
        public Dictionary<string, string> Sequences { get; } = new();
        public Dictionary<string, Config> Configs { get; } = new();
        public List<Row> Rows { get; } = new();
        public Dictionary<int, string> GeneticCodeErrors { get; } = new();
    }

    private static readonly Lazy<Oracle> Data = new(ReadOracle);

    private static Oracle ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.orffinder_cli_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var oracle = new Oracle();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            switch (f[0])
            {
                case "S":
                    oracle.Sequences[f[1]] = f[2];
                    break;
                case "C":
                    oracle.Configs[f[1]] = new Config(f[1], int.Parse(f[2]), int.Parse(f[3]), int.Parse(f[4]), f[5] == "1", f[6], f[7],
                        int.Parse(f[8]), int.Parse(f[9]));
                    break;
                case "R":
                    oracle.Rows.Add(new Row(oracle.Configs[f[1]], f[2], f[3], f.Length > 4 ? f[4] : ""));
                    break;
                case "E":
                    oracle.GeneticCodeErrors[int.Parse(f[1])] = f[2];
                    break;
                default:
                    throw new InvalidOperationException("Unknown oracle line: " + line);
            }
        }

        return oracle;
    }

    // ORFfinder -g/-s/-ml/-n/-strand/-b/-e as OrfSearchOptions; -c has no counterpart (the binary ignores it).
    private static OrfSearchOptions OptionsOf(Config c) => OrfSearchOptions.OrfFinderCompatible with
    {
        GeneticCode = GeneticCode.GetByTableNumber(c.Table),
        StartCodons = c.StartMode switch { 0 => OrfStartCodons.AtgOnly, 1 => OrfStartCodons.GeneticCodeStarts, _ => OrfStartCodons.AnySenseCodon },
        IgnoreNestedOrfs = c.IgnoreNested,
        Strands = c.Strand switch { "plus" => OrfStrands.Plus, "minus" => OrfStrands.Minus, _ => OrfStrands.Both },
        RegionStart = c.Begin,
        RegionEnd = c.End,
    };

    private static string Flag(OrfInfo o) => (o.IsStartPartial, !o.IsTerminated) switch
    {
        (true, true) => "B",
        (true, false) => "5",
        (false, true) => "3",
        _ => "0",
    };

    private static bool FlagMatches(string expected, string actual) => expected switch
    {
        "?" => true,
        "P" => actual != "0",
        _ => expected == actual,
    };

    // ORFfinder protein: the table's translation of the codons without the stop; a 5'-complete ORF starts with Met whatever its first
    // codon (alternative start, or any sense codon for -s 2), a 5'-partial ORF keeps the translation of its first codon.
    private static string Protein(OrfInfo o, GeneticCode code)
    {
        int codons = o.Sequence.Length / 3 - (o.IsTerminated ? 1 : 0);
        var sb = new StringBuilder(codons);
        for (int i = 0; i < codons; i++)
            sb.Append(code.Translate(o.Sequence.Substring(3 * i, 3)));
        if (sb.Length > 0 && !o.IsStartPartial)
            sb[0] = 'M';
        return sb.ToString();
    }

    private static string Hash(IEnumerable<string> proteins)
        => Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes(string.Join("|", proteins)))).Substring(0, 8).ToLowerInvariant();

    [Test]
    public void FindOpenReadingFrames_EqualsNcbiOrfFinderBinary_OnEveryOracleRow()
    {
        var oracle = Data.Value;
        Assert.That(oracle.Rows, Has.Count.EqualTo(54820), "oracle row count");
        Assert.That(oracle.Configs, Has.Count.EqualTo(260), "oracle option sets");

        var mismatches = new List<string>();
        int orfs = 0, errorRows = 0, dnaChecked = 0;
        foreach (var row in oracle.Rows)
        {
            string sequence = oracle.Sequences[row.SeqId];
            var options = OptionsOf(row.Config);
            string what = $"[{row.Config.Id} g{row.Config.Table} s{row.Config.StartMode} ml{row.Config.MinLength} n{(row.Config.IgnoreNested ? 1 : 0)} " +
                          $"{row.Config.Strand} c{row.Config.Circular} b{row.Config.Begin} e{row.Config.End}] {row.SeqId}";

            if (row.Orfs.StartsWith("ERR:", StringComparison.Ordinal))
            {
                errorRows++;
                if (!ThrowsArgumentOutOfRange(() => GenomicAnalyzer.FindOpenReadingFrames(sequence, row.Config.MinLength, options)))
                    mismatches.Add($"{what}: ORFfinder rejected the region ({row.Orfs}) but no ArgumentOutOfRangeException");
                continue;
            }

            var actual = GenomicAnalyzer.FindOpenReadingFrames(sequence, row.Config.MinLength, options).ToList();
            if (++dnaChecked % 25 == 0 && sequence.All(ch => "ACGT".Contains(ch)) && sequence.Length > 0)
            {
                // the DnaSequence overload is the same engine
                var viaDna = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(sequence), row.Config.MinLength, options).ToList();
                if (!viaDna.Select(o => (o.Sequence, o.Position, o.IsReverseComplement, o.IsTerminated)).SequenceEqual(
                        actual.Select(o => (o.Sequence, o.Position, o.IsReverseComplement, o.IsTerminated))))
                    mismatches.Add($"{what}: DnaSequence overload differs from the string overload");
            }

            var expected = row.Orfs == "-" ? [] : row.Orfs.Split(';').Select(x => x.Split(',')).ToList();
            orfs += expected.Count;

            string Describe(IEnumerable<string> x) => string.Join(";", x);
            string[] actualSpans = actual.Select(o =>
            {
                var (s, e) = o.GetForwardStrandRange();
                return $"{(o.IsReverseComplement ? '-' : '+')},{s},{e}";
            }).ToArray();
            string[] expectedSpans = expected.Select(x => $"{x[0]},{x[1]},{x[2]}").ToArray();
            if (!actualSpans.SequenceEqual(expectedSpans))
            {
                mismatches.Add($"{what}\n  expected {Describe(expectedSpans)}\n  actual   {Describe(actualSpans)}");
                continue;
            }

            for (int i = 0; i < actual.Count; i++)
            {
                if (!FlagMatches(expected[i][3], Flag(actual[i])))
                    mismatches.Add($"{what}: ORF {expectedSpans[i]} partial flag expected {expected[i][3]} actual {Flag(actual[i])}");
            }

            if (row.ProteinHash != "-" && Hash(actual.Select(o => Protein(o, options.GeneticCode!))) != row.ProteinHash)
                mismatches.Add($"{what}: protein sequences differ ({Describe(expectedSpans)})");
        }

                Assert.That(orfs, Is.EqualTo(232_162), "ORFs compared (plus 360 rows of rejected regions)");
        Assert.That(errorRows, Is.EqualTo(360));
        Assert.That(mismatches, Is.Empty, mismatches.Count + " mismatching rows:\n" + string.Join("\n", mismatches.Take(5)));
    }

    private static bool ThrowsArgumentOutOfRange(Action action)
    {
        try
        {
            action();
            return false;
        }
        catch (ArgumentOutOfRangeException)
        {
            return true;
        }
    }

    // Every option dimension of the binary is present in the oracle: the matrix of the hand-off package, as measured.
    [Test]
    public void Oracle_CoversEveryOptionDimension()
    {
        var oracle = Data.Value;
        var configs = oracle.Configs.Values.ToList();
        Assert.Multiple(() =>
        {
            Assert.That(configs.Select(c => c.StartMode).Distinct().Order(), Is.EqualTo(new[] { 0, 1, 2 }), "-s");
            Assert.That(configs.Select(c => c.Strand).Distinct().Order(), Is.EqualTo(new[] { "both", "minus", "plus" }), "-strand");
            Assert.That(configs.Select(c => c.IgnoreNested).Distinct().Count(), Is.EqualTo(2), "-n");
            Assert.That(configs.Select(c => c.Circular).Distinct().Order(), Is.EqualTo(new[] { "0", "01", "1" }), "-c");
            Assert.That(configs.Select(c => c.MinLength).Distinct(), Is.SupersetOf(new[] { -5, 0, 1, 10, 29, 30, 31, 74, 75, 76, 1000 }), "-ml");
            Assert.That(configs.Where(c => c.Begin != 1 || c.End != 0).Select(c => (c.Begin, c.End)).Distinct().Count(), Is.GreaterThanOrEqualTo(20), "-b/-e");
            // -g: ORFfinder knows exactly the NCBI gc.prt tables of GeneticCode.SupportedTableNumbers; 0, 7, 8, 17-20 are rejected.
            Assert.That(configs.Select(c => c.Table).Distinct().Order(), Is.EqualTo(GeneticCode.SupportedTableNumbers.Order()), "-g");
            Assert.That(oracle.GeneticCodeErrors.Keys.Order(), Is.EqualTo(new[] { 0, 7, 8, 17, 18, 19, 20 }));
            Assert.That(oracle.GeneticCodeErrors.Values, Has.All.StartWith("Unable to find genetic code number"));
            foreach (int bad in oracle.GeneticCodeErrors.Keys)
                Assert.That(GeneticCode.SupportedTableNumbers, Does.Not.Contain(bad), $"table {bad}");
        });
    }

    // -c true: the binary's output is byte-identical to -c false on every circular genome (36 pairs), region and no-stop genome
    // (USAGE.txt: "Under development"); the oracle merges the two as "01". Therefore OrfSearchOptions has no circular option.
    [Test]
    public void CircularOption_HasNoObservableEffectInTheBinary()
    {
        var oracle = Data.Value;
        var both = oracle.Configs.Values.Where(c => c.Circular == "01").ToList();
        Assert.That(both, Has.Count.GreaterThanOrEqualTo(20), "option sets measured with -c true and -c false");
        var onlyTrue = oracle.Configs.Values.Where(c => c.Circular == "1").ToList();
        Assert.That(onlyTrue.Count + both.Count, Is.GreaterThan(0));
        // Wrap-around ORFs (cut by the origin) do not exist in the output: no oracle ORF of a circular genome crosses the origin, i.e. every
        // ORF is a plain forward-coordinate span lo <= hi inside the sequence.
        foreach (var row in oracle.Rows.Where(r => r.SeqId.StartsWith("circular:", StringComparison.Ordinal) && r.Config.Circular != "0" && r.Orfs != "-"))
        {
            int length = oracle.Sequences[row.SeqId].Length;
            foreach (var x in row.Orfs.Split(';').Select(x => x.Split(',')))
                Assert.That(int.Parse(x[1]) <= int.Parse(x[2]) && int.Parse(x[2]) < length, Is.True, row.SeqId);
        }
    }

    #region Worked examples (values read from the oracle)

    private static string Spans(IEnumerable<OrfInfo> orfs) => string.Join(";", orfs.Select(o =>
    {
        var (s, e) = o.GetForwardStrandRange();
        return $"{(o.IsReverseComplement ? '-' : '+')},{s},{e}";
    }));

    private static OrfSearchOptions Ncbi => OrfSearchOptions.OrfFinderCompatible;

    // designed sequence strand_equal_span_s2_12 (48 nt), -s 2 -ml 30: the plus ORF 7..45 and the minus ORF 7..45 have identical spans;
    // with -n true ORFfinder keeps the minus-strand ORF and drops the plus-strand one (same for 3..44).
    [Test]
    public void IgnoreNestedOrfs_IdenticalSpans_KeepTheMinusStrandOrf()
    {
        const string seq = "TAATTAGCTAAAGCTGTCCACCATGCTGTCCTTCCAGGTAGCTAATTA";
        var any = Ncbi with { StartCodons = OrfStartCodons.AnySenseCodon };

        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, any)), Is.EqualTo("+,3,44;+,7,45;-,3,44;-,2,40;-,7,45"));
        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, any with { IgnoreNestedOrfs = true })), Is.EqualTo("-,3,44;-,2,40;-,7,45"));
    }

    // designed sequence xnest_edge_a4_e4 (120 nt), defaults of -s 1 -ml 30: the ORF 4..108 (frame 2) lies inside the ORF 0..113
    // (frame 1) and the minus-strand ORFs 2..67 and 40..93 lie inside it too: -n true keeps the two outer plus-strand ORFs.
    [Test]
    public void IgnoreNestedOrfs_DropsOrfsOfOtherFramesAndStrands()
    {
        const string seq = "ATGGATGGAGACGAATCTTCTCGGAGGAGACGGCCCAGCATCAACCCGGAACCGCCGTTCTCAAGCAAGCCCGGGCGCTGCTACAACACTCCAAACCGAAAAATCCTAAATTAACCCCCC";
        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi)), Is.EqualTo("+,0,113;+,4,108;+,77,118;-,2,67;-,40,93"));
        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { IgnoreNestedOrfs = true })), Is.EqualTo("+,0,113;+,77,118"));
        // -strand plus: the minus-strand ORFs are not part of the result, so the plus ORF 4..108 is still nested in 0..113.
        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { IgnoreNestedOrfs = true, Strands = OrfStrands.Plus })),
            Is.EqualTo("+,0,113;+,77,118"));
    }

    // The container must pass the length filter itself and be a reported ORF of the visible strands: with -strand minus a plus ORF
    // cannot nest anything.
    [Test]
    public void IgnoreNestedOrfs_OnlyOrfsOfTheSelectedStrandsContain()
    {
        const string seq = "ATGGATGGAGACGAATCTTCTCGGAGGAGACGGCCCAGCATCAACCCGGAACCGCCGTTCTCAAGCAAGCCCGGGCGCTGCTACAACACTCCAAACCGAAAAATCCTAAATTAACCCCCC";
        var minusOnly = Ncbi with { Strands = OrfStrands.Minus };
        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, minusOnly)), Is.EqualTo("-,2,67;-,40,93"));
        Assert.That(Spans(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, minusOnly with { IgnoreNestedOrfs = true })), Is.EqualTo("-,2,67;-,40,93"));
    }

    // NCBI USAGE: "-ml: Value less than 30 is automatically changed by 30"; measured: the length is used in whole codons.
    // ORF = ATG + 9 codons + TAA: 30 coding nt (33 incl. stop).
    [TestCase(-5, true)]
    [TestCase(0, true)]
    [TestCase(29, true)]
    [TestCase(30, true)]
    [TestCase(31, true)]   // 31 -> 30 whole codons: an ORF of 30 coding nt is reported
    [TestCase(32, true)]
    [TestCase(33, false)]
    public void MinimumLength_ClampsToThirtyAndCountsWholeCodons(int ml, bool reported)
    {
        string seq = "ATG" + string.Concat(Enumerable.Repeat("GCC", 9)) + "TAA";
        var found = GenomicAnalyzer.FindOpenReadingFrames(seq, ml, Ncbi with { Strands = OrfStrands.Plus });
        Assert.That(found.Count(), Is.EqualTo(reported ? 1 : 0), $"-ml {ml}");
    }

    [Test]
    public void MinimumLength_BelowThirty_IsRaisedToThirty_ShorterOrfsAreNeverReported()
    {
        string orf27 = "ATG" + string.Concat(Enumerable.Repeat("GCC", 8)) + "TAA"; // 27 coding nt
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(orf27, 1, Ncbi with { Strands = OrfStrands.Plus }), Is.Empty, "-ml 1 acts as -ml 30");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(orf27, 1, Ncbi with { Strands = OrfStrands.Plus, OrfFinderMinimumLength = false }).Count(), Is.EqualTo(1),
            "without the ORFfinder rule -ml 1 admits it");
    }

    // 5'-partial = begins within the first 3 nt of the scanned strand / region (also for an ORF that starts with ATG); 3'-partial = no stop.
    [Test]
    public void StartPartial_IsStartWithinTheFirstThreeNucleotides()
    {
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + "TAA";
        foreach (int lead in new[] { 0, 1, 2, 3, 4 })
        {
            var o = GenomicAnalyzer.FindOpenReadingFrames(new string('C', lead) + orf, 30, Ncbi with { Strands = OrfStrands.Plus }).Single();
            Assert.That(o.IsStartPartial, Is.EqualTo(lead < 3), $"lead {lead}");
            Assert.That(o.IsTerminated, Is.True);
        }

        var open = GenomicAnalyzer.FindOpenReadingFrames("CCC" + orf[..^3], 30, Ncbi with { Strands = OrfStrands.Plus }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(open.IsTerminated, Is.False, "no stop: 3'-partial");
            Assert.That(open.IsStartPartial, Is.False);
        });
    }

    // -b/-e: 1-based inclusive, the region is cut out and scanned alone. region.fa-style: an ORF crossing the region end is open there.
    [Test]
    public void Region_IsCutOutThenScanned_CoordinatesStayInInputSpace()
    {
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + "TAA"; // 42 nt
        string seq = new string('C', 20) + orf + new string('C', 20);              // ORF at 20..61

        var whole = GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { Strands = OrfStrands.Plus }).Single();
        Assert.That(whole.GetForwardStrandRange(), Is.EqualTo((20, 61)));

        // -b 21 -e 62 is exactly the ORF: complete, but starts at region offset 0 -> 5'-partial (ORFfinder rule).
        var exact = GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { Strands = OrfStrands.Plus, RegionStart = 21, RegionEnd = 62 }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(exact.GetForwardStrandRange(), Is.EqualTo((20, 61)));
            Assert.That(exact.Position, Is.EqualTo(0), "Position is relative to the region");
            Assert.That(exact.RegionStart, Is.EqualTo(20));
            Assert.That(exact.IsStartPartial, Is.True);
            Assert.That(exact.IsTerminated, Is.True);
        });

        // A region ending inside the ORF: the ORF is open (3'-partial) at the cut, ending at the last complete codon.
        var cut = GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { Strands = OrfStrands.Plus, RegionEnd = 50 }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(cut.IsTerminated, Is.False);
            Assert.That(cut.GetForwardStrandRange(), Is.EqualTo((20, 49)));
        });

        // A region starting inside the ORF removes the start codon: nothing left to report.
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { Strands = OrfStrands.Plus, RegionStart = 24 }), Is.Empty);
    }

    // Reverse-strand frames are anchored at the region end: the reverse complement of the region is scanned.
    [Test]
    public void Region_ReverseStrandCoordinatesMapBackToTheInput()
    {
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + "TAA";
        string rc = DnaSequence.GetReverseComplementString(orf);
        string seq = new string('C', 25) + rc + new string('C', 25);               // minus ORF at forward 25..66
        var o = GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { Strands = OrfStrands.Minus, RegionStart = 11, RegionEnd = 80 }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(o.IsReverseComplement, Is.True);
            Assert.That(o.GetForwardStrandRange(), Is.EqualTo((25, 66)));
            Assert.That(o.GetForwardStrandRange(seq.Length), Is.EqualTo((25, 66)));
            Assert.That(o.Sequence, Is.EqualTo(orf));
            Assert.That(o.Position, Is.EqualTo(80 - 66 - 1), "position in the reverse complement of the region 11..80");
        });
    }

    [Test]
    public void Region_Conventions_MatchTheBinary()
    {
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + "TAA";
        string seq = new string('C', 20) + orf + new string('C', 20);
        var plus = Ncbi with { Strands = OrfStrands.Plus };
        int whole = GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus).Count();

        // -b 0 is -b 1; -e 0 and -e beyond the sequence mean "to the end" (measured: b0_e0, b0_e100, b1_e301, b1_e1000).
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionStart = 0 }).Count(), Is.EqualTo(whole));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionEnd = seq.Length + 1 }).Count(), Is.EqualTo(whole));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionEnd = 1000 }).Count(), Is.EqualTo(whole));
        // ORFfinder errors per sequence ("Interval Start value is expected to be smaller than Stop", "Requested interval start -5 is out of range").
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionStart = 50, RegionEnd = 40 }),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionStart = -5, RegionEnd = 100 }),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionEnd = -1 }),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        // A one-nucleotide region (-b 5 -e 5) is valid and has no ORF.
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(seq, 30, plus with { RegionStart = 5, RegionEnd = 5 }), Is.Empty);
    }

    // -g 27 / 28 / 31: the binary accepts them and the codons the table lists as stops end ORFs (measured on 488 designed rows each).
    [Test]
    public void DualCodingTables_AreAcceptedWithTheOrfFinderProfile()
    {
        foreach (int table in new[] { 27, 28, 31 })
        {
            var code = GeneticCode.GetByTableNumber(table);
            string stop = code.StopCodons.Order().First().Replace('U', 'T');
            string seq = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + stop;
            var o = GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { GeneticCode = code, Strands = OrfStrands.Plus }).Single();
            Assert.That(o.IsTerminated, Is.True, $"table {table}: {stop} ends the ORF");
            Assert.That(o.Sequence, Is.EqualTo(seq));
            Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(seq, 30, Ncbi with { GeneticCode = code, DualCodingStopsTerminate = false }),
                NUnit.Framework.Throws.TypeOf<ArgumentException>(), "default options still reject the dual-coding tables");
        }
    }

    // OrfFinderCompatible is the profile the oracle test uses; its members are the measured conventions.
    [Test]
    public void OrfFinderCompatible_ProfileMembers()
    {
        var p = OrfSearchOptions.OrfFinderCompatible;
        Assert.Multiple(() =>
        {
            Assert.That(p.StartCodons, Is.EqualTo(OrfStartCodons.GeneticCodeStarts), "-s 1");
            Assert.That(p.NestedOrfs, Is.False, "longest ORF per stop");
            Assert.That(p.IncludeOpenEnded, Is.True, "3'-partial ORFs");
            Assert.That(p.MinLengthExcludesStop, Is.True);
            Assert.That(p.OrfFinderMinimumLength, Is.True);
            Assert.That(p.DualCodingStopsTerminate, Is.True);
            Assert.That(p.IgnoreNestedOrfs, Is.False, "-n false");
            Assert.That(p.Strands, Is.EqualTo(OrfStrands.Both));
            Assert.That(p.RegionStart, Is.EqualTo(1));
            Assert.That(p.RegionEnd, Is.EqualTo(0));
        });
    }

    [Test]
    public void SearchBothStrands_IsAViewOfStrands()
    {
        Assert.That(new OrfSearchOptions { SearchBothStrands = false }.Strands, Is.EqualTo(OrfStrands.Plus));
        Assert.That(new OrfSearchOptions { SearchBothStrands = true }.Strands, Is.EqualTo(OrfStrands.Both));
        Assert.That((new OrfSearchOptions { Strands = OrfStrands.Minus }).SearchBothStrands, Is.False);
        Assert.That(OrfSearchOptions.Default.SearchBothStrands, Is.True);
        Assert.That(OrfSearchOptions.Default, Is.EqualTo(new OrfSearchOptions()), "record equality still value based");
    }

    [Test]
    public void FindNestedSpans_Rules()
    {
        var spans = new[]
        {
            new OrfScanner.OrfSpan(0, 99, false),   // container
            new OrfScanner.OrfSpan(0, 50, false),   // shares the start: nested (non-strict)
            new OrfScanner.OrfSpan(40, 99, true),   // shares the end: nested
            new OrfScanner.OrfSpan(10, 20, true),   // inside
            new OrfScanner.OrfSpan(90, 120, true),  // overlaps only: kept
            new OrfScanner.OrfSpan(200, 230, false),
            new OrfScanner.OrfSpan(200, 230, true), // identical to the previous one on the other strand: the plus one is nested
        };
        Assert.That(OrfScanner.FindNestedSpans(spans), Is.EqualTo(new[] { false, true, true, true, false, true, false }));
        Assert.That(OrfScanner.FindNestedSpans([]), Is.Empty);
    }

    [Test]
    public void OrfInfo_GetForwardStrandRange_WithoutRegionInformation_NeedsTheLength()
    {
        var legacy = new OrfInfo("ATGTAA", 2, 3, false);
        Assert.That(legacy.GetForwardStrandRange(10), Is.EqualTo((2, 7)));
        Assert.That(() => legacy.GetForwardStrandRange(), NUnit.Framework.Throws.TypeOf<InvalidOperationException>());
    }

    #endregion
}
