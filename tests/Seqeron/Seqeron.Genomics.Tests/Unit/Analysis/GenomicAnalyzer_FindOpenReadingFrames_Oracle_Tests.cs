// GENOMIC-ORF-001 (B09 F21-F26) — FindOpenReadingFrames against independent oracles.
// Evidence: docs/Evidence/GENOMIC-ORF-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-ORF-001.md
// Oracle: scripts/analysis/generate_orf_oracle.py (TestData/Analysis/orf_oracle.tsv), 1746 (sequence, options) rows, each asserted
//         equal by three independent implementations at generation time:
//           (1) definition-based Python finder (Biopython 1.88 translate / CodonTable / reverse_complement),
//           (2) line-by-line port of NCBI C++ Toolkit COrf::FindOrfs (the NCBI ORFfinder engine, orf.cpp),
//           (3) EMBOSS getorf 6.6.0 (-find 3, longest ORF per stop, alternative starts, open ends) for 377 rows.
//         Coordinates in the oracle are 0-based inclusive forward-strand spans including the stop codon.

using System.Diagnostics;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests
{
    private sealed record OracleRow(int Table, string StartMode, bool Nested, bool Open, bool ExcludeStop, int MinLength,
        bool Both, string Sequence, string Orfs);

    private static List<OracleRow> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.orf_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<OracleRow>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new OracleRow(int.Parse(f[0]), f[1], f[2] == "1", f[3] == "1", f[4] == "1", int.Parse(f[5]), f[6] == "1",
                f[7] == "-" ? "" : f[7], f[8]));
        }
        return rows;
    }

    private static OrfSearchOptions OptionsOf(OracleRow r) => new()
    {
        GeneticCode = GeneticCode.GetByTableNumber(r.Table),
        StartCodons = r.StartMode switch
        {
            "atg" => OrfStartCodons.AtgOnly,
            "alt" => OrfStartCodons.GeneticCodeStarts,
            _ => OrfStartCodons.AnySenseCodon,
        },
        NestedOrfs = r.Nested,
        IncludeOpenEnded = r.Open,
        MinLengthExcludesStop = r.ExcludeStop,
        SearchBothStrands = r.Both,
    };

    // strand,frame,start,end,terminated — forward-strand 0-based inclusive span incl. the stop codon.
    private static string Format(IEnumerable<OrfInfo> orfs, int length)
    {
        string s = string.Join(";", orfs.Select(o =>
        {
            var (start, end) = o.GetForwardStrandRange(length);
            return $"{(o.IsReverseComplement ? '-' : '+')},{o.Frame},{start},{end},{(o.IsTerminated ? 1 : 0)}";
        }));
        return s.Length == 0 ? "-" : s;
    }

    [Test]
    public void FindOpenReadingFrames_MatchesBiopythonNcbiCOrfAndGetorfOracle()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(1746), "oracle row count");

        var mismatches = new List<string>();
        foreach (var row in rows)
        {
            var dna = new DnaSequence(row.Sequence);
            string actual = Format(GenomicAnalyzer.FindOpenReadingFrames(dna, row.MinLength, OptionsOf(row)), dna.Length);
            if (actual != row.Orfs)
                mismatches.Add($"[table {row.Table} {row.StartMode} nested={row.Nested} open={row.Open} exclStop={row.ExcludeStop} " +
                               $"min={row.MinLength} both={row.Both}] {row.Sequence}\n  expected {row.Orfs}\n  actual   {actual}");
        }

        Assert.That(mismatches, Is.Empty, string.Join("\n", mismatches.Take(5)));
    }

    // Every oracle row also reproduces through the original two-argument API when its options are the defaults.
    [Test]
    public void FindOpenReadingFrames_DefaultOverload_EqualsDefaultOptions()
    {
        foreach (var row in ReadOracle().Where(r => r is { Table: 1, StartMode: "atg", Nested: true, Open: false, ExcludeStop: false, Both: true }))
        {
            var dna = new DnaSequence(row.Sequence);
            Assert.That(Format(GenomicAnalyzer.FindOpenReadingFrames(dna, row.MinLength), dna.Length), Is.EqualTo(row.Orfs), row.Sequence);
        }
    }

    // Reverse-strand coordinates: the forward-strand span is the one whose reverse complement is the reported ORF
    // (Biopython-independent definition, asserted on random sequences for every ORF).
    [Test]
    public void GetForwardStrandRange_ReverseComplementOfSpan_EqualsOrfSequence()
    {
        var rnd = new Random(20261118);
        for (int trial = 0; trial < 200; trial++)
        {
            string seq = new string(Enumerable.Range(0, rnd.Next(0, 400)).Select(_ => "ACGTATGTAA"[rnd.Next(10)]).ToArray());
            var dna = new DnaSequence(seq);
            foreach (var orf in GenomicAnalyzer.FindOpenReadingFrames(dna, 6, new OrfSearchOptions { IncludeOpenEnded = true }))
            {
                var (start, end) = orf.GetForwardStrandRange(seq.Length);
                string span = seq.Substring(start, end - start + 1);
                string expected = orf.IsReverseComplement ? DnaSequence.GetReverseComplementString(span) : span;
                Assert.That(expected, Is.EqualTo(orf.Sequence), $"{seq} @ {start}-{end}");
            }
        }
    }

    [Test]
    public void GetForwardStrandRange_WorkedExamples()
    {
        // TTAGGGGGGCAT: reverse complement ATGCCCCCCTAA is the ORF at RC position 0 -> forward span 0..11.
        var rev = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence("TTAGGGGGGCAT"), 1).Single();
        Assert.That(rev.GetForwardStrandRange(12), Is.EqualTo((0, 11)));

        // CC + TTAGGGCAT + GG: the reverse complement CC ATGCCCTAA GG holds the ORF at RC position 2 -> forward span 2..10.
        string seq = "CC" + "TTAGGGCAT" + "GG";
        var orf = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(seq), 9).Single(o => o.IsReverseComplement);
        Assert.Multiple(() =>
        {
            Assert.That(orf.Sequence, Is.EqualTo("ATGCCCTAA"));
            Assert.That(orf.Position, Is.EqualTo(2), "RC = CC ATGCCCTAA GG");
            Assert.That(orf.GetForwardStrandRange(seq.Length), Is.EqualTo((2, 10)));
            Assert.That(() => orf.GetForwardStrandRange(5), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    #region Stop codons come from the genetic code table

    // Biopython 1.88: table 2 (vertebrate mitochondrial) stops = TAA TAG AGA AGG (TGA = Trp); table 1 stops = TAA TAG TGA.
    [Test]
    public void FindOpenReadingFrames_TableTwo_TgaReadsThrough_AgaTerminates()
    {
        var table2 = new OrfSearchOptions { GeneticCode = GeneticCode.GetByTableNumber(2) };

        var tga = new DnaSequence("ATGTGAAAATAA");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(tga, 1).Single(o => !o.IsReverseComplement).Sequence, Is.EqualTo("ATGTGA"));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(tga, 1, table2).Single(o => !o.IsReverseComplement).Sequence,
            Is.EqualTo("ATGTGAAAATAA"), "TGA is Trp in table 2");

        var aga = new DnaSequence("ATGCCCAGAGGG");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(aga, 1), Is.Empty, "no stop in table 1");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(aga, 1, table2).Single().Sequence, Is.EqualTo("ATGCCCAGA"), "AGA is a stop in table 2");
    }

    [TestCase(27)]
    [TestCase(28)]
    [TestCase(31)]
    public void FindOpenReadingFrames_DualCodingStopTables_ThrowEagerly(int table)
    {
        var options = new OrfSearchOptions { GeneticCode = GeneticCode.GetByTableNumber(table) };
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence("ATGTAA"), 1, options),
            NUnit.Framework.Throws.TypeOf<ArgumentException>(), "tables 27/28/31: every stop codon also codes an amino acid");
    }

    [Test]
    public void FindOpenReadingFrames_NullArguments_Throw()
    {
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(null!, 1, OrfSearchOptions.Default), NUnit.Framework.Throws.TypeOf<ArgumentNullException>());
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence("ATG"), 1, null!), NUnit.Framework.Throws.TypeOf<ArgumentNullException>());
    }

    #endregion

    #region Options: starts, nesting, open ends, length convention

    [Test]
    public void StartCodons_AtgOnly_AlternativeStarts_AnySenseCodon()
    {
        var alt = new OrfSearchOptions { StartCodons = OrfStartCodons.GeneticCodeStarts };
        // Table 1 (NCBI gc.prt): starts TTG, CTG, ATG.
        foreach (var start in new[] { "TTG", "CTG" })
        {
            var dna = new DnaSequence(start + "AAATAA");
            Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 1).Where(o => !o.IsReverseComplement), Is.Empty, $"{start}: ATG only");
            Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 1, alt).Single(o => !o.IsReverseComplement).Sequence,
                Is.EqualTo(start + "AAATAA"));
        }

        // GTG is a start in table 11 (bacterial), not in table 1.
        var gtg = new DnaSequence("GTGAAATAA");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(gtg, 1, alt), Is.Empty);
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(gtg, 1,
            alt with { GeneticCode = GeneticCode.GetByTableNumber(11) }).Single().Sequence, Is.EqualTo("GTGAAATAA"));

        // Any sense codon: stop-to-stop regions; nested = every sense codon starts one.
        var any = new DnaSequence("CCCAAATAAGGG");
        var stopToStop = new OrfSearchOptions { StartCodons = OrfStartCodons.AnySenseCodon, NestedOrfs = false };
        Assert.That(Format(GenomicAnalyzer.FindOpenReadingFrames(any, 1, stopToStop), any.Length), Is.EqualTo("+,1,0,8,1"));
        Assert.That(Format(GenomicAnalyzer.FindOpenReadingFrames(any, 1, stopToStop with { NestedOrfs = true }), any.Length),
            Is.EqualTo("+,1,0,8,1;+,1,3,8,1"));
    }

    [Test]
    public void NestedOrfs_False_ReportsOnlyTheMostUpstreamStartPerStop()
    {
        var dna = new DnaSequence("ATGGGGATGCCCTAA");
        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, 1, new OrfSearchOptions { NestedOrfs = false }).ToList();
        Assert.That(orfs, Has.Count.EqualTo(1), "NCBI ORFfinder reports each ORF from the most upstream start");
        Assert.That(orfs[0].Position, Is.EqualTo(0));
        Assert.That(orfs[0].Sequence, Is.EqualTo("ATGGGGATGCCCTAA"));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 1).Count(o => !o.IsReverseComplement), Is.EqualTo(2), "Rosalind default: both");
    }

    [Test]
    public void IncludeOpenEnded_ReportsThreePrimePartialOrfsToTheLastCompleteCodon()
    {
        var dna = new DnaSequence("ATGAAACCCGG");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 1), Is.Empty, "no stop: Rosalind does not report it");

        var open = GenomicAnalyzer.FindOpenReadingFrames(dna, 1, new OrfSearchOptions { IncludeOpenEnded = true }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(open.Sequence, Is.EqualTo("ATGAAACCC"), "ends at the last complete codon, the trailing GG is dropped");
            Assert.That(open.IsTerminated, Is.False);
            Assert.That(open.Length % 3, Is.EqualTo(0));
            Assert.That(open.GetForwardStrandRange(dna.Length), Is.EqualTo((0, 8)));
        });

        // A terminated ORF reports IsTerminated = true.
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence("ATGTAA"), 1).Single(o => !o.IsReverseComplement).IsTerminated, Is.True);
    }

    // NCBI COrf: "length in base pairs not including the stop" (orf.cpp FindForwardOrfs); EMBOSS getorf -minsize likewise.
    [Test]
    public void MinLengthExcludesStop_FollowsNcbiOrfFinderConvention()
    {
        var dna = new DnaSequence("ATGAAATAA"); // 9 nt with stop, 6 nt without
        var ncbi = new OrfSearchOptions { MinLengthExcludesStop = true };

        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 9).Count(o => !o.IsReverseComplement), Is.EqualTo(1), "default: Length 9 >= 9");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 10).Count(o => !o.IsReverseComplement), Is.EqualTo(0));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 6, ncbi).Count(o => !o.IsReverseComplement), Is.EqualTo(1), "6 nt without stop >= 6");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 7, ncbi).Count(o => !o.IsReverseComplement), Is.EqualTo(0), "6 nt without stop < 7");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 9, ncbi).Count(o => !o.IsReverseComplement), Is.EqualTo(0));
    }

    [Test]
    public void SearchBothStrands_False_ReportsForwardOnly()
    {
        var dna = new DnaSequence("TTAGGGGGGCAT");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 1).Count(), Is.EqualTo(1));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 1, new OrfSearchOptions { SearchBothStrands = false }), Is.Empty);
    }

    #endregion

    #region Differential tests against the other two ORF engines

    // Translator.FindOrfs (B02, EMBOSS getorf -find 1 model): first START per STOP, starts of the genetic code, open ORFs reported,
    // length in amino acids (excluding the stop), strand coordinates, frames +-1..3, End = last base incl. the stop.
    [Test]
    public void Differential_TranslatorFindOrfs_EqualsOutermostStartsOfTheGeneticCode()
    {
        var rnd = new Random(7);
        var options = new OrfSearchOptions { StartCodons = OrfStartCodons.GeneticCodeStarts, NestedOrfs = false, IncludeOpenEnded = true, MinLengthExcludesStop = true };
        foreach (int table in new[] { 1, 2, 4, 11 })
        {
            var code = GeneticCode.GetByTableNumber(table);
            for (int trial = 0; trial < 60; trial++)
            {
                string seq = RandomCodingRichSequence(rnd, rnd.Next(0, 300));
                var dna = new DnaSequence(seq);
                int minAa = rnd.Next(1, 6);

                var expected = Translator.FindOrfs(dna, code, minAa, true)
                    .Select(o => (Reverse: o.Frame < 0, Frame: Math.Abs(o.Frame), o.StartPosition, o.EndPosition))
                    .OrderBy(t => t.Reverse).ThenBy(t => t.Frame).ThenBy(t => t.StartPosition).ToList();
                var actual = GenomicAnalyzer.FindOpenReadingFrames(dna, 3 * minAa, options with { GeneticCode = code })
                    .Select(o => (Reverse: o.IsReverseComplement, o.Frame, StartPosition: o.Position, EndPosition: o.Position + o.Length - 1))
                    .OrderBy(t => t.Reverse).ThenBy(t => t.Frame).ThenBy(t => t.StartPosition).ToList();

                Assert.That(actual, Is.EqualTo(expected), $"table {table}, min {minAa} aa: {seq}");
            }
        }
    }

    // GenomeAnnotator.FindOrfs (ATG/GTG/TTG starts, every start, terminated ORFs, aa length excluding the stop, forward coordinates,
    // half-open end): its ATG-initiated ORFs are exactly FindOpenReadingFrames' ORFs.
    [Test]
    public void Differential_GenomeAnnotatorFindOrfs_AtgInitiatedOrfsAreTheSameSet()
    {
        var rnd = new Random(11);
        for (int trial = 0; trial < 150; trial++)
        {
            string seq = RandomCodingRichSequence(rnd, rnd.Next(0, 360));
            var dna = new DnaSequence(seq);
            int minAa = rnd.Next(0, 5);

            var expected = GenomeAnnotator.FindOrfs(seq, minAa, true)
                .Where(o => o.Sequence.StartsWith("ATG", StringComparison.Ordinal))
                .Select(o => (o.IsReverseComplement, o.Start, o.End, o.Sequence)).OrderBy(t => t.IsReverseComplement).ThenBy(t => t.Start).ThenBy(t => t.End).ToList();
            var actual = GenomicAnalyzer.FindOpenReadingFrames(dna, 3 * minAa + 3)
                .Select(o =>
                {
                    var (s, e) = o.GetForwardStrandRange(seq.Length);
                    return (o.IsReverseComplement, Start: s, End: e + 1, o.Sequence);
                }).OrderBy(t => t.IsReverseComplement).ThenBy(t => t.Start).ThenBy(t => t.End).ToList();

            Assert.That(actual, Is.EqualTo(expected), $"min {minAa} aa: {seq}");
        }
    }

    private static string RandomCodingRichSequence(Random rnd, int length)
    {
        string[] pieces = { "ATG", "TTG", "CTG", "GTG", "ATA", "TAA", "TAG", "TGA", "AGA", "AAA", "GGC", "CCC" };
        var sb = new System.Text.StringBuilder();
        while (sb.Length < length)
        {
            if (rnd.NextDouble() < 0.6) sb.Append(pieces[rnd.Next(pieces.Length)]);
            else sb.Append("ACGT"[rnd.Next(4)]);
        }

        return sb.ToString(0, Math.Min(length, sb.Length));
    }

    #endregion

    #region Complexity

    // Before the single-pass scanner every ATG rescanned to the next stop: ATG repeated n times with one stop at the end was O(n^2).
    [Test]
    public void FindOpenReadingFrames_LongStopFreeAtgRun_IsLinear()
    {
        const int repeats = 300_000;
        var dna = new DnaSequence(string.Concat(Enumerable.Repeat("ATG", repeats)) + "TAA");

        var sw = Stopwatch.StartNew();
        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, 899_000).Where(o => !o.IsReverseComplement).ToList();
        sw.Stop();

        // The ORF of the k-th ATG (k = 0..299999) has length 900003 - 3k; >= 899000 gives 3k <= 1003, i.e. k = 0..334: 335 ORFs.
        Assert.That(orfs, Has.Count.EqualTo(335));
        Assert.That(orfs[0].Position, Is.EqualTo(0));
        Assert.That(orfs.All(o => o.Length >= 899_000), Is.True);
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(30), "O(n) scan; the previous per-ATG rescan took minutes");
    }

    [Test]
    public void FindOpenReadingFrames_FiveMegabaseRandomSequence_ScansAllSixFramesQuickly()
    {
        var rnd = new Random(3);
        var chars = new char[5_000_000];
        for (int i = 0; i < chars.Length; i++) chars[i] = "ACGT"[rnd.Next(4)];
        var dna = new DnaSequence(new string(chars));

        var sw = Stopwatch.StartNew();
        int count = GenomicAnalyzer.FindOpenReadingFrames(dna, 300).Count();
        sw.Stop();

        Assert.That(count, Is.GreaterThan(0));
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(60));
    }

    #endregion
}
