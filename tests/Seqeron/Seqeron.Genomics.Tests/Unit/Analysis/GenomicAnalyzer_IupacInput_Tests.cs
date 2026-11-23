// GENOMIC-TANDEM-001 / GENOMIC-REPEAT-001 / GENOMIC-ORF-001 (B09 A1-21) — string overloads of GenomicAnalyzer.FindTandemRepeats,
// FindPalindromes and FindOpenReadingFrames for sequences with N / IUPAC symbols (DnaSequence is A/C/G/T only).
// Evidence: docs/Evidence/GENOMIC-TANDEM-001-Evidence.md, GENOMIC-REPEAT-001, GENOMIC-ORF-001; TestSpecs of the three units.
// Oracles (the existing generators in --ambiguous mode, scripts/analysis/):
//   tandem_iupac_oracle.tsv     892 rows: three independent methods (brute-force runs, boolean runs, regex) + the MISA [acgt] unit rule
//   palindrome_iupac_oracle.tsv 861 rows: Biopython reverse_complement on windows over A/C/G/T only (REVP / einverted alphabet)
//   orf_ambiguous_oracle.tsv    2844 rows: Biopython definition (ambiguous stop = all concrete codons stop; START unambiguous) and
//                               EMBOSS getorf 6.6.0 for the 184 A/C/G/T/N rows it resolves exactly

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_IupacInput_Tests
{
    private static IEnumerable<string[]> ReadRows(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Analysis." + resource;
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length > 0) yield return line.Split('\t');
        }
    }

    #region Tandem repeats

    private static string Format(IEnumerable<TandemRepeat> items)
    {
        string s = string.Join(";", items.Select(t => $"{t.Unit}:{t.Position}:{t.Repetitions}"));
        return s.Length == 0 ? "-" : s;
    }

    [Test]
    public void FindTandemRepeats_String_MatchesAmbiguousOracle()
    {
        var rows = ReadRows("tandem_iupac_oracle.tsv").ToList();
        Assert.That(rows, Has.Count.EqualTo(892));
        int nonEmpty = 0;
        foreach (var f in rows)
        {
            int maxUnit = int.Parse(f[3]);
            var actual = (maxUnit == 0
                ? GenomicAnalyzer.FindTandemRepeats(f[0], int.Parse(f[1]), int.Parse(f[2]))
                : GenomicAnalyzer.FindTandemRepeats(f[0], int.Parse(f[1]), int.Parse(f[2]), maxUnit)).ToList();
            Assert.That(actual, Has.Count.EqualTo(int.Parse(f[4])), $"{f[0]} {f[1]} {f[2]} {f[3]}");
            Assert.That(Format(actual), Is.EqualTo(f[5]), $"{f[0]} {f[1]} {f[2]} {f[3]}");
            if (actual.Count > 0) nonEmpty++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(200));
    }

    [Test]
    public void FindTandemRepeats_String_WorkedExamples_AmbiguityNeverFormsOrBridgesARepeat()
    {
        // N next to a repeat does not disturb it; the N run itself and units containing N are not repeats.
        Assert.That(Format(GenomicAnalyzer.FindTandemRepeats("ACACNACAC", 1, 2)), Is.EqualTo("AC:0:2;AC:5:2"));
        Assert.That(GenomicAnalyzer.FindTandemRepeats("NNNNNNNN", 1, 2), Is.Empty);
        Assert.That(GenomicAnalyzer.FindTandemRepeats("ACNACNACN", 1, 2), Is.Empty);
        Assert.That(GenomicAnalyzer.FindTandemRepeats("ACGNACGNACGN", 1, 2), Is.Empty);
        // Case-insensitive; RNA U is not A/C/G/T (canonical microsatellite alphabet), so no repeat.
        Assert.That(Format(GenomicAnalyzer.FindTandemRepeats("acacacnGt", 2, 2)), Is.EqualTo("AC:0:3"));
        Assert.That(GenomicAnalyzer.FindTandemRepeats("AUAUAUAU", 2, 2), Is.Empty);
    }

    [Test]
    public void FindTandemRepeats_String_AcgtInput_EqualsTheDnaSequenceOverload_BothEngines()
    {
        var rnd = new Random(2026);
        foreach (int n in new[] { 40, 150, 700 }) // 700 -> run-enumeration engine (more than 256 admissible unit lengths)
        {
            for (int trial = 0; trial < 15; trial++)
            {
                string s = new string(Enumerable.Range(0, n).Select(_ => "ACGTAC"[rnd.Next(6)]).ToArray());
                s = s.Insert(rnd.Next(n), string.Concat(Enumerable.Repeat("ACG", 4)));
                var a = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(s), 1, 2).Select(t => (t.Unit, t.Position, t.Repetitions));
                var b = GenomicAnalyzer.FindTandemRepeats(s, 1, 2).Select(t => (t.Unit, t.Position, t.Repetitions));
                Assert.That(b, Is.EqualTo(a));
            }
        }
    }

    [Test]
    public void FindTandemRepeats_String_WideUnitRange_RunEngineEqualsTheCanonicalPerUnitScan_WithAmbiguity()
    {
        // n = 700 with an unbounded range: > 256 admissible unit lengths, so the run engine answers; the canonical engine
        // (RepeatFinder.FindMicrosatellites(string)) scans every unit length. They must agree on strings with N / IUPAC symbols.
        var rnd = new Random(77);
        for (int trial = 0; trial < 12; trial++)
        {
            var chars = Enumerable.Range(0, 700).Select(_ => "ACGTACGTNNRY"[rnd.Next(12)]).ToArray();
            string s = new string(chars);
            s = s.Insert(rnd.Next(600), string.Concat(Enumerable.Repeat("ACGTN", 3)))
                 .Insert(rnd.Next(600), string.Concat(Enumerable.Repeat("TTAGC", 4)));
            var viaFacade = GenomicAnalyzer.FindTandemRepeats(s, 1, 2).Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
            var canonical = RepeatFinder.FindMicrosatellites(s, 1, s.Length / 2, 2)
                .Select(m => (m.RepeatUnit, m.Position, m.RepeatCount)).ToList();
            Assert.That(viaFacade, Is.EqualTo(canonical));
            Assert.That(viaFacade, Has.Count.GreaterThan(0));
        }
    }

    [Test]
    public void FindTandemRepeats_String_Contract()
    {
        Assert.That(() => GenomicAnalyzer.FindTandemRepeats((string)null!, 1, 2), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(GenomicAnalyzer.FindTandemRepeats("", 1, 2), Is.Empty);
        Assert.That(() => GenomicAnalyzer.FindTandemRepeats("ACAC", 0, 2), NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        Assert.That(() => GenomicAnalyzer.FindTandemRepeats("ACAC", 1, 1), NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        Assert.That(() => GenomicAnalyzer.FindTandemRepeats("ACAC", 3, 2, 2), NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
    }

    #endregion

    #region Palindromes

    [Test]
    public void FindPalindromes_String_MatchesAmbiguousOracle()
    {
        var rows = ReadRows("palindrome_iupac_oracle.tsv").ToList();
        Assert.That(rows, Has.Count.EqualTo(861));
        int withHits = 0;
        foreach (var f in rows)
        {
            var expected = f[3] == "-"
                ? Array.Empty<(int, int)>()
                : f[3].Split(',').Select(h => { var a = h.Split(':'); return (int.Parse(a[0]), int.Parse(a[1])); }).ToArray();
            var actual = GenomicAnalyzer.FindPalindromes(f[0], int.Parse(f[1]), int.Parse(f[2])).Select(p => (p.Position, p.Length)).ToArray();
            Assert.That(actual, Is.EqualTo(expected), $"{f[0]} [{f[1]},{f[2]}]");
            if (actual.Length > 0) withHits++;
        }

        Assert.That(withHits, Is.GreaterThan(200));
    }

    [Test]
    public void FindPalindromes_String_WorkedExamples_AmbiguousWindowsNeverPair()
    {
        // EcoRI around an N: the two sites are found, windows touching N are not (NNNN / SSSS equal their symbolic reverse complement).
        var eco = GenomicAnalyzer.FindPalindromes("GAATTCNGAATTC", 6, 6).Select(p => (p.Position, p.Sequence)).ToArray();
        Assert.That(eco, Is.EqualTo(new[] { (0, "GAATTC"), (7, "GAATTC") }));
        Assert.That(GenomicAnalyzer.FindPalindromes("NNNN"), Is.Empty);
        Assert.That(GenomicAnalyzer.FindPalindromes("SSSS"), Is.Empty);
        Assert.That(GenomicAnalyzer.FindPalindromes("gaattc").Select(p => p.Sequence), Is.EqualTo(new[] { "GAATTC", "AATT" }).Or.EqualTo(new[] { "GAATTC", "AATT" }.Reverse()),
            "case-insensitive; windows are reported in position order, the 4-mer AATT starts at 1");
    }

    [Test]
    public void FindPalindromes_String_AcgtInput_EqualsTheDnaSequenceOverload_AndNullContract()
    {
        var rnd = new Random(5);
        for (int t = 0; t < 40; t++)
        {
            string s = new string(Enumerable.Range(0, rnd.Next(0, 120)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var a = GenomicAnalyzer.FindPalindromes(new DnaSequence(s), 4, 10).Select(p => (p.Position, p.Sequence));
            var b = GenomicAnalyzer.FindPalindromes(s, 4, 10).Select(p => (p.Position, p.Sequence));
            Assert.That(b, Is.EqualTo(a));
        }

        Assert.That(() => GenomicAnalyzer.FindPalindromes((string)null!), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(GenomicAnalyzer.FindPalindromes(""), Is.Empty);
        Assert.That(() => GenomicAnalyzer.FindPalindromes("ACGT", 5, 8), NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
    }

    #endregion

    #region ORFs

    private sealed record OrfRow(int Table, string StartMode, bool Nested, bool Open, bool ExcludeStop, int MinLength, bool Both,
        string Sequence, string Orfs);

    private static OrfSearchOptions OptionsOf(OrfRow r) => new()
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
    public void FindOpenReadingFrames_String_MatchesAmbiguousOracle()
    {
        var rows = ReadRows("orf_ambiguous_oracle.tsv")
            .Select(f => new OrfRow(int.Parse(f[0]), f[1], f[2] == "1", f[3] == "1", f[4] == "1", int.Parse(f[5]), f[6] == "1",
                f[7] == "-" ? "" : f[7], f[8])).ToList();
        Assert.That(rows, Has.Count.EqualTo(2844));

        var mismatches = new List<string>();
        foreach (var row in rows)
        {
            string actual = Format(GenomicAnalyzer.FindOpenReadingFrames(row.Sequence, row.MinLength, OptionsOf(row)), row.Sequence.Length);
            if (actual != row.Orfs)
                mismatches.Add($"[table {row.Table} {row.StartMode} nested={row.Nested} open={row.Open} exclStop={row.ExcludeStop} " +
                               $"min={row.MinLength} both={row.Both}] {row.Sequence}\n  expected {row.Orfs}\n  actual   {actual}");
        }

        Assert.That(mismatches, Is.Empty, string.Join("\n", mismatches.Take(5)));
        Assert.That(rows.Count(r => r.Orfs != "-"), Is.GreaterThan(2000));
    }

    [Test]
    public void FindOpenReadingFrames_String_WorkedExamples()
    {
        // NNN is a sense codon the ORF reads through; TAR (TAA/TAG) is a stop; an N codon never starts an ORF.
        var through = GenomicAnalyzer.FindOpenReadingFrames("ATGNNNTAA", 1).Single();
        Assert.That((through.Position, through.Length, through.Sequence), Is.EqualTo((0, 9, "ATGNNNTAA")));
        var ends = GenomicAnalyzer.FindOpenReadingFrames("ATGTARCCCTAA", 1, new OrfSearchOptions { SearchBothStrands = false }).Single();
        Assert.That(ends.Sequence, Is.EqualTo("ATGTAR"), "TAR = TAA|TAG, every concrete codon is a stop");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames("NNNATGTAA", 1, new OrfSearchOptions { SearchBothStrands = false }).Single().Position, Is.EqualTo(3));
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames("NTGCCCTAA", 1, new OrfSearchOptions { SearchBothStrands = false }), Is.Empty,
            "an ambiguous codon is never a start");
        // Lower case and RNA are read as upper-case DNA (U -> T); the reverse strand is the IUPAC reverse complement.
        var rna = GenomicAnalyzer.FindOpenReadingFrames("augcccuaa", 1, new OrfSearchOptions { SearchBothStrands = false }).Single();
        Assert.That(rna.Sequence, Is.EqualTo("ATGCCCTAA"));
        var rev = GenomicAnalyzer.FindOpenReadingFrames("TTANNNCAT", 1).Single(o => o.IsReverseComplement);
        Assert.That(rev.Sequence, Is.EqualTo("ATGNNNTAA"));
        Assert.That(rev.GetForwardStrandRange(9), Is.EqualTo((0, 8)));
    }

    [Test]
    public void FindOpenReadingFrames_String_AcgtInput_EqualsTheDnaSequenceOverload_AndContracts()
    {
        var rnd = new Random(31);
        for (int t = 0; t < 60; t++)
        {
            string s = new string(Enumerable.Range(0, rnd.Next(0, 300)).Select(_ => "ACGTATGTAA"[rnd.Next(10)]).ToArray());
            var opts = new OrfSearchOptions { IncludeOpenEnded = t % 2 == 0, NestedOrfs = t % 3 != 0, StartCodons = (OrfStartCodons)(t % 3) };
            string a = Format(GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(s), 6, opts), s.Length);
            string b = Format(GenomicAnalyzer.FindOpenReadingFrames(s, 6, opts), s.Length);
            Assert.That(b, Is.EqualTo(a), s);
            Assert.That(Format(GenomicAnalyzer.FindOpenReadingFrames(s, 6), s.Length),
                Is.EqualTo(Format(GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(s), 6), s.Length)));
        }

        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames((string)null!, 1, OrfSearchOptions.Default), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames("ATG", 1, null!), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames("", 1), Is.Empty);
        var bad = Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindOpenReadingFrames("ATGXTAA", 1));
        Assert.That(bad!.Message, Does.Contain("'X'").And.Contain("position 3"));
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames("ATG-TAA", 1), NUnit.Framework.Throws.ArgumentException);
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames("ATGTAA", 1, new OrfSearchOptions { GeneticCode = GeneticCode.GetByTableNumber(27) }),
            NUnit.Framework.Throws.ArgumentException, "dual-coding stop tables are rejected eagerly");
    }

    #endregion
}
