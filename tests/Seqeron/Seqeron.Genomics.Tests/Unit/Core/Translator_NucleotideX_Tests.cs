namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Opt-in Biopython-compatible nucleotide <c>X</c> in translation (review 2026-09, B02 F33).
/// Test Units: TRANS-PROT-001 (Translator.Translate string overload), TRANS-CODON-001
/// (GeneticCode.Translate / IsStartCodon / IsStopCodon overloads with allowNucleotideX).
///
/// Oracle: Biopython 1.88 (<c>Bio.Seq.translate</c>, <c>Bio.Data.CodonTable.ambiguous_dna_by_id</c>,
/// <c>IUPACData.ambiguous_dna_values['X'] == 'GATC'</c>), embedded as
/// TestData/GeneticCode/biopython_nucleotide_x.tsv (generator: gen_biopython_nucleotide_x.py next to it).
/// Biopython resolves an X-codon like the N-spelled codon, except that a possible-stop codon
/// (stop and amino-acid expansions) raises TranslationError because X is not in its ambiguous
/// nucleotide alphabet GATCRYWSMKHBVDN.
/// </summary>
[TestFixture]
public class Translator_NucleotideX_Tests
{
    private const string XAlphabet = "ACGTRYSWKMBDHVNX";

    private sealed record CodonRow(int TableId, string Translation, HashSet<string> Starts, HashSet<string> Stops);

    private sealed record RandomRow(int TableId, bool ToStop, string Input, string Expected);

    private static (List<CodonRow> Codons, List<RandomRow> Random) LoadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.GeneticCode.biopython_nucleotide_x.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var codons = new List<CodonRow>();
        var random = new List<RandomRow>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line.StartsWith('#')) continue;
            var f = line.Split('\t');
            if (f[0] == "C")
            {
                static HashSet<string> Set(string s) => s == "-" ? new HashSet<string>() : s.Split(',').ToHashSet();
                codons.Add(new CodonRow(int.Parse(f[1]), f[2], Set(f[3]), Set(f[4])));
            }
            else
            {
                random.Add(new RandomRow(int.Parse(f[1]), f[2] == "1", f[3], f[4]));
            }
        }
        return (codons, random);
    }

    private static string[] XCodons()
    {
        var list = new List<string>();
        foreach (char a in XAlphabet)
            foreach (char b in XAlphabet)
                foreach (char c in XAlphabet)
                {
                    var codon = new string(new[] { a, b, c });
                    if (codon.Contains('X')) list.Add(codon);
                }
        return list.ToArray();
    }

    /// <summary>
    /// All 721 X-containing codons over ACGTRYSWKMBDHVNX × 27 NCBI tables (19 467 codons), in DNA,
    /// RNA and lower-case spelling, translate exactly as Biopython 1.88; '!' (TranslationError) ↔
    /// ArgumentException. Reference counts: 12 407 resolve (= the N-spelled codon), 7 060 throw.
    /// </summary>
    [Test]
    public void Translate_AllXCodons_AllTables_MatchBiopython()
    {
        var codons = XCodons();
        Assert.That(codons, Has.Length.EqualTo(721));
        var oracle = LoadOracle().Codons;
        Assert.That(oracle.Select(r => r.TableId), Is.EqualTo(GeneticCode.SupportedTableNumbers));

        int resolved = 0, rejected = 0;
        foreach (var row in oracle)
        {
            var code = GeneticCode.GetByTableNumber(row.TableId);
            Assert.That(row.Translation, Has.Length.EqualTo(codons.Length));
            for (int i = 0; i < codons.Length; i++)
            {
                char expected = row.Translation[i];
                foreach (var spelling in new[] { codons[i], codons[i].Replace('T', 'U'), codons[i].ToLowerInvariant() })
                {
                    if (expected == '!')
                    {
                        Assert.Throws<ArgumentException>(() => code.Translate(spelling, allowNucleotideX: true),
                            $"Table {row.TableId} {spelling}");
                    }
                    else
                    {
                        Assert.That(code.Translate(spelling, allowNucleotideX: true), Is.EqualTo(expected),
                            $"Table {row.TableId} {spelling}");
                        // Biopython X ≡ N whenever the codon resolves.
                        Assert.That(expected, Is.EqualTo(code.Translate(codons[i].Replace('X', 'N'))));
                    }
                }
                if (expected == '!') rejected++; else resolved++;
            }
        }
        Assert.That(resolved, Is.EqualTo(12407));
        Assert.That(rejected, Is.EqualTo(7060));
    }

    /// <summary>
    /// IsStartCodon / IsStopCodon with allowNucleotideX match Biopython 1.88 ambiguous_dna_by_id
    /// start_codons / stop_codons for every X-containing codon × 27 tables (X-starts: ATX in 2/4/5/11/32,
    /// XTG in 4/11/24/32/33; no X-containing stop codon in any table).
    /// </summary>
    [Test]
    public void IsStartStopCodon_AllXCodons_AllTables_MatchBiopython()
    {
        var codons = XCodons();
        foreach (var row in LoadOracle().Codons)
        {
            var code = GeneticCode.GetByTableNumber(row.TableId);
            foreach (var c in codons)
            {
                Assert.That(code.IsStartCodon(c, allowNucleotideX: true), Is.EqualTo(row.Starts.Contains(c)),
                    $"Table {row.TableId} start {c}");
                Assert.That(code.IsStopCodon(c, allowNucleotideX: true), Is.EqualTo(row.Stops.Contains(c)),
                    $"Table {row.TableId} stop {c}");
                Assert.That(code.IsStartCodon(c.Replace('T', 'U').ToLowerInvariant(), allowNucleotideX: true),
                    Is.EqualTo(row.Starts.Contains(c)), $"Table {row.TableId} start rna {c}");
                // default overloads keep rejecting X
                Assert.That(code.IsStartCodon(c), Is.False);
                Assert.That(code.IsStopCodon(c), Is.False);
            }
        }
        Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("ATX", allowNucleotideX: true), Is.True);
        Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("xtg", allowNucleotideX: true), Is.True);
        Assert.That(GeneticCode.Standard.IsStartCodon("XTG", allowNucleotideX: true), Is.False);
    }

    /// <summary>
    /// 600 random strings (DNA or RNA spelling, mixed case, with X; lengths 0–40; all 27 tables;
    /// to_stop on/off) translate as Biopython 1.88 <c>Seq(s).translate(table, to_stop)</c>; '!' =
    /// Biopython TranslationError / ValueError ↔ ArgumentException.
    /// </summary>
    [Test]
    public void Translate_RandomStringsWithX_MatchBiopython()
    {
        var rows = LoadOracle().Random;
        Assert.That(rows, Has.Count.EqualTo(600));
        foreach (var row in rows)
        {
            var code = GeneticCode.GetByTableNumber(row.TableId);
            if (row.Expected == "!")
            {
                Assert.Throws<ArgumentException>(
                    () => Translator.Translate(row.Input, code, 0, row.ToStop, allowNucleotideX: true),
                    $"Table {row.TableId} to_stop={row.ToStop} '{row.Input}'");
            }
            else
            {
                Assert.That(Translator.Translate(row.Input, code, 0, row.ToStop, allowNucleotideX: true).Sequence,
                    Is.EqualTo(row.Expected), $"Table {row.TableId} to_stop={row.ToStop} '{row.Input}'");
            }
        }
        Assert.That(rows.Count(r => r.Expected != "!" && r.Input.Contains('X', StringComparison.OrdinalIgnoreCase)),
            Is.EqualTo(206), "X-containing strings that Biopython translates");
    }

    /// <summary>
    /// Spot values (Biopython 1.88): GCX → A, gcx → A, GGX → G, YTX → X (F/L), ATX → X (I/M, table 1),
    /// TAX / XXX / ATGXXX → TranslationError; "ATGGCX" → "MA"; to_stop ends before a later XXX.
    /// </summary>
    [Test]
    public void Translate_NucleotideX_SpotValues_MatchBiopython()
    {
        var std = GeneticCode.Standard;
        Assert.Multiple(() =>
        {
            Assert.That(std.Translate("GCX", allowNucleotideX: true), Is.EqualTo('A'));
            Assert.That(std.Translate("gcx", allowNucleotideX: true), Is.EqualTo('A'));
            Assert.That(std.Translate("GCU", allowNucleotideX: true), Is.EqualTo('A'));
            Assert.That(std.Translate("GGX", allowNucleotideX: true), Is.EqualTo('G'));
            Assert.That(std.Translate("YTX", allowNucleotideX: true), Is.EqualTo('X'));
            Assert.That(std.Translate("ATX", allowNucleotideX: true), Is.EqualTo('X'));
            Assert.Throws<ArgumentException>(() => std.Translate("TAX", allowNucleotideX: true));
            Assert.Throws<ArgumentException>(() => std.Translate("XXX", allowNucleotideX: true));
            Assert.Throws<ArgumentException>(() => std.Translate("GCZ", allowNucleotideX: true));
            Assert.Throws<ArgumentException>(() => std.Translate("GC", allowNucleotideX: true));
            Assert.That(Translator.Translate("ATGGCX", null, 0, false, allowNucleotideX: true).Sequence, Is.EqualTo("MA"));
            Assert.That(Translator.Translate("AUGGCXA", null, 0, false, allowNucleotideX: true).Sequence, Is.EqualTo("MA"));
            Assert.Throws<ArgumentException>(() => Translator.Translate("ATGXXX", null, 0, false, allowNucleotideX: true));
            Assert.That(Translator.Translate("ATGTAAXXX", null, 0, true, allowNucleotideX: true).Sequence, Is.EqualTo("M"));
            Assert.That(Translator.Translate("ATGGCC", null, 0, false, allowNucleotideX: false).Sequence, Is.EqualTo("MA"));
        });
    }

    /// <summary>
    /// Default behaviour is unchanged: without the opt-in, X is not an IUPAC nucleotide code and is
    /// rejected by GeneticCode.Translate and every Translator.Translate string overload.
    /// </summary>
    [Test]
    public void Translate_NucleotideX_RejectedByDefault()
    {
        var std = GeneticCode.Standard;
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => std.Translate("GCX"));
            Assert.Throws<ArgumentException>(() => std.Translate("GCX", allowNucleotideX: false));
            Assert.Throws<ArgumentException>(() => Translator.Translate("ATGGCX"));
            Assert.Throws<ArgumentException>(() => Translator.Translate("ATGGCX", null, 0, false, allowNucleotideX: false));
            Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("ATX"), Is.False);
            Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("ATX", allowNucleotideX: false), Is.False);
        });
    }
}
