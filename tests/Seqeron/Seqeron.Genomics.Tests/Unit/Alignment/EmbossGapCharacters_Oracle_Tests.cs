using System.Globalization;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A4-3 (B09 F153): the EMBOSS programs read a PLAIN input sequence through the EMBOSS reader, which drops every character that is not a letter or <c>*</c> (gap characters <c>- . ~ ?</c>, <c>#</c>, <c>+</c>, <c>_</c>, digits,
/// blanks, punctuation, non-ASCII letters) before any length check, coordinate or score: <c>ajseqread.c seqAppend</c> + <c>ajseqtype.c ajSeqTypeCheckIn</c> + <c>ajstr.c ajStrRemoveGap</c>. One shared cleaner
/// (<see cref="EmbossSequenceReader"/>) serves every port; this fixture locks it against the installed 6.6.0 programs: &gt;= 200 rows per program whose inputs contain such characters
/// (TestData/Alignment/emboss_gap_oracle.tsv, scripts/alignment/generate_emboss_gap_oracle.py; the per-character measurement is scripts/alignment/measure_emboss_reader_chars.py).
/// </summary>
[TestFixture]
public class EmbossGapCharacters_Oracle_Tests
{
    private sealed record Row(int Line, string Program, string Kind, string Params, string[] Seqs, string Expected)
    {
        public override string ToString() => $"line {Line} {Program} {Kind} [{Params}] {string.Join(" / ", Seqs.Where(s => s.Length > 0))}";
    }

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_gap_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0)
                continue;
            var f = text.Split('\t');
            rows.Add(new Row(line, f[1], f[2], f[3], [f[4], f[5], f[6]], f[7]));
        }

        return rows;
    }

    private static string Score(double s) => Math.Round(s, 3).ToString("0.0##", CultureInfo.InvariantCulture);

    private static string Pair(AlignmentResult r) => $"{Score(r.ExactScore)}|{r.AlignedSequence1}|{r.AlignedSequence2}";

    private static string Unescape(string s) => s.Replace("\\n", "\n", StringComparison.Ordinal).Replace("\\t", "\t", StringComparison.Ordinal).Replace("\\r", "\r", StringComparison.Ordinal).Replace("\\\\", "\\", StringComparison.Ordinal);

    private static string Hits(IEnumerable<string> hits)
    {
        string s = string.Join(" ", hits);
        return s.Length == 0 ? "-" : s;
    }

    private static string Produce(Row r)
    {
        string a = r.Seqs[0], b = r.Seqs[1], c = r.Seqs[2];
        bool nuc = r.Kind == "n";
        switch (r.Program)
        {
            case "needle":
                return Pair(SequenceAligner.NeedleAlign(a, b));
            case "water":
                return Pair(SequenceAligner.WaterAlign(a, b));
            case "stretcher":
                {
                    var s = SequenceAligner.StretcherAlign(a, b);
                    return string.Create(CultureInfo.InvariantCulture, $"{s.Score}|{s.AlignedSequence1}|{s.AlignedSequence2}");
                }
            case "matcher":
                {
                    int k = int.Parse(r.Params, CultureInfo.InvariantCulture);
                    var list = nuc
                        ? SequenceAligner.LocalAlignAffineTopN(a, b, new ScoringMatrix(5, -4, -10, -4), k, 0, TopNTieOrder.Sim)
                        : SequenceAligner.LocalAlignAffineTopN(a, b, new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -4), k, 0, TopNTieOrder.Sim);
                    if (list.Count < k)
                        return "ABORT";
                    return string.Join(';', list.Select(x => string.Create(CultureInfo.InvariantCulture, $"{x.Score}|{x.StartPosition1}|{x.StartPosition2}|{x.EndPosition1}|{x.EndPosition2}|{x.AlignedSequence1}|{x.AlignedSequence2}")));
                }
            case "merger":
                return SequenceAligner.Merger(a, b).MergedSequence;
            case "megamerger":
                {
                    var m = GenomicAnalyzer.MegaMerge(a, b, int.Parse(r.Params, CultureInfo.InvariantCulture));
                    return m.MergedSequence.Length == 0 ? "-" : m.MergedSequence;
                }
            case "seqmatchall":
                {
                    var seqs = r.Seqs.Where(s => s.Length > 0).ToList();
                    var matches = GenomicAnalyzer.SeqMatchAll(seqs, int.Parse(r.Params, CultureInfo.InvariantCulture))
                        .OrderBy(m => m.Sequence1).ThenBy(m => m.Sequence2).ThenBy(m => m.Position1).ThenBy(m => m.Position2).ThenBy(m => m.Length);
                    return Hits(matches.Select(m => $"{m.Sequence1}:{m.Sequence2}:{m.Position1}:{m.Position2}:{m.Length}"));
                }
            case "supermatcher":
            case "needleall":
                {
                    var queries = new[] { new EmbossSequence("q1", a) };
                    var targets = new[] { b, c }.Where(s => s.Length > 0).Select((s, i) => new EmbossSequence($"t{i + 1}", s)).ToList();
                    var result = r.Program == "supermatcher"
                        ? SequenceAligner.Supermatcher(queries, targets, wordLength: 4)
                        : SequenceAligner.NeedleAll(queries, targets);
                    // needleall / supermatcher print the alignments in the order of the program's loops; the oracle lists them as printed
                    string s = string.Join(';', result.Alignments.Select(x => $"{x.First.Row}|{x.Second.Row}"));
                    return s.Length == 0 ? "-" : s;
                }
            case "fuzznuc":
                {
                    var p = r.Params.Split('|');
                    string Strand(PatternHit h) => h.IsReverseStrand ? "-" : "+";
                    var hits = GenomicAnalyzer.FindPatternHits(a, p[0], int.Parse(p[1], CultureInfo.InvariantCulture), p[2] == "1", default, sequenceAmbiguity: true)
                        .Select(h => (Start: h.Position + 1, End: h.Position + h.Length, Strand: Strand(h)))
                        .OrderBy(h => h.Start).ThenBy(h => h.End).ThenBy(h => h.Strand, StringComparer.Ordinal);
                    return $"L={EmbossSequenceReader.Clean(a).Length};{Hits(hits.Select(h => $"{h.Start}:{h.End}:{h.Strand}"))}";
                }
            case "fuzzpro":
                {
                    var hits = GenomicAnalyzer.FindProteinPatternHits(a, r.Params).OrderBy(h => h.Position).ThenBy(h => h.Length);
                    return $"L={EmbossSequenceReader.Clean(a).Length};{Hits(hits.Select(h => $"{h.Position + 1}:{h.Position + h.Length}"))}";
                }
            case "fuzztran":
                {
                    var p = r.Params.Split('|');
                    var hits = GenomicAnalyzer.FindTranslatedPatternHits(a, p[0], p[1])
                        .OrderBy(h => h.Frame).ThenBy(h => h.Position).ThenBy(h => h.Length);
                    return $"L={EmbossSequenceReader.Clean(a).Length};{Hits(hits.Select(h => $"{h.Frame}:{h.EmbossStart}:{h.EmbossEnd}"))}";
                }
            case "etandem":
                {
                    var p = r.Params.Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
                    var hits = GenomicAnalyzer.FindTandemRepeatsEtandem(a, p[0], p[1], p[2], p[3] == 1, p[4] == 1);
                    return Hits2(hits.Select(h => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount, CFormat.Fixed(h.Identity, 1), h.Consensus)));
                }
            case "equicktandem":
                {
                    var p = r.Params.Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
                    var hits = GenomicAnalyzer.FindTandemRepeatsEquicktandem(a, p[0], p[1]);
                    return Hits2(hits.Select(h => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount)));
                }
            case "dotmatcher":
                {
                    var p = r.Params.Split(',');
                    var segs = DotMatcher.Compute(a, b, int.Parse(p[1], CultureInfo.InvariantCulture), int.Parse(p[2], CultureInfo.InvariantCulture), SubstitutionMatrix.Get(p[0]));
                    return Hits(segs.Select(s => $"{s.Start2}:{s.Start1}:{s.End2}:{s.End1}"));
                }
            case "consambig":
                {
                    string cons = MotifFinder.GenerateEmbossAmbiguousConsensus(r.Seqs.Where(s => s.Length > 0));
                    return cons.Length == 0 ? "-" : cons;
                }
            case "showalign":
                {
                    var rows = r.Seqs.Where(s => s.Length > 0).ToList();
                    var names = Enumerable.Range(0, rows.Count).Select(i => $"r{i}").ToList();
                    return MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Width = 30 });
                }
            case "distmat":
                {
                    int method = int.Parse(r.Params, CultureInfo.InvariantCulture);
                    var rows = r.Seqs.Where(s => s.Length > 0).ToList();
                    var d = MultipleAlignmentStatistics.CalculateDistanceMatrix(rows, new MsaDistanceOptions { Method = (MsaDistanceMethod)method });
                    var cells = new List<string>();
                    for (int i = 0; i < rows.Count - 1; i++)
                        for (int j = i + 1; j < rows.Count; j++)
                            cells.Add(PrintedDistance(d[i, j]));
                    return string.Join(",", cells);
                }
            default:
                throw new InvalidOperationException(r.Program);
        }
    }

    /// <summary>distmat's print: D = D * (float)100, then <c>%.2f</c> (glibc: nan / inf).</summary>
    private static string PrintedDistance(double distance)
    {
        float v = (float)distance * 100f;
        if (float.IsNaN(v))
            return "nan";
        if (float.IsInfinity(v))
            return v > 0 ? "inf" : "-inf";
        return v.ToString("F2", CultureInfo.InvariantCulture);
    }

    private static string Hits2(IEnumerable<string> rows)
    {
        string s = string.Join("|", rows);
        return s.Length == 0 ? "-" : s;
    }

    private static string Expected(Row r) => r.Program switch
    {
        "showalign" => Unescape(r.Expected),
        "distmat" => r.Expected.Replace("-nan", "nan", StringComparison.Ordinal),      // glibc prints the sign of a NaN; the value is the same
        _ => r.Expected,
    };

    private static void Check(string program, int minimum)
    {
        var rows = Oracle.Value.Where(r => r.Program == program).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(minimum), $"oracle rows of {program}");
        var failures = new List<string>();
        int dirty = 0;
        foreach (var row in rows)
        {
            try
            {
                string got = Produce(row);
                if (got != Expected(row))
                    failures.Add($"{row}\n   expected {Expected(row)}\n   actual   {got}");
                if (row.Seqs.Any(s => s.Any(ch => !char.IsAsciiLetter(ch))))
                    dirty++;
                // the cleaned inputs give the same result as the dirty ones (the reader is the only difference between them)
                var cleaned = row with { Seqs = [.. row.Seqs.Select(s => program is "consambig" or "showalign" or "distmat" ? s : EmbossSequenceReader.Clean(s))] };
                if (program is not ("consambig" or "showalign" or "distmat") && Produce(cleaned) != got)
                    failures.Add($"{row}: the result of the cleaned input differs from the result of the input");
            }
            catch (Exception e) when (e is not AssertionException)
            {
                failures.Add($"{row}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(dirty, Is.GreaterThanOrEqualTo(minimum), "at least `minimum` rows contain a gap character or another non-sequence character (the generator inserts at least one into every sequence of a plain program)");
    }

    [Test] public void Needle_DropsGapCharactersLikeEmboss_OnAllRows() => Check("needle", 200);
    [Test] public void Water_DropsGapCharactersLikeEmboss_OnAllRows() => Check("water", 200);
    [Test] public void Stretcher_DropsGapCharactersLikeEmboss_OnAllRows() => Check("stretcher", 200);
    [Test] public void Matcher_DropsGapCharactersLikeEmboss_OnAllRows() => Check("matcher", 200);
    [Test] public void Merger_DropsGapCharactersLikeEmboss_OnAllRows() => Check("merger", 200);
    [Test] public void Megamerger_DropsGapCharactersLikeEmboss_OnAllRows() => Check("megamerger", 200);
    [Test] public void Seqmatchall_DropsGapCharactersLikeEmboss_OnAllRows() => Check("seqmatchall", 200);
    [Test] public void Supermatcher_DropsGapCharactersLikeEmboss_OnAllRows() => Check("supermatcher", 200);
    [Test] public void Needleall_DropsGapCharactersLikeEmboss_OnAllRows() => Check("needleall", 200);
    [Test] public void Fuzznuc_DropsGapCharactersLikeEmboss_OnAllRows() => Check("fuzznuc", 200);
    [Test] public void Fuzzpro_DropsGapCharactersLikeEmboss_OnAllRows() => Check("fuzzpro", 200);
    [Test] public void Fuzztran_DropsGapCharactersLikeEmboss_OnAllRows() => Check("fuzztran", 200);
    [Test] public void Etandem_DropsGapCharactersLikeEmboss_OnAllRows() => Check("etandem", 200);
    [Test] public void Equicktandem_DropsGapCharactersLikeEmboss_OnAllRows() => Check("equicktandem", 200);
    [Test] public void Dotmatcher_DropsGapCharactersLikeEmboss_OnAllRows() => Check("dotmatcher", 200);
    [Test] public void Consambig_ReadsAlignmentsLikeEmboss_OnAllRows() => Check("consambig", 200);
    [Test] public void Showalign_ReadsAlignmentsLikeEmboss_OnAllRows() => Check("showalign", 200);
    [Test] public void Distmat_ReadsAlignmentsLikeEmboss_OnAllRows() => Check("distmat", 200);

    [Test]
    public void NeedleScoresFiftyEightWhereTheProgramDoes_TheAuditExample()
    {
        // audit A4-3: needle 6.6.0 `AC-GTTAGC.TT*` vs `ACGTTAGCTT` = score 58, Length 11 (it was 43 with the gap characters kept)
        var r = SequenceAligner.NeedleAlign("AC-GTTAGC.TT*", "ACGTTAGCTT");
        Assert.That(r.ExactScore, Is.EqualTo(58.0));
        Assert.That(r.AlignedSequence1, Has.Length.EqualTo(11));
        Assert.That(r.AlignedSequence1, Is.EqualTo("ACGTTAGCTT*"));
        Assert.That(SequenceAligner.WaterAlign("AC-GTTAGC.TT*", "ACGTTAGCTT").ExactScore, Is.EqualTo(58.0));
    }

    [Test]
    public void Clean_KeepsLettersAndStar_AndDropsEveryOtherCharacter()
    {
        Assert.That(EmbossSequenceReader.Clean("AC-GT.T~A?C#G+T_1 2\t3é*Ω"), Is.EqualTo("ACGTTACGT*"));
        string clean = "ACGT";
        Assert.That(EmbossSequenceReader.Clean(clean), Is.SameAs(clean));          // nothing to remove: the argument itself
        Assert.That(EmbossSequenceReader.Clean("---"), Is.Empty);
        for (int c = 0; c < 256; c++)
            Assert.That(EmbossSequenceReader.Clean(((char)c).ToString()), Is.EqualTo(char.IsAsciiLetter((char)c) || c == '*' ? ((char)c).ToString() : ""), $"U+{c:X4}");
        Assert.That(EmbossSequenceReader.CleanGapped("AC-GT.T~A?C#G+T_1 2*é"), Is.EqualTo("AC-GT.T~A?C?G?T*"));
        string gapped = "AC-GT.~?*";
        Assert.That(EmbossSequenceReader.CleanGapped(gapped), Is.SameAs(gapped));
        Assert.Throws<ArgumentNullException>(() => EmbossSequenceReader.Clean((string)null!));
    }

    [Test]
    public void IsNucleotide_IsTheAjSeqIsNucRule_OfEveryPort()
    {
        Assert.That(EmbossSequenceReader.IsNucleotide("ACGTUBDHKMNRSVWXYacgtu?.~-"), Is.True);
        Assert.That(EmbossSequenceReader.IsNucleotide(""), Is.True);
        Assert.That(EmbossSequenceReader.IsNucleotide("ACGE"), Is.False);
        Assert.That(EmbossSequenceReader.IsNucleotide("ACG*"), Is.False);
        Assert.That(EmbossSequenceReader.IsNucleotide("ACGſ"), Is.False);      // U+017F (long s) upper-cases to S under ToUpperInvariant but is not a nucleotide letter
    }
}
