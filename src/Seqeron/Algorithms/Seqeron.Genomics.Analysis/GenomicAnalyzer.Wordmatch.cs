using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>One exact match of <see cref="GenomicAnalyzer.Wordmatch"/>: 0-based starts, the match is extended to the right beyond its word.</summary>
/// <param name="Pattern">Index of the pattern sequence (of the <c>-asequence</c> set).</param>
/// <param name="PatternStart">0-based start in the pattern sequence.</param>
/// <param name="QueryStart">0-based start in the query sequence.</param>
/// <param name="Length">Match length (at least the word size).</param>
public sealed record WordmatchHit(int Pattern, int PatternStart, int QueryStart, int Length);

/// <summary>The statistics of one pattern word (<c>EmbPWordRK</c>).</summary>
/// <param name="Word">The word as first read (case kept).</param>
/// <param name="Sequences">The number of pattern sequences (names) that contain the word.</param>
/// <param name="Matches">The number of matches of the word.</param>
/// <param name="TotalLength">The summed length of those matches.</param>
public sealed record WordmatchWordStatistics(string Word, int Sequences, long Matches, long TotalLength);

/// <summary>The result of one query sequence and of the whole run of <see cref="GenomicAnalyzer.Wordmatch"/>.</summary>
/// <param name="Matches">Per query sequence, per pattern sequence, the matches in the order the program pushes them.</param>
/// <param name="PatternCount">The number of distinct words of the pattern set (<c>Number of patterns/words found</c>).</param>
/// <param name="AllMatches">The number of matches.</param>
/// <param name="SumLength">The summed match length.</param>
/// <param name="MatchesPerSequence">The matches of each pattern sequence.</param>
/// <param name="Words">The words in the order of the program's sorted pattern table (by Rabin-Karp hash).</param>
public sealed record WordmatchResult(IReadOnlyList<IReadOnlyList<IReadOnlyList<WordmatchHit>>> Matches, int PatternCount, long AllMatches, long SumLength, IReadOnlyList<long> MatchesPerSequence, IReadOnlyList<WordmatchWordStatistics> Words);

public static partial class GenomicAnalyzer
{
    #region wordmatch (exact matches of a pattern set against large sequences, Rabin-Karp)

    private const ulong RabinKarpModulus = 1073741789UL;
    private const ulong RabinKarpRadix = 256UL;

    private sealed class WordSeqLocs(int sequence)
    {
        public int Sequence { get; } = sequence;

        public List<int> Locs { get; } = [];

        public long Matches { get; set; }
    }

    private sealed class WordRecord(string word)
    {
        public string Word { get; } = word;

        public Dictionary<string, WordSeqLocs> SeqLocs { get; } = new(StringComparer.Ordinal);

        public ulong Hash { get; set; }

        public long Matches { get; set; }

        public long Length { get; set; }
    }

    /// <summary>
    /// EMBOSS <c>wordmatch</c> (6.6.0; <c>emboss/wordmatch.c</c>, <c>nucleus/embword.c</c> <c>embWordGetTable</c> / <c>embWordRabinKarpInit</c> / <c>embWordRabinKarpSearch</c>; B09 A2-28d, F198): the table of all words of
    /// <paramref name="wordSize"/> symbols of the pattern sequences (the words around the unknown symbol of a sequence's type are left out as in <see cref="SeqMatchAll"/>; a pattern sequence is identified by its NAME: two sequences of one
    /// name share the locations of a word), searched in each query sequence by Rabin-Karp (radix 256, modulus 1073741789, upper-cased symbols). A window whose hash is in the table is verified symbol by symbol and the match is extended to the
    /// right over the pattern sequence; a pattern sequence is not examined again before the end of its last match (the window must end after it: <c>lastlocation</c>), and the symbol after a stop by mismatch counts as scanned.
    /// Words with equal hashes (distinct words of five or more symbols can collide) are not distinguished: the program's binary search returns an arbitrary one of them (this port the first).
    /// </summary>
    /// <param name="patterns">The pattern sequences (<c>-asequence</c>): name and text.</param>
    /// <param name="queries">The query texts (<c>-bsequence</c>).</param>
    /// <param name="wordSize">The word size (<c>-wordsize</c>, default 4, at least 2).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="wordSize"/> &lt; 2.</exception>
    public static WordmatchResult Wordmatch(IReadOnlyList<(string Name, string Text)> patterns, IReadOnlyList<string> queries, int wordSize = 4, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(patterns);
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 2);
        var table = new Dictionary<string, WordRecord>(StringComparer.Ordinal);
        for (int s = 0; s < patterns.Count; s++)
        {
            AlignmentCancellation.Poll(cancellationToken);
            string text = patterns[s].Text;
            if (text.Length < wordSize)
                continue;
            string upper = text.ToUpperInvariant();
            char skip = EmbossSequenceReader.IsNucleotide(text) ? 'N' : 'X';
            var allowed = WordTableStarts(upper, wordSize, skip);
            for (int p = 0; p + wordSize <= text.Length; p++)
            {
                if (!allowed[p])
                    continue;
                string key = upper.Substring(p, wordSize);
                if (!table.TryGetValue(key, out var record))
                    table[key] = record = new WordRecord(text.Substring(p, wordSize));
                if (!record.SeqLocs.TryGetValue(patterns[s].Name, out var locs))
                    record.SeqLocs[patterns[s].Name] = locs = new WordSeqLocs(s);
                locs.Locs.Add(p);
            }
        }

        var words = new List<WordRecord>(table.Values);
        foreach (var w in words)
        {
            ulong hash = 0;
            foreach (char c in w.Word)
                hash = (RabinKarpRadix * hash + char.ToUpperInvariant(c)) % RabinKarpModulus;
            w.Hash = hash;
        }

        var sorted = words.OrderBy(w => w.Hash).ToList();          // qsort by hash (stable here; equal hashes are in no defined order in the program)
        var byHash = new Dictionary<ulong, WordRecord>();
        foreach (var w in sorted)
            byHash.TryAdd(w.Hash, w);
        ulong rm = 1;
        for (int i = 1; i <= wordSize - 1; i++)
            rm = RabinKarpRadix * rm % RabinKarpModulus;

        var perQuery = new List<IReadOnlyList<IReadOnlyList<WordmatchHit>>>();
        long all = 0;
        if (table.Count > 0)
        {
            foreach (string query in queries)
            {
                AlignmentCancellation.Poll(cancellationToken);
                var lists = new List<WordmatchHit>[patterns.Count];
                for (int i = 0; i < lists.Length; i++)
                    lists[i] = [];
                var lastLocation = new int[patterns.Count];
                all += SearchRabinKarp(query, patterns, byHash, wordSize, rm, lastLocation, lists);
                perQuery.Add(lists);
            }
        }

        var perSequence = new long[patterns.Count];
        long sumLength = 0;
        foreach (var w in sorted)
        {
            sumLength += w.Length;
            foreach (var loc in w.SeqLocs.Values)
                perSequence[loc.Sequence] += loc.Matches;
        }

        return new WordmatchResult(perQuery, table.Count, all, sumLength, perSequence, [.. sorted.Select(w => new WordmatchWordStatistics(w.Word, w.SeqLocs.Count, w.Matches, w.Length))]);
    }

    // the C string reads NUL beyond its end (a location of a same-named sequence can lie past the end of the sequence that created the word)
    private static char At(string s, int i) => i < s.Length ? s[i] : '\0';

    // embWordRabinKarpSearch
    private static long SearchRabinKarp(string text, IReadOnlyList<(string Name, string Text)> patterns, Dictionary<ulong, WordRecord> byHash, int plen, ulong rm, int[] lastLocation, List<WordmatchHit>[] lists)
    {
        long matches = 0;
        int tlen = text.Length;
        ulong textHash = 0;
        for (int i = 0; i < plen; i++)
            textHash = (RabinKarpRadix * textHash + (i < tlen ? char.ToUpperInvariant(text[i]) : 0UL)) % RabinKarpModulus;
        for (int i = plen; i <= tlen; i++)
        {
            if (byHash.TryGetValue(textHash, out var word))
            {
                int seq2start = i - plen;
                foreach (var loc in word.SeqLocs.Values)
                {
                    int sidx = loc.Sequence;
                    if (lastLocation[sidx] >= i)
                        continue;
                    string seq = patterns[sidx].Text;
                    int maxloc = 0;
                    foreach (int pos in loc.Locs)
                    {
                        int matchlen = 0;
                        int ii = seq2start;
                        while (matchlen < plen)
                        {
                            if (char.ToUpperInvariant(At(seq, pos + matchlen)) != char.ToUpperInvariant(text[ii++]))
                                break;
                            matchlen++;
                        }

                        if (matchlen < plen)
                            continue;
                        while (ii < tlen && pos + matchlen < seq.Length)
                        {
                            bool equal = char.ToUpperInvariant(seq[pos + matchlen]) == char.ToUpperInvariant(text[ii]);
                            ii++;
                            if (!equal)
                                break;
                            matchlen++;
                        }

                        matches++;
                        lists[sidx].Add(new WordmatchHit(sidx, pos, seq2start, matchlen));
                        if (ii > maxloc)
                            maxloc = ii;
                        word.Length += matchlen;
                        word.Matches++;
                        loc.Matches++;
                    }

                    if (maxloc > 0)
                        lastLocation[sidx] = maxloc;
                }
            }

            if (i < tlen)
                textHash = ((textHash + char.ToUpperInvariant(text[i - plen]) * (RabinKarpModulus - rm)) * RabinKarpRadix + char.ToUpperInvariant(text[i])) % RabinKarpModulus;
        }

        return matches;
    }

    #endregion
}

/// <summary>What a <c>wordmatch</c> run needs besides the sequences.</summary>
/// <param name="Fasta1">The pattern set (<c>-asequence</c>): every FASTA record is a pattern sequence.</param>
/// <param name="Fasta2">The sequences to scan (<c>-bsequence</c>): every FASTA record is scanned.</param>
/// <param name="WordSize"><c>-wordsize</c> (default 4, at least 2).</param>
/// <param name="DumpAlign"><c>-dumpalign</c> (default true): the matches as alignments (<c>match</c> format).</param>
/// <param name="DumpFeature"><c>-dumpfeat</c> (default true): the matches as GFF3 feature files.</param>
/// <param name="CommandLine">The command-line entries printed in the alignment header.</param>
/// <param name="RunDate">The run date; null = now.</param>
/// <param name="Usa1">The first input as typed (printed in the log: <c>a.fa</c>).</param>
/// <param name="Usa2">The second input as typed.</param>
/// <param name="AlignFile">The alignment file name (<c>Report_file:</c>).</param>
/// <param name="AlignFormat">The <c>-aformat</c> (default <c>match</c>).</param>
/// <param name="CancellationToken">Cancellation.</param>
public sealed record EmbossWordmatchRequest(string Fasta1, string Fasta2, int WordSize = 4, bool DumpAlign = true, bool DumpFeature = true, IReadOnlyList<string>? CommandLine = null, DateTime? RunDate = null,
    string Usa1 = "a.fa", string Usa2 = "b.fa", string AlignFile = "stdout", string AlignFormat = "match", CancellationToken CancellationToken = default);

/// <summary>The files <c>wordmatch</c> writes.</summary>
/// <param name="Alignments">The <c>-outfile</c> alignment text (empty without <c>-dumpalign</c>).</param>
/// <param name="Log">The <c>-logfile</c> statistics.</param>
/// <param name="FeatureFile1">The <c>-aoutfeat</c> file: the features of the pattern sequences (GFF3).</param>
/// <param name="FeatureFile2">The <c>-boutfeat</c> file: the features of the scanned sequences.</param>
public sealed record EmbossWordmatchOutput(string Alignments, string Log, string FeatureFile1, string FeatureFile2);

/// <summary>The files of the EMBOSS 6.6.0 program <c>wordmatch</c> (<c>wordmatch.c</c>; B09 A2-28d, F198): <see cref="GenomicAnalyzer.Wordmatch"/> + the alignment writer + the log + the GFF3 feature files.</summary>
public static class EmbossWordmatchPrograms
{
    // the alignments are ungapped exact matches: every identical pair counts as similar (a matrix with a positive diagonal; the program's own matrix would give N / X / * pairs a score of 0 or less)
    private static readonly SubstitutionMatrix IdentityMatrix = BuildIdentityMatrix();

    private static SubstitutionMatrix BuildIdentityMatrix()
    {
        const string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*";
        var scores = new int[alphabet.Length, alphabet.Length];
        for (int i = 0; i < alphabet.Length; i++)
            for (int j = 0; j < alphabet.Length; j++)
                scores[i, j] = i == j ? 1 : -1;
        return new SubstitutionMatrix("IDENTITY", alphabet, scores);
    }

    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    /// <summary>Runs <c>wordmatch</c>.</summary>
    /// <exception cref="ArgumentException">An input has no FASTA record or the alignment format is unknown.</exception>
    public static EmbossWordmatchOutput Run(EmbossWordmatchRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Fasta1);
        ArgumentNullException.ThrowIfNull(r.Fasta2);
        var records1 = EmbossFastaReader.Read(r.Fasta1);
        var records2 = EmbossFastaReader.Read(r.Fasta2);
        if (records1.Count == 0 || records2.Count == 0)
            throw new ArgumentException("An input has no FASTA record.", nameof(r));
        if (!EmbossAlignmentWriter.TryParseFormat(r.AlignFormat, out var format))
            throw new ArgumentException($"Unknown alignment format '{r.AlignFormat}'.", nameof(r));
        bool nuc = EmbossSequenceReader.IsNucleotide(EmbossSequenceReader.Clean(records1[0].Text));
        var type2 = nuc ? EmbossSequenceType.Nucleotide : EmbossSequenceType.StopProtein;
        var set = records1.Select(x => (x.Name, Text: EmbossInputSequence.Read(x.Text, null, EmbossSequenceType.Any).Residues)).ToList();
        var queries = records2.Select(x => (x.Name, Text: EmbossInputSequence.Read(x.Text, null, type2).Residues)).ToList();
        var result = GenomicAnalyzer.Wordmatch(set, [.. queries.Select(q => q.Text)], r.WordSize, r.CancellationToken);
        DateTime date = r.RunDate ?? DateTime.Now;

        var log = new StringBuilder();
        log.Append("Small sequence/file for constructing target patterns: ").Append(r.Usa1).Append('\n');
        log.Append("Large sequence/file to be scanned for patterns: ").Append(r.Usa2).Append('\n');
        log.Append("Number of sequences in the patterns file: ").Append(Int(set.Count)).Append('\n');
        log.Append("Pattern/word length: ").Append(Int(r.WordSize)).Append('\n');
        string matrixName = "";          // no Matrix line in the match format
        if (format != EmbossAlignFormat.Match)
            matrixName = nuc ? "EDNAFULL" : "EBLOSUM62";
        var alignOptions = new EmbossAlignOptions { Program = "wordmatch", CommandLine = r.CommandLine ?? [], RunDate = date, ReportFile = r.AlignFile, IntegerScores = true, MatrixName = matrixName, Matrix = IdentityMatrix };         // ajAlignSetScoreI; the header names the matrix (not in the match format)
        var writer = new EmbossAlignmentWriter(format, alignOptions);
        var tables1 = new List<EmbossFeature>?[set.Count];
        var features2 = new StringBuilder();
        bool featureOutput = r.DumpFeature;
        string featureType = nuc ? "SO:0000110" : "polypeptide_region";            // ajFeatNew("misc_feature") on a protein table is a polypeptide_region without strand
        char featureStrand = nuc ? '+' : '\0';
        if (result.PatternCount > 0)
        {
            log.Append("Number of patterns/words found: ").Append(Int(result.PatternCount)).Append('\n');
            for (int q = 0; q < queries.Count; q++)
            {
                if (!r.DumpAlign && !r.DumpFeature)
                    break;                                              // check mode: statistics only
                for (int i = 0; i < set.Count; i++)
                {
                    var hits = result.Matches[q][i];
                    if (hits.Count == 0)
                        continue;
                    if (r.DumpAlign)
                    {
                        var pairs = hits.Select(h => new EmbossPairAlignment(
                            new EmbossAlignedSequence(set[i].Name, set[i].Text.Substring(h.PatternStart, h.Length), h.PatternStart, set[i].Text.Length - h.PatternStart - h.Length, Description: records1[i].Description),
                            new EmbossAlignedSequence(queries[q].Name, queries[q].Text.Substring(h.QueryStart, h.Length), h.QueryStart, queries[q].Text.Length - h.QueryStart - h.Length, Description: records2[q].Description), h.Length)).ToList();
                        writer.Write(pairs);
                    }

                    if (r.DumpFeature)
                    {
                        var list1 = tables1[i] ??= [];
                        var list2 = new List<EmbossFeature>();
                        foreach (var h in hits)
                        {
                            list1.Add(new EmbossFeature(featureType, "wordmatch", h.PatternStart + 1, h.PatternStart + h.Length, 1.0, featureStrand, 0, [new("note", queries[q].Name)]));
                            list2.Add(new EmbossFeature(featureType, "wordmatch", h.QueryStart + 1, h.QueryStart + h.Length, 1.0, featureStrand, 0, [new("note", set[i].Name)]));
                        }

                        features2.Append(FeatureText(queries[q].Name, queries[q].Text, nuc, list2, date));
                    }
                }
            }

            log.Append("Number of sequences in the file scanned for patterns: ").Append(Int(queries.Count)).Append('\n');
            log.Append("Number of all matches: ").Append(Int(result.AllMatches)).Append(" (wordmatch finds exact matches only)\n");
            var features1 = new StringBuilder();
            if (result.AllMatches > 0)
            {
                log.Append("Sum of match lengths: ").Append(Int(result.SumLength)).Append('\n');
                log.Append("Average match length: ").Append(CFormat.Fixed(result.SumLength * 1.0 / result.AllMatches, 2)).Append('\n');
                log.Append("\nDistribution of the matches among pattern sequences:\n");
                log.Append("----------------------------------------------------\n");
                for (int i = 0; i < set.Count; i++)
                {
                    if (result.MatchesPerSequence[i] > 0)
                        log.Append(set[i].Name.PadRight(42)).Append(": ").Append(Int(result.MatchesPerSequence[i]).PadLeft(8)).Append('\n');
                    if (tables1[i] is { } t)
                        features1.Append(FeatureText(set[i].Name, set[i].Text, nuc, t, date));
                }

                log.Append("\nPattern statistics:\n");
                log.Append("-------------------\n");
                log.Append("Pattern ").Append(r.WordSize > 7 ? new string(' ', r.WordSize - 7) : "").Append("  #pat-sequences  #all-matches  avg-match-length\n");
                foreach (var w in result.Words)
                    if (w.Matches > 0)
                        log.Append(w.Word.PadRight(7)).Append(": ").Append(Int(w.Sequences).PadLeft(12)).Append("  ").Append(Int(w.Matches).PadLeft(12)).Append(' ').Append(CFormat.Fixed(w.TotalLength * 1.0 / w.Matches, 2).PadLeft(17)).Append('\n');
            }

            return new EmbossWordmatchOutput(r.DumpAlign ? writer.Close() : "", log.ToString(), featureOutput ? features1.ToString() : "", featureOutput ? features2.ToString() : "");
        }

        return new EmbossWordmatchOutput(r.DumpAlign ? writer.Close() : "", log.ToString(), "", "");
    }

    // the GFF3 text of one feature table (ajFeattableWrite): the table of a sequence with its matches
    private static string FeatureText(string name, string residues, bool nuc, List<EmbossFeature> features, DateTime date)
    {
        var seq = new EmbossReportSequence(name, name, "", "", residues, 1, residues.Length, 0, false, nuc);
        var table = new EmbossReport(seq, features, !nuc, nuc ? "N" : "P", TableLength: residues.Length, WriteWhenEmpty: true);
        return EmbossReportWriter.Write(table, EmbossReportFormat.Gff, date, new EmbossReportOptions { Program = "wordmatch", RunDate = date, Format = EmbossReportFormat.Gff });
    }
}
