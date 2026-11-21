using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>One region of difference between two sequences (<c>PosODiff</c> of <c>diffseq.c</c>; 1-based, a gap has length 0 and <c>Start</c> after the gap, <c>End</c> before it).</summary>
/// <param name="Start1">Start of the difference in the first sequence.</param>
/// <param name="Start2">Start in the second sequence.</param>
/// <param name="End1">End in the first sequence.</param>
/// <param name="End2">End in the second sequence.</param>
/// <param name="Len1">Length in the first sequence (0 = a gap: something is inserted in the second sequence).</param>
/// <param name="Len2">Length in the second sequence (0 = something is inserted in the first sequence).</param>
public sealed record DiffseqDifference(int Start1, int Start2, int End1, int End2, int Len1, int Len2)
{
    /// <summary>The first sequence was reverse-complemented (<c>-sreverse1</c>).</summary>
    public bool Rev1 { get; init; }

    /// <summary>The second sequence was reverse-complemented.</summary>
    public bool Rev2 { get; init; }
}

/// <summary>The differences found by <see cref="GenomicAnalyzer.Diffseq(string, string, int, bool, CancellationToken)"/>.</summary>
/// <param name="Differences">The differences in order of position.</param>
/// <param name="OverlapStart1">Start of the overlap in the first sequence (1-based; -1 when no region of identity was found).</param>
/// <param name="OverlapEnd1">End of the overlap in the first sequence.</param>
/// <param name="OverlapStart2">Start of the overlap in the second sequence.</param>
/// <param name="OverlapEnd2">End of the overlap in the second sequence.</param>
/// <param name="Snps">Differences of one base against one base.</param>
/// <param name="Transitions">Of those, the transitions (purine to purine, pyrimidine to pyrimidine).</param>
/// <param name="Transversions">Of those, the transversions.</param>
public sealed record DiffseqResult(IReadOnlyList<DiffseqDifference> Differences, int OverlapStart1, int OverlapEnd1, int OverlapStart2, int OverlapEnd2, int Snps, int Transitions, int Transversions);

public static partial class GenomicAnalyzer
{
    #region diffseq (differences between two nearly identical sequences)

    /// <summary>
    /// EMBOSS <c>diffseq</c> (6.6.0; <c>emboss/diffseq.c</c>, <c>nucleus/embword.c</c> <c>embWordGetTable</c> / <c>embWordBuildMatchTable</c> / <c>embWordMatchMin</c>, <c>nucleus/embprop.c</c> <c>embPropTransition</c> / <c>embPropTransversion</c>;
    /// B09 A5-8c3, F196): the exact word matches of <paramref name="wordSize"/> symbols of the two sequences (the same engine and minimal non-overlapping set as <see cref="MegaMerge(string, string, int, bool, string, string, CancellationToken)"/>) are the
    /// regions of identity; the stretches between them are the differences (<c>diffseq_DiffList</c>). With <paramref name="globalDifferences"/> the differences before the first and after the last match are reported too, where a start (end) that
    /// matches by itself is not counted (the small matches the words missed).
    /// </summary>
    /// <param name="sequence1">The first sequence (letters; case as read).</param>
    /// <param name="sequence2">The second sequence.</param>
    /// <param name="wordSize">Word size (<c>-wordsize</c>, default 10, at least 2).</param>
    /// <param name="globalDifferences">Also report the differences at the ends (<c>-globaldifferences</c>).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>The differences; <see cref="DiffseqResult.OverlapStart1"/> is -1 when no match was found (no difference is reported then).</returns>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="wordSize"/> &lt; 2.</exception>
    public static DiffseqResult Diffseq(string sequence1, string sequence2, int wordSize = 10, bool globalDifferences = false, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return DiffseqCore(EmbossSequenceReader.Clean(sequence1), EmbossSequenceReader.Clean(sequence2), wordSize, globalDifferences, false, false, cancellationToken) ?? new DiffseqResult([], -1, 0, 0, 0, 0, 0, 0);
    }

    // null = the program builds no match list (the first sequence has no word table): it writes an empty report without header text or tail
    internal static DiffseqResult? DiffseqCore(string sequence1, string sequence2, int wordSize, bool global, bool rev1, bool rev2, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 2);
        AlignmentCancellation.Poll(cancellationToken);
        var matcher = new PairWordMatcher(sequence1, wordSize);
        if (!matcher.HasTable)
            return null;
        var matches = sequence2.Length < wordSize
            ? []
            : matcher.Matches(sequence2).Select(m => new MegaMatch(m.PositionInText, m.PositionInQuery, m.Length)).ToList();
        AlignmentCancellation.Poll(cancellationToken);
        MinimalMatchSet(matches, wordSize, cancellationToken);

        int len1 = sequence1.Length, len2 = sequence2.Length;
        var diffs = new List<DiffseqDifference>();
        int misstart1 = -1, misstart2 = -1, over1start = -1, over2start = 0;
        foreach (var p in matches)
        {
            int misend1 = p.Start1, misend2 = p.Start2;
            if (misstart1 == -1)
            {
                over1start = p.Start1 + 1;
                over2start = p.Start2 + 1;
                if (global && (p.Start1 != 0 || p.Start2 != 0))
                {
                    int start1 = 1, start2 = 1, l1 = misend1, l2 = misend2;
                    for (int i = 0; l1 > 0 && l2 > 0 && i < len1 && i < len2 && char.ToLowerInvariant(sequence1[i]) == char.ToLowerInvariant(sequence2[i]); i++)
                    {
                        l1--;
                        l2--;
                        start1++;
                        start2++;
                    }

                    diffs.Add(new DiffseqDifference(start1, start2, misend1, misend2, l1, l2) { Rev1 = rev1, Rev2 = rev2 });
                }
            }
            else
            {
                diffs.Add(new DiffseqDifference(misstart1, misstart2, misend1, misend2, misend1 - misstart1 + 1, misend2 - misstart2 + 1) { Rev1 = rev1, Rev2 = rev2 });
            }

            misstart1 = p.Start1 + p.Length + 1;
            misstart2 = p.Start2 + p.Length + 1;
        }

        int over1end = misstart1 - 1, over2end = misstart2 - 1;
        if (global && (misstart1 <= len1 || misstart2 <= len2))
        {
            int end1 = len1, end2 = len2, l1 = end1 - misstart1 + 1, l2 = end2 - misstart2 + 1;
            for (int i = len1, j = len2; l1 != 0 && l2 != 0 && i > 0 && j > 0 && char.ToLowerInvariant(sequence1[i - 1]) == char.ToLowerInvariant(sequence2[j - 1]); i--, j--)
            {
                l1--;
                l2--;
                end1--;
                end2--;
            }

            diffs.Add(new DiffseqDifference(misstart1, misstart2, end1, end2, l1, l2) { Rev1 = rev1, Rev2 = rev2 });
        }

        int snps = 0, transitions = 0, transversions = 0;
        foreach (var d in diffs)
            if (d.Len1 == 1 && d.Len2 == 1)
            {
                snps++;
                string t1 = EmbossText.SubstringFromEnd(sequence1, d.Start1 - 1, d.End1 - 1), t2 = EmbossText.SubstringFromEnd(sequence2, d.Start2 - 1, d.End2 - 1);
                char b1 = t1.Length > 0 ? t1[0] : '\0', b2 = t2.Length > 0 ? t2[0] : '\0';
                if (EmbossPurinePyrimidine.Transition(b1, b2))
                    transitions++;
                if (EmbossPurinePyrimidine.Transversion(b1, b2))
                    transversions++;
            }

        return new DiffseqResult(diffs, over1start, over1end, over2start, over2end, snps, transitions, transversions);
    }

    #endregion
}


/// <summary><c>embPropTransition</c> / <c>embPropTransversion</c> (<c>nucleus/embprop.c</c>): the substitution classes of two nucleotide symbols.</summary>
internal static class EmbossPurinePyrimidine
{
    // embPropPurine: A G R (and lower case); embPropPyrimidine: C T U Y
    private static bool Purine(char c) => char.ToUpperInvariant(c) is 'A' or 'G' or 'R';

    private static bool Pyrimidine(char c) => char.ToUpperInvariant(c) is 'C' or 'T' or 'U' or 'Y';

    internal static bool Transition(char base1, char base2)
    {
        bool bu1 = Purine(base1), bu2 = Purine(base2), by1 = Pyrimidine(base1), by2 = Pyrimidine(base2);
        if (!bu1 && !by1)
            return false;
        if (!bu2 && !by2)
            return false;
        char l1 = char.ToLowerInvariant(base1), l2 = char.ToLowerInvariant(base2);
        if (l1 == l2)
            return false;
        if (l1 == 't' && l2 == 'u')
            return false;
        if (l1 == 'u' && l2 == 't')
            return false;
        if (bu1 && l2 == 'r')
            return false;
        if (bu2 && l1 == 'r')
            return false;
        if (by1 && l2 == 'y')
            return false;
        if (by2 && l1 == 'y')
            return false;
        return bu1 == bu2;
    }

    internal static bool Transversion(char base1, char base2)
    {
        bool bu1 = Purine(base1), bu2 = Purine(base2), by1 = Pyrimidine(base1), by2 = Pyrimidine(base2);
        if (!bu1 && !by1)
            return false;
        if (!bu2 && !by2)
            return false;
        return bu1 != bu2;
    }
}

/// <summary>What <see cref="EmbossDiffseqPrograms.Run"/> needs besides the sequences (the arguments of a <c>diffseq</c> run).</summary>
/// <param name="Fasta1">The first input file (<c>-asequence</c>): the first FASTA record is read.</param>
/// <param name="Fasta2">The second input file (<c>-bsequence</c>).</param>
/// <param name="WordSize">The word size (<c>-wordsize</c>, default 10).</param>
/// <param name="GlobalDifferences"><c>-globaldifferences</c>.</param>
/// <param name="ReportQualifiers">The report qualifiers as typed (<c>-rformat table -raccshow ...</c>); null = the program's default <c>diffseq</c> format.</param>
/// <param name="CommandLine">The command-line entries printed in the report header.</param>
/// <param name="RunDate">The run date (<c>Rundate</c>, <c>#!Date</c>); null = now.</param>
/// <param name="FastaFile1">The name of the first input file (the USA is <c>fasta::file:name</c>).</param>
/// <param name="FastaFile2">The name of the second input file.</param>
/// <param name="ReportFile">The report file name (<c>Report_file:</c>).</param>
/// <param name="FeatureFile1">The feature file of the first sequence (<c>-aoutfeat</c>; default <c>name.diffgff</c>).</param>
/// <param name="FeatureFile2">The feature file of the second sequence (<c>-boutfeat</c>).</param>
/// <param name="CancellationToken">Cancellation of the search.</param>
/// <param name="FeatureInputFiles">The feature files the inputs name with <c>-ufo1</c> / <c>-fformat1</c> / <c>-fopenfile1</c> (and the second input's), by file name (A5-8c3e, B09 F241); a file that is not here does not exist.</param>
public sealed record EmbossDiffseqRequest(string Fasta1, string Fasta2, int WordSize = 10, bool GlobalDifferences = false, string? ReportQualifiers = null, IReadOnlyList<string>? CommandLine = null,
    DateTime? RunDate = null, string FastaFile1 = "a.fa", string FastaFile2 = "b.fa", string ReportFile = "stdout", string? FeatureFile1 = null, string? FeatureFile2 = null, CancellationToken CancellationToken = default,
    IReadOnlyDictionary<string, string>? FeatureInputFiles = null);

/// <summary>The files <c>diffseq</c> writes.</summary>
/// <param name="Report">The <c>-outfile</c> report (7 report tags, default format <c>diffseq</c>).</param>
/// <param name="FeatureFile1"><c>-aoutfeat</c>: the GFF3 feature file of the first sequence (empty when the program writes nothing).</param>
/// <param name="FeatureFile2"><c>-boutfeat</c>: the GFF3 feature file of the second sequence.</param>
public sealed record EmbossDiffseqOutput(string Report, string FeatureFile1, string FeatureFile2)
{
    /// <summary>The name of the file the first feature table is written to (the value of <c>-aoutfeat</c> without a <c>format::</c> prefix, or <c>-ofopenfile</c>); null when the program writes no feature file.</summary>
    public string? FeatureFileName1 { get; init; }

    /// <summary>The name of the file the second feature table is written to.</summary>
    public string? FeatureFileName2 { get; init; }

    /// <summary>The differences found (positions in the coordinates of the sequences as the program holds them after <c>-sbegin</c> / <c>-send</c> / <c>-sreverse</c>); null when the sequences share no word of the word size.</summary>
    public DiffseqResult? Result { get; init; }
}

/// <summary>
/// The files of the EMBOSS 6.6.0 program <c>diffseq</c> (<c>diffseq.c</c>; B09 A5-8c3, F196): <see cref="GenomicAnalyzer.Diffseq(string, string, int, bool, CancellationToken)"/> for the differences, the seven-tag report
/// (<c>reportWriteDiffseq</c> and every other <c>-rformat</c>) written by the one <see cref="EmbossReportWriter"/>, and the GFF3 feature files of both sequences (<c>diffseq_WordMatchListConvDiffToFeat</c>).
/// </summary>
public static class EmbossDiffseqPrograms
{
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static EmbossFeatureTag Note(string text) => new("note", text);

    /// <summary>Runs <c>diffseq</c>.</summary>
    /// <exception cref="ArgumentException">An input has no FASTA record, the second sequence is of another kind (nucleotide / protein) than the first, or the sequence / report qualifiers are refused.</exception>
    /// <exception cref="EmbossReportException">The report format is not written for the sequence type.</exception>
    public static EmbossDiffseqOutput Run(EmbossDiffseqRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Fasta1);
        ArgumentNullException.ThrowIfNull(r.Fasta2);
        var (inputs, rq) = EmbossSequenceQualifiers.ParseInputsAndReport("diffseq", r.ReportQualifiers, true);
        var context = new EmbossFeatureReadContext();
        var entry1 = ReadEntry(r.Fasta1, inputs[0], context, r.FeatureInputFiles);
        var entry2 = ReadEntry(r.Fasta2, inputs[1], context, r.FeatureInputFiles);
        var rec1 = entry1.Record;
        var rec2 = entry2.Record;
        var view1 = EmbossInputSequence.Read(rec1.Text, inputs[0], EmbossSequenceType.Any);
        var view2 = EmbossInputSequence.Read(rec2.Text, inputs[1], view1.IsNucleotide ? EmbossSequenceType.Nucleotide : EmbossSequenceType.StopProtein);
        // the feature tables follow the sequence: ajSeqSetRange reverses the table of a reversed sequence (ilen = the length of the whole sequence), ajSeqTrim cuts it (ajFeattableTrimOff)
        var table1 = PrepareFeatures(entry1.Features, view1, entry1.Format);
        var table2 = PrepareFeatures(entry2.Features, view2, entry2.Format);
        string s1 = view1.Residues, s2 = view2.Residues;
        int len1 = s1.Length, len2 = s2.Length;
        bool nuc = view1.IsNucleotide;

        string formatText = rq.Format ?? "diffseq";
        if (!EmbossReportFormats.TryParse(formatText, out var format))
            throw new EmbossReportException($"Unknown report format '{formatText}'");
        // seqMakeUsa after ajSeqTrim: the range is cut (begin = end = 0) and only a reversed sequence keeps a suffix, "[::r]"
        string usa1 = entry1.Format + "::" + r.FastaFile1 + ":" + rec1.Name + (view1.Reversed ? "[::r]" : ""), usa2 = entry2.Format + "::" + r.FastaFile2 + ":" + rec2.Name + (view2.Reversed ? "[::r]" : "");
        string name2 = rq.ShowUsa ? usa2 : rec2.Name;
        var out1 = OpenFeatureOutput(r.FeatureFile1 ?? rec1.Name + ".diffgff", rq.FeatureOutputs?[0]);
        var out2 = OpenFeatureOutput(r.FeatureFile2 ?? rec2.Name + ".diffgff", rq.FeatureOutputs?[1]);

        var result = GenomicAnalyzer.DiffseqCore(s1, s2, r.WordSize, r.GlobalDifferences, view1.Reversed, view2.Reversed, r.CancellationToken);
        var reportFeatures = new List<EmbossFeature>();
        var features1 = new List<EmbossFeature>();
        var features2 = new List<EmbossFeature>();
        string header = "", tail = "";
        if (result is not null)
        {
            var sb = new System.Text.StringBuilder();
            sb.Append("Compare: ").Append(name2).Append("     from: 1   to: ").Append(Int(len2)).Append("\n\n");
            if (result.OverlapStart1 != -1)
            {
                sb.Append(rq.ShowUsa ? usa1 : rec1.Name).Append(" overlap starts at ").Append(Int(result.OverlapStart1)).Append('\n');
                sb.Append(name2).Append(" overlap starts at ").Append(Int(result.OverlapStart2)).Append("\n\n");
            }

            header = sb.ToString();
            foreach (var d in result.Differences)
            {
                bool rev1 = d.Rev1, rev2 = d.Rev2;
                int start, end;
                bool between = false;
                if (d.Len1 > 0)
                {
                    (start, end) = rev1 ? (len1 - d.End1 + 1, len1 - d.Start1 + 1) : (d.Start1, d.End1);
                }
                else
                {
                    start = d.Start1 - 1;
                    end = d.Start1;
                    between = true;
                }

                // diffseq_Features: the features of each sequence that overlap the difference, first those of the first sequence then those of the second, before the program's own notes
                var tags = new List<EmbossFeatureTag>();
                tags.AddRange(FeatureNotes("first_feature", table1, d.Start1, d.End1));
                tags.AddRange(FeatureNotes("second_feature", table2, d.Start2, d.End2));
                tags.Add(Note("*name " + name2));
                tags.Add(Note("*length " + Int(d.Len2)));
                if (rev2)
                    tags.Add(Note("*strand -"));
                if (d.Len2 > 0)
                {
                    tags.Add(Note("*start " + Int(rev2 ? len2 - d.End2 + 1 : d.Start2)));
                    tags.Add(Note("*end " + Int(rev2 ? len2 - d.Start2 + 1 : d.End2)));
                    tags.Add(Note("*sequence " + EmbossText.SubstringFromEnd(s2, d.Start2 - 1, d.End2 - 1)));
                }
                else
                {
                    tags.Add(Note("*start " + Int(rev2 ? len2 - d.Start2 + 1 : d.End2)));
                    tags.Add(Note("*end " + Int(rev2 ? len2 - d.Start2 + 1 : d.End2)));
                }

                char featureStrand = d.Len1 > 0 && rev1 ? '-' : '\0';
                // ajFeatNewII in a protein table: a polypeptide_region (SO:0000839) that is not itself a protein feature
                reportFeatures.Add(new EmbossFeature(nuc ? "SO:0000110" : "SO:0000839", "diffseq", start, end, 0, featureStrand, 0, tags, between, nuc ? null : false));
                AddConflictFeatures(d, s1, s2, len1, len2, rec1.Name, rec2.Name, nuc, features1, features2);
            }

            var tail2 = new System.Text.StringBuilder();
            if (result.OverlapStart1 != -1)
            {
                tail2.Append("Overlap_end: ").Append(Int(result.OverlapEnd1)).Append(" in ").Append(rq.ShowUsa ? usa1 : rec1.Name).Append('\n');
                tail2.Append("Overlap_end: ").Append(Int(result.OverlapEnd2)).Append(" in ").Append(name2).Append('\n');
                if (nuc)
                {
                    tail2.Append('\n');
                    tail2.Append("SNP_count: ").Append(Int(result.Snps)).Append('\n');
                    tail2.Append("Transitions: ").Append(Int(result.Transitions)).Append('\n');
                    tail2.Append("Transversions: ").Append(Int(result.Transversions)).Append('\n');
                }
            }
            else
            {
                tail2.Append("No regions of alignment found.\n");
            }

            tail = tail2.ToString();
        }

        DateTime runDate = r.RunDate ?? DateTime.Now;
        var featureState = new EmbossFeatureWriteState();
        string FeatureText(string name, string usa, string residues, bool reversed, List<EmbossFeature> list, int format)
        {
            var seq = new EmbossReportSequence(name, usa, "", "", residues, 1, residues.Length, 0, reversed, nuc);
            var t = new EmbossReport(seq, list, !nuc, nuc ? "N" : "P", TableLength: Math.Max(residues.Length, EmbossReportWriter.UnsignedRegionEnd(list)), WriteWhenEmpty: true);
            return EmbossReportWriter.WriteFeatureTable(t, format, runDate, "diffseq", featureState) ?? "";
        }

        // the feature files are written before the report (their writers share a static string with the writers of the report)
        string featureText1 = "", featureText2 = "";
        if (result is not null)
        {
            featureText1 = FeatureText(rec1.Name, usa1, s1, view1.Reversed, features1, out1.Format);
            featureText2 = FeatureText(rec2.Name, usa2, s2, view2.Reversed, features2, out2.Format);
        }

        var options = new EmbossReportOptions
        {
            Program = "diffseq",
            CommandLine = r.CommandLine ?? [],
            RunDate = runDate,
            FeatureLeftover = featureState.Leftover,
            ReportFile = r.ReportFile,
            Format = format,
            FormatText = formatText,
            Tags = [new("int", "start", "start"), new("int", "end", "end"), new("int", "length", "length"), new("str", "name", "name"), new("str", "sequence", "sequence"), new("str", "first_feature", "first_feature"), new("str", "second_feature", "second_feature")],
            Header = header,
            Tail = tail,
            ShowScore = rq.ShowScore ?? true,
            ShowStrand = rq.ShowStrand ?? true,
            ShowAccession = rq.ShowAccession,
            ShowDescription = rq.ShowDescription,
            ShowUsa = rq.ShowUsa,
            MaxAll = rq.MaxAll,
            MaxSeq = rq.MaxSeq,
            Precision = 3,
            TotalSequences = 2,
            TotalLength = len1 + len2,
            AdditionalFiles = [new(out1.File, "Feature file for first sequence"), new(out2.File, "Feature file for second sequence")],
        };
        var sequence = new EmbossReportSequence(rec1.Name, usa1, rec1.Accession, rec1.Description, s1, 1, len1, view1.Offset, view1.Reversed, nuc);
        var table = new EmbossReport(sequence, reportFeatures, !nuc, nuc ? "N" : "P", TableLength: Math.Max(len1, EmbossReportWriter.UnsignedRegionEnd(reportFeatures)), Tags: options.Tags, WriteWhenEmpty: true);
        string report = EmbossReportWriter.WriteFile(options, [table]);

        if (result is null)
            return new EmbossDiffseqOutput(report, "", "") { FeatureFileName1 = out1.File, FeatureFileName2 = out2.File };
        return new EmbossDiffseqOutput(report, featureText1, featureText2)
        {
            Result = result,
            FeatureFileName1 = out1.File,
            FeatureFileName2 = out2.File,
        };
    }

    // acdSetFeatout + ajFeattabOutOpen + featoutUfoProcess: the format and the file name of one feature output from its value on the command line (the default is name.diffgff) and -offormat / -ofopenfile
    internal static (int Format, string File) OpenFeatureOutput(string ufo, EmbossFeatureOutQualifiers? q, string defaultFormat = "")
    {
        string fmt = q?.Format ?? "";
        string formatstr, file;
        if (!string.IsNullOrEmpty(q?.OpenFile))
        {
            formatstr = fmt.Length > 0 ? fmt : defaultFormat;      // the ACD attribute offormat of the output (wordmatch: gff; diffseq has none)
            file = q.OpenFile;
        }
        else
        {
            var m = System.Text.RegularExpressions.Regex.Match(ufo, "^([A-Za-z0-9]*):+(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant);
            string rest = ufo;
            if (m.Success)
            {
                formatstr = m.Groups[1].Value.Length > 0 ? m.Groups[1].Value : "unknown";
                rest = m.Groups[2].Value;
            }
            else
                formatstr = fmt.Length > 0 ? fmt : "gff";
            file = System.Text.RegularExpressions.Regex.IsMatch(rest, "^([^:]+)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant) ? rest : "unknown." + formatstr;
        }

        int index = EmbossReportWriter.FeatureFormatIndex(formatstr);
        if (index < 0)
            index = formatstr.Length == 0 ? 0 : 1;               // featFindOutFormatC: a name that is not a format is written as gff3, an empty one stays "unknown"
        return (index, file);
    }

    // the first entry of the input text in the format of -sformat (or the one the program finds); the plain FASTA reader keeps its own messages
    private static EmbossEntry ReadEntry(string text, EmbossSequenceQualifiers qualifiers, EmbossFeatureReadContext context, IReadOnlyDictionary<string, string>? featureFiles)
    {
        try
        {
            // seqinUfoLocal: with -ufo the features of the sequence file are not read; with none found (or none named) the feature file of -ufo / -fformat / -fopenfile is
            bool ufoGiven = !string.IsNullOrEmpty(qualifiers.FeaturesUfo);
            var entry = EmbossEntryReader.Read(text, qualifiers.Format, context, readFeatures: !ufoGiven);
            if (entry.Features is null)
            {
                var table = EmbossFeatureUfo.Read(qualifiers.FeaturesUfo ?? "", qualifiers.FeaturesFormat, qualifiers.FeaturesFile, entry.Record.Name,
                    name => featureFiles is not null && featureFiles.TryGetValue(name, out var t) ? t : null, context);
                if (table is not null)
                {
                    table.Length = entry.Record.Text.Length;
                    entry = entry with { Features = table };
                }
            }

            return entry;
        }
        catch (ArgumentException e) when (e.Message.StartsWith("Unable to read sequence", StringComparison.Ordinal) && qualifiers.Format is null)
        {
            throw new ArgumentException("An input has no FASTA record.", nameof(text), e);
        }
    }

    private static EmbossSeqFeatureTable? PrepareFeatures(EmbossSeqFeatureTable? table, EmbossSequenceView view, string format)
    {
        if (table is null)
            return null;
        // ajSeqSetNuc / ajSeqSetProt set the type of the table of the sequence when it has none - but a sequence the format reader made nucleotide (seqReadEmbl; the readers of swiss, gff, nbrf, gcg, msf do the same) has its type already, and the call returns at once
        if (format != "embl")
            table.SetType(view.IsNucleotide ? "N" : "P");
        if (view.Reversed)
        {
            // ajSeqReverseDo: ajFeattableSetReverse + ajFeattableReverse (a table that is not known to be nucleotide, as those of the GFF readers, is left as it is)
            table.Rev = !table.Rev;
            table.Reversed = false;
            table.ReverseTable();
        }

        table.TrimOff(view.Offset, view.Residues.Length);
        return table.NewFtable();                                              // ajSeqGetFeatCopy
    }

    // diffseq_Features + diffseq_AddTags: "*first_feature TYPE start-end tag='value' ..." for every feature (but source) that overlaps [start, end]; the translation is left out
    private static IEnumerable<EmbossFeatureTag> FeatureNotes(string kind, EmbossSeqFeatureTable? table, int start, int end)
    {
        if (table is null)
            yield break;
        foreach (var f in table.Features)
        {
            if ((uint)start > (uint)f.GetEnd() || (uint)end < (uint)f.GetStart())        // ajuint parameters and results: a negative position (the first difference of -globaldifferences with no match starts at -1) is a huge number
                continue;
            string typeName = EmbossFeatureVocab.TypeExternal((f.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types, f.Type);       // ajFeatGetType
            if (f.Type.Length == 0 || typeName == "source")
                continue;
            var sb = new System.Text.StringBuilder("*").Append(kind).Append(' ').Append(typeName).Append(' ').Append(Int(f.GetStart())).Append('-').Append(Int(f.GetEnd()));
            foreach (var (name, value) in f.Tags)
            {
                if (name != "translation")
                    sb.Append(' ').Append(name).Append("='").Append(value).Append('\'');
            }

            yield return Note(sb.ToString());
        }
    }

    // diffseq_WordMatchListConvDiffToFeat: a "conflict" feature in each sequence that has a stretch of the difference
    private static void AddConflictFeatures(DiffseqDifference d, string s1, string s2, int len1, int len2, string name1, string name2, bool nuc, List<EmbossFeature> tab1, List<EmbossFeature> tab2)
    {
        if (d.Len1 != 0)
        {
            char forward = nuc ? '+' : '\0';          // ajFeatNew(..., '+') on a protein table has no strand
            var (start, end, strand) = d.Rev1 ? (len1 - d.Start1 + 1, len1 - d.End1 + 1, '-') : (d.Start1, d.End1, forward);
            string note = name2;
            if (d.Len1 == 1 && d.Len2 == 1)
                note = "SNP in " + name2;
            else if (d.Len2 == 0)
                note = $"Insertion of {Int(d.Len1)} bases in {name1}";
            string replace = d.Len2 > 0 ? EmbossText.SubstringFromEnd(s2, d.Start2 - 1, d.End2 - 1) : "";
            tab1.Add(new EmbossFeature("conflict", "diffseq", start, end, 1.0, strand, 0, ConflictTags(note, replace, d.Len2 > 0 ? EmbossText.SubstringFromEnd(s1, d.Start1 - 1, d.End1 - 1) : "", nuc, true)));
        }

        if (d.Len2 != 0)
        {
            char forward = nuc ? '+' : '\0';
            var (start, end, strand) = d.Rev2 ? (len2 - d.Start2 + 1, len2 - d.End2 + 1, '-') : (d.Start2, d.End2, forward);
            string note = name1;
            if (d.Len2 == 1 && d.Len1 == 1)
                note = "SNP in " + name1;
            else if (d.Len1 == 0)
                note = $"Insertion of {Int(d.Len2)} bases in {name2}";
            string replace = d.Len1 > 0 ? EmbossText.SubstringFromEnd(s1, d.Start1 - 1, d.End1 - 1) : "";
            tab2.Add(new EmbossFeature("conflict", "diffseq", start, end, 1.0, strand, 0, ConflictTags(note, replace, EmbossText.SubstringFromEnd(s2, d.Start2 - 1, d.End2 - 1), nuc, false)));
        }
    }

    // nucleotide table: note, replace; protein table: the note is replaced by "source -> replacement" (or MISSING)
    private static List<EmbossFeatureTag> ConflictTags(string note, string replace, string source, bool nuc, bool first)
    {
        _ = first;
        if (nuc)
            return [new("note", note), new("replace", replace)];
        return [new("note", replace.Length > 0 ? source + " -> " + replace : "MISSING")];
    }
}
