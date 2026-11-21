namespace Seqeron.Genomics.Infrastructure;

/// <summary>An input sequence file as EMBOSS 6.6 reads it: the first entry of a FASTA, GenBank or EMBL file with its name, accession, description, residues and, for GenBank / EMBL, the feature table.</summary>
internal sealed record EmbossEntry(string Format, EmbossFastaRecord Record, EmbossSeqFeatureTable? Features);

/// <summary>
/// The sequence readers of EMBOSS 6.6 for GenBank and EMBL entries (<c>ajseqread.c</c> <c>seqReadGenbank</c>, <c>seqReadEmbl</c>, <c>seqSetName</c>, <c>seqAccSave</c>, the feature lines handed to <see cref="EmbossSeqFeatureTable"/>) and the format finder
/// of the program for these three formats (a file that starts with <c>&gt;</c> is FASTA, <c>LOCUS</c> GenBank, <c>ID</c> EMBL). A5-8c3b, B09 F218. Only the first entry of the text is read, as the programs read one sequence.
/// </summary>
internal static class EmbossEntryReader
{
    /// <summary>The format of the text as the program finds it (<c>fasta</c>, <c>genbank</c>, <c>embl</c>), or null.</summary>
    public static string? DetectFormat(string text)
    {
        foreach (string raw in text.Split('\n'))
        {
            string line = raw.TrimEnd('\r');
            if (string.IsNullOrWhiteSpace(line))
                continue;
            if (line.StartsWith('>'))
                return "fasta";
            if (line.StartsWith("LOCUS ", StringComparison.Ordinal))
                return "genbank";
            if (line.StartsWith("ID   ", StringComparison.Ordinal) && !IsSwissId(line))
                return "embl";
            return null;
        }

        return null;
    }

    private static bool IsSwissId(string line)
        => line.Contains(" PRT; ", StringComparison.Ordinal) || line.Contains(" Unreviewed; ", StringComparison.Ordinal) || line.Contains(" Reviewed; ", StringComparison.Ordinal) || line.Contains(" Preliminary; ", StringComparison.Ordinal);

    /// <summary>
    /// Reads the first entry. <paramref name="format"/> is <c>-sformat</c> (null = found from the text; <c>fasta</c>, <c>genbank</c> / <c>gb</c>, <c>embl</c> / <c>em</c>).
    /// </summary>
    /// <exception cref="ArgumentException">The text is not an entry of the format ("Unable to read sequence"), or the format is not read by this port.</exception>
    public static EmbossEntry Read(string text, string? format = null) => Read(text, format, null);

    /// <summary><paramref name="context"/> is the state of the feature library shared by the entries of one program run.</summary>
    internal static EmbossEntry Read(string text, string? format, EmbossFeatureReadContext? context)
    {
        ArgumentNullException.ThrowIfNull(text);
        string f = (format ?? DetectFormat(text) ?? "fasta").ToLowerInvariant();
        f = f switch { "gb" => "genbank", "em" => "embl", _ => f };
        switch (f)
        {
            case "fasta":
                {
                    var records = EmbossFastaReader.Read(text);
                    if (records.Count == 0)
                        throw new ArgumentException("Unable to read sequence: no FASTA record.", nameof(text));
                    return new EmbossEntry("fasta", records[0], null);
                }
            case "genbank":
                return ReadGenbank(text, context);
            case "embl":
                return ReadEmbl(text, context);
            default:
                throw new ArgumentException($"The sequence format '{format}' is not read by this port (fasta, genbank and embl are).", nameof(format));
        }
    }

    private static List<string> Lines(string text) => [.. text.Split('\n').Select(l => l.TrimEnd('\r'))];

    // ajStrTokenAssignC(" ") + step over the first word + ajStrTokenNextParseC("\n\r"): the text after the first word and its blanks
    private static string AfterFirstWord(string line)
    {
        int i = 0;
        while (i < line.Length && line[i] == ' ')
            i++;
        while (i < line.Length && line[i] != ' ')
            i++;
        while (i < line.Length && line[i] == ' ')
            i++;
        return line[i..];
    }

    // seqPrefixGenbank: the record type of a line (the first characters)
    private enum Gb { Unknown, Accession, Base, Comment, Definition, Features, Keywords, Locus, Origin, Reference, Source, Version, End, More }

    private static Gb Prefix(string line)
    {
        if (line.StartsWith("ACCESSION", StringComparison.Ordinal)) return Gb.Accession;
        if (line.StartsWith("BASE COUNT", StringComparison.Ordinal)) return Gb.Base;
        if (line.StartsWith("COMMENT", StringComparison.Ordinal)) return Gb.Comment;
        if (line.StartsWith("DEFINITION", StringComparison.Ordinal)) return Gb.Definition;
        if (line.StartsWith("FEATURES", StringComparison.Ordinal)) return Gb.Features;
        if (line.StartsWith("KEYWORDS", StringComparison.Ordinal)) return Gb.Keywords;
        if (line.StartsWith("LOCUS", StringComparison.Ordinal)) return Gb.Locus;
        if (line.StartsWith("ORIGIN", StringComparison.Ordinal)) return Gb.Origin;
        if (line.StartsWith("REFERENCE", StringComparison.Ordinal)) return Gb.Reference;
        if (line.StartsWith("SEQVERSION", StringComparison.Ordinal)) return Gb.Version;
        if (line.StartsWith("SOURCE", StringComparison.Ordinal)) return Gb.Source;
        if (line.StartsWith("//", StringComparison.Ordinal)) return Gb.End;
        if (line.StartsWith("  ", StringComparison.Ordinal)) return Gb.More;
        return Gb.Unknown;
    }

    // seqPrefixGenbankMore: only a line that starts with NINE blanks is GB_MORE_MORE (the program compares from its SECOND character, so the keywords AUTHORS / TITLE / JOURNAL / ORGANISM never match)
    private static bool IsMoreMore(string line) => line.StartsWith("         ", StringComparison.Ordinal);

    // seqReadGenbank: the loop over the header records with its skips: a record handled WITHOUT "done" is followed by the reading of one more line, so the line after the lines of a SOURCE block
    // (ORGANISM and its lineage: all of them GB_MORE) is skipped - a FEATURES line there is lost together with its features (measured on 6.6.0)
    private static EmbossEntry ReadGenbank(string text, EmbossFeatureReadContext? context)
    {
        var lines = Lines(text);
        int i = 0;
        while (i < lines.Count && string.IsNullOrWhiteSpace(lines[i]))
            i++;
        if (i >= lines.Count || !lines[i].StartsWith("LOCUS ", StringComparison.Ordinal))
            throw new ArgumentException("Unable to read sequence: no LOCUS line.", nameof(text));
        string[] locus = lines[i].Split([' ', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
        if (locus.Length == 9 || (locus.Length > 3 && locus[3] == "aa"))
            throw new ArgumentException("GenPept / RefSeq protein entries are not read by this port.", nameof(text));
        string name = locus.Length > 1 ? EmbossFastaReader.SetName(locus[1]) : "";
        string acc = "", desc = "";
        var featureLines = new List<string>();
        bool ok = true;
        string line = "";
        Gb prefix = Gb.Unknown;
        void Next()
        {
            i++;
            ok = i < lines.Count;
            line = ok ? lines[i] : "";
            prefix = ok ? Prefix(line) : Gb.Unknown;
        }

        Next();
        while (ok && prefix != Gb.End && prefix != Gb.Origin && prefix != Gb.Base)
        {
            bool done = false;
            switch (prefix)
            {
                case Gb.Definition:
                    {
                        var definition = new System.Text.StringBuilder(AfterFirstWord(line));
                        Next();
                        done = true;
                        while (ok && prefix == Gb.More)
                        {
                            definition.Append(' ').Append(line.TrimStart(' '));
                            Next();
                        }

                        desc = definition.ToString();
                        break;
                    }

                case Gb.Accession:
                    {
                        string[] tokens = line.Split([' ', ';', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
                        if (acc.Length == 0 && tokens.Length > 1)
                            acc = tokens[1];
                        break;
                    }

                case Gb.Source:
                    Next();
                    done = true;
                    while (ok && prefix == Gb.More)
                    {
                        done = false;
                        Next();
                    }

                    break;
                case Gb.Features:
                    featureLines.Add(line);
                    Next();
                    done = true;
                    while (ok && prefix == Gb.More)
                    {
                        featureLines.Add(line);
                        Next();
                    }

                    break;
                case Gb.Reference:
                case Gb.Comment:
                    Next();
                    done = true;
                    while (ok && IsMoreMore(line))
                        Next();
                    break;
                case Gb.Keywords:
                    Next();
                    done = true;
                    while (ok && prefix == Gb.More)
                        Next();
                    break;
            }

            if (!done)
                Next();
        }

        var seq = new System.Text.StringBuilder();
        if (ok && prefix != Gb.End)
            Next();
        while (ok && prefix != Gb.End)
        {
            if (prefix != Gb.Origin && prefix != Gb.Base)
                seq.Append(line);
            Next();
        }

        var features = EmbossSeqFeatureTable.Read(featureLines, context);
        string residues = EmbossSequenceReader.Clean(seq.ToString());
        features.Length = residues.Length;
        features.SeqId = name;
        return new EmbossEntry("genbank", new EmbossFastaRecord(name, acc, desc, residues), features);
    }

    private static EmbossEntry ReadEmbl(string text, EmbossFeatureReadContext? context)
    {
        var lines = Lines(text);
        int i = 0;
        while (i < lines.Count && string.IsNullOrWhiteSpace(lines[i]))
            i++;
        if (i >= lines.Count || !lines[i].StartsWith("ID   ", StringComparison.Ordinal) || IsSwissId(lines[i]))
            throw new ArgumentException("Unable to read sequence: no EMBL ID line.", nameof(text));
        string[] id = lines[i].Split([' ', ';', '\t', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
        string name = id.Length > 1 ? EmbossFastaReader.SetName(id[1]) : "";
        string acc = "", desc = "";
        var featureLines = new List<string>();
        i++;
        for (; i < lines.Count && !lines[i].StartsWith("SQ", StringComparison.Ordinal) && !lines[i].StartsWith("//", StringComparison.Ordinal); i++)
        {
            string l = lines[i];
            if (l.StartsWith("AC", StringComparison.Ordinal))
            {
                string[] tokens = l.Split([' ', ';', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
                if (acc.Length == 0 && tokens.Length > 1)
                    acc = tokens[1];
            }
            else if (l.StartsWith("DE", StringComparison.Ordinal))
            {
                string part = AfterFirstWord(l);
                desc = desc.Length > 0 ? desc + " " + part : part;
            }
            else if (l.StartsWith("FT   ", StringComparison.Ordinal))
            {
                featureLines.Add(l);
            }
        }

        var seq = new System.Text.StringBuilder();
        if (i < lines.Count && lines[i].StartsWith("SQ", StringComparison.Ordinal))
            i++;
        for (; i < lines.Count && !lines[i].StartsWith("//", StringComparison.Ordinal); i++)
            seq.Append(lines[i]);
        var features = EmbossSeqFeatureTable.Read(featureLines, context);
        string residues = EmbossSequenceReader.Clean(seq.ToString());
        features.Length = residues.Length;
        features.SeqId = name;
        return new EmbossEntry("embl", new EmbossFastaRecord(name, acc, desc, residues), features);
    }
}
