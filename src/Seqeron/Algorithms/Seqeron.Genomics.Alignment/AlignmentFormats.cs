using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The canonical alignment file layer for the multiple-alignment text formats of Biopython 1.88 <c>Bio.Align</c> (<c>fasta</c>, <c>a2m</c>, <c>clustal</c>, <c>phylip</c>, <c>nexus</c>, <c>stockholm</c>): writers byte for byte as
/// <c>Bio.Align.write</c> and readers that give what <c>Bio.Align.parse</c> gives (records with ids, descriptions and printed rows, alignment and column annotations). The EMBOSS formats of
/// <see cref="EmbossAlignmentWriter"/> and the ClustalW program layouts of <see cref="ClustalWAlignmentWriter"/> are different dialects (different widths, headers and name rules) and stay separate files of this layer.
/// </summary>
public static partial class AlignmentFormats
{
    /// <summary>Parses a Biopython format name (case-insensitive).</summary>
    public static bool TryParseFormat(string? name, out AlignmentFileFormat format)
    {
        format = default;
        switch (name?.Trim().ToLowerInvariant())
        {
            case "fasta": format = AlignmentFileFormat.Fasta; return true;
            case "a2m": format = AlignmentFileFormat.A2m; return true;
            case "clustal": format = AlignmentFileFormat.Clustal; return true;
            case "phylip": format = AlignmentFileFormat.Phylip; return true;
            case "nexus": format = AlignmentFileFormat.Nexus; return true;
            case "stockholm": format = AlignmentFileFormat.Stockholm; return true;
            default: return false;
        }
    }

    /// <summary>The names of the supported formats.</summary>
    public static string FormatName(AlignmentFileFormat format) => format.ToString().ToLowerInvariant();

    /// <summary>Writes the alignments as <c>Bio.Align.write(alignments, handle, format)</c> does (the Clustal header once, the NEXUS and A2M formats take exactly one alignment).</summary>
    public static string Write(IReadOnlyList<BioAlignment> alignments, AlignmentFileFormat format, AlignmentWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        options ??= new AlignmentWriteOptions();
        var sb = new StringBuilder();
        if (format is AlignmentFileFormat.A2m or AlignmentFileFormat.Nexus && alignments.Count != 1)
            throw new ArgumentException($"Alignment files in the {FormatName(format)} format can contain a single alignment only.", nameof(alignments));
        if (format == AlignmentFileFormat.Clustal)
            sb.Append(options.ClustalProgram).Append(' ').Append(options.ClustalVersion).Append(" multiple sequence alignment\n\n\n");
        foreach (var a in alignments)
        {
            sb.Append(format switch
            {
                AlignmentFileFormat.Fasta => FormatFasta(a),
                AlignmentFileFormat.A2m => FormatA2m(a),
                AlignmentFileFormat.Clustal => FormatClustal(a),
                AlignmentFileFormat.Phylip => FormatPhylip(a),
                AlignmentFileFormat.Nexus => FormatNexus(a, options.NexusInterleave),
                _ => FormatStockholm(a),
            });
        }

        return sb.ToString();
    }

    /// <summary>Writes one alignment (<c>alignment.format(format)</c>: no file header).</summary>
    public static string Format(BioAlignment alignment, AlignmentFileFormat format, AlignmentWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        options ??= new AlignmentWriteOptions();
        return format switch
        {
            AlignmentFileFormat.Fasta => FormatFasta(alignment),
            AlignmentFileFormat.A2m => FormatA2m(alignment),
            AlignmentFileFormat.Clustal => FormatClustal(alignment),
            AlignmentFileFormat.Phylip => FormatPhylip(alignment),
            AlignmentFileFormat.Nexus => FormatNexus(alignment, options.NexusInterleave),
            _ => FormatStockholm(alignment),
        };
    }

    private static int Columns(BioAlignment a) => a.Records.Count == 0 ? 0 : a.Records[0].Row.Length;

    private static void RequireRectangular(BioAlignment a)
    {
        if (a.Records.Count == 0)
            throw new ArgumentException("Must have at least one sequence");
        int cols = a.Records[0].Row.Length;
        if (a.Records.Any(r => r.Row.Length != cols))
            throw new ArgumentException("All rows of an alignment must have the same number of columns.");
        if (cols == 0)
            throw new ArgumentException("Non-empty sequences are required");
    }

    // ----------------------------------------------------------------- fasta ----

    private static string FormatFasta(BioAlignment a)
    {
        var lines = new List<string>();
        foreach (var r in a.Records)
        {
            lines.Add(r.Description.Length > 0 ? $">{r.Id} {r.Description}" : $">{r.Id}");
            lines.Add(r.Row);
        }

        return string.Join('\n', lines) + "\n";
    }

    // ------------------------------------------------------------------- a2m ----

    private static string FormatA2m(BioAlignment alignment)
    {
        string state = alignment.ColumnAnnotations.TryGetPropertyValue("state", out var node) && node is System.Text.Json.Nodes.JsonValue value && value.TryGetValue<string>(out var text)
            ? text
            : throw new ArgumentException("The a2m format needs the column annotation \"state\" (a string with one 'D' (match) or 'I' (insert) per column; Biopython raises KeyError: 'state' here).", "alignment");
        var lines = new List<string>();
        foreach (var r in alignment.Records)
        {
            lines.Add(r.Description.Length > 0 ? $">{r.Id} {r.Description}" : $">{r.Id}");
            var s = new StringBuilder();
            int n = Math.Min(r.Row.Length, state.Length);
            for (int i = 0; i < n; i++)
            {
                char c = r.Row[i];
                if (state[i] == 'D')
                    s.Append(char.ToUpperInvariant(c));
                else if (state[i] == 'I')
                    s.Append(c == '-' ? '.' : char.ToLowerInvariant(c));
            }

            lines.Add(s.ToString());
        }

        return string.Join('\n', lines) + "\n";
    }

    // --------------------------------------------------------------- clustal ----

    private static string FormatClustal(BioAlignment a)
    {
        RequireRectangular(a);
        int length = Columns(a);
        string? consensus = a.ColumnAnnotations.TryGetPropertyValue("clustal_consensus", out var c) ? c?.GetValue<string>() : null;
        var names = a.Records.Select(r =>
        {
            string name = r.Id.Length > 30 ? r.Id[..30] : r.Id;
            return name.Replace(' ', '_').PadRight(36);
        }).ToList();
        var sb = new StringBuilder();
        int start = 0;
        while (start != length)
        {
            int stop = Math.Min(start + 50, length);
            for (int i = 0; i < names.Count; i++)
                sb.Append(names[i]).Append(a.Records[i].Row, start, stop - start).Append('\n');
            if (consensus is not null)
                sb.Append(' ', 36).Append(SafeSubstring(consensus, start, stop)).Append('\n');
            sb.Append('\n');
            start = stop;
        }

        sb.Append('\n');
        return sb.ToString();
    }

    // Python slice [start:stop] of a string that may be shorter
    private static string SafeSubstring(string s, int start, int stop)
    {
        if (start >= s.Length)
            return "";
        return s.Substring(start, Math.Min(stop, s.Length) - start);
    }

    // ---------------------------------------------------------------- phylip ----

    private const int PhylipIdWidth = 10;

    private static string FormatPhylip(BioAlignment a)
    {
        RequireRectangular(a);
        var names = new List<string>();
        foreach (var r in a.Records)
        {
            string name = r.Id.Trim();
            foreach (char ch in "[](),")
                name = name.Replace(ch.ToString(), "", StringComparison.Ordinal);
            foreach (char ch in ":;")
                name = name.Replace(ch, '|');
            names.Add(name.Length > PhylipIdWidth ? name[..PhylipIdWidth] : name);
        }

        var sb = new StringBuilder();
        sb.Append(a.Records.Count.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(Columns(a).ToString(CultureInfo.InvariantCulture)).Append('\n');
        for (int i = 0; i < names.Count; i++)
            sb.Append(names[i].PadRight(PhylipIdWidth)).Append(a.Records[i].Row).Append('\n');
        return sb.ToString();
    }

    // ----------------------------------------------------------------- nexus ----

    private const string NexusPunctuation = "()[]{}\\,;:=*\\'\"`+-<>";
    private const string NexusWhitespace = " \t\n";

    private static string NexusSafeName(string name)
    {
        string safe = name.Replace("'", "''", StringComparison.Ordinal);
        if (safe.Any(c => NexusPunctuation.Contains(c, StringComparison.Ordinal) || NexusWhitespace.Contains(c, StringComparison.Ordinal)))
            safe = "'" + safe + "'";
        return safe;
    }

    private static string NexusUniqueLabel(ICollection<string> previous, string label)
    {
        while (previous.Contains(label))
        {
            var split = label.Split('.');
            string last = split[^1];
            if (last.StartsWith("copy", StringComparison.Ordinal))
            {
                int copyNum = 1;
                if (last != "copy")
                    copyNum = int.Parse(last[4..], CultureInfo.InvariantCulture) + 1;
                label = $"{string.Join('.', split[..^1])}.copy{copyNum.ToString(CultureInfo.InvariantCulture)}";
            }
            else
            {
                label = string.Concat(label, ".copy");
            }
        }

        return label;
    }

    private static string FormatNexus(BioAlignment a, bool? interleaveOption)
    {
        RequireRectangular(a);
        int columns = Columns(a);
        var types = a.Records.Select(r => r.Annotations.TryGetPropertyValue("molecule_type", out var v) ? v?.GetValue<string>() : null).ToList();
        string datatype;
        if (types.All(t => !string.IsNullOrEmpty(t) && t.Contains("DNA", StringComparison.Ordinal)))
            datatype = "dna";
        else if (types.All(t => !string.IsNullOrEmpty(t) && t.Contains("RNA", StringComparison.Ordinal)))
            datatype = "rna";
        else if (types.All(t => !string.IsNullOrEmpty(t) && t.Contains("protein", StringComparison.Ordinal)))
            datatype = "protein";
        else
            throw new ArgumentException("Need the molecule type to be defined");
        var labels = new List<string>();
        foreach (var r in a.Records)
        {
            string residues = r.Row.Replace("-", "", StringComparison.Ordinal);
            if (datatype == "dna" && residues.Contains('U', StringComparison.Ordinal))
                throw new ArgumentException($"{r.Id} contains U, but DNA alignment");
            if (datatype == "rna" && residues.Contains('T', StringComparison.Ordinal))
                throw new ArgumentException($"{r.Id} contains T, but RNA alignment");
            if (r.Id.Length == 0)
                throw new ArgumentException("New sequence must have a name");
            labels.Add(NexusUniqueLabel(labels, r.Id));
        }

        bool interleave = interleaveOption ?? columns > 1000;
        int blocksize = interleave ? 70 : columns;
        var sb = new StringBuilder();
        sb.Append("#NEXUS\nbegin data;\n");
        sb.Append(CultureInfo.InvariantCulture, $"dimensions ntax={a.Records.Count} nchar={columns};\n");
        sb.Append("format datatype=").Append(datatype).Append(" missing=? gap=-");
        if (interleave)
            sb.Append(" interleave");
        sb.Append(";\nmatrix\n");
        int nameLength = labels.Max(l => NexusSafeName(l).Length);
        if (interleave)
        {
            for (int seek = 0; seek < columns; seek += blocksize)
            {
                for (int i = 0; i < labels.Count; i++)
                    sb.Append(NexusSafeName(labels[i]).PadRight(nameLength + 1)).Append(a.Records[i].Row, seek, Math.Min(blocksize, columns - seek)).Append('\n');
                sb.Append('\n');
            }
        }
        else
        {
            for (int i = 0; i < labels.Count; i++)
            {
                if (blocksize < columns)
                    sb.Append(NexusSafeName(labels[i])).Append('\n');
                else
                    sb.Append(NexusSafeName(labels[i]).PadRight(nameLength + 1));
                for (int seek = 0; seek < columns; seek += blocksize)
                    sb.Append(a.Records[i].Row, seek, Math.Min(blocksize, columns - seek)).Append('\n');
            }
        }

        sb.Append(";\nend;\n");
        return sb.ToString();
    }

    // ------------------------------------------------------------- stockholm ----

    private static readonly (string Key, string Code)[] StockholmGf =
    [
        ("identifier", "ID"), ("accession", "AC"), ("definition", "DE"), ("author", "AU"), ("source of seed", "SE"), ("source of structure", "SS"), ("gathering method", "GA"),
        ("trusted cutoff", "TC"), ("noise cutoff", "NC"), ("build method", "BM"), ("search method", "SM"), ("type", "TP"), ("previous identifier", "PI"), ("comment", "CC"),
        ("clan", "CL"), ("wikipedia", "WK"), ("calibration method", "CB"), ("**", "**"),
    ];

    private static readonly Dictionary<string, string> StockholmGs = new(StringComparer.Ordinal)
    {
        ["accession"] = "AC", ["organism"] = "OS", ["organism classification"] = "OC", ["look"] = "LO",
    };

    private static readonly (string Key, string Code)[] StockholmGr =
    [
        ("secondary structure", "SS"), ("posterior probability", "PP"), ("Catalytic Site Atlas", "CSA"), ("surface accessibility", "SA"), ("transmembrane", "TM"), ("ligand binding", "LI"),
        ("active site", "AS"), ("active site - Pfam predicted", "pAS"), ("active site - from SwissProt", "sAS"), ("intron", "IN"),
    ];

    private static readonly Dictionary<string, string> StockholmGrByKey = StockholmGr.ToDictionary(x => x.Key, x => x.Code, StringComparer.Ordinal);

    private static readonly Dictionary<string, string> StockholmGfByCode = StockholmGf.ToDictionary(x => x.Code, x => x.Key, StringComparer.Ordinal);
    private static readonly Dictionary<string, string> StockholmGsByCode = StockholmGs.ToDictionary(x => x.Value, x => x.Key, StringComparer.Ordinal);
    private static readonly Dictionary<string, string> StockholmGrByCode = StockholmGr.ToDictionary(x => x.Code, x => x.Key, StringComparer.Ordinal);
    // gc_mapping reversed (long name -> code): the named ones, "<gr code>_cons" -> "consensus <gr name>", and the Rfam keywords (underscores become blanks)
    private static readonly Dictionary<string, string> StockholmGcByKey = BuildGcMapping();

    private static readonly Dictionary<string, string> StockholmGcByCode = StockholmGcByKey.ToDictionary(x => x.Value, x => x.Key, StringComparer.Ordinal);

    private static Dictionary<string, string> BuildGcMapping()
    {
        var gc = new Dictionary<string, string>(StringComparer.Ordinal)
        {
            ["reference coordinate annotation"] = "RF", ["consensus sequence"] = "seq_cons", ["consensus score"] = "scorecons", ["consensus score 70"] = "scorecons_70",
            ["consensus score 80"] = "scorecons_80", ["consensus score 90"] = "scorecons_90", ["model mask"] = "MM",
        };
        foreach (var (key, code) in StockholmGr)
            gc["consensus " + key] = code + "_cons";
        foreach (string keyword in new[]
        {
            "RNA_elements", "RNA_structural_element", "RNA_structural_elements", "RNA_ligand_AdoCbl", "RNA_ligand_AqCbl", "RNA_ligand_FMN", "RNA_ligand_Guanidinium", "RNA_ligand_SAM",
            "RNA_ligand_THF_1", "RNA_ligand_THF_2", "RNA_ligand_TPP", "RNA_ligand_preQ1", "RNA_motif_k_turn", "Repeat_unit", "2L3J_B_SS", "CORE", "PK", "PK_SS", "cons",
        })
        {
            gc[keyword.Replace('_', ' ')] = keyword;
        }

        return gc;
    }

    private static string FormatStockholm(BioAlignment a)
    {
        RequireRectangular(a);
        int columns = Columns(a);
        var ann = a.Annotations;
        var lines = new List<string> { "# STOCKHOLM 1.0\n" };
        foreach (var (key, feature) in StockholmGf)
        {
            if (key == "comment")
                continue;
            if (!ann.TryGetPropertyValue(key, out var value) || value is null)
                continue;
            if (key is "author" or "wikipedia")
            {
                foreach (var item in value.AsArray())
                    lines.Add($"#=GF {feature}   {item!.GetValue<string>()}\n");
            }
            else
            {
                lines.Add($"#=GF {feature}   {JsonText(value)}\n");
            }
        }

        if (ann.TryGetPropertyValue("nested domains", out var nested) && nested is not null)
        {
            foreach (var d in nested.AsArray())
            {
                var o = d!.AsObject();
                if (o.TryGetPropertyValue("accession", out var acc) && acc is not null)
                    lines.Add($"#=GF NE   {acc.GetValue<string>()}\n");
                if (o.TryGetPropertyValue("location", out var loc) && loc is not null)
                    lines.Add($"#=GF NL   {loc.GetValue<string>()}\n");
            }
        }

        if (ann.TryGetPropertyValue("references", out var refs) && refs is not null)
        {
            foreach (var rr in refs.AsArray())
            {
                var o = rr!.AsObject();
                o.TryGetPropertyValue("comment", out var comment);
                lines.Add(PythonTextWrap.Fill(comment is null ? null : comment.GetValue<string>(), "#=GF RC   "));
                lines.Add($"#=GF RN   [{JsonText(o["number"])}]\n");
                lines.Add($"#=GF RM   {JsonText(o["medline"])}\n");
                lines.Add(PythonTextWrap.Fill(o["title"]!.GetValue<string>(), "#=GF RT   "));
                lines.Add($"#=GF RA   {JsonText(o["author"])}\n");
                lines.Add($"#=GF RL   {JsonText(o["location"])}\n");
            }
        }

        if (ann.TryGetPropertyValue("database references", out var dbr) && dbr is not null)
        {
            foreach (var d in dbr.AsArray())
            {
                var o = d!.AsObject();
                lines.Add($"#=GF DR   {JsonText(o["reference"])}\n");
                if (o.TryGetPropertyValue("comment", out var c) && c is not null)
                    lines.Add($"#=GF DC   {JsonText(c)}\n");
            }
        }

        if (ann.TryGetPropertyValue("comment", out var cmt) && cmt is not null)
            lines.Add(PythonTextWrap.Fill(cmt.GetValue<string>(), "#=GF CC   "));
        foreach (var kv in ann)
        {
            if (StockholmGf.Any(g => g.Key == kv.Key) || kv.Key is "nested domains" or "references" or "database references")
                continue;
            throw new ArgumentException($"Unknown annotation {kv.Key} found in alignment.annotations");
        }

        lines.Add($"#=GF SQ   {a.Records.Count.ToString(CultureInfo.InvariantCulture)}\n");
        int width = a.Records.Max(r => r.Id.Length);
        int start = Math.Max(width, 20) + 12;
        foreach (var r in a.Records)
        {
            string name = r.Id.PadRight(width);
            foreach (var kv in r.Annotations)
            {
                string feature = StockholmGs.TryGetValue(kv.Key, out var f) ? f : kv.Key;
                lines.Add($"#=GS {name}  {feature} {JsonText(kv.Value)}\n");
            }

            if (r.Description.Length > 0)
                lines.Add($"#=GS {name}  DE {r.Description}\n");
            foreach (string value in r.DbXrefs)
                lines.Add($"#=GS {name}  DR {value}\n");
        }

        foreach (var r in a.Records)
        {
            string aligned = r.Row;
            lines.Add(r.Id.PadRight(start) + aligned + "\n");
            string grName = r.Id.PadRight(width);
            foreach (var kv in r.LetterAnnotations)
            {
                string feature = StockholmGrByKey.TryGetValue(kv.Key, out var f) ? f : kv.Key;
                string value = kv.Value!.GetValue<string>();
                var values = new char[aligned.Length];
                Array.Fill(values, '.');
                int j = 0;
                for (int i = 0; i < aligned.Length; i++)
                {
                    if (aligned[i] is not ('.' or '-'))
                        values[i] = value[j++];
                }

                lines.Add($"#=GR {grName}  {feature} ".PadRight(start) + new string(values) + "\n");
            }
        }

        foreach (var kv in a.ColumnAnnotations)
        {
            string feature = StockholmGcByKey.TryGetValue(kv.Key, out var f) ? f : kv.Key;
            lines.Add($"#=GC {feature} ".PadRight(start) + kv.Value!.GetValue<string>() + "\n");
        }

        lines.Add("//\n");
        _ = columns;
        return string.Concat(lines);
    }

    // the Python str() of an annotation value: strings as they are, numbers in invariant text
    private static string JsonText(JsonNode? node)
    {
        if (node is null)
            return "None";
        if (node is JsonValue v)
        {
            if (v.TryGetValue<string>(out var s))
                return s;
            return v.ToJsonString();
        }

        return node.ToJsonString();
    }
}

/// <summary>Python 3.13 <c>textwrap.fill(text, width=79, break_long_words=False, initial_indent=prefix, subsequent_indent=prefix) + "\n"</c> (the Stockholm long-text lines).</summary>
internal static class PythonTextWrap
{
    private const string Whitespace = "\t\n\x0b\x0c\r ";
    private const int Width = 79;

    private static readonly Regex WordSep = BuildWordSep();

    private static Regex BuildWordSep()
    {
        const string ws = @"[\t\n\x0b\x0c\r ]";
        const string nws = @"[^\t\n\x0b\x0c\r ]";
        const string wp = @"[\w!""'&.,?]";
        const string lt = @"[^\d\W]";
        string pattern = $@"({ws}+|(?<={wp})-{{2,}}(?=\w)|{nws}+?(?:-(?:(?<={lt}{{2}}-)|(?<={lt}-{lt}-))(?={lt}-?{lt})|(?={ws}|\Z)|(?<={wp})(?=-{{2,}}\w)))";
        return new Regex(pattern, RegexOptions.CultureInvariant);
    }

    /// <summary>The wrapped text, every line ending with a newline; null gives the empty string (as <c>_format_long_text</c>).</summary>
    internal static string Fill(string? text, string prefix)
    {
        if (text is null)
            return "";
        var sb = new StringBuilder(text.Length);
        foreach (char c in text)
            sb.Append(Whitespace.Contains(c, StringComparison.Ordinal) ? ' ' : c);
        var chunks = WordSep.Split(sb.ToString()).Where(c => c.Length > 0).ToList();
        chunks.Reverse();
        var lines = new List<string>();
        while (chunks.Count > 0)
        {
            var cur = new List<string>();
            int curLen = 0;
            int width = Width - prefix.Length;
            if (chunks[^1].Trim().Length == 0 && lines.Count > 0)
                chunks.RemoveAt(chunks.Count - 1);
            while (chunks.Count > 0)
            {
                int l = chunks[^1].Length;
                if (curLen + l <= width)
                {
                    cur.Add(chunks[^1]);
                    chunks.RemoveAt(chunks.Count - 1);
                    curLen += l;
                }
                else
                {
                    break;
                }
            }

            if (chunks.Count > 0 && chunks[^1].Length > width && cur.Count == 0)
            {
                cur.Add(chunks[^1]);
                chunks.RemoveAt(chunks.Count - 1);
            }

            if (cur.Count > 0 && cur[^1].Trim().Length == 0)
                cur.RemoveAt(cur.Count - 1);
            if (cur.Count > 0)
                lines.Add(prefix + string.Concat(cur));
        }

        return string.Join('\n', lines) + "\n";
    }
}
