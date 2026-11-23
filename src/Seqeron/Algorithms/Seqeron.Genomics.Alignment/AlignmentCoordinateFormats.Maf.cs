using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    private static readonly string[] MafTrackKeys = ["name", "description", "frames", "mafDot", "visibility", "speciesOrder"];

    private static string MetaString(JsonNode? n)
    {
        if (n is null)
            return "";
        return n.GetValueKind() == System.Text.Json.JsonValueKind.String ? n.GetValue<string>() : PyStr(n);
    }

    /// <summary><c>Bio/Align/maf.py</c> <c>AlignmentWriter.write_header</c> / <c>_write_trackline</c>.</summary>
    private static string MafHeader(JsonObject? metadata)
    {
        metadata ??= new JsonObject { ["MAF Version"] = "1" };
        var sb = new StringBuilder();
        if (MafTrackKeys.Any(metadata.ContainsKey))
        {
            sb.Append("track");
            foreach (var (key, node) in metadata)
            {
                string value;
                switch (key)
                {
                    case "name" or "description" or "frames":
                        value = MetaString(node);
                        break;
                    case "mafDot":
                        value = MetaString(node);
                        if (value is not ("on" or "off"))
                            throw new ArgumentException($"mafDot value must be 'on' or 'off' (received '{value}')");
                        break;
                    case "visibility":
                        value = MetaString(node);
                        if (value is not ("dense" or "pack" or "full"))
                            throw new ArgumentException($"visibility value must be 'dense', 'pack', or 'full' (received '{value}')");
                        break;
                    case "speciesOrder":
                        value = string.Join(' ', node!.AsArray().Select(x => x!.GetValue<string>()));
                        break;
                    default:
                        continue;
                }

                if (value.Contains(' '))
                    value = "\"" + value + "\"";
                sb.Append(' ').Append(key).Append('=').Append(value);
            }

            sb.Append('\n');
        }

        sb.Append("##maf");
        foreach (var (key0, node) in metadata)
        {
            if (MafTrackKeys.Contains(key0) || key0 == "Comments")
                continue;
            string key;
            switch (key0)
            {
                case "MAF Version":
                    if (MetaString(node) != "1")
                        throw new ArgumentException("MAF version must be 1");
                    key = "version";
                    break;
                case "Scoring":
                    key = "scoring";
                    break;
                case "Program":
                    key = "program";
                    break;
                default:
                    throw new ArgumentException($"Unexpected key '{key0}' for header");
            }

            sb.Append(' ').Append(key).Append('=').Append(MetaString(node));
        }

        sb.Append('\n');
        if (metadata.TryGetPropertyValue("Comments", out var comments) && comments is JsonArray ca)
        {
            foreach (var c in ca)
                sb.Append("# ").Append(c!.GetValue<string>()).Append('\n');
            sb.Append('\n');
        }

        return sb.ToString();
    }

    /// <summary><c>Alignment._get_row</c>: the printed row of a sequence of the alignment (<c>alignment[index]</c>).</summary>
    internal static string GetRow(BioCoordAlignment a, int index)
    {
        int n = a.Coordinates.Length;
        int cols = a.Coordinates[0].Length;
        var steps = new long[n][];
        for (int i = 0; i < n; i++)
        {
            steps[i] = new long[cols - 1];
            for (int c = 0; c < cols - 1; c++)
                steps[i][c] = a.Coordinates[i][c + 1] - a.Coordinates[i][c];
        }

        var aligned = new bool[cols - 1];
        for (int c = 0; c < cols - 1; c++)
            aligned[c] = Enumerable.Range(0, n).Count(i => steps[i][c] != 0) > 1;
        long[] coordinates = a.Coordinates[index];
        var sequence = a.Sequences[index];
        for (int i = 0; i < n; i++)
        {
            int pos = 0, neg = 0;
            for (int c = 0; c < cols - 1; c++)
            {
                if (!aligned[c])
                    continue;
                if (steps[i][c] > 0)
                    pos++;
                else if (steps[i][c] < 0)
                    neg++;
            }

            if (pos < neg)
            {
                for (int c = 0; c < cols - 1; c++)
                    steps[i][c] = -steps[i][c];
                if (i == index)
                {
                    sequence = ReverseComplement(sequence);
                    long len = sequence.Length;
                    coordinates = [.. coordinates.Select(x => len - x)];
                }
            }
        }

        var line = new StringBuilder();
        long k = coordinates[0];
        for (int c = 0; c < cols - 1; c++)
        {
            long step = steps[index][c];
            long gap = Enumerable.Range(0, n).Max(i => steps[i][c]);
            if (step > 0)
            {
                long j = k + step;
                line.Append(sequence.Slice(k, j) ?? throw new ArgumentException("The sequence contents are undefined (UndefinedSequenceError)."));
                k = j;
            }
            else if (step < 0)
            {
                k += step;
            }
            else
            {
                line.Append('-', (int)gap);
            }
        }

        return line.ToString();
    }

    private static string RecordId(BioCoordSequence s, string fallback) => s.Id ?? fallback;

    /// <summary><c>Bio/Align/maf.py</c> <c>AlignmentWriter.format_alignment</c>.</summary>
    private static string FormatMaf(BioCoordAlignment alignment)
    {
        var annotations = alignment.Attributes.TryGetPropertyValue("annotations", out var an) && an is JsonObject ao ? ao : [];
        var sb = new StringBuilder();
        if (TryAttr(alignment, "score", out var score))
            sb.Append("a score=").Append(CFormat.Fixed(score.AsValue().TryGetValue<long>(out var sl) ? sl : score.GetValue<double>(), 6));
        else
            sb.Append('a');
        if (annotations.TryGetPropertyValue("pass", out var pass) && pass is not null)
            sb.Append(" pass=").Append(PyStr(pass));
        sb.Append('\n');
        int nameWidth = 0, startWidth = 0, sizeWidth = 0, lengthWidth = 0;
        int n = alignment.Sequences.Count;
        var empties = annotations.TryGetPropertyValue("empty", out var em) && em is JsonArray ea ? ea : [];
        for (int i = 0; i < n; i++)
        {
            var record = alignment.Sequences[i];
            var co = alignment.Coordinates[i];
            string name = RecordId(record, "sequence_" + i);
            long start = co[0], end = co[^1], length = record.Length, size;
            if (start < end)
            {
                size = end - start;
            }
            else
            {
                size = start - end;
                start = length - start;
            }

            nameWidth = Math.Max(nameWidth, name.Length);
            startWidth = Math.Max(startWidth, Str(start).Length);
            sizeWidth = Math.Max(sizeWidth, Str(size).Length);
            lengthWidth = Math.Max(lengthWidth, Str(length).Length);
        }

        for (int i = 0; i < empties.Count; i++)
        {
            var e = (JsonObject)empties[i]!;
            var rec = e["record"]!;
            string name = rec["id"]?.GetValue<string>() ?? "sequence_" + (i + n);
            nameWidth = Math.Max(nameWidth, name.Length);
            long start = e["segment"]![0]!.GetValue<long>(), end = e["segment"]![1]!.GetValue<long>(), length = rec["length"]!.GetValue<long>(), size;
            if (start <= end)
            {
                size = end - start;
            }
            else
            {
                size = start - end;
                start = length - start;
            }

            startWidth = Math.Max(startWidth, Str(start).Length);
            sizeWidth = Math.Max(sizeWidth, Str(size).Length);
            lengthWidth = Math.Max(lengthWidth, Str(length).Length);
        }

        int qualityWidth = nameWidth + startWidth + sizeWidth + lengthWidth + 5;
        for (int i = 0; i < n; i++)
        {
            var record = alignment.Sequences[i];
            var co = alignment.Coordinates[i];
            string recordId = RecordId(record, "sequence_" + i);
            long start = co[0], end = co[^1], length = record.Length, size;
            string strand;
            if (start < end)
            {
                size = end - start;
                strand = "+";
            }
            else
            {
                size = start - end;
                start = length - start;
                strand = "-";
            }

            string text = GetRow(alignment, i);
            sb.Append("s ").Append(recordId.PadRight(nameWidth)).Append(' ').Append(Str(start).PadLeft(startWidth)).Append(' ').Append(Str(size).PadLeft(sizeWidth)).Append(' ').Append(strand).Append(' ')
              .Append(Str(length).PadLeft(lengthWidth)).Append(' ').Append(text).Append('\n');
            var ra = record.Annotations;
            if (ra.TryGetPropertyValue("quality", out var quality) && quality is not null)
            {
                string q = quality.GetValue<string>();
                var gq = new StringBuilder();
                int j = 0;
                foreach (char letter in text)
                {
                    if (letter == '-')
                        gq.Append('-');
                    else
                        gq.Append(q[j++]);
                }

                sb.Append("q ").Append(recordId.PadRight(qualityWidth)).Append(' ').Append(gq).Append('\n');
            }

            if (ra.ContainsKey("leftStatus") && ra.ContainsKey("leftCount") && ra.ContainsKey("rightStatus") && ra.ContainsKey("rightCount"))
                sb.Append("i ").Append(recordId.PadRight(nameWidth)).Append(' ').Append(PyStr(ra["leftStatus"]!)).Append(' ').Append(PyStr(ra["leftCount"]!)).Append(' ').Append(PyStr(ra["rightStatus"]!)).Append(' ').Append(PyStr(ra["rightCount"]!)).Append('\n');
        }

        for (int i = 0; i < empties.Count; i++)
        {
            var e = (JsonObject)empties[i]!;
            var rec = e["record"]!;
            string name = (rec["id"]?.GetValue<string>() ?? "sequence_" + (i + n)).PadRight(nameWidth);
            long start = e["segment"]![0]!.GetValue<long>(), end = e["segment"]![1]!.GetValue<long>(), length = rec["length"]!.GetValue<long>(), size;
            string strand;
            if (start <= end)
            {
                size = end - start;
                strand = "+";
            }
            else
            {
                size = start - end;
                start = length - start;
                strand = "-";
            }

            sb.Append("e ").Append(name).Append(' ').Append(Str(start).PadLeft(startWidth)).Append(' ').Append(Str(size).PadLeft(sizeWidth)).Append(' ').Append(strand).Append(' ')
              .Append(Str(length).PadLeft(lengthWidth)).Append(' ').Append(e["status"]!.GetValue<string>()).Append('\n');
        }

        sb.Append('\n');
        return sb.ToString();
    }

    /// <summary><c>shlex.split</c> of a track line (double and single quotes, backslash escapes outside single quotes).</summary>
    private static List<string> ShlexSplit(string line)
    {
        var words = new List<string>();
        var cur = new StringBuilder();
        bool inWord = false;
        int i = 0;
        while (i < line.Length)
        {
            char c = line[i++];
            if (c is ' ' or '\t' or '\n' or '\r')
            {
                if (inWord)
                {
                    words.Add(cur.ToString());
                    cur.Clear();
                    inWord = false;
                }

                continue;
            }

            inWord = true;
            if (c == '"' || c == '\'')
            {
                int end = line.IndexOf(c, i);
                if (end < 0)
                    throw new ArgumentException("No closing quotation");
                cur.Append(line, i, end - i);
                i = end + 1;
            }
            else if (c == '\\' && i < line.Length)
            {
                cur.Append(line[i++]);
            }
            else
            {
                cur.Append(c);
            }
        }

        if (inWord)
            words.Add(cur.ToString());
        return words;
    }

    private static (string Key, string Value) SplitKeyValue(string word)
    {
        var parts = word.Split('=');
        if (parts.Length != 2)
            throw new ArgumentException("not enough values to unpack");
        return (parts[0], parts[1]);
    }

    /// <summary>Bio/Align/maf.py <c>AlignmentIterator</c>.</summary>
    private static List<BioCoordAlignment> ParseMaf(string text)
    {
        var all = Lines(text).ToList();
        int p = 0;
        string line = p < all.Count ? all[p++] : "";
        if (line.StartsWith("track ", StringComparison.Ordinal))
        {
            foreach (string word in ShlexSplit(line).Skip(1))
            {
                var (key, value) = SplitKeyValue(word);
                switch (key)
                {
                    case "name" or "description" or "frames" or "speciesOrder":
                        break;
                    case "mafDot":
                        if (value is not ("on" or "off"))
                            throw new ArgumentException($"Variable mafDot in track line has unexpected value '{value}'");
                        break;
                    case "visibility":
                        if (value is not ("dense" or "pack" or "full"))
                            throw new ArgumentException($"Variable visibility in track line has unexpected value '{value}'");
                        break;
                    default:
                        throw new ArgumentException($"Unexpected variable '{key}' in track line");
                }
            }

            line = p < all.Count ? all[p++] : "";
        }

        var hw = Words(line);
        if (hw.Length == 0 || hw[0] != "##maf")
            throw new ArgumentException("header line does not start with ##maf");
        string? version = null;
        foreach (string word in hw.Skip(1))
        {
            var (key, value) = SplitKeyValue(word);
            if (key == "version")
                version = value;
            else if (key is not ("scoring" or "program"))
                throw new ArgumentException($"Unexpected variable '{key}' in header line");
        }

        if (version != "1")
            throw new ArgumentException("MAF version must be 1");
        string? aline = null;
        while (p < all.Count)
        {
            string l = all[p++];
            if (l.Trim().Length == 0)
                continue;
            if (!l.StartsWith('#'))
            {
                if (!l.StartsWith('a'))
                    throw new ArgumentException("AssertionError");
                aline = l;
                break;
            }
        }

        var result = new List<BioCoordAlignment>();
        while (aline is not null)
        {
            var records = new List<BioCoordSequence>();
            var starts = new List<long>();
            var sizes = new List<long>();
            var strands = new List<string>();
            var aligned = new List<string>();
            var annotations = new JsonObject();
            double? score = null;
            foreach (string word in Words(aline[1..]))
            {
                var (key, value) = SplitKeyValue(word);
                if (key == "score")
                {
                    score = ParseDouble(value);
                }
                else if (key == "pass")
                {
                    long v = PyInt(value);
                    if (v <= 0)
                        throw new ArgumentException($"pass value must be positive (found {v})");
                    annotations["pass"] = v;
                }
                else
                {
                    throw new ArgumentException($"Unknown annotation variable '{key}'");
                }
            }

            string? src = null;
            JsonObject? lastRecordAnnotations = null;
            var empties = new JsonArray();
            string? next = null;
            while (p < all.Count)
            {
                string l = all[p++];
                if (l.StartsWith('#'))
                    continue;
                if (l.StartsWith('a'))
                {
                    next = l;
                    break;
                }

                if (l.StartsWith("s ", StringComparison.Ordinal))
                {
                    var w = Words(l.Trim());
                    if (w.Length != 7)
                        throw new ArgumentException("Error parsing alignment - 's' line must have 7 fields");
                    src = w[1];
                    long start = PyInt(w[2]), size = PyInt(w[3]), srcSize = PyInt(w[5]);
                    string text2 = w[6].Replace('.', '-').Replace('=', '-').Replace('_', '-');
                    aligned.Add(text2);
                    lastRecordAnnotations = [];
                    records.Add(new BioCoordSequence { Id = src, Length = srcSize, Parts = null, Annotations = lastRecordAnnotations });
                    starts.Add(start);
                    sizes.Add(size);
                    strands.Add(w[4]);
                }
                else if (l.StartsWith("i ", StringComparison.Ordinal))
                {
                    var w = Words(l.Trim());
                    if (w.Length != 6 || w[1] != src || lastRecordAnnotations is null)
                        throw new ArgumentException("AssertionError");
                    string[] statusChars = ["C", "I", "N", "n", "M", "T"];
                    if (!statusChars.Contains(w[2]) || !statusChars.Contains(w[4]))
                        throw new ArgumentException("AssertionError");
                    lastRecordAnnotations["leftStatus"] = w[2];
                    lastRecordAnnotations["leftCount"] = PyInt(w[3]);
                    lastRecordAnnotations["rightStatus"] = w[4];
                    lastRecordAnnotations["rightCount"] = PyInt(w[5]);
                }
                else if (l.StartsWith('e'))
                {
                    var w = Words(l[1..]);
                    if (w.Length != 6)
                        throw new ArgumentException("AssertionError");
                    long start = PyInt(w[1]), size = PyInt(w[2]), srcSize = PyInt(w[4]);
                    if (w[5] is not ("C" or "I" or "M" or "n"))
                        throw new ArgumentException("AssertionError");
                    long end = start + size;
                    long[] segment = w[3] == "+" ? [start, end] : [srcSize - start, srcSize - end];
                    empties.Add(new JsonObject
                    {
                        ["record"] = new JsonObject { ["id"] = w[0], ["length"] = srcSize, ["parts"] = null },
                        ["segment"] = new JsonArray(segment[0], segment[1]),
                        ["status"] = w[5],
                    });
                }
                else if (l.StartsWith("q ", StringComparison.Ordinal))
                {
                    var w = Words(l.Trim());
                    if (w.Length != 3 || w[1] != src || lastRecordAnnotations is null)
                        throw new ArgumentException("AssertionError");
                    lastRecordAnnotations["quality"] = w[2].Replace("-", "", StringComparison.Ordinal);
                }
                else if (l.Trim().Length == 0)
                {
                    continue;
                }
                else
                {
                    throw new ArgumentException($"Error parsing alignment - unexpected line:\n{l}");
                }
            }

            if (empties.Count > 0)
                annotations["empty"] = empties;
            var (seqs, co) = ParsePrintedAlignment(aligned);
            for (int i = 0; i < records.Count; i++)
            {
                long srcSize = records[i].Length;
                if (seqs[i].Length != sizes[i])
                    throw new ArgumentException($"sequence size is incorrect (found {seqs[i].Length}, expected {sizes[i]})");
                var parts = new SortedDictionary<long, string> { [starts[i]] = seqs[i] };
                records[i] = new BioCoordSequence { Id = records[i].Id, Length = srcSize, Parts = parts, Annotations = records[i].Annotations };
            }

            for (int i = 0; i < records.Count; i++)
            {
                if (strands[i] == "-")
                {
                    long last = co[i][^1], first = co[i][0];
                    for (int c = 0; c < co[i].Length; c++)
                        co[i][c] = last - first - co[i][c];
                    records[i] = ReverseComplement(records[i]);
                }

                long definedStart = records[i].Parts!.First().Key;
                for (int c = 0; c < co[i].Length; c++)
                    co[i][c] += definedStart;
            }

            var attrs = new JsonObject { ["annotations"] = annotations };
            if (score is not null)
                attrs["score"] = score.Value;
            result.Add(new BioCoordAlignment { Sequences = records, Coordinates = co, Attributes = attrs });
            aline = next;
        }

        return result;
    }

    /// <summary><c>Alignment.parse_printed_alignment</c>: the ungapped sequences and the block coordinates of printed rows (a new coordinate column whenever the gap pattern of the columns changes).</summary>
    internal static (List<string> Sequences, long[][] Coordinates) ParsePrintedAlignment(IReadOnlyList<string> rows)
    {
        int n = rows.Count;
        int cols = rows.Count == 0 ? 0 : rows.Max(r => r.Length);
        var seqs = rows.Select(r => r.Replace("-", "", StringComparison.Ordinal)).ToList();
        var coords = Enumerable.Range(0, n).Select(_ => new List<long> { 0 }).ToList();
        var current = new long[n];
        bool[]? pattern = null;
        long runLength = 0;
        void Flush()
        {
            if (pattern is null)
                return;
            for (int i = 0; i < n; i++)
            {
                if (pattern[i])
                    current[i] += runLength;
                coords[i].Add(current[i]);
            }
        }

        for (int c = 0; c < cols; c++)
        {
            var pat = new bool[n];
            for (int i = 0; i < n; i++)
                pat[i] = c < rows[i].Length && rows[i][c] != '-';
            if (pattern is not null && pat.SequenceEqual(pattern))
            {
                runLength++;
                continue;
            }

            Flush();
            pattern = pat;
            runLength = 1;
        }

        Flush();
        if (cols == 0)
        {
            // PrintedAlignmentParser: no column at all gives the two coordinates 0, 0
            foreach (var c in coords)
                c.Add(0);
        }

        return (seqs, [.. coords.Select(l => l.ToArray())]);
    }
}
