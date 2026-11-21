using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    private static string PyValueStr(JsonNode? n)
    {
        if (n is null)
            return "None";
        if (n is JsonArray arr)
            return "[" + string.Join(", ", arr.Select(x => x is JsonValue v && v.TryGetValue<string>(out var s) ? "'" + s + "'" : PyValueStr(x))) + "]";
        if (n is JsonValue val && val.TryGetValue<bool>(out var b))
            return b ? "True" : "False";
        return PyStr(n);
    }

    /// <summary><c>Bio/Align/sam.py</c> <c>AlignmentWriter.write_header</c>.</summary>
    private static string SamHeader(CoordinateWriteOptions options)
    {
        var metadata = options.Metadata ?? [];
        var sb = new StringBuilder();
        if (metadata.TryGetPropertyValue("HD", out var hd) && hd is JsonObject values)
        {
            if (!values.TryGetPropertyValue("VN", out var vn))
                throw new ArgumentException("KeyError: 'VN'");
            var f = new List<string> { "@HD", "VN:" + PyValueStr(vn) };
            foreach (var (key, value) in values)
            {
                if (key == "VN")
                    continue;
                f.Add(key + ":" + PyValueStr(value));
            }

            sb.Append(string.Join('\t', f)).Append('\n');
        }

        foreach (var record in options.Targets ?? [])
        {
            var f = new List<string> { "@SQ", "SN:" + (record.Id ?? throw new ArgumentException("TypeError: sequence item: expected str instance, NoneType found")) };
            if (record.Length < 0)
                throw new ArgumentException("TypeError: object of type 'NoneType' has no len()");
            f.Add("LN:" + Str(record.Length));
            foreach (var (key, value) in record.Annotations)
            {
                switch (key)
                {
                    case "alternate_locus": f.Add("AH:" + PyValueStr(value)); break;
                    case "names":
                        f.Add("AN:" + string.Join(",", (value?.AsArray() ?? throw new ArgumentException("TypeError")).Select(x => PyValueStr(x is JsonValue v && v.TryGetValue<string>(out var s) ? JsonValue.Create(s) : x))));
                        break;
                    case "assembly": f.Add("AS:" + PyValueStr(value)); break;
                    case "MD5": f.Add("M5:" + PyValueStr(value)); break;
                    case "species": f.Add("SP:" + PyValueStr(value)); break;
                    case "topology":
                        string t = PyValueStr(value);
                        if (t is not ("linear" or "circular"))
                            throw new ArgumentException("AssertionError");
                        f.Add("PP:" + t);
                        break;
                    case "URI": f.Add("UR:" + PyValueStr(value)); break;
                    default: f.Add((key.Length > 2 ? key[..2] : key) + ":" + PyValueStr(value)); break;
                }
            }

            if (record.Description != "<unknown description>")
                f.Add("DS:" + record.Description);
            sb.Append(string.Join('\t', f)).Append('\n');
        }

        foreach (var (tag, rowsNode) in metadata)
        {
            if (tag == "HD")
                continue;
            foreach (var row in rowsNode?.AsArray() ?? throw new ArgumentException("TypeError"))
            {
                var f = new List<string> { "@" + tag };
                foreach (var (key, value) in row!.AsObject())
                    f.Add(key + ":" + PyValueStr(value));
                sb.Append(string.Join('\t', f)).Append('\n');
            }
        }

        return sb.ToString();
    }

    private static string RevCompDna(string s)
    {
        var sb = new StringBuilder(s.Length);
        for (int i = s.Length - 1; i >= 0; i--)
            sb.Append(ReverseComplement(s[i].ToString()));
        return sb.ToString();
    }

    /// <summary>The query / target residues of a record as Python's <c>bytes(seq)</c>; null = the contents are (partly) undefined.</summary>
    private static string? Contents(BioCoordSequence s) => s.Length < 0 ? null : s.Slice(0, s.Length);

    /// <summary><c>Bio/Align/sam.py</c> <c>AlignmentWriter.format_alignment</c>.</summary>
    private static string FormatSam(BioCoordAlignment alignment, CoordinateWriteOptions options)
    {
        if (alignment.NoCoordinates)
            throw new ArgumentException("AttributeError: 'NoneType' object has no attribute 'transpose'");
        if (alignment.Sequences.Count != 2)
            throw new ArgumentException("ValueError: not enough values to unpack");
        var target = alignment.Sequences[0];
        var query = alignment.Sequences[1];
        int m = alignment.Coordinates[0].Length;
        if (m == 0)
            throw new ArgumentException("IndexError");
        var ts = (long[])alignment.Coordinates[0].Clone();
        var qs = (long[])alignment.Coordinates[1].Clone();
        long? hardLeft = null, hardRight = null;
        if (query.Annotations.TryGetPropertyValue("hard_clip_left", out var hl) && hl is not null)
            hardLeft = Int(hl);
        if (query.Annotations.TryGetPropertyValue("hard_clip_right", out var hr) && hr is not null)
            hardRight = Int(hr);
        string qual = "*";
        if (query.LetterAnnotations.TryGetPropertyValue("phred_quality", out var ph) && ph is JsonArray phArr)
            qual = new string(phArr.Select(v => (char)(Int(v!) + 33)).ToArray());
        string qName = query.Id ?? throw new ArgumentException("TypeError: sequence item 0: expected str instance, NoneType found");
        string rname = target.Id ?? throw new ArgumentException("TypeError: sequence item 2: expected str instance, NoneType found");
        long qSize = query.Length;
        if (ts[0] > ts[^1])
        {
            Array.Reverse(ts);
            Array.Reverse(qs);
        }

        long flag;
        string? qContents = Contents(query);
        if (qs[0] < qs[^1])
        {
            flag = 0;
        }
        else
        {
            flag = 16;
            if (qContents is not null)
                qContents = RevCompDna(qContents);
            for (int i = 0; i < qs.Length; i++)
                qs[i] = qSize - qs[i];
            (hardLeft, hardRight) = (hardRight, hardLeft);
        }

        if (TryAttr(alignment, "flag", out var fl))
            flag |= Int(fl);
        string queryText = qContents ?? "*";
        long tStart = ts[0], qStart = qs[0];
        long pos = tStart + 1;
        var cigar = new StringBuilder();
        if (hardLeft is not null)
            cigar.Append(Str(hardLeft.Value)).Append('H');
        if (qStart > 0)
            cigar.Append(Str(qStart)).Append('S');
        string? operations = TryAttr(alignment, "operations", out var opn) ? opn.GetValue<string>() : null;
        if (m == 1)
            throw new ArgumentException("UnboundLocalError: cannot access local variable 'qEnd'");
        long qEnd = 0;
        for (int i = 1; i < m; i++)
        {
            long tEnd = ts[i];
            qEnd = qs[i];
            long tCount = tEnd - tStart, qCount = qEnd - qStart;
            if (operations is null)
            {
                if (tCount == 0)
                {
                    cigar.Append(Str(qCount)).Append('I');
                    qStart = qEnd;
                }
                else if (qCount == 0)
                {
                    cigar.Append(Str(tCount)).Append('D');
                    tStart = tEnd;
                }
                else
                {
                    if (tCount != qCount)
                        throw new ArgumentException("Unequal step sizes in alignment");
                    cigar.Append(Str(tCount)).Append('M');
                    tStart = tEnd;
                    qStart = qEnd;
                }
            }
            else
            {
                if (i - 1 >= operations.Length)
                    break; // zip stops with the shorter sequence
                char operation = operations[i - 1];
                if (tCount == 0)
                {
                    if (operation != 'I')
                        throw new ArgumentException("AssertionError");
                    cigar.Append(Str(qCount)).Append('I');
                    qStart = qEnd;
                }
                else if (qCount == 0)
                {
                    if (operation == 'N')
                        cigar.Append(Str(tCount)).Append('N');
                    else if (operation == 'D')
                        cigar.Append(Str(tCount)).Append('D');
                    else
                        throw new ArgumentException($"Unexpected operation {(int)operation}");
                    tStart = tEnd;
                }
                else
                {
                    if (tCount != qCount)
                        throw new ArgumentException("Unequal step sizes in alignment");
                    if (operation != 'M')
                        throw new ArgumentException("AssertionError");
                    cigar.Append(Str(tCount)).Append('M');
                    tStart = tEnd;
                    qStart = qEnd;
                }
            }
        }

        if (operations is not null && operations.Length < m - 1)
        {
            // qEnd is the value of the last pair that was zipped
            qEnd = operations.Length == 0 ? throw new ArgumentException("UnboundLocalError: qEnd") : qs[operations.Length];
        }

        if (qEnd < qSize)
            cigar.Append(Str(qSize - qEnd)).Append('S');
        if (hardRight is not null)
            cigar.Append(Str(hardRight.Value)).Append('H');
        string mapq = TryAttr(alignment, "mapq", out var mq) ? PyStr(mq) : "255";
        string rnext = "*";
        if (TryAttr(alignment, "rnext", out var rn))
        {
            rnext = PyStr(rn);
            if (rnext == rname)
                rnext = "=";
        }

        string pnext = TryAttr(alignment, "pnext", out var pn) ? Str(Int(pn) + 1) : "0";
        var fields = new List<string> { qName, Str(flag), rname, Str(pos), mapq, cigar.ToString(), rnext, pnext, "0", queryText, qual };
        if (options.Md)
        {
            if (queryText == "*")
                throw new ArgumentException("requested MD tag with undefined sequence");
            long tS = ts[0], qS = qs[0];
            long number = 0;
            var md = new StringBuilder();
            for (int i = 1; i < m; i++)
            {
                long tEnd = ts[i], qE = qs[i];
                long tCount = tEnd - tS, qCount = qE - qS;
                if (operations is not null && i - 1 >= operations.Length)
                    break;
                if (tCount == 0)
                {
                    qS = qE;
                }
                else if (qCount == 0)
                {
                    if (operations is null || operations[i - 1] != 'N')
                    {
                        if (number != 0)
                        {
                            md.Append(Str(number));
                            number = 0;
                        }

                        md.Append('^').Append(target.Slice(tS, tEnd) ?? throw new ArgumentException("UndefinedSequenceError"));
                    }

                    tS = tEnd;
                }
                else
                {
                    if (tCount != qCount)
                        throw new ArgumentException("Unequal step sizes in alignment");
                    string tt = target.Slice(tS, tEnd) ?? throw new ArgumentException("UndefinedSequenceError");
                    for (long k = 0; k < tt.Length && qS + k < queryText.Length; k++)
                    {
                        char tc = tt[(int)k], qc = queryText[(int)(qS + k)];
                        if (tc == qc)
                        {
                            number++;
                        }
                        else
                        {
                            md.Append(Str(number)).Append(tc);
                            number = 0;
                        }
                    }

                    tS = tEnd;
                    qS = qE;
                }
            }

            if (number != 0)
                md.Append(Str(number));
            fields.Add("MD:Z:" + md);
        }

        if (TryAttr(alignment, "score", out var sc))
        {
            var v = sc.AsValue();
            double score = v.TryGetValue<long>(out var sl) ? sl : v.GetValue<double>();
            fields.Add("AS:i:" + CFormat.Fixed(score, 0));
        }

        if (alignment.Attributes.TryGetPropertyValue("annotations", out var anns) && anns is JsonObject annotations)
        {
            foreach (var (key, value) in annotations)
            {
                string datatype, text;
                if (value is JsonValue jv && jv.TryGetValue<long>(out var iv))
                {
                    datatype = "i";
                    text = Str(iv);
                }
                else if (value is JsonValue jf && !jf.TryGetValue<string>(out _) && jf.TryGetValue<double>(out var dv))
                {
                    datatype = "f";
                    text = PyFloat(dv);
                }
                else if (value is JsonValue js && js.TryGetValue<string>(out var sv))
                {
                    datatype = sv.Length == 1 ? "A" : "Z";
                    text = sv;
                }
                else if (value is JsonObject jo && jo.ContainsKey("__bytes__"))
                {
                    datatype = "H";
                    text = string.Concat(jo["__bytes__"]!.AsArray().Select(x => Str(Int(x!))));
                }
                else if (value is JsonObject)
                {
                    throw new ArgumentException("TypeError: isinstance() arg 2 must be a type");
                }
                else
                {
                    throw new ArgumentException("UnboundLocalError: datatype");
                }

                fields.Add(key + ":" + datatype + ":" + text);
            }
        }

        return string.Join('\t', fields) + "\n";
    }

    private static JsonValue JInt(long v) => JsonValue.Create(v)!;

    private static readonly char[] SamBases = ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'];

    /// <summary><c>Bio/Align/sam.py</c> <c>AlignmentIterator</c>.</summary>
    private static CoordinateAlignmentFile ParseSam(string text)
    {
        var lines = Lines(text.Replace("\r\n", "\n").Replace('\r', '\n')).ToList();
        var metadata = new JsonObject();
        var targets = new List<BioCoordSequence>();
        int li = 0;
        for (; li < lines.Count; li++)
        {
            string line = lines[li];
            if (!line.StartsWith('@'))
                break;
            var fields = line[1..].Trim().Split('\t');
            string tag = fields[0];
            if (tag == "SQ")
            {
                var annotations = new JsonObject();
                string? description = null, rname = null;
                long? length = null;
                foreach (var field in fields.Skip(1))
                {
                    int c = field.IndexOf(':');
                    if (c < 0)
                        throw new ArgumentException("not enough values to unpack");
                    string key = field[..c], value = field[(c + 1)..];
                    if (key.Length != 2)
                        throw new ArgumentException("AssertionError");
                    switch (key)
                    {
                        case "SN": rname = value; break;
                        case "LN": length = PyInt(value); break;
                        case "AH": annotations["alternate_locus"] = value; break;
                        case "AN": annotations["names"] = new JsonArray(value.Split(',').Select(x => (JsonNode)JsonValue.Create(x)!).ToArray()); break;
                        case "AS": annotations["assembly"] = value; break;
                        case "DS": description = value; break;
                        case "M5": annotations["MD5"] = value; break;
                        case "SP": annotations["species"] = value; break;
                        case "TP":
                            if (value is not ("linear" or "circular"))
                                throw new ArgumentException("AssertionError");
                            annotations["topology"] = value;
                            break;
                        case "UR": annotations["URI"] = value; break;
                        default: annotations[key] = value; break;
                    }
                }

                if (rname is null || length is null)
                    throw new ArgumentException("UnboundLocalError");
                targets.Add(new BioCoordSequence { Id = rname, Length = length.Value, Parts = UndefinedParts(length.Value), Description = description ?? "", Annotations = annotations });
            }
            else
            {
                var values = new JsonObject();
                foreach (var field in fields.Skip(1))
                {
                    int c = field.IndexOf(':');
                    if (c < 0)
                        throw new ArgumentException("not enough values to unpack");
                    string key = field[..c];
                    if (key.Length != 2)
                        throw new ArgumentException("AssertionError");
                    values[key] = field[(c + 1)..];
                }

                if (tag == "HD")
                {
                    metadata[tag] = values;
                }
                else
                {
                    if (!metadata.ContainsKey(tag))
                        metadata[tag] = new JsonArray();
                    metadata[tag]!.AsArray().Add(values);
                }
            }
        }

        var targetIndex = new Dictionary<string, int>();
        for (int i = 0; i < targets.Count; i++)
            targetIndex[targets[i].Id!] = i;
        var alignments = new List<BioCoordAlignment>();
        for (; li < lines.Count; li++)
            alignments.Add(ParseSamLine(lines[li], targets, targetIndex));
        return new CoordinateAlignmentFile { Alignments = alignments, Metadata = metadata, Targets = targets };
    }

    private static BioCoordAlignment ParseSamLine(string line, List<BioCoordSequence> targets, Dictionary<string, int> targetIndices)
    {
        var fields = Words(line);
        if (fields.Length < 11)
            throw new ArgumentException($"line has {fields.Length} columns; expected at least 11");
        string qname = fields[0];
        long flag = PyInt(fields[1]);
        string rname = fields[2];
        long targetPos = PyInt(fields[3]) - 1;
        long mapq = PyInt(fields[4]);
        string cigar = fields[5];
        string rnext = fields[6];
        long pnext = PyInt(fields[7]) - 1;
        long tlen = PyInt(fields[8]);
        string query = fields[9];
        string qual = fields[10];
        string? md = null;
        long? score = null;
        var annotations = new JsonObject();
        foreach (var field in fields.Skip(11))
        {
            var p = field.Split(':', 3);
            if (p.Length != 3)
                throw new ArgumentException("not enough values to unpack");
            string tag = p[0], datatype = p[1], value = p[2];
            if (tag == "AS")
            {
                if (datatype != "i")
                    throw new ArgumentException("AssertionError");
                score = PyInt(value);
            }
            else if (tag == "MD")
            {
                if (datatype != "Z")
                    throw new ArgumentException("AssertionError");
                md = value;
            }
            else
            {
                JsonNode v;
                switch (datatype)
                {
                    case "i": v = JInt(PyInt(value)); break;
                    case "f": v = JsonValue.Create(PyFloatParse(value)); break;
                    case "H":
                        var bytes = new JsonArray();
                        for (int i = 0; i < value.Length; i += 2)
                            bytes.Add(JInt(PyInt(value.Substring(i, Math.Min(2, value.Length - i)))));
                        v = new JsonObject { ["__bytes__"] = bytes };
                        break;
                    case "B":
                        if (value.Length == 0)
                            throw new ArgumentException("IndexError");
                        char letter = value[0];
                        var items = value[1..].Split(',');
                        var arr = new JsonArray();
                        string dtype;
                        if ("cCsSiI".Contains(letter))
                        {
                            dtype = "int";
                            foreach (var it in items)
                                arr.Add(JInt(PyInt(it)));
                        }
                        else if (letter == 'f')
                        {
                            dtype = "float";
                            foreach (var it in items)
                                arr.Add(JsonValue.Create(PyFloatParse(it)));
                        }
                        else
                        {
                            throw new ArgumentException($"Unknown number type '{letter}' in tag '{field}'");
                        }

                        v = new JsonObject { ["__array__"] = arr, ["dtype"] = dtype };
                        break;
                    default: v = JsonValue.Create(value)!; break;
                }

                annotations[tag] = v;
            }
        }

        bool minus = (flag & 0x10) != 0;
        long? hardLeft = null, hardRight = null;
        bool storeOperations = false;
        bool unmapped = (flag & 0x4) != 0;
        BioCoordSequence? target = null;
        long[]? tcoords = null, qcoords = null;
        long queryPos = 0;
        var operations = new StringBuilder();
        if (unmapped)
        {
            queryPos = long.MinValue; // unbound
        }
        else if (md is null)
        {
            var tl = new List<long> { targetPos };
            var ql = new List<long> { 0 };
            var number = new StringBuilder();
            foreach (char letter in cigar)
            {
                long length;
                switch (letter)
                {
                    case 'M':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        queryPos += length;
                        break;
                    case '=':
                    case 'X':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        queryPos += length;
                        storeOperations = true;
                        break;
                    case 'I':
                        length = PyInt(number.ToString());
                        queryPos += length;
                        break;
                    case 'S':
                        length = PyInt(number.ToString());
                        if (queryPos == 0)
                            ql[0] += length;
                        queryPos += length;
                        number.Clear();
                        continue;
                    case 'D':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        break;
                    case 'N':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        storeOperations = true;
                        break;
                    case 'H':
                        if (queryPos == 0)
                            hardLeft = PyInt(number.ToString());
                        else
                            hardRight = PyInt(number.ToString());
                        number.Clear();
                        continue;
                    case 'P':
                        throw new ArgumentException("padding operator is not yet implemented");
                    default:
                        number.Append(letter);
                        continue;
                }

                tl.Add(targetPos);
                ql.Add(queryPos);
                operations.Append(letter);
                number.Clear();
            }

            tcoords = [.. tl];
            qcoords = [.. ql];
            if (!targetIndices.TryGetValue(rname, out int index))
            {
                if (targets.Count > 0)
                    throw new ArgumentException($"Found target {rname} missing from header");
                target = new BioCoordSequence { Id = rname, Length = -1, Parts = null, Description = "" };
            }
            else
            {
                target = targets[index];
            }
        }
        else
        {
            var tl = new List<long> { targetPos };
            var ql = new List<long> { 0 };
            string seq = query;
            var tgt = new StringBuilder();
            var starts = new List<long> { targetPos };
            long size = 0;
            var sizes = new List<long>();
            var number = new StringBuilder();
            foreach (char letter in cigar)
            {
                long length;
                switch (letter)
                {
                    case 'M':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        queryPos += length;
                        tgt.Append(Take(seq, length));
                        seq = Drop(seq, length);
                        size += length;
                        break;
                    case '=':
                    case 'X':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        queryPos += length;
                        tgt.Append(Take(seq, length));
                        seq = Drop(seq, length);
                        size += length;
                        storeOperations = true;
                        break;
                    case 'I':
                        length = PyInt(number.ToString());
                        queryPos += length;
                        seq = Drop(seq, length);
                        break;
                    case 'S':
                        length = PyInt(number.ToString());
                        if (queryPos == 0)
                            ql[0] += length;
                        queryPos += length;
                        seq = Drop(seq, length);
                        number.Clear();
                        continue;
                    case 'D':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        size += length;
                        starts.Add(targetPos);
                        sizes.Add(size);
                        size = 0;
                        break;
                    case 'N':
                        length = PyInt(number.ToString());
                        targetPos += length;
                        starts.Add(targetPos);
                        sizes.Add(size);
                        size = 0;
                        storeOperations = true;
                        break;
                    case 'H':
                        if (queryPos == 0)
                            hardLeft = PyInt(number.ToString());
                        else
                            hardRight = PyInt(number.ToString());
                        number.Clear();
                        continue;
                    case 'P':
                        throw new ArgumentException("padding operator is not yet implemented");
                    default:
                        number.Append(letter);
                        continue;
                }

                tl.Add(targetPos);
                ql.Add(queryPos);
                operations.Append(letter);
                number.Clear();
            }

            sizes.Add(size);
            tcoords = [.. tl];
            qcoords = [.. ql];
            seq = tgt.ToString();
            var outTarget = new StringBuilder();
            number.Clear();
            int mi = 0;
            while (mi < md.Length)
            {
                char letter = md[mi++];
                if (SamBases.Contains(letter))
                {
                    if (number.Length > 0)
                    {
                        long n = PyInt(number.ToString());
                        outTarget.Append(Take(seq, n));
                        seq = Drop(seq, n);
                        number.Clear();
                    }

                    outTarget.Append(letter);
                    seq = Drop(seq, 1);
                }
                else if (letter == '^')
                {
                    if (number.Length > 0)
                    {
                        long n = PyInt(number.ToString());
                        outTarget.Append(Take(seq, n));
                        seq = Drop(seq, n);
                        number.Clear();
                    }

                    bool broke = false;
                    while (mi < md.Length)
                    {
                        char l2 = md[mi++];
                        if (!SamBases.Contains(l2))
                        {
                            number.Clear().Append(l2);
                            broke = true;
                            break;
                        }

                        outTarget.Append(l2);
                    }

                    if (!broke)
                        break;
                }
                else
                {
                    number.Append(letter);
                }
            }

            if (number.Length > 0)
            {
                long n = PyInt(number.ToString());
                outTarget.Append(Take(seq, n));
            }

            seq = outTarget.ToString();
            if (!targetIndices.TryGetValue(rname, out int rIdx))
                throw new ArgumentException("KeyError: " + rname);
            var rt = targets[rIdx];
            long tlength = rt.Length;
            var data = new Dictionary<long, string>();
            int idx = 0;
            for (int i = 0; i < Math.Min(starts.Count, sizes.Count); i++)
            {
                data[starts[i]] = Take(Drop(seq, idx), sizes[i]);
                idx += (int)sizes[i];
            }

            target = new BioCoordSequence { Id = rt.Id, Length = tlength, Parts = SeqParts(data, tlength), Description = rt.Description, Annotations = (JsonObject)rt.Annotations.DeepClone() };
        }

        if (!unmapped && minus)
        {
            for (int i = 0; i < qcoords!.Length; i++)
                qcoords[i] = queryPos - qcoords[i];
        }

        BioCoordSequence querySeq;
        if (query == "*")
        {
            if (queryPos == long.MinValue)
                throw new ArgumentException("UnboundLocalError: query_pos");
            querySeq = new BioCoordSequence { Id = qname, Length = queryPos, Parts = UndefinedParts(queryPos) };
        }
        else
        {
            string sequence = query;
            if (!unmapped)
            {
                if (query.Length != queryPos)
                    throw new ArgumentException("AssertionError");
                if (minus)
                    sequence = RevCompDna(sequence);
            }

            querySeq = new BioCoordSequence { Id = qname, Length = sequence.Length, Parts = SeqParts([new KeyValuePair<long, string>(0, sequence)], sequence.Length) };
        }

        if (minus)
            (hardLeft, hardRight) = (hardRight, hardLeft);
        var qAnn = new JsonObject();
        if (hardLeft is not null)
            qAnn["hard_clip_left"] = JInt(hardLeft.Value);
        if (hardRight is not null)
            qAnn["hard_clip_right"] = JInt(hardRight.Value);
        var letterAnn = new JsonObject();
        if (qual != "*")
            letterAnn["phred_quality"] = new JsonArray(qual.Select(c => (JsonNode)JInt(c - 33)).ToArray());
        querySeq = new BioCoordSequence { Id = querySeq.Id, Length = querySeq.Length, Parts = querySeq.Parts, Description = "", Annotations = qAnn, LetterAnnotations = letterAnn };
        var attrs = new JsonObject { ["flag"] = JInt(flag) };
        if (mapq != 255)
            attrs["mapq"] = JInt(mapq);
        if (rnext == "=")
            attrs["rnext"] = rname;
        else if (rnext != "*")
            attrs["rnext"] = rnext;
        if (pnext >= 0)
            attrs["pnext"] = JInt(pnext);
        if (tlen != 0)
            attrs["tlen"] = JInt(tlen);
        if (score is not null)
            attrs["score"] = JInt(score.Value);
        if (annotations.Count > 0)
            attrs["annotations"] = annotations;
        if (hardLeft is not null)
            attrs["hard_clip_left"] = JInt(hardLeft.Value);
        if (hardRight is not null)
            attrs["hard_clip_right"] = JInt(hardRight.Value);
        if (storeOperations)
            attrs["operations"] = operations.ToString();
        if (unmapped)
        {
            return new BioCoordAlignment
            {
                Sequences = [new BioCoordSequence { Id = null, Length = -1 }, querySeq],
                Coordinates = [],
                NoCoordinates = true,
                Attributes = attrs,
            };
        }

        return new BioCoordAlignment { Sequences = [target!, querySeq], Coordinates = [tcoords!, qcoords!], Attributes = attrs };
    }

    private static string Take(string s, long n)
    {
        if (n <= 0)
            return "";
        return n >= s.Length ? s : s[..(int)n];
    }

    private static string Drop(string s, long n)
    {
        if (n <= 0)
            return s;
        return n >= s.Length ? "" : s[(int)n..];
    }

    /// <summary>Python's <c>float(text)</c> for the forms that occur in SAM tags (decimal, exponent, <c>nan</c>, <c>inf</c>).</summary>
    private static double PyFloatParse(string s)
    {
        string t = s.Trim().ToLowerInvariant();
        switch (t)
        {
            case "nan": case "+nan": case "-nan": return double.NaN;
            case "inf": case "+inf": case "infinity": case "+infinity": return double.PositiveInfinity;
            case "-inf": case "-infinity": return double.NegativeInfinity;
        }

        if (t.Length == 0 || t.Contains('_') || !double.TryParse(t, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double d))
            throw new ArgumentException($"could not convert string to float: '{s}'");
        return d;
    }
}
