using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    // ------------------------------------------------------------------ HHR

    /// <summary>Python's <c>"%-Ns" % text</c> after <c>" " * start</c>: <c>start</c> spaces, then the text left-justified in <c>width</c> columns.</summary>
    private static string HhrPad(long start, long width, string text)
    {
        string s = new string(' ', (int)start) + text;
        int w = (int)Math.Max(width, 0);
        return text.Length >= w ? s : s + new string(' ', w - text.Length);
    }

    private static string[] SixWords(string line)
    {
        var w = Words(line);
        if (w.Length != 6)
            throw new ArgumentException("ValueError: unpacking");
        return w;
    }

    private static long HhrTotal(string total)
    {
        if (!total.StartsWith('(') || !total.EndsWith(')'))
            throw new ArgumentException("AssertionError");
        return PyInt(total[1..^1]);
    }

    /// <summary><c>Bio/Align/hhr.py</c> <c>AlignmentIterator</c>.</summary>
    private static CoordinateAlignmentFile ParseHhr(string text)
    {
        var stream = new LineStream(text);
        var metadata = new JsonObject();
        string? queryName = null;
        string? l;
        while ((l = stream.Next()) is not null)
        {
            string line = l.Trim();
            if (line.Length == 0)
                break;
            var (key, value) = SplitWs1(line);
            switch (key)
            {
                case "Query": queryName = value; break;
                case "Match_columns": metadata[key] = PyInt(value); break;
                case "No_of_seqs":
                    var p = value.Split(" out of ");
                    if (p.Length != 2)
                        throw new ArgumentException("ValueError: unpacking");
                    metadata[key] = new JsonArray(JsonValue.Create(PyInt(p[0])), JsonValue.Create(PyInt(p[1])));
                    break;
                case "Neff":
                case "Template_Neff":
                    metadata[key] = PyFloatParse(value);
                    break;
                case "Searched_HMMs": metadata[key] = PyInt(value); break;
                case "Date": metadata["Rundate"] = value; break;
                case "Command": metadata["Command line"] = value; break;
                default: throw new ArgumentException($"Unknown key '{key}'");
            }
        }

        string head = stream.ReadLine();
        if (head.Length == 0)
            throw new ArgumentException("Truncated file.");
        string[] expected = ["No", "Hit", "Prob", "E-value", "P-value", "Score", "SS", "Cols", "Query", "HMM", "Template", "HMM"];
        if (!Words(head).SequenceEqual(expected))
            throw new ArgumentException("AssertionError");
        long counter = 0;
        while ((l = stream.Next()) is not null)
        {
            if (l.Trim().Length == 0)
                break;
            counter++;
            var (word, _) = SplitWs1(l);
            if (PyInt(word) != counter)
                throw new ArgumentException("AssertionError");
        }

        long total = counter;
        long seen = 0;
        var alignments = new List<BioCoordAlignment>();
        while (true)
        {
            var (alignment, stop) = ReadHhrAlignment(stream, metadata, queryName, total, ref seen);
            if (alignment is null)
            {
                if (stop && seen != total)
                    throw new ArgumentException($"Expected {total} alignments, found {seen}");
                break;
            }

            alignments.Add(alignment);
        }

        return new CoordinateAlignmentFile { Alignments = alignments, Metadata = metadata };
    }

    private static (BioCoordAlignment? Alignment, bool AtEnd) ReadHhrAlignment(LineStream stream, JsonObject metadata, string? queryName, long length, ref long counterField)
    {
        long queryStart = -1, targetStart = -1, queryLength = 0, targetLength = 0;
        bool haveQueryStart = false, haveTargetStart = false;
        var querySequence = new StringBuilder();
        var queryConsensus = new StringBuilder();
        var querySsPred = new StringBuilder();
        var targetSequence = new StringBuilder();
        var targetConsensus = new StringBuilder();
        var targetSsPred = new StringBuilder();
        var targetSsDssp = new StringBuilder();
        var columnScore = new StringBuilder();
        var confidence = new StringBuilder();
        string consensus = "";
        string hmmName = "", hmmDescription = "", targetName = "";
        var alignmentAnnotations = new JsonObject();
        string? l;

        BioCoordAlignment? Create()
        {
            int n = targetSequence.Length;
            if (querySequence.Length != n)
                throw new ArgumentException("AssertionError");
            if (n == 0)
                return null;
            var (seqData, co) = ParsePrintedAlignment([targetSequence.ToString(), querySequence.ToString()]);
            for (int c = 0; c < co[0].Length; c++)
            {
                co[0][c] += targetStart;
                co[1][c] += queryStart;
            }

            var query = new BioCoordSequence
            {
                Id = queryName ?? throw new ArgumentException("AttributeError: query_name"),
                Length = queryLength,
                Parts = SeqParts([new KeyValuePair<long, string>(queryStart, seqData[1])], queryLength),
                Description = BiopythonUnknownDescription,
            };
            var tLetters = new JsonObject();
            string tfmt(string s) => HhrPad(targetStart, targetLength - targetStart, s);
            tLetters["Consensus"] = tfmt(targetConsensus.ToString().Replace("-", ""));
            tLetters["ss_pred"] = tfmt(targetSsPred.ToString().Replace("-", ""));
            tLetters["ss_dssp"] = tfmt(targetSsDssp.ToString().Replace("-", ""));
            var conf = new StringBuilder();
            string tseq = targetSequence.ToString(), cs = confidence.ToString();
            for (int i = 0; i < Math.Min(tseq.Length, cs.Length); i++)
            {
                if (tseq[i] != '-')
                    conf.Append(cs[i]);
            }

            string alignmentConfidence = tfmt(conf.ToString());
            var qLetters = new JsonObject();
            string qfmt(string s) => HhrPad(queryStart, queryLength - queryStart, s);
            if (queryConsensus.Length > 0)
                qLetters["Consensus"] = qfmt(queryConsensus.ToString().Replace("-", ""));
            if (querySsPred.Length > 0)
                qLetters["ss_pred"] = qfmt(querySsPred.ToString().Replace("-", ""));
            foreach (var v in tLetters.Select(kv => kv.Value))
            {
                if (v!.GetValue<string>().Length != targetLength)
                    throw new ArgumentException("Unexpected letter annotations length");
            }

            foreach (var v in qLetters.Select(kv => kv.Value))
            {
                if (v!.GetValue<string>().Length != queryLength)
                    throw new ArgumentException("Unexpected letter annotations length");
            }

            var target = new BioCoordSequence
            {
                Id = targetName,
                Length = targetLength,
                Parts = SeqParts([new KeyValuePair<long, string>(targetStart, seqData[0])], targetLength),
                Description = BiopythonUnknownDescription,
                Annotations = new JsonObject { ["hmm_name"] = hmmName, ["hmm_description"] = hmmDescription },
                LetterAnnotations = tLetters,
            };
            query = new BioCoordSequence { Id = query.Id, Length = query.Length, Parts = query.Parts, Description = query.Description, LetterAnnotations = qLetters };
            var attrs = new JsonObject
            {
                ["annotations"] = alignmentAnnotations.DeepClone(),
                ["column_annotations"] = new JsonObject { ["column score"] = columnScore.ToString(), ["Confidence"] = alignmentConfidence },
            };
            return new BioCoordAlignment { Sequences = [target, query], Coordinates = co, Attributes = attrs };
        }

        while ((l = stream.Next()) is not null)
        {
            string line = l.TrimEnd();
            if (line.Length == 0)
                continue;
            if (line.StartsWith('>'))
            {
                var (n1, n2) = SplitWs1(line[1..]);
                hmmName = n1;
                hmmDescription = n2;
                var words = Words(stream.ReadLine());
                alignmentAnnotations = [];
                foreach (var word in words)
                {
                    var kv = word.Split('=');
                    if (kv.Length != 2)
                        throw new ArgumentException("ValueError: unpacking");
                    if (kv[0] == "Aligned_cols")
                        continue;
                    string v = kv[1];
                    if (kv[0] == "Identities")
                        v = v.TrimEnd('%');
                    alignmentAnnotations[kv[0]] = PyFloatParse(v);
                }
            }
            else if (line == "Done!")
            {
                if (stream.ReadLine().Length > 0)
                    throw new ArgumentException("Found additional data after 'Done!'; corrupt file?");
            }
            else if (line.StartsWith(' '))
            {
                columnScore.Append(line.Trim());
            }
            else if (line.StartsWith("No ", StringComparison.Ordinal))
            {
                long counter = counterField;
                counterField++;
                var w = Words(line);
                if (w.Length != 2)
                    throw new ArgumentException("ValueError: unpacking");
                if (PyInt(w[1]) != counterField)
                    throw new ArgumentException("AssertionError");
                if (counterField > length)
                    throw new ArgumentException($"Expected {length} alignments, found {counterField}");
                if (counter > 0)
                    return (Create(), false);
            }
            else if (line.StartsWith("Confidence", StringComparison.Ordinal))
            {
                int k = consensus.Length;
                confidence.Append(k == 0 || k >= line.Length ? line : line[^k..]);
            }
            else if (line.StartsWith("Q ss_pred ", StringComparison.Ordinal))
            {
                querySsPred.Append(RSplitLast(line));
            }
            else if (line.StartsWith("Q Consensus ", StringComparison.Ordinal))
            {
                var w = SixWords(line);
                _ = PyInt(w[2]);
                _ = PyInt(w[4]);
                HhrTotal(w[5]);
                consensus = w[3];
                queryConsensus.Append(consensus);
            }
            else if (line.StartsWith("Q ", StringComparison.Ordinal))
            {
                var w = SixWords(line);
                if (queryName is null || !queryName.StartsWith(w[1], StringComparison.Ordinal))
                    throw new ArgumentException("AssertionError");
                long start = PyInt(w[2]) - 1;
                _ = PyInt(w[4]);
                queryLength = HhrTotal(w[5]);
                if (!metadata.TryGetPropertyValue("Match_columns", out var mc) || queryLength != Int(mc!))
                    throw new ArgumentException("AssertionError");
                if (!haveQueryStart)
                {
                    queryStart = start;
                    haveQueryStart = true;
                }

                querySequence.Append(w[3]);
            }
            else if (line.StartsWith("T ss_pred ", StringComparison.Ordinal))
            {
                targetSsPred.Append(RSplitLast(line));
            }
            else if (line.StartsWith("T ss_dssp ", StringComparison.Ordinal))
            {
                targetSsDssp.Append(RSplitLast(line));
            }
            else if (line.StartsWith("T Consensus ", StringComparison.Ordinal))
            {
                var w = SixWords(line);
                _ = PyInt(w[2]);
                _ = PyInt(w[4]);
                HhrTotal(w[5]);
                consensus = w[3];
                targetConsensus.Append(consensus);
            }
            else if (line.StartsWith("T ", StringComparison.Ordinal))
            {
                var w = SixWords(line);
                targetName = w[1];
                long start = PyInt(w[2]) - 1;
                _ = PyInt(w[4]);
                targetLength = HhrTotal(w[5]);
                if (!haveTargetStart)
                {
                    targetStart = start;
                    haveTargetStart = true;
                }

                targetSequence.Append(w[3]);
            }
            else
            {
                throw new ArgumentException($"Failed to parse line '{(line.Length > 30 ? line[..30] : line)}...'");
            }
        }

        var alignment = Create();
        if (alignment is null && counterField != length)
            throw new ArgumentException($"Expected {length} alignments, found {counterField}");
        return (alignment, true);
    }

    /// <summary><c>line.rsplit(None, 1)[1]</c>.</summary>
    private static string RSplitLast(string line)
    {
        string t = line.TrimEnd();
        int i = t.Length;
        while (i > 0 && !char.IsWhiteSpace(t[i - 1]))
            i--;
        if (i == 0 || t[..i].Trim().Length == 0)
            throw new ArgumentException("not enough values to unpack (expected 2, got 1)");
        return t[i..];
    }

    // ------------------------------------------------------------------ tabular

    private sealed class TabularState
    {
        public string[]? Fields;
        public string? QueryId;
        public bool HasQuery;
        public string? QueryDescription;
        public long? QuerySize;
        public string FinalPrefix = "";
        public JsonObject Metadata = [];
    }

    private static readonly string[] BlastPrograms = ["BLASTN", "BLASTP", "BLASTX", "TBLASTN", "TBLASTX", "DELTABLAST", "PSIBLAST", "RPSBLAST", "RPSTBLASTN"];

    private static void TabularHeader(LineStream stream, string line, TabularState st)
    {
        var metadata = new JsonObject();
        bool blast = false;
        string program = "", version = "";
        var words = line[2..].Split((char[]?)null, 2, StringSplitOptions.RemoveEmptyEntries);
        if (words.Length == 2 && BlastPrograms.Contains(words[0]))
        {
            blast = true;
            program = words[0];
            version = line[2..].TrimStart()[words[0].Length..].TrimStart();
        }

        if (!blast)
        {
            metadata["Command line"] = line[2..];
            string l2 = stream.ReadLine();
            if (!l2.StartsWith("# ", StringComparison.Ordinal))
                throw new ArgumentException("AssertionError");
            var (p, v) = SplitWs1(l2[2..].TrimEnd());
            metadata["Program"] = p;
            metadata["Version"] = v;
            st.FinalPrefix = "# FASTA processed ";
        }
        else
        {
            metadata["Program"] = program;
            metadata["Version"] = version;
            st.FinalPrefix = "# BLAST processed ";
        }

        string? l;
        while ((l = stream.Next()) is not null)
        {
            string t = l.Trim();
            if (!t.StartsWith("# ", StringComparison.Ordinal))
                throw new ArgumentException("AssertionError");
            var pr = t[2..].Split(": ");
            if (pr.Length != 2)
            {
                const string suffix = " hits found";
                if (!t.EndsWith(suffix, StringComparison.Ordinal))
                    throw new ArgumentException("AssertionError");
                _ = PyInt(t[2..^suffix.Length]);
                break;
            }

            string prefix = pr[0], value = pr[1];
            switch (prefix)
            {
                case "Query":
                    string queryLine;
                    if (metadata["Program"]!.GetValue<string>() == "FASTA")
                    {
                        int dash = value.LastIndexOf(" - ", StringComparison.Ordinal);
                        if (dash < 0)
                            throw new ArgumentException("ValueError: unpacking");
                        queryLine = value[..dash];
                        var sz = Words(value[(dash + 3)..]);
                        if (sz.Length != 2)
                            throw new ArgumentException("ValueError: unpacking");
                        st.QuerySize = PyInt(sz[0]);
                        if (sz[1] is not ("nt" or "aa"))
                            throw new ArgumentException("AssertionError");
                    }
                    else
                    {
                        queryLine = value;
                        st.QuerySize = null;
                    }

                    st.HasQuery = true;
                    var qw = queryLine.Split((char[]?)null, 2, StringSplitOptions.RemoveEmptyEntries);
                    if (qw.Length == 2)
                    {
                        st.QueryId = qw[0];
                        st.QueryDescription = queryLine.TrimStart()[qw[0].Length..].TrimStart();
                    }
                    else
                    {
                        st.QueryId = queryLine.Trim();
                        st.QueryDescription = null;
                    }

                    break;
                case "Database": metadata["Database"] = value; break;
                case "Fields": st.Fields = value.Split(", "); break;
                case "RID": metadata["RID"] = value; break;
            }
        }

        st.Metadata = metadata;
    }

    private static CoordinateAlignmentFile ParseTabular(string text)
    {
        var stream = new LineStream(text);
        string first = stream.ReadLine();
        if (first.Length == 0)
            throw new ArgumentException("Empty file.");
        if (!first.StartsWith("# ", StringComparison.Ordinal))
            throw new ArgumentException("Missing header.");
        var st = new TabularState();
        TabularHeader(stream, first.TrimEnd(), st);
        var alignments = new List<BioCoordAlignment>();
        while (true)
        {
            var a = ReadTabularAlignment(stream, st);
            if (a is null)
                break;
            alignments.Add(a);
        }

        return new CoordinateAlignmentFile { Alignments = alignments, Metadata = st.Metadata };
    }

    private static BioCoordAlignment? ReadTabularAlignment(LineStream stream, TabularState st)
    {
        string? line = null;
        string? l;
        while ((l = stream.Next()) is not null)
        {
            line = l.TrimEnd();
            if (line.StartsWith("# ", StringComparison.Ordinal))
            {
                if (line.StartsWith(st.FinalPrefix, StringComparison.Ordinal) && line.EndsWith(" queries", StringComparison.Ordinal))
                {
                    if (st.Fields is null || !st.HasQuery)
                        throw new ArgumentException("AttributeError: _fields");
                    return null;
                }

                TabularHeader(stream, line, st);
            }
            else
            {
                break;
            }
        }

        if (line is null)
            throw new ArgumentException("UnboundLocalError: line");
        if (st.Fields is null)
            throw new ArgumentException("AttributeError: _fields");
        long? alignmentLength = null;
        string? btop = null;
        string? cigar = null;
        long? score = null;
        string? queryId = null, targetId = null;
        long? queryStart = null, queryEnd = null, targetStart = null, targetEnd = null;
        string? querySequence = null, targetSequence = null;
        long? targetLength = null;
        long[][]? coordinates = null;
        long? querySize = st.QuerySize;
        var columns = line.Split('\t');
        if (columns.Length != st.Fields.Length)
            throw new ArgumentException("AssertionError");
        var annotations = new JsonObject();
        var queryAnnotations = new JsonObject();
        var targetAnnotations = new JsonObject();
        for (int i = 0; i < columns.Length; i++)
        {
            string column = columns[i], field = st.Fields[i];
            switch (field)
            {
                case "query id":
                    queryId = column;
                    if (st.QueryId is not null && queryId != st.QueryId)
                        throw new ArgumentException("AssertionError");
                    break;
                case "subject id": targetId = column; break;
                case "% identity": annotations[field] = PyFloatParse(column); break;
                case "alignment length": alignmentLength = PyInt(column); break;
                case "mismatches":
                case "gap opens":
                    annotations[field] = PyInt(column);
                    break;
                case "q. start": queryStart = PyInt(column); break;
                case "q. end": queryEnd = PyInt(column); break;
                case "s. start": targetStart = PyInt(column); break;
                case "s. end": targetEnd = PyInt(column); break;
                case "evalue": annotations["evalue"] = PyFloatParse(column); break;
                case "bit score": annotations["bit score"] = PyFloatParse(column); break;
                case "BTOP":
                    btop = column;
                    coordinates = ParseBtop(column);
                    break;
                case "aln_code":
                    cigar = column;
                    coordinates = ParseTabularCigar(column);
                    break;
                case "query gi": queryAnnotations["gi"] = column; break;
                case "query acc.": queryAnnotations["acc."] = column; break;
                case "query acc.ver":
                    queryAnnotations["acc.ver"] = column;
                    queryId ??= column;
                    break;
                case "query length":
                    long ql = PyInt(column);
                    if (querySize is null)
                        querySize = ql;
                    else if (querySize != ql)
                        throw new ArgumentException("AssertionError");
                    break;
                case "subject ids": targetAnnotations["ids"] = column; break;
                case "subject gi": targetAnnotations["gi"] = column; break;
                case "subject gis": targetAnnotations["gis"] = column; break;
                case "subject acc.": targetAnnotations["acc."] = column; break;
                case "subject accs.": targetAnnotations["accs."] = column; break;
                case "subject tax ids": targetAnnotations["tax ids"] = column; break;
                case "subject sci names": targetAnnotations["sci names"] = column; break;
                case "subject com names": targetAnnotations["com names"] = column; break;
                case "subject blast names": targetAnnotations["blast names"] = column; break;
                case "subject super kingdoms": targetAnnotations["super kingdoms"] = column; break;
                case "subject title": targetAnnotations["title"] = column; break;
                case "subject titles": targetAnnotations["titles"] = column; break;
                case "subject strand": targetAnnotations["strand"] = column; break;
                case "% subject coverage": targetAnnotations["% coverage"] = PyFloatParse(column); break;
                case "subject acc.ver":
                    targetAnnotations["acc.ver"] = column;
                    targetId ??= column;
                    break;
                case "subject length": targetLength = PyInt(column); break;
                case "query seq": querySequence = column; break;
                case "subject seq": targetSequence = column; break;
                case "score": score = PyInt(column); break;
                case "identical": annotations[field] = PyInt(column); break;
                case "positives":
                case "gaps":
                    annotations[field] = PyInt(column);
                    break;
                case "% positives":
                case "% hsp coverage":
                    annotations[field] = PyFloatParse(column);
                    break;
                case "query/sbjct frames": annotations[field] = column; break;
                case "query frame": queryAnnotations["frame"] = column; break;
                case "sbjct frame": targetAnnotations["frame"] = column; break;
                default: throw new ArgumentException($"Unexpected field '{field}'");
            }
        }

        _ = btop;
        _ = cigar;
        string program = st.Metadata["Program"]!.GetValue<string>();
        if (coordinates is null && alignmentLength is not null)
            annotations["alignment length"] = alignmentLength.Value;
        if (queryStart is not null && queryEnd is not null)
        {
            if (queryStart < queryEnd)
                queryStart--;
            else
                queryEnd--;
        }

        if (targetStart is not null && targetEnd is not null)
        {
            if (targetStart < targetEnd)
                targetStart--;
            else
                targetEnd--;
        }

        if (coordinates is null || program is "BLASTX" or "TBLASTX")
        {
            if (queryStart is not null)
                queryAnnotations["start"] = queryStart.Value;
            if (queryEnd is not null)
                queryAnnotations["end"] = queryEnd.Value;
        }
        else
        {
            if (queryStart is null || queryEnd is null)
                throw new ArgumentException("TypeError: '<' not supported between instances of 'NoneType' and 'NoneType'");
            if (queryStart < queryEnd)
            {
                for (int c = 0; c < coordinates[1].Length; c++)
                    coordinates[1][c] += queryStart.Value;
            }
            else
            {
                for (int c = 0; c < coordinates[1].Length; c++)
                    coordinates[1][c] = queryStart.Value - coordinates[1][c];
            }
        }

        if (coordinates is null || program is "TBLASTN" or "TBLASTX")
        {
            if (targetStart is not null)
                targetAnnotations["start"] = targetStart.Value;
            if (targetEnd is not null)
                targetAnnotations["end"] = targetEnd.Value;
        }
        else
        {
            if (targetStart is null)
                throw new ArgumentException("TypeError");
            for (int c = 0; c < coordinates[0].Length; c++)
                coordinates[0][c] += targetStart.Value;
        }

        BioCoordSequence querySeq;
        if (querySequence is null)
        {
            querySeq = querySize is null
                ? new BioCoordSequence { Id = queryId, Length = -1 }
                : new BioCoordSequence { Id = queryId, Length = querySize.Value, Parts = UndefinedParts(querySize.Value) };
        }
        else
        {
            querySequence = querySequence.Replace("-", "");
            if (program == "TBLASTN")
            {
                if (queryStart is null || queryEnd is null || querySequence.Length != queryEnd - queryStart)
                    throw new ArgumentException("AssertionError");
                if (querySize is null)
                    throw new ArgumentException("length must not be None if data is a dictionary");
                querySeq = new BioCoordSequence { Id = queryId, Length = querySize.Value, Parts = SeqParts([new KeyValuePair<long, string>(queryStart.Value, querySequence)], querySize.Value) };
            }
            else if (program == "TBLASTX")
            {
                queryAnnotations["start"] = queryStart is null ? null : JsonValue.Create(queryStart.Value);
                queryAnnotations["end"] = queryEnd is null ? null : JsonValue.Create(queryEnd.Value);
                querySeq = new BioCoordSequence { Id = queryId, Length = querySequence.Length, Parts = SeqParts([new KeyValuePair<long, string>(0, querySequence)], querySequence.Length) };
            }
            else
            {
                throw new ArgumentException($"Unknown program {program}");
            }
        }

        querySeq = new BioCoordSequence
        {
            Id = querySeq.Id,
            Length = querySeq.Length,
            Parts = querySeq.Parts,
            Description = st.QueryDescription ?? BiopythonUnknownDescription,
            Annotations = queryAnnotations,
        };
        BioCoordSequence targetSeq;
        if (program is "TBLASTN" or "TBLASTX")
        {
            targetAnnotations["length"] = targetLength is null ? null : JsonValue.Create(targetLength.Value);
            if (targetSequence is null)
            {
                targetSeq = new BioCoordSequence { Id = targetId, Length = -1 };
            }
            else
            {
                targetSequence = targetSequence.Replace("-", "");
                targetSeq = new BioCoordSequence { Id = targetId, Length = targetSequence.Length, Parts = SeqParts([new KeyValuePair<long, string>(0, targetSequence)], targetSequence.Length) };
            }
        }
        else if (targetSequence is null)
        {
            targetSeq = targetEnd is null
                ? new BioCoordSequence { Id = targetId, Length = -1 }
                : new BioCoordSequence { Id = targetId, Length = targetEnd.Value, Parts = UndefinedParts(targetEnd.Value) };
        }
        else
        {
            targetSequence = targetSequence.Replace("-", "");
            if (targetStart is null || targetEnd is null)
                throw new ArgumentException("UnboundLocalError: target_seq");
            if (targetSequence.Length != targetEnd - targetStart)
                throw new ArgumentException("AssertionError");
            targetSeq = new BioCoordSequence { Id = targetId, Length = targetEnd.Value, Parts = SeqParts([new KeyValuePair<long, string>(targetStart.Value, targetSequence)], targetEnd.Value) };
        }

        targetSeq = new BioCoordSequence { Id = targetSeq.Id, Length = targetSeq.Length, Parts = targetSeq.Parts, Description = BiopythonUnknownDescription, Annotations = targetAnnotations };
        querySeq = new BioCoordSequence { Id = querySeq.Id, Length = querySeq.Length, Parts = querySeq.Parts, Description = querySeq.Description, Annotations = querySeq.Annotations };
        var attrs = new JsonObject { ["annotations"] = annotations };
        if (score is not null)
            attrs["score"] = score.Value;
        var records = new List<BioCoordSequence> { targetSeq, querySeq };
        if (coordinates is null)
            return AlignmentWithoutCoordinates(records, attrs);
        return new BioCoordAlignment { Sequences = records, Coordinates = coordinates, Attributes = attrs };
    }

    private enum BtopState
    {
        None,
        Match,
        QueryGap,
        TargetGap,
    }

    /// <summary><c>AlignmentIterator.parse_btop</c>.</summary>
    private static long[][] ParseBtop(string btop)
    {
        var t = new List<long> { 0 };
        var q = new List<long> { 0 };
        var state = BtopState.None;
        foreach (Match m in Regex.Matches(btop, "([A-Z-*]{2}|[0-9]+)"))
        {
            string token = m.Value;
            if (token.StartsWith('-'))
            {
                if (state != BtopState.QueryGap)
                {
                    t.Add(t[^1]);
                    q.Add(q[^1]);
                    state = BtopState.QueryGap;
                }

                t[^1]++;
            }
            else if (token.EndsWith('-'))
            {
                if (state != BtopState.TargetGap)
                {
                    t.Add(t[^1]);
                    q.Add(q[^1]);
                    state = BtopState.TargetGap;
                }

                q[^1]++;
            }
            else
            {
                long length = long.TryParse(token, NumberStyles.None, CultureInfo.InvariantCulture, out long v) ? v : 1;
                if (state == BtopState.Match)
                {
                    t[^1] += length;
                    q[^1] += length;
                }
                else
                {
                    t.Add(t[^1] + length);
                    q.Add(q[^1] + length);
                    state = BtopState.Match;
                }
            }
        }

        return [[.. t], [.. q]];
    }

    /// <summary><c>AlignmentIterator.parse_cigar</c> (the aln_code of FASTA's tabular output).</summary>
    private static long[][] ParseTabularCigar(string cigar)
    {
        var t = new List<long> { 0 };
        var q = new List<long> { 0 };
        long tc = 0, qc = 0;
        var tokens = Regex.Matches(cigar, "(M|D|I|[0-9]+)").Select(m => m.Value).ToList();
        for (int i = 0; i + 1 < tokens.Count; i += 2)
        {
            long length = PyInt(tokens[i]);
            switch (tokens[i + 1])
            {
                case "M":
                    tc += length;
                    qc += length;
                    break;
                case "I": tc += length; break;
                case "D": qc += length; break;
            }

            t.Add(tc);
            q.Add(qc);
        }

        return [[.. t], [.. q]];
    }
}
