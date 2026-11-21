using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary><c>Bio/Align/exonerate.py</c> <c>AlignmentWriter.write_header</c>.</summary>
    private static string ExonerateHeader(CoordinateWriteOptions options)
    {
        if (options.ExonerateFormat is not ("vulgar" or "cigar"))
            throw new ArgumentException($"argument fmt should be 'vulgar' or 'cigar' (received {options.ExonerateFormat})");
        string commandLine = "", hostname = "";
        if (options.Metadata is not null)
        {
            if (options.Metadata.TryGetPropertyValue("Command line", out var c) && c is not null)
                commandLine = PyValueStr(c);
            if (options.Metadata.TryGetPropertyValue("Hostname", out var h) && h is not null)
                hostname = PyValueStr(h);
        }

        return $"Command line: [{commandLine}]\nHostname: [{hostname}]\n";
    }

    private static string? MoleculeType(BioCoordSequence s) => s.Annotations.TryGetPropertyValue("molecule_type", out var n) && n is JsonValue v && v.TryGetValue<string>(out var t) ? t : null;

    /// <summary><c>Bio/Align/exonerate.py</c> <c>_format_alignment_cigar</c> / <c>_format_alignment_vulgar</c>.</summary>
    private static string FormatExonerate(BioCoordAlignment alignment, CoordinateWriteOptions options)
    {
        if (options.ExonerateFormat is not ("vulgar" or "cigar"))
            throw new ArgumentException($"argument fmt should be 'vulgar' or 'cigar' (received {options.ExonerateFormat})");
        bool vulgar = options.ExonerateFormat == "vulgar";
        var co = alignment.Coordinates;
        if (alignment.NoCoordinates || co.Length != 2 || alignment.Sequences.Count != 2)
            throw new ArgumentException("Expected an alignment of two sequences (target and query).");
        var target = alignment.Sequences[0];
        var query = alignment.Sequences[1];
        int m = co[0].Length;
        long targetStart = co[0][0], targetEnd = co[0][^1], queryStart = co[1][0], queryEnd = co[1][^1];
        var tsteps = new long[m - 1];
        var qsteps = new long[m - 1];
        for (int i = 0; i < m - 1; i++)
        {
            tsteps[i] = co[0][i + 1] - co[0][i];
            qsteps[i] = co[1][i + 1] - co[1][i];
        }

        string queryId = query.Id ?? "None", targetId = target.Id ?? "None";
        string? tmt = MoleculeType(target), qmt = MoleculeType(query);
        string targetStrand, queryStrand;
        if (tmt == "protein")
        {
            targetStrand = ".";
        }
        else if (targetStart <= targetEnd)
        {
            targetStrand = "+";
        }
        else
        {
            targetStrand = "-";
            for (int i = 0; i < tsteps.Length; i++)
                tsteps[i] = -tsteps[i];
        }

        if (qmt == "protein")
        {
            queryStrand = ".";
        }
        else if (queryStart <= queryEnd)
        {
            queryStrand = "+";
        }
        else
        {
            queryStrand = "-";
            for (int i = 0; i < qsteps.Length; i++)
                qsteps[i] = -qsteps[i];
        }

        if (!TryAttr(alignment, "score", out var scoreNode))
            throw new ArgumentException("AttributeError: score");
        string score = G(scoreNode);
        var words = new List<string> { vulgar ? "vulgar:" : "cigar:", queryId, Str(queryStart), Str(queryEnd), queryStrand, targetId, Str(targetStart), Str(targetEnd), targetStrand, score };
        string? operations = TryAttr(alignment, "operations", out var opn) ? opn.GetValue<string>() : null;
        bool tProt = tmt == "protein", qProt = qmt == "protein";
        int steps = m - 1;
        if (operations is null)
        {
            for (int i = 0; i < steps; i++)
            {
                long ts = tsteps[i], qs = qsteps[i];
                if (vulgar)
                {
                    string op;
                    if (ts == qs)
                        op = "M";
                    else if (qs == 0)
                        op = "G";
                    else if (ts == 0)
                        op = "G";
                    else if ((qProt && !tProt) || (!qProt && tProt))
                        op = "M";
                    else
                        throw new ArgumentException("Both target and query step are zero");
                    words.Add(op);
                    words.Add(Str(qs));
                    words.Add(Str(ts));
                }
                else
                {
                    var (op, step) = CigarStep(ts, qs, tProt, qProt, tmt, qmt);
                    words.Add(op);
                    words.Add(Str(step));
                }
            }
        }
        else if (vulgar)
        {
            int n = Math.Min(operations.Length, steps);
            int i = 0;
            while (i < n)
            {
                long ts = tsteps[i], qs = qsteps[i];
                string operation = operations[i].ToString();
                switch (operation)
                {
                    case "M":
                        if (ts != qs)
                            CheckCodonSteps(ts, qs, tProt, qProt);
                        break;
                    case "5":
                    case "3":
                        if (!(ts == 2 || qs == 2))
                            throw new ArgumentException("AssertionError");
                        break;
                    case "N":
                        operation = "I";
                        if (!(qs == 0 || ts == 0))
                            throw new ArgumentException("AssertionError");
                        break;
                    case "C":
                        if (ts != qs)
                            throw new ArgumentException("AssertionError");
                        break;
                    case "D":
                        if (qs != 0)
                            throw new ArgumentException("AssertionError");
                        operation = "G";
                        break;
                    case "I":
                        if (ts != 0)
                            throw new ArgumentException("AssertionError");
                        operation = "G";
                        break;
                    case "U":
                        if (ts == 0)
                        {
                            if (qs <= 0)
                                throw new ArgumentException("AssertionError");
                            i++;
                            if (i >= steps)
                                throw new ArgumentException("IndexError");
                            ts = tsteps[i];
                            if (qsteps[i] != 0)
                                throw new ArgumentException("AssertionError");
                        }

                        if (qs == 0)
                        {
                            if (ts <= 0)
                                throw new ArgumentException("AssertionError");
                            i++;
                            if (i >= steps)
                                throw new ArgumentException("IndexError");
                            if (tsteps[i] != 0)
                                throw new ArgumentException("AssertionError");
                            qs = qsteps[i];
                        }

                        if (i >= operations.Length)
                            throw new ArgumentException("IndexError");
                        operation = operations[i].ToString();
                        if (operation != "U")
                            throw new ArgumentException("AssertionError");
                        operation = "N";
                        break;
                    case "S":
                    case "F":
                        break;
                    default:
                        throw new ArgumentException($"Unknown operation {operation}");
                }

                words.Add(operation);
                words.Add(Str(qs));
                words.Add(Str(ts));
                i++;
            }
        }
        else
        {
            int n = Math.Min(operations.Length, steps);
            for (int i = 0; i < n; i++)
            {
                long ts = tsteps[i], qs = qsteps[i];
                string operation = operations[i].ToString();
                long step;
                switch (operation)
                {
                    case "M":
                        if (ts == qs)
                        {
                            step = ts;
                        }
                        else if (ts == 3 * qs)
                        {
                            step = ts;
                            if (!qProt || tProt)
                                throw new ArgumentException("AssertionError");
                        }
                        else if (qs == 3 * ts)
                        {
                            step = qs;
                            if (qProt || !tProt)
                                throw new ArgumentException("AssertionError");
                        }
                        else
                        {
                            throw new ArgumentException("Unexpected steps target %d, query %s for operation 'M'");
                        }

                        break;
                    case "5":
                    case "3":
                        if (qs == 0)
                        {
                            step = ts;
                            operation = "D";
                        }
                        else if (ts == 0)
                        {
                            step = qs;
                            operation = "I";
                        }
                        else
                        {
                            if (qs != ts)
                                throw new ArgumentException("AssertionError");
                            step = ts;
                            operation = "M";
                        }

                        break;
                    case "N":
                        if (qs == 0)
                        {
                            step = ts;
                            operation = "D";
                        }
                        else if (ts == 0)
                        {
                            step = qs;
                            operation = "I";
                        }
                        else
                        {
                            throw new ArgumentException($"Unexpected intron with steps target {ts}, query {qs}");
                        }

                        break;
                    case "C":
                        if (ts != qs)
                            throw new ArgumentException("AssertionError");
                        step = ts;
                        operation = "M";
                        break;
                    case "D":
                        if (qs != 0)
                            throw new ArgumentException("AssertionError");
                        step = ts;
                        break;
                    case "I":
                        if (ts != 0)
                            throw new ArgumentException("AssertionError");
                        step = qs;
                        break;
                    case "U":
                    case "S":
                        if (ts > 0)
                        {
                            words.Add("D");
                            words.Add(Str(ts));
                        }

                        if (qs > 0)
                        {
                            words.Add("I");
                            words.Add(Str(qs));
                        }

                        continue;
                    case "F":
                        if (ts == 0)
                        {
                            step = qs;
                            operation = "I";
                        }
                        else if (qs == 0)
                        {
                            step = ts;
                            operation = "D";
                        }
                        else
                        {
                            throw new ArgumentException("Expected target step or query step to be 0");
                        }

                        break;
                    default:
                        throw new ArgumentException($"Unknown operation {operation}");
                }

                words.Add(operation);
                words.Add(Str(step));
            }
        }

        return string.Join(' ', words) + "\n";
    }

    private static void CheckCodonSteps(long ts, long qs, bool tProt, bool qProt)
    {
        if (ts == 3 * qs)
        {
            if (!qProt || tProt)
                throw new ArgumentException("AssertionError");
        }
        else if (qs == 3 * ts)
        {
            if (qProt || !tProt)
                throw new ArgumentException("AssertionError");
        }
        else
        {
            throw new ArgumentException($"Unexpected steps target {ts}, query {qs} for operation 'M'");
        }
    }

    private static (string Operation, long Step) CigarStep(long ts, long qs, bool tProt, bool qProt, string? tmt, string? qmt)
    {
        if (ts == qs)
            return ("M", ts);
        if (qs == 0)
            return ("D", ts);
        if (ts == 0)
            return ("I", qs);
        if (!tProt && qProt)
            return ("M", ts);
        if (tProt && !qProt)
            return ("M", qs);
        throw new ArgumentException($"Unexpected step target {ts}, query {qs} for molecule type {tmt ?? "None"}, {qmt ?? "None"}");
    }

    /// <summary><c>Bio/Align/exonerate.py</c> <c>AlignmentIterator</c>.</summary>
    private static CoordinateAlignmentFile ParseExonerate(string text)
    {
        var lines = Lines(text.Replace("\r\n", "\n").Replace('\r', '\n')).ToList();
        int li = 0;
        string ReadLine() => li < lines.Count ? lines[li++] : "";
        var metadata = new JsonObject { ["Program"] = "exonerate" };
        string line = ReadLine();
        const string prefix1 = "Command line: ";
        if (!line.StartsWith(prefix1, StringComparison.Ordinal))
            throw new ArgumentException("AssertionError");
        string commandline = line[prefix1.Length..].Trim();
        if (!commandline.StartsWith('[') || !commandline.EndsWith(']'))
            throw new ArgumentException("AssertionError");
        metadata["Command line"] = commandline.Length >= 2 ? commandline[1..^1] : "";
        line = ReadLine();
        const string prefix2 = "Hostname: ";
        if (!line.StartsWith(prefix2, StringComparison.Ordinal))
            throw new ArgumentException("AssertionError");
        string hostname = line[prefix2.Length..].Trim();
        if (!hostname.StartsWith('[') || !hostname.EndsWith(']'))
            throw new ArgumentException("AssertionError");
        metadata["Hostname"] = hostname.Length >= 2 ? hostname[1..^1] : "";
        var alignments = new List<BioCoordAlignment>();
        while (true)
        {
            if (li >= lines.Count)
                throw new ArgumentException("Failed to find 'completed exonerate analysis'; truncated file?");
            string l = lines[li++].Trim();
            if (l == "-- completed exonerate analysis")
            {
                if (ReadLine().Length > 0)
                    throw new ArgumentException("Found additional data after 'completed exonerate analysis'; corrupt file?");
                break;
            }

            if (l.StartsWith("vulgar: ", StringComparison.Ordinal))
                alignments.Add(ParseVulgar(Words(l[8..])));
            else if (l.StartsWith("cigar: ", StringComparison.Ordinal))
                alignments.Add(ParseCigar(Words(l[7..])));
            else
                throw new ArgumentException("UnboundLocalError: alignment");
        }

        return new CoordinateAlignmentFile { Alignments = alignments, Metadata = metadata };
    }

    private static (long TargetLength, long QueryLength, string? TargetType, string? QueryType) ExonerateStrands(long[] tco, long[] qco, string targetStrand, long targetStart, long targetEnd, string queryStrand, long queryStart, long queryEnd, bool vulgar)
    {
        long targetLength, queryLength;
        string? tmt = null, qmt = null;
        if (targetStrand == "+")
        {
            for (int i = 0; i < tco.Length; i++)
                tco[i] += targetStart;
            targetLength = targetEnd;
        }
        else if (targetStrand == "-")
        {
            for (int i = 0; i < tco.Length; i++)
                tco[i] = targetStart - tco[i];
            targetLength = targetStart;
        }
        else if (targetStrand == ".")
        {
            if (!vulgar && queryStrand != ".")
            {
                for (int i = 0; i < tco.Length; i++)
                    tco[i] = FloorDiv(tco[i] + 2, 3);
            }

            for (int i = 0; i < tco.Length; i++)
                tco[i] += targetStart;
            tmt = "protein";
            targetLength = targetEnd;
        }
        else
        {
            throw new ArgumentException("UnboundLocalError: target_molecule_type");
        }

        if (queryStrand == "+")
        {
            for (int i = 0; i < qco.Length; i++)
                qco[i] += queryStart;
            queryLength = queryEnd;
        }
        else if (queryStrand == "-")
        {
            for (int i = 0; i < qco.Length; i++)
                qco[i] = queryStart - qco[i];
            queryLength = queryStart;
        }
        else if (queryStrand == ".")
        {
            if (!vulgar && targetStrand != ".")
            {
                for (int i = 0; i < qco.Length; i++)
                    qco[i] = -FloorDiv(qco[i], -3);
            }

            for (int i = 0; i < qco.Length; i++)
                qco[i] += queryStart;
            qmt = "protein";
            queryLength = queryEnd;
        }
        else
        {
            throw new ArgumentException("UnboundLocalError: query_molecule_type");
        }

        return (targetLength, queryLength, tmt, qmt);
    }

    private static long FloorDiv(long a, long b)
    {
        long q = a / b;
        if ((a % b != 0) && ((a < 0) != (b < 0)))
            q--;
        return q;
    }

    private static BioCoordAlignment ExonerateAlignment(long[] tco, long[] qco, string targetId, string queryId, long targetLength, long queryLength, string? tmt, string? qmt, double score, string? operations)
    {
        if (targetLength < 0 || queryLength < 0)
            throw new ArgumentException("length must not be negative.");
        var tann = new JsonObject();
        var qann = new JsonObject();
        if (tmt is not null)
            tann["molecule_type"] = tmt;
        if (qmt is not null)
            qann["molecule_type"] = qmt;
        var attrs = new JsonObject();
        if (operations is not null)
            attrs["operations"] = operations;
        attrs["score"] = score;
        return new BioCoordAlignment
        {
            Sequences =
            [
                new BioCoordSequence { Id = targetId, Length = targetLength, Parts = UndefinedParts(targetLength), Annotations = tann },
                new BioCoordSequence { Id = queryId, Length = queryLength, Parts = UndefinedParts(queryLength), Annotations = qann },
            ],
            Coordinates = [tco, qco],
            Attributes = attrs,
        };
    }

    private static void NeedWords(string[] words, int n)
    {
        if (words.Length < n)
            throw new ArgumentException("IndexError: list index out of range");
    }

    private static BioCoordAlignment ParseCigar(string[] words)
    {
        NeedWords(words, 9);
        string queryId = words[0];
        long queryStart = PyInt(words[1]), queryEnd = PyInt(words[2]);
        string queryStrand = words[3], targetId = words[4];
        long targetStart = PyInt(words[5]), targetEnd = PyInt(words[6]);
        string targetStrand = words[7];
        double score = PyFloatParse(words[8]);
        if ((words.Length - 9) % 2 != 0)
            throw new ArgumentException("uninitialised coordinates (np.empty): odd operation list");
        int n = (words.Length - 8) / 2;
        var tco = new long[n + 1];
        var qco = new long[n + 1];
        long qs = 0, ts = 0;
        for (int i = 0; 10 + 2 * i < words.Length; i++)
        {
            string operation = words[9 + 2 * i];
            long step = PyInt(words[10 + 2 * i]);
            switch (operation)
            {
                case "M": ts += step; qs += step; break;
                case "I": qs += queryStrand == "." && targetStrand != "." ? step * 3 : step; break;
                case "D": ts += targetStrand == "." && queryStrand != "." ? step * 3 : step; break;
                default: throw new ArgumentException($"Unknown operation {operation} in cigar string");
            }

            tco[i + 1] = ts;
            qco[i + 1] = qs;
        }

        var (tl, ql, tmt, qmt) = ExonerateStrands(tco, qco, targetStrand, targetStart, targetEnd, queryStrand, queryStart, queryEnd, false);
        return ExonerateAlignment(tco, qco, targetId, queryId, tl, ql, tmt, qmt, score, null);
    }

    private static BioCoordAlignment ParseVulgar(string[] words)
    {
        NeedWords(words, 9);
        string queryId = words[0];
        long queryStart = PyInt(words[1]), queryEnd = PyInt(words[2]);
        string queryStrand = words[3], targetId = words[4];
        long targetStart = PyInt(words[5]), targetEnd = PyInt(words[6]);
        string targetStrand = words[7];
        double score = PyFloatParse(words[8]);
        var ops = new List<string>();
        for (int i = 9; i < words.Length; i += 3)
            ops.Add(words[i]);
        int n = (words.Length - 8) / 3 + ops.Count(o => o == "N");
        var tco = new long[n + 1];
        var qco = new long[n + 1];
        var operations = new char[n];
        long qs = 0, ts = 0;
        int k = 0;
        for (int j = 0; 11 + 3 * j < words.Length; j++)
        {
            string operation = words[9 + 3 * j];
            long queryStep = PyInt(words[10 + 3 * j]), targetStep = PyInt(words[11 + 3 * j]);
            switch (operation)
            {
                case "M": break;
                case "5":
                case "3":
                    if (!(targetStep == 2 || queryStep == 2))
                        throw new ArgumentException("AssertionError");
                    break;
                case "I": operation = "N"; break;
                case "C":
                    if (targetStep % 3 != 0 || queryStep % 3 != 0)
                        throw new ArgumentException("AssertionError");
                    break;
                case "G":
                    if (queryStep == 0)
                        operation = "D";
                    else if (targetStep == 0)
                        operation = "I";
                    else
                        throw new ArgumentException($"Unexpected gap operation with steps {queryStep}, {targetStep} in vulgar line");
                    break;
                case "N":
                    operation = "U";
                    if (targetStep > 0)
                    {
                        ts += targetStep;
                        Put(tco, qco, operations, k, ts, qs, operation);
                        k++;
                    }

                    if (queryStep > 0)
                    {
                        qs += queryStep;
                        Put(tco, qco, operations, k, ts, qs, operation);
                        k++;
                    }

                    continue;
                case "S":
                case "F": break;
                default: throw new ArgumentException($"Unknown operation {operation} in vulgar string");
            }

            ts += targetStep;
            qs += queryStep;
            Put(tco, qco, operations, k, ts, qs, operation);
            k++;
        }

        var (tl, ql, tmt, qmt) = ExonerateStrands(tco, qco, targetStrand, targetStart, targetEnd, queryStrand, queryStart, queryEnd, true);
        return ExonerateAlignment(tco, qco, targetId, queryId, tl, ql, tmt, qmt, score, new string(operations));
    }

    private static void Put(long[] tco, long[] qco, char[] operations, int k, long ts, long qs, string operation)
    {
        if (k + 1 >= tco.Length)
            throw new ArgumentException("IndexError");
        tco[k + 1] = ts;
        qco[k + 1] = qs;
        operations[k] = operation[0];
    }
}
