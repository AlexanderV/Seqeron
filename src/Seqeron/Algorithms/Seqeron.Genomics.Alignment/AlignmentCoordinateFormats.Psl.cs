using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    private const string PslHeader = "psLayout version 3\n\nmatch\tmis- \trep. \tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        \tT   \tT    \tT  \tblock\tblockSizes \tqStarts\t tStarts\n     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     \tsize\tstart\tend\tname     \tsize\tstart\tend\tcount\n---------------------------------------------------------------------------------------------------------------------------------------------------------------\n";

    /// <summary><c>Bio/Align/psl.py</c> <c>AlignmentWriter.format_alignment</c>.</summary>
    private static string FormatPsl(BioCoordAlignment alignment, CoordinateWriteOptions options)
    {
        var coordinates = Copy(alignment.Coordinates);
        if (coordinates.Length == 0 || coordinates[0].Length == 0)
            return "";
        var (target, query) = Pair(alignment);
        string qName = Id(query, "None"), tName = Id(target, "None");
        long tSize = target.Length, qSize = query.Length;
        int n = coordinates[0].Length;
        long tCount = 0, qCount = 0;
        for (int k = 1; k < n; k++)
        {
            long dt = Math.Abs(coordinates[0][k] - coordinates[0][k - 1]), dq = Math.Abs(coordinates[1][k] - coordinates[1][k - 1]);
            if (Math.Min(dt, dq) > 0)
            {
                tCount += dt;
                qCount += dq;
            }
        }

        bool dnax;
        string strand;
        if (tCount == qCount)
        {
            dnax = false;
            if (coordinates[0][0] > coordinates[0][^1])
                ReverseColumns(coordinates);
            if (coordinates[1][0] > coordinates[1][^1])
            {
                strand = "-";
                query = ReverseComplement(query);
                for (int i = 0; i < coordinates[1].Length; i++)
                    coordinates[1][i] = qSize - coordinates[1][i];
            }
            else
            {
                strand = "+";
            }
        }
        else if (tCount == 3 * qCount)
        {
            dnax = true;
            if (coordinates[0][0] > coordinates[0][^1])
            {
                strand = "-";
                target = ReverseComplement(target);
                for (int i = 0; i < coordinates[0].Length; i++)
                    coordinates[0][i] = tSize - coordinates[0][i];
            }
            else
            {
                strand = "+";
            }
        }
        else
        {
            throw new ArgumentException($"inconsistent step sizes {tCount}, {qCount}");
        }

        int? wildcard = null;
        if (options.Wildcard is not null)
            wildcard = options.Mask == "upper" ? char.ToLowerInvariant(options.Wildcard[0]) : char.ToUpperInvariant(options.Wildcard[0]);
        string? mask = options.Mask;
        long matches = 0, misMatches = 0, repMatches = 0, nCount = 0, qNumInsert = 0, qBaseInsert = 0, tNumInsert = 0, tBaseInsert = 0;
        var blockSizes = new List<long>();
        var qStarts = new List<long>();
        var tStarts = new List<long>();
        long tStart = coordinates[0][0], qStart = coordinates[1][0];
        long lastTCount = 0, lastQCount = 0;
        for (int k = 1; k < n; k++)
        {
            long tEnd = coordinates[0][k], qEnd = coordinates[1][k];
            if (tStart == tEnd)
            {
                if (qStart > 0 && qEnd < qSize)
                {
                    qNumInsert++;
                    qBaseInsert += qEnd - qStart;
                }

                qStart = qEnd;
            }
            else if (qStart == qEnd)
            {
                if (tStart > 0 && tEnd < tSize)
                {
                    tNumInsert++;
                    tBaseInsert += tEnd - tStart;
                }

                tStart = tEnd;
            }
            else
            {
                lastTCount = tEnd - tStart;
                lastQCount = qEnd - qStart;
                tStarts.Add(tStart);
                qStarts.Add(qStart);
                blockSizes.Add(lastQCount);
                string? tSeq = target.Slice(tStart, tEnd);
                string? qSeq = query.Slice(qStart, qEnd);
                if (tSeq is null || qSeq is null)
                {
                    matches += lastQCount;
                }
                else if (mask is "lower" or "upper")
                {
                    for (int i = 0; i < tSeq.Length && i < qSeq.Length; i++)
                    {
                        char c1 = tSeq[i];
                        char u1 = mask == "lower" ? char.ToUpperInvariant(c1) : char.ToLowerInvariant(c1);
                        char u2 = mask == "lower" ? char.ToUpperInvariant(qSeq[i]) : char.ToLowerInvariant(qSeq[i]);
                        if (wildcard is int w && (u1 == w || u2 == w))
                            nCount++;
                        else if (u1 == u2)
                        {
                            if (u1 == c1)
                                matches++;
                            else
                                repMatches++;
                        }
                        else
                            misMatches++;
                    }
                }
                else
                {
                    for (int i = 0; i < tSeq.Length && i < qSeq.Length; i++)
                    {
                        char u1 = char.ToUpperInvariant(tSeq[i]), u2 = char.ToUpperInvariant(qSeq[i]);
                        if (wildcard is int w && (u1 == w || u2 == w))
                            nCount++;
                        else if (u1 == u2)
                            matches++;
                        else
                            misMatches++;
                    }
                }

                tStart = tEnd;
                qStart = qEnd;
            }
        }

        if (TryAttr(alignment, "matches", out var am)) matches = Int(am);
        if (TryAttr(alignment, "misMatches", out var amm)) misMatches = Int(amm);
        if (TryAttr(alignment, "repMatches", out var arm)) repMatches = Int(arm);
        if (TryAttr(alignment, "nCount", out var anc)) nCount = Int(anc);
        if (tStarts.Count == 0)
            throw new ArgumentException("The alignment has no aligned block.");
        tStart = tStarts[0];
        qStart = qStarts[0];
        long tEndAll = tStarts[^1] + lastTCount;
        long qEndAll = qStarts[^1] + lastQCount;
        if (strand == "-")
        {
            if (dnax)
                (tStart, tEndAll) = (tSize - tEndAll, tSize - tStart);
            else
                (qStart, qEndAll) = (qSize - qEndAll, qSize - qStart);
        }

        string sizes = string.Join(',', blockSizes.Select(Str)) + ",";
        string qs = string.Join(',', qStarts.Select(Str)) + ",";
        string ts = string.Join(',', tStarts.Select(Str)) + ",";
        if (dnax)
            strand = "+" + strand;
        var words = new[]
        {
            Str(matches), Str(misMatches), Str(repMatches), Str(nCount), Str(qNumInsert), Str(qBaseInsert), Str(tNumInsert), Str(tBaseInsert), strand, qName, Str(qSize), Str(qStart), Str(qEndAll),
            tName, Str(tSize), Str(tStart), Str(tEndAll), Str(blockSizes.Count), sizes, qs, ts,
        };
        return string.Join('\t', words) + "\n";
    }

    /// <summary><c>Bio/Align/psl.py</c> <c>AlignmentIterator</c> (the header, PSL and PSLX lines; the translated-protein PSLX with its CDS feature is not read).</summary>
    private static List<BioCoordAlignment> ParsePsl(string text)
    {
        var result = new List<BioCoordAlignment>();
        var all = Lines(text).ToList();
        int p = 0;
        if (all.Count == 0)
            throw new ArgumentException("line has 0 columns; expected 21 or 23"); // the first readline() of an empty stream is '' and is parsed as a line
        if (all[0].StartsWith("psLayout ", StringComparison.Ordinal))
        {
            var hw = Words(all[0]);
            if (hw.Length < 2 || hw[1] != "version")
                throw new ArgumentException($"Unexpected word '{(hw.Length > 1 ? hw[1] : "")}' in header line");
            p = 4;
            string last = p < all.Count ? all[p] : "";
            p++;
            if (last.TrimStart('-').Trim() != "")
                throw new ArgumentException("End of header not found");
        }

        for (; p < all.Count; p++)
        {
            var words = Words(all[p]);
            bool pslx;
            if (words.Length == 23)
                pslx = true;
            else if (words.Length == 21)
                pslx = false;
            else
                throw new ArgumentException($"line has {words.Length} columns; expected 21 or 23");
            string strand = words[8], qName = words[9];
            long qSize = PyInt(words[10]);
            string tName = words[13];
            long tSize = PyInt(words[14]);
            long blockCount = PyInt(words[17]);
            var blockSizes = IntList(words[18]);
            var qStartsL = IntList(words[19]);
            var tStartsL = IntList(words[20]);
            if (blockSizes.Count != blockCount)
                throw new ArgumentException($"Inconsistent number of blocks ({blockSizes.Count} found, expected {blockCount})");
            if (qStartsL.Count != blockCount)
                throw new ArgumentException($"Inconsistent number of query start positions ({qStartsL.Count} found, expected {blockCount})");
            if (tStartsL.Count != blockCount)
                throw new ArgumentException($"Inconsistent number of target start positions ({tStartsL.Count} found, expected {blockCount})");
            bool protein = strand is "++" or "+-";
            long qPosition = qStartsL[0], tPosition = tStartsL[0];
            var ts = new List<long> { tPosition };
            var qs = new List<long> { qPosition };
            for (int i = 0; i < blockSizes.Count; i++)
            {
                long tBlock = protein ? 3 * blockSizes[i] : blockSizes[i];
                if (tStartsL[i] != tPosition)
                {
                    ts.Add(tStartsL[i]);
                    qs.Add(qPosition);
                    tPosition = tStartsL[i];
                }

                if (qStartsL[i] != qPosition)
                {
                    ts.Add(tPosition);
                    qs.Add(qStartsL[i]);
                    qPosition = qStartsL[i];
                }

                tPosition += tBlock;
                qPosition += blockSizes[i];
                ts.Add(tPosition);
                qs.Add(qPosition);
            }

            var co = new[] { ts.ToArray(), qs.ToArray() };
            long qNumInsert = 0, qBaseInsert = 0, tNumInsert = 0, tBaseInsert = 0;
            long tS = co[0][0], qS = co[1][0];
            for (int k = 1; k < co[0].Length; k++)
            {
                long tE = co[0][k], qE = co[1][k];
                long tCount = tE - tS, qCount = qE - qS;
                if (tCount == 0)
                {
                    if (qS > 0 && qE < qSize)
                    {
                        qNumInsert++;
                        qBaseInsert += qCount;
                    }

                    qS = qE;
                }
                else if (qCount == 0)
                {
                    if (tS > 0 && tE < tSize)
                    {
                        tNumInsert++;
                        tBaseInsert += tCount;
                    }

                    tS = tE;
                }
                else
                {
                    tS = tE;
                    qS = qE;
                }
            }

            if (qNumInsert != PyInt(words[4]))
                throw new ArgumentException($"Inconsistent qNumInsert found ({words[4]}, expected {qNumInsert})");
            if (qBaseInsert != PyInt(words[5]))
                throw new ArgumentException($"Inconsistent qBaseInsert found ({words[5]}, expected {qBaseInsert})");
            if (tNumInsert != PyInt(words[6]))
                throw new ArgumentException($"Inconsistent tNumInsert found ({words[6]}, expected {tNumInsert})");
            if (tBaseInsert != PyInt(words[7]))
                throw new ArgumentException($"Inconsistent tBaseInsert found ({words[7]}, expected {tBaseInsert})");
            long qStart = PyInt(words[11]), qEnd = PyInt(words[12]), tStart = PyInt(words[15]), tEnd = PyInt(words[16]);
            if (strand == "-")
            {
                (qStart, qEnd) = (qEnd, qStart);
                for (int i = 0; i < co[1].Length; i++) co[1][i] = qSize - co[1][i];
            }
            else if (strand == "+-")
            {
                (tStart, tEnd) = (tEnd, tStart);
                for (int i = 0; i < co[0].Length; i++) co[0][i] = tSize - co[0][i];
            }

            if (tStart != co[0][0])
                throw new ArgumentException($"Inconsistent tStart found ({tStart}, expected {co[0][0]})");
            if (tEnd != co[0][^1])
                throw new ArgumentException($"Inconsistent tEnd found ({tEnd}, expected {co[0][^1]})");
            if (qStart != co[1][0])
                throw new ArgumentException($"Inconsistent qStart found ({qStart}, expected {co[1][0]})");
            if (qEnd != co[1][^1])
                throw new ArgumentException($"Inconsistent qEnd found ({qEnd}, expected {co[1][^1]})");
            BioCoordSequence target, query;
            if (pslx)
            {
                if (protein)
                    throw new NotSupportedException("A PSLX line of a protein aligned to translated DNA (strand ++ / +-) carries a CDS feature with the translation; it is not read by this library.");
                var qSeqs = words[21].TrimEnd(',').Split(',');
                var tSeqs = words[22].TrimEnd(',').Split(',');
                var qParts = new SortedDictionary<long, string>();
                var tParts = new SortedDictionary<long, string>();
                for (int i = 0; i < Math.Min(qStartsL.Count, qSeqs.Length); i++) qParts[qStartsL[i]] = qSeqs[i];
                for (int i = 0; i < Math.Min(tStartsL.Count, tSeqs.Length); i++) tParts[tStartsL[i]] = tSeqs[i];
                target = new BioCoordSequence { Id = tName, Length = tSize, Parts = tParts };
                query = new BioCoordSequence { Id = qName, Length = qSize, Parts = qParts };
                if (strand == "-")
                    query = ReverseComplement(query);
            }
            else
            {
                target = Unknown(tName, tSize);
                query = Unknown(qName, qSize);
            }

            var attrs = new JsonObject { ["matches"] = PyInt(words[0]), ["misMatches"] = PyInt(words[1]), ["repMatches"] = PyInt(words[2]), ["nCount"] = PyInt(words[3]) };
            result.Add(new BioCoordAlignment { Sequences = [target, query], Coordinates = co, Attributes = attrs });
        }

        return result;
    }
}
