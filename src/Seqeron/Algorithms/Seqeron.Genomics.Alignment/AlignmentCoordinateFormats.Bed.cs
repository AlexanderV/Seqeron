using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary><c>Bio/Align/bed.py</c> <c>AlignmentWriter.format_alignment</c>.</summary>
    private static string FormatBed(BioCoordAlignment alignment, int bedN)
    {
        if (bedN < 3 || bedN > 12)
            throw new ArgumentException("bedN must be between 3 and 12");
        var coordinates = Copy(alignment.Coordinates);
        if (coordinates.Length == 0 || coordinates[0].Length == 0)
            return "";
        var (target, query) = Pair(alignment);
        string chrom = Id(target, "target");
        if (coordinates[0][0] > coordinates[0][^1])
            ReverseColumns(coordinates);
        string strand = coordinates[1][0] > coordinates[1][^1] ? "-" : "+";
        var blockSizes = new List<long>();
        var blockStarts = new List<long>();
        long tStart = coordinates[0][0], qStart = coordinates[1][0], blockSize = 0;
        for (int k = 1; k < coordinates[0].Length; k++)
        {
            long tEnd = coordinates[0][k], qEnd = coordinates[1][k];
            if (tStart == tEnd)
            {
                qStart = qEnd;
            }
            else if (qStart == qEnd)
            {
                tStart = tEnd;
            }
            else
            {
                blockSize = tEnd - tStart;
                blockStarts.Add(tStart);
                blockSizes.Add(blockSize);
                tStart = tEnd;
                qStart = qEnd;
            }
        }

        long chromStart, chromEnd;
        if (blockStarts.Count > 0)
        {
            chromStart = blockStarts[0];
            chromEnd = blockStarts[^1] + blockSize;
        }
        else
        {
            chromStart = chromEnd = tStart;
        }

        var fields = new List<string> { chrom, Str(chromStart), Str(chromEnd) };
        string Done() => string.Join('\t', fields) + "\n";
        if (bedN == 3)
            return Done();
        fields.Add(Id(query, "query"));
        if (bedN == 4)
            return Done();
        fields.Add(TryAttr(alignment, "score", out var score) ? G(score) : "0");
        if (bedN == 5)
            return Done();
        fields.Add(strand);
        if (bedN == 6)
            return Done();
        fields.Add(TryAttr(alignment, "thickStart", out var ts) ? PyStr(ts) : Str(chromStart));
        if (bedN == 7)
            return Done();
        fields.Add(TryAttr(alignment, "thickEnd", out var te) ? PyStr(te) : Str(chromEnd));
        if (bedN == 8)
            return Done();
        fields.Add(TryAttr(alignment, "itemRgb", out var rgb) ? PyStr(rgb) : "0");
        if (bedN == 9)
            return Done();
        fields.Add(Str(blockSizes.Count));
        if (bedN == 10)
            return Done();
        fields.Add(string.Join(',', blockSizes.Select(Str)) + ",");
        if (bedN == 11)
            return Done();
        fields.Add(string.Join(',', blockStarts.Select(s => Str(s - chromStart))) + ",");
        return Done();
    }

    /// <summary><c>Bio/Align/bed.py</c> <c>AlignmentIterator</c>.</summary>
    private static List<BioCoordAlignment> ParseBed(string text)
    {
        var result = new List<BioCoordAlignment>();
        foreach (string line in Lines(text))
        {
            var words = Words(line);
            int bedN = words.Length;
            if (bedN < 3 || bedN > 12)
                throw new ArgumentException($"expected between 3 and 12 columns, found {bedN}");
            string chrom = words[0];
            long chromStart = PyInt(words[1]), chromEnd = PyInt(words[2]);
            string? name = bedN > 3 ? words[3] : null;
            string strand = bedN > 5 ? words[5] : "+";
            long[][] coordinates;
            long qSize;
            if (bedN > 9)
            {
                if (bedN < 12)
                    throw new ArgumentException("list index out of range"); // Biopython reads words[10] / words[11] of a line with 10 or 11 columns (IndexError)
                long blockCount = PyInt(words[9]);
                var blockSizes = IntList(words[10]);
                var blockStarts = IntList(words[11]);
                if (blockSizes.Count != blockCount)
                    throw new ArgumentException($"Inconsistent number of block sizes ({blockSizes.Count} found, expected {blockCount})");
                if (blockStarts.Count != blockCount)
                    throw new ArgumentException($"Inconsistent number of block start positions ({blockStarts.Count} found, expected {blockCount})");
                long tPosition = 0, qPosition = 0;
                var ts = new List<long> { tPosition };
                var qs = new List<long> { qPosition };
                for (int i = 0; i < blockSizes.Count; i++)
                {
                    if (blockStarts[i] != tPosition)
                    {
                        ts.Add(blockStarts[i]);
                        qs.Add(qPosition);
                        tPosition = blockStarts[i];
                    }

                    tPosition += blockSizes[i];
                    qPosition += blockSizes[i];
                    ts.Add(tPosition);
                    qs.Add(qPosition);
                }

                coordinates = [[.. ts], [.. qs]];
                qSize = blockSizes.Sum();
            }
            else
            {
                long blockSize = chromEnd - chromStart;
                coordinates = [[0, blockSize], [0, blockSize]];
                qSize = blockSize;
            }

            for (int c = 0; c < coordinates[0].Length; c++)
                coordinates[0][c] += chromStart;
            if (strand == "-")
            {
                for (int c = 0; c < coordinates[1].Length; c++)
                    coordinates[1][c] = qSize - coordinates[1][c];
            }

            if (chromStart != coordinates[0][0])
                throw new ArgumentException($"Inconsistent chromStart found ({chromStart}, expected {coordinates[0][0]})");
            if (chromEnd != coordinates[0][^1])
                throw new ArgumentException($"Inconsistent chromEnd found ({chromEnd}, expected {coordinates[0][^1]})");
            var attrs = new JsonObject();
            if (bedN > 4)
            {
                attrs["score"] = double.TryParse(words[4], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double sc) && !words[4].Contains('_') ? sc : words[4];
                if (bedN > 6)
                    attrs["thickStart"] = PyInt(words[6]);
                if (bedN > 7)
                    attrs["thickEnd"] = PyInt(words[7]);
                if (bedN > 8)
                    attrs["itemRgb"] = words[8];
            }

            result.Add(new BioCoordAlignment { Sequences = [Unknown(chrom, long.MaxValue), Unknown(name, qSize)], Coordinates = coordinates, Attributes = attrs });
        }

        return result;
    }
}
