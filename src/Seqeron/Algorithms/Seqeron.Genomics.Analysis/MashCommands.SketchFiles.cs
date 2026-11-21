using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The sketch file commands of Mash 2.3: <c>mash sketch -o</c>, <c>mash info</c> (<c>CommandInfo.cpp</c>) and <c>mash paste</c> (<c>CommandPaste.cpp</c>) on the <c>.msh</c> file of <see cref="MashSketchFile"/> (A5-10b, B09 F193).</summary>
public static partial class MashCommands
{
    /// <summary>
    /// <c>mash sketch</c>: the content of the sketch file of the inputs (one sketch per file, or per sequence with <paramref name="individual"/>; a read set with <paramref name="reads"/> also keeps the multiplicities of the
    /// hashes). <see cref="SketchFile"/> writes it byte for byte as the program writes it.
    /// </summary>
    public static MashSketchFileData Sketch(IReadOnlyList<MashInput> inputs, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, CancellationToken cancellationToken = default, MashReadsOptions? reads = null, string? idOverride = null, string? commentOverride = null, MashAlphabet? alphabet = null)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        alphabet ??= MashAlphabet.Nucleotide();
        canonical &= !alphabet.ImpliesNoncanonical;                 // -a and -z imply -n
        RequireSketchParameters(k, sketchSize);
        // -r: all the files are ONE read set (initFromReads): one sketch named after the first file; the files are read in turn, one record of each (paired reads), a file leaves the rotation at its end
        if (reads is not null && inputs.Count > 1)
            inputs = [new MashInput(inputs[0].Name, Interleave(inputs, k))];
        var list = Sketches(inputs, nameof(inputs), k, sketchSize, seed, canonical, individual, cancellationToken, reads, false, alphabet);
        // -I / -C: the id and the comment of the FIRST sketch (the program warns when there are several files and they are not read sets)
        if (list.Count > 0 && (idOverride is not null || commentOverride is not null))
            list[0] = list[0] with { Name = idOverride ?? list[0].Name, Comment = commentOverride ?? list[0].Comment };
        var references = list.Select(n => new MashSketchFileReference(n.Name, n.Comment, (ulong)n.Sketch.Length, n.Sketch.Hashes, reads is not null ? n.Sketch.Counts : null)).ToList();
        return new MashSketchFileData(k, sketchSize, seed, alphabet.Characters, !canonical, alphabet.PreserveCase, !individual, alphabet.Uses64BitHash(k), references);
    }

    /// <summary>
    /// A FASTA text as Mash reads it (kseq): the first word of a header is the id, the rest the comment. The comment that names the sketch of a FILE (<see cref="MashRecord.StaleComment"/>) of a record WITHOUT a comment is that
    /// of the nearest earlier record of the file that had one (kseq clears the length of the comment but not its text, and Mash reads the C string: measured on mash 2.3; a record before any comment has none); the comment of
    /// a per-sequence sketch (<c>-i</c>) is the record's own. Whitespace inside sequence lines is dropped.
    /// </summary>
    public static MashInput ParseFasta(string name, string fasta)
    {
        ArgumentNullException.ThrowIfNull(name);
        ArgumentNullException.ThrowIfNull(fasta);
        var records = new List<MashRecord>();
        string id = "", comment = "", stale = "";
        var seq = new StringBuilder();
        bool open = false;
        void Close()
        {
            if (open)
                records.Add(new MashRecord(id, comment, seq.ToString()) { StaleComment = stale });
            seq.Clear();
        }

        foreach (string raw in fasta.Split('\n'))
        {
            string line = raw.TrimEnd('\r');
            if (line.StartsWith('>'))
            {
                Close();
                open = true;
                string header = line[1..].Trim();
                int sp = header.IndexOfAny([' ', '\t']);
                id = sp < 0 ? header : header[..sp];
                comment = sp >= 0 ? header[(sp + 1)..].Trim() : "";
                if (sp >= 0)
                    stale = comment;
            }
            else if (open)
            {
                foreach (char c in line)
                    if (!char.IsWhiteSpace(c))
                        seq.Append(c);
            }
        }

        Close();
        return new MashInput(name, records);
    }

    // Sketch.cpp sketchFile: a rotation over the open files; a read shorter than k is skipped WITHOUT passing on to the next file, a processed read does
    private static List<MashRecord> Interleave(IReadOnlyList<MashInput> inputs, int k)
    {
        var queues = inputs.Select(i => new Queue<MashRecord>(i.Records)).ToList();
        var order = new List<MashRecord>();
        int it = 0;
        while (queues.Count > 0)
        {
            if (queues[it].Count == 0)
            {
                queues.RemoveAt(it);
                if (it >= queues.Count)
                    it = 0;
                continue;
            }

            var record = queues[it].Dequeue();
            if (record.Sequence.Length < k)
                continue;
            order.Add(record);
            it = (it + 1) % queues.Count;
        }

        return order;
    }

    /// <summary>The bytes of the <c>.msh</c> file <c>mash sketch -o x</c> writes (<see cref="MashSketchFile.Write"/> of <see cref="Sketch"/>).</summary>
    public static byte[] SketchFile(IReadOnlyList<MashInput> inputs, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, CancellationToken cancellationToken = default, MashReadsOptions? reads = null, string? idOverride = null, string? commentOverride = null, MashAlphabet? alphabet = null)
        => MashSketchFile.Write(Sketch(inputs, k, sketchSize, seed, canonical, individual, cancellationToken, reads, idOverride, commentOverride, alphabet));

    /// <summary>
    /// <c>mash paste out.msh a.msh b.msh ...</c> (<c>CommandPaste.cpp</c>): one sketch file of the sketches of several. The parameters (k, sketch size, seed, alphabet, strand mode, case) are those of the first file; a later file with
    /// another alphabet, seed, k or strand mode, or a smaller sketch size is skipped (the program warns on the error stream), a larger sketch is reduced to the first file's size; the hash counts of read-set sketches are not
    /// kept (the pasted set is not in reads mode). <see cref="MashSketchFile.Write"/> writes the result as the program does.
    /// </summary>
    /// <exception cref="ArgumentException">No file.</exception>
    public static MashSketchFileData Paste(IReadOnlyList<MashSketchFileData> files)
    {
        ArgumentNullException.ThrowIfNull(files);
        if (files.Count == 0)
            throw new ArgumentException("At least one sketch file is required.", nameof(files));
        var first = files[0];
        var references = new List<MashSketchFileReference>();
        foreach (var f in files)
        {
            if (f.Alphabet != first.Alphabet || f.HashSeed != first.HashSeed || f.KmerSize != first.KmerSize || f.MinHashesPerWindow < first.MinHashesPerWindow || f.Noncanonical != first.Noncanonical)
                continue;
            foreach (var r in f.References)
                references.Add(new MashSketchFileReference(r.Name, r.Comment, r.Length, r.Hashes.Count > first.MinHashesPerWindow ? r.Hashes.Take(first.MinHashesPerWindow).ToArray() : r.Hashes));
        }

        return new MashSketchFileData(first.KmerSize, first.MinHashesPerWindow, first.HashSeed, first.Alphabet, first.Noncanonical, first.PreserveCase, first.Concatenated, first.Use64, references, first.WindowSize, first.Error);
    }

    /// <summary>
    /// The seed of <c>-S value</c> as the program stores it (measured on mash 2.3): a value of 4294967295 (0xFFFFFFFF) or more is stored as 0 (the option parser accepts 4294967296 too; 5000000000 is refused), any other
    /// value up to 4294967294 is the seed.
    /// </summary>
    public static uint SeedFromOption(long value) => value is >= 0 and < uint.MaxValue ? (uint)value : 0u;

    /// <summary>The <c>mash info</c> output options.</summary>
    /// <param name="HeaderOnly"><c>-H</c>: only the header (incompatible with <c>-d</c>, <c>-t</c>, <c>-c</c>).</param>
    /// <param name="Tabular"><c>-t</c>: tab separated, no header (incompatible with <c>-d</c>, <c>-H</c>, <c>-c</c>).</param>
    /// <param name="Counts"><c>-c</c>: the hash count histogram of each sketch (incompatible with <c>-d</c>, <c>-H</c>, <c>-t</c>).</param>
    /// <param name="Dump"><c>-d</c>: the sketches as JSON (incompatible with <c>-H</c>, <c>-t</c>, <c>-c</c>).</param>
    public sealed record MashInfoOptions(bool HeaderOnly = false, bool Tabular = false, bool Counts = false, bool Dump = false);

    /// <summary>
    /// <c>mash info [-H] [-t] [-c] [-d] file.msh</c> (<c>CommandInfo.cpp</c>): the header (hash function and seed, k with the hash width, the alphabet with "(canonical)" / "(case-sensitive)", the target sketch size and the number of
    /// sketches) and a column table of the sketches (<c>printColumns</c> with indent 2, spacing 2, "-" for an empty field; the width of the terminal is never reached in a pipe), or the tab separated table, the count histogram or
    /// the JSON dump. An incompatible pair of options or a missing count list is reported as the program does (an <see cref="ArgumentException"/> with the program's message).
    /// </summary>
    public static string Info(MashSketchFileData sketch, MashInfoOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sketch);
        var o = options ?? new MashInfoOptions();
        if (o.HeaderOnly && o.Tabular)
            throw new ArgumentException("The options -H and -t are incompatible.");
        if (o.HeaderOnly && o.Counts)
            throw new ArgumentException("The options -H and -c are incompatible.");
        if (o.Tabular && o.Counts)
            throw new ArgumentException("The options -t and -c are incompatible.");
        if (o.Dump && o.Tabular)
            throw new ArgumentException("The options -d and -t are incompatible.");
        if (o.Dump && o.HeaderOnly)
            throw new ArgumentException("The options -d and -H are incompatible.");
        if (o.Dump && o.Counts)
            throw new ArgumentException("The options -d and -c are incompatible.");
        var sb = new StringBuilder();
        if (o.Counts)
            return InfoCounts(sketch);
        if (o.Dump)
            return InfoJson(sketch);
        if (o.Tabular)
        {
            sb.Append("#Hashes\tLength\tID\tComment\n");
        }
        else
        {
            sb.Append("Header:\n");
            sb.Append("  Hash function (seed):          MurmurHash3_x64_128 (").Append(sketch.HashSeed.ToString(CultureInfo.InvariantCulture)).Append(")\n");
            sb.Append("  K-mer size:                    ").Append(sketch.KmerSize.ToString(CultureInfo.InvariantCulture)).Append(" (").Append(sketch.Use64 ? "64" : "32").Append("-bit hashes)\n");
            sb.Append("  Alphabet:                      ").Append(sketch.Alphabet).Append(sketch.Noncanonical ? "" : " (canonical)").Append(sketch.PreserveCase ? " (case-sensitive)" : "").Append('\n');
            sb.Append("  Target min-hashes per sketch:  ").Append(sketch.MinHashesPerWindow.ToString(CultureInfo.InvariantCulture)).Append('\n');
            sb.Append("  Sketches:                      ").Append(sketch.References.Count.ToString(CultureInfo.InvariantCulture)).Append('\n');
        }

        if (o.HeaderOnly)
            return sb.ToString();
        if (o.Tabular)
        {
            foreach (var r in sketch.References)
                sb.Append(r.Hashes.Count.ToString(CultureInfo.InvariantCulture)).Append('\t').Append(r.Length.ToString(CultureInfo.InvariantCulture)).Append('\t').Append(r.Name).Append('\t').Append(r.Comment).Append('\n');
            return sb.ToString();
        }

        sb.Append("\nSketches:\n");
        var columns = new List<string[]> { new[] { "[Hashes]", "[Length]", "[ID]", "[Comment]" } };
        foreach (var r in sketch.References)
            columns.Add([r.Hashes.Count.ToString(CultureInfo.InvariantCulture), r.Length.ToString(CultureInfo.InvariantCulture), r.Name, r.Comment]);
        sb.Append(PrintColumns(columns, 2, 2, sketch.References.Count == 0));
        return sb.ToString();
    }

    // printColumns of Command.cpp without the line wrapping (the terminal width of a non-terminal run is far beyond the lines): rows of padded columns, an empty field is "-", every row is followed by a blank line
    // narrow: the terminal width the program reads is uninitialised memory when the standard input is no terminal; measured on mash 2.3 it behaves as a width of 0 for the table of a file without sketches (every column on its own line)
    // and as an unlimited width otherwise
    private static string PrintColumns(List<string[]> rows, int indent, int spacing, bool narrow = false)
    {
        int columns = rows[0].Length;
        var widths = new int[columns];
        foreach (var row in rows)
            for (int j = 0; j < columns; j++)
                widths[j] = Math.Max(widths[j], Math.Max(row[j].Length, 1));
        var sb = new StringBuilder();
        foreach (var row in rows)
        {
            int offset = 0, target = indent;
            for (int j = 0; j < columns; j++)
            {
                sb.Append(' ', Math.Max(target - offset, 0));
                sb.Append(row[j].Length == 0 ? "-" : row[j]);
                offset = target + row[j].Length;
                if (narrow)
                {
                    if (j < columns - 1)
                        sb.Append('\n');
                    offset = 0;
                }
                else
                {
                    target += widths[j] + spacing;
                }
            }

            sb.Append("\n\n");
        }

        return sb.ToString();
    }

    private static string InfoCounts(MashSketchFileData sketch)
    {
        if (sketch.References.Count == 0)
            throw new ArgumentException("Sketch file contains no sketches");
        if (!sketch.References.Any(r => r.Counts is { Count: > 0 }))
            throw new ArgumentException("Sketch file does not have hash counts. Re-sketch with counting enabled to use this feature.");
        var sb = new StringBuilder("#Sketch\tBin\tFrequency\n");
        foreach (var r in sketch.References)
        {
            var histogram = new SortedDictionary<uint, ulong>();
            foreach (uint c in r.Counts ?? [])
                histogram[c] = histogram.TryGetValue(c, out ulong have) ? have + 1 : 1;
            foreach (var (bin, frequency) in histogram)
                sb.Append(r.Name).Append('\t').Append(bin.ToString(CultureInfo.InvariantCulture)).Append('\t').Append(frequency.ToString(CultureInfo.InvariantCulture)).Append('\n');
        }

        return sb.ToString();
    }

    private static string InfoJson(MashSketchFileData sketch)
    {
        var sb = new StringBuilder("{\n");
        sb.Append("\t\"kmer\" : ").Append(sketch.KmerSize.ToString(CultureInfo.InvariantCulture)).Append(",\n");
        sb.Append("\t\"alphabet\" : \"").Append(sketch.Alphabet).Append("\",\n");
        sb.Append("\t\"preserveCase\" : ").Append(sketch.PreserveCase ? "true" : "false").Append(",\n");
        sb.Append("\t\"canonical\" : ").Append(sketch.Noncanonical ? "false" : "true").Append(",\n");
        sb.Append("\t\"sketchSize\" : ").Append(sketch.MinHashesPerWindow.ToString(CultureInfo.InvariantCulture)).Append(",\n");
        sb.Append("\t\"hashType\" : \"MurmurHash3_x64_128\",\n");
        sb.Append("\t\"hashBits\" : ").Append(sketch.Use64 ? 64 : 32).Append(",\n");
        sb.Append("\t\"hashSeed\" : ").Append(sketch.HashSeed.ToString(CultureInfo.InvariantCulture)).Append(",\n");
        sb.Append(" \t\"sketches\" :\n\t[\n");
        for (int i = 0; i < sketch.References.Count; i++)
        {
            var r = sketch.References[i];
            sb.Append("\t\t{\n");
            sb.Append("\t\t\t\"name\" : \"").Append(r.Name).Append("\",\n");
            sb.Append("\t\t\t\"length\" : ").Append(r.Length.ToString(CultureInfo.InvariantCulture)).Append(",\n");
            sb.Append("\t\t\t\"comment\" : \"").Append(r.Comment).Append("\",\n");
            sb.Append("\t\t\t\"hashes\" :\n\t\t\t[\n");
            for (int j = 0; j < r.Hashes.Count; j++)
            {
                ulong h = sketch.Use64 ? r.Hashes[j] : (uint)r.Hashes[j];
                sb.Append("\t\t\t\t").Append(h.ToString(CultureInfo.InvariantCulture));
                if (j < r.Hashes.Count - 1)
                    sb.Append(',');
                sb.Append('\n');
            }

            sb.Append("\t\t\t]\n");
            sb.Append(i < sketch.References.Count - 1 ? "\t\t},\n" : "\t\t}\n");
        }

        sb.Append("\t]\n}\n");
        return sb.ToString();
    }
}
