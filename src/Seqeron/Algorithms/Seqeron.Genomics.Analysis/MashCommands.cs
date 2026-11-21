using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A named sequence record of a <see cref="MashInput"/> (FASTA: the first word of the header is the id, the rest the comment).</summary>
public readonly record struct MashRecord(string Id, string Comment, string Sequence)
{
    /// <summary>The comment string Mash reads for the record when it names the sketch of a file: kseq leaves the text of the previous comment in its buffer when a record has none (<see cref="MashCommands.ParseFasta"/>); null = <see cref="Comment"/>.</summary>
    public string? StaleComment { get; init; }
}

/// <summary>
/// One input file of a Mash command: its name (what the program prints, the path as given) and its records, or - for a sketch file (<c>x.msh</c>, <see cref="MashSketchFile.Read"/>) - its content in <paramref name="SketchFile"/>
/// (<paramref name="Records"/> is then ignored): its sketches are used as they are, with their own names. As in the program, the k-mer size, seed, strand mode and sketch size of the FIRST reference sketch file are inherited by
/// the command (a <c>-k</c> / <c>-s</c> / <c>-S</c> / <c>-n</c> given too is ignored here: the program refuses it); a later sketch file with another alphabet, seed, k, strand mode or a smaller sketch size is skipped (the
/// program warns on the error stream and goes on), a larger sketch is reduced to the current size.
/// </summary>
public sealed record MashInput(string Name, IReadOnlyList<MashRecord> Records, MashSketchFileData? SketchFile = null);

/// <summary>
/// The text commands of Mash 2.3 on top of the sketching and comparison of <see cref="KmerAnalyzer"/> (audit A5-10, B09 F174): <c>mash bounds</c>, <c>mash dist</c> (with <c>-t</c>, <c>-v</c>, <c>-d</c>, <c>-i</c>) and
/// <c>mash triangle</c> (Phylip matrix, edge list <c>-E</c>, <c>-v</c>, <c>-d</c>), written as the program writes them (<c>CommandBounds.cpp</c>, <c>CommandDistance.cpp</c>, <c>CommandTriangle.cpp</c>, marbl/Mash v2.3 at raw.githubusercontent.com).
/// </summary>
public static partial class MashCommands
{
    // C++ `cout << double`: %g with six significant digits ("inf" / "nan" as glibc prints them)
    /// <summary>
    /// The value of a Mash option of the type <c>Size</c> (<c>-b</c>, <c>-g</c>: <c>Command.cpp</c>): a whole number optionally followed by one of <c>kKmMgGtT</c> (x 1000 ... 10^12); the number is read with <c>stof</c> (a float), so large values
    /// are rounded to 24 bits of mantissa.
    /// </summary>
    /// <exception cref="ArgumentException">The text is no whole positive number with an optional unit.</exception>
    public static long ParseSize(string argument, string option = "b")
    {
        ArgumentNullException.ThrowIfNull(argument);
        if (argument.Length == 0)
            return 0;
        long factor = 1;
        char suffix = argument[^1];
        string number = argument;
        if (suffix is < '0' or > '9')
        {
            factor = suffix switch
            {
                'k' or 'K' => 1000L,
                'm' or 'M' => 1_000_000L,
                'g' or 'G' => 1_000_000_000L,
                't' or 'T' => 1_000_000_000_000L,
                _ => throw new ArgumentException($"Unrecognized unit (\"{suffix}\") in argument to -{option}. If specified, unit must be one of [kKmMgGtT]."),
            };
            number = argument[..^1];
        }

        if (!float.TryParse(number.TrimStart(), System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out float value) || value <= 0 || Math.Floor(value) != value)
            throw new ArgumentException($"Argument to -{option} must be a whole number, optionally followed by one of [kKmMgGtT].");
        return checked((long)((double)value * factor));
    }

    private static string Num(double v)
    {
        if (double.IsNaN(v))
            return "nan";
        if (double.IsInfinity(v))
            return v > 0 ? "inf" : "-inf";
        return CFormat.General(v, 6);
    }

    /// <summary>
    /// <c>mash bounds -k K -p P</c>: the table of error bounds of the Mash distance and of the screen distance for nine sketch sizes and eight distances (<c>CommandBounds.cpp</c>): for the binomial model of the number of shared
    /// hashes the smallest <c>x</c> whose cumulative probability exceeds <c>(1 - p) / 2</c> gives the lower Jaccard bound, converted back to a distance and printed as the error (bound minus distance, <c>inf</c> when
    /// no hash is certain to be shared).
    /// </summary>
    /// <param name="k">K-mer size (1 - 32, default 21).</param>
    /// <param name="probability">Confidence of the bounds (0 - 1, default 0.99).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> or <paramref name="probability"/> outside the program's range.</exception>
    public static string Bounds(int k = 21, double probability = 0.99)
    {
        if (k is < 1 or > 32)
            throw new ArgumentOutOfRangeException(nameof(k), k, "k must be in 1..32 (mash bounds -k).");
        if (!(probability > 0 && probability < 1))
            throw new ArgumentOutOfRangeException(nameof(probability), probability, "p must be in (0, 1) (mash bounds -p).");
        int[] sketchSizes = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000];
        double[] dists = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4];
        double q2 = (1.0 - probability) / 2.0;
        var sb = new StringBuilder();
        sb.Append('\n').Append("Parameters (run with -h for details):\n").Append("   k:   ").Append(k).Append('\n').Append("   p:   ").Append(Num(probability)).Append("\n\n");
        for (int screen = 0; screen < 2; screen++)
        {
            sb.Append(screen == 1 ? "\tScreen distance\n" : "\tMash distance\n").Append("Sketch");
            foreach (double d in dists)
                sb.Append('\t').Append(Num(d));
            sb.Append('\n');
            foreach (int s in sketchSizes)
            {
                sb.Append(s);
                foreach (double dist in dists)
                {
                    double m2j = screen == 1 ? Math.Pow(1.0 - dist, k) : 1.0 / (2.0 * Math.Exp(k * dist) - 1.0);
                    // the smallest x < s whose cumulative probability exceeds q2 (s when none does): the CDF is monotone, so bisect instead of the program's linear scan
                    int lo = 0, hi = s;
                    while (lo < hi)
                    {
                        int mid = lo + (hi - lo) / 2;
                        if (StatisticsHelper.BinomialCdf(mid, s, m2j) > q2)
                            hi = mid;
                        else
                            lo = mid + 1;
                    }

                    int x = lo;
                    double je = (double)x / s;
                    double j2m = screen == 1 ? 1.0 - Math.Pow(je, 1.0 / k) : -1.0 / k * Math.Log(2.0 * je / (1.0 + je));
                    sb.Append('\t').Append(Num(j2m - dist));
                }

                sb.Append('\n');
            }

            sb.Append('\n');
        }

        return sb.ToString();
    }

    internal sealed record Named(string Name, string Comment, MinHashSketch Sketch);

    /// <summary>The read-set options of <c>mash sketch</c> / <c>dist</c> / <c>triangle</c> (<c>-r -m -c -g -b</c>).</summary>
    /// <param name="MinCopies"><c>-m</c>: minimum copies of each k-mer (the noise filter, default 1).</param>
    /// <param name="TargetCoverage"><c>-c</c>: sketching concludes when this coverage is reached (0 = read everything).</param>
    /// <param name="GenomeSize"><c>-g</c>: the genome size used as the sketch length (0 = estimated).</param>
    /// <param name="BloomFilterBytes"><c>-b</c>: the size of a Bloom filter that filters out the k-mers seen once (0 = none; cannot be combined with <c>-m</c>; see <see cref="MashBloomFilter"/>).</param>
    public sealed record MashReadsOptions(int MinCopies = 1, double TargetCoverage = 0.0, long GenomeSize = 0, long BloomFilterBytes = 0);

    private static string FileComment(IReadOnlyList<MashRecord> usable, int count)
    {
        // Sketch.cpp sketchFile: the comment of a file's sketch is "name comment" of its first record; with more records "[N seqs] ... [...]"
        string comment = usable[0].Id + " " + (usable[0].StaleComment ?? usable[0].Comment ?? "");
        if (count > 1)
            comment = "[" + count.ToString(System.Globalization.CultureInfo.InvariantCulture) + " seqs] " + comment + " [...]";
        return comment;
    }

    /// <summary>
    /// The sketches of the inputs, as <c>Sketch.cpp</c> <c>sketchFile</c> / <c>sketchFileOrSequence</c> make them: one per file (all its records of at least k symbols together; a read set with <paramref name="reads"/>) or, with <paramref name="individual"/>, one per sequence
    /// of at least k symbols, named by its id. A file without such a record stops the program (exit status 1: "WARNING: All fasta records in input files were shorter than the k-mer size (k)." when records were skipped, else "ERROR: Did not
    /// find fasta records in "input files"." - the program prints the words "input files", not the file name); here it is an <see cref="ArgumentException"/> with that text.
    /// </summary>
    private static List<Named> Sketches(IReadOnlyList<MashInput> inputs, string paramName, int k, int sketchSize, uint seed, bool canonical, bool individual, CancellationToken cancellationToken, MashReadsOptions? reads = null, bool contain = false,
        MashAlphabet? alphabet = null)
    {
        alphabet ??= MashAlphabet.Nucleotide();
        if (reads is not null && individual)
            throw new ArgumentException("Reads mode (-r) is incompatible with -i.", paramName);
        var list = new List<Named>();
        for (int index = 0; index < inputs.Count; index++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            var input = inputs[index] ?? throw new ArgumentNullException(paramName, $"Input {index} is null.");
            if (input.Name is null)
                throw new ArgumentException($"Input {index} has no Name.", paramName);
            if (input.SketchFile is not null)
            {
                list.AddRange(LoadSketchFile(input.SketchFile, k, sketchSize, seed, canonical, contain, alphabet));
                continue;
            }

            if (input.Records is null)
                throw new ArgumentException($"Input {index} ('{input.Name}') has no Records list.", paramName);
            for (int r = 0; r < input.Records.Count; r++)
            {
                if (input.Records[r].Sequence is null)
                    throw new ArgumentException($"Record {r} of input {index} ('{input.Name}') has no Sequence.", paramName);
                if (individual && input.Records[r].Id is null)
                    throw new ArgumentException($"Record {r} of input {index} ('{input.Name}') has no Id (the name of its sketch with -i).", paramName);
            }

            if (individual)
            {
                foreach (var r in input.Records.Where(r => r.Sequence.Length >= k))              // Sketch.cpp: a sequence shorter than k is not sketched
                    list.Add(new Named(r.Id, r.Comment ?? "", KmerAnalyzer.CreateMinHashSketch([r.Sequence], alphabet, k, sketchSize, canonical, seed, cancellationToken)));
                continue;
            }

            var usable = input.Records.Where(r => r.Sequence.Length >= k).ToList();
            if (usable.Count == 0)
                throw new ArgumentException(NoUsableRecordMessage(input, k), paramName);
            if (reads is null)
            {
                list.Add(new Named(input.Name, FileComment(usable, usable.Count), KmerAnalyzer.CreateMinHashSketch(usable.Select(r => r.Sequence), alphabet, k, sketchSize, canonical, seed, cancellationToken)));
                continue;
            }

            var sketch = KmerAnalyzer.CreateReadsMinHashSketch([.. usable.Select(r => r.Sequence)], k, sketchSize, canonical, seed, reads.MinCopies, reads.TargetCoverage, reads.GenomeSize, out int used, cancellationToken, reads.BloomFilterBytes, alphabet);
            list.Add(new Named(input.Name, FileComment(usable, used), sketch));
        }

        return list;
    }

    // initFromFiles with a sketch file: the sketches of a file that agrees with the current parameters (a larger sketch is reduced to the bottom sketchSize hashes), else none (skipped with a warning)
    private static List<Named> LoadSketchFile(MashSketchFileData file, int k, int sketchSize, uint seed, bool canonical, bool contain = false, MashAlphabet? alphabet = null)
    {
        var list = new List<Named>();
        alphabet ??= MashAlphabet.Nucleotide();
        // Sketch.cpp initFromFiles compares the alphabet STRINGS (the file's one folded by the file's own case flag, the current one by the current flag)
        var fileAlphabet = MashAlphabet.FromSketchFile(file.Alphabet, file.PreserveCase);
        if (fileAlphabet.Characters != alphabet.Characters || file.HashSeed != seed || file.KmerSize != k || (!contain && file.MinHashesPerWindow < sketchSize) || file.Noncanonical == canonical)
            return list;
        foreach (var r in file.References)
        {
            var hashes = r.Hashes.Count > sketchSize ? r.Hashes.Take(sketchSize).ToArray() : r.Hashes;
            var counts = r.Counts is not null && r.Counts.Count > sketchSize ? r.Counts.Take(sketchSize).ToArray() : r.Counts;
            list.Add(new Named(r.Name, r.Comment, new MinHashSketch(k, sketchSize, canonical, seed, file.Use64, (long)r.Length, hashes, counts) { AlphabetSize = fileAlphabet.Size }));
        }

        return list;
    }

    // the parameters of the first sketch file are the command's (CommandDistance / CommandTriangle / CommandScreen: initFromFiles without enforced parameters): k, size, seed, strand mode, alphabet and its case flag
    private static void InheritFromSketchFile(IReadOnlyList<MashInput> first, ref int k, ref int sketchSize, ref uint seed, ref bool canonical, ref MashAlphabet? alphabet)
    {
        if (first.Count > 0 && first[0]?.SketchFile is { } f)
        {
            k = f.KmerSize;
            sketchSize = f.MinHashesPerWindow;
            seed = f.HashSeed;
            canonical = !f.Noncanonical;
            alphabet = MashAlphabet.FromSketchFile(f.Alphabet, f.PreserveCase);
        }
    }

    private static void RequireSketchParameters(int k, int sketchSize)
    {
        if (!KmerAnalyzer.IsMashKmerSize(k))
            throw new ArgumentOutOfRangeException(nameof(k), k, $"K must be in 1..{KmerAnalyzer.MaxMashKmerSize} (mash -k).");
        ArgumentOutOfRangeException.ThrowIfLessThan(sketchSize, 1);
    }

    private static bool Passes(MashComparison c, double maxPValue, double maxDistance) => c.Distance <= maxDistance && c.PValue <= maxPValue;

    /// <summary>
    /// <c>mash dist [-t] [-v P] [-d D] [-i] reference... query...</c>: for every query (in input order) and every reference, one line <c>reference TAB query TAB distance TAB p-value TAB shared/total</c> of the pairs whose p-value is at most
    /// <paramref name="maxPValue"/> and distance at most <paramref name="maxDistance"/>; with <paramref name="table"/> a header <c>#query TAB reference...</c> and a row per query with the distances (a failing pair leaves its field blank).
    /// </summary>
    /// <param name="references">The reference files.</param>
    /// <param name="queries">The query files.</param>
    /// <param name="k">K-mer size (<c>-k</c>, 1 - 32).</param>
    /// <param name="sketchSize">Sketch size (<c>-s</c>).</param>
    /// <param name="seed">Hash seed (<c>-S</c>).</param>
    /// <param name="canonical">False for <c>-n</c> (strand preserved).</param>
    /// <param name="individual"><c>-i</c>: one sketch per sequence (named by its id) instead of per file.</param>
    /// <param name="table"><c>-t</c>: table output.</param>
    /// <param name="maxPValue"><c>-v</c> (0 - 1).</param>
    /// <param name="maxDistance"><c>-d</c> (0 - 1).</param>
    /// <param name="cancellationToken">Polled per input, per sketched window block and per compared pair.</param>
    /// <param name="reads">The inputs are read sets (<c>-r -m -c -g</c>); null = genomes.</param>
    /// <param name="showComments"><c>-C</c>: <c>name:comment</c> (the comment of a file sketch is the first record's, <c>[N seqs] ... [...]</c> for several).</param>
    /// <exception cref="ArgumentException">A file without a record of at least k symbols (the program's exit status 1 with "WARNING: All fasta records in input files were shorter than the k-mer size (k)." or "ERROR: Did not find fasta records in &quot;input files&quot;."), a null input / record field, or no reference sequence at all with <paramref name="individual"/>.</exception>
    public static string Dist(IReadOnlyList<MashInput> references, IReadOnlyList<MashInput> queries, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, bool table = false, double maxPValue = 1.0, double maxDistance = 1.0, CancellationToken cancellationToken = default, MashReadsOptions? reads = null, bool showComments = false, MashAlphabet? alphabet = null)
    {
        ArgumentNullException.ThrowIfNull(references);
        ArgumentNullException.ThrowIfNull(queries);
        RequireFraction(maxPValue, nameof(maxPValue), "-v");
        RequireFraction(maxDistance, nameof(maxDistance), "-d");
        InheritFromSketchFile(references, ref k, ref sketchSize, ref seed, ref canonical, ref alphabet);
        canonical &= !(alphabet?.ImpliesNoncanonical ?? false);
        RequireSketchParameters(k, sketchSize);
        var refs = Sketches(references, nameof(references), k, sketchSize, seed, canonical, individual, cancellationToken, reads, false, alphabet);
        var qs = Sketches(queries, nameof(queries), k, sketchSize, seed, canonical, individual, cancellationToken, reads, false, alphabet);
        if (refs.Count == 0)                                                                       // mash dist -i with no reference sequence of k symbols: division by zero in the program (SIGFPE)
            throw new ArgumentException($"No reference sequence is at least k = {k} symbols long (mash dist -i has nothing to compare with and crashes).", nameof(references));
        var sb = new StringBuilder();
        if (table)
        {
            sb.Append("#query");
            foreach (var r in refs)
                sb.Append('\t').Append(r.Name);
            sb.Append('\n');
        }

        foreach (var q in qs)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (table)
                sb.Append(q.Name);
            foreach (var r in refs)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var c = KmerAnalyzer.CompareMinHashSketches(r.Sketch, q.Sketch);
                bool pass = Passes(c, maxPValue, maxDistance);
                if (table)
                {
                    sb.Append('\t');
                    if (pass)
                        sb.Append(Num(c.Distance));
                }
                else if (pass)
                {
                    sb.Append(r.Name).Append(showComments ? ":" + r.Comment : "").Append('\t').Append(q.Name).Append(showComments ? ":" + q.Comment : "").Append('\t').Append(Num(c.Distance)).Append('\t').Append(Num(c.PValue)).Append('\t').Append(c.SharedHashes).Append('/').Append(c.Denominator).Append('\n');
                }
            }

            if (table)
                sb.Append('\n');
        }

        return sb.ToString();
    }

    /// <summary>
    /// <c>mash triangle [-E] [-v P] [-d D] [-i] sequence...</c>: the lower-triangular distance matrix in relaxed Phylip format (the number of sequences, then one row per sequence: its name and the distances to the sequences before it), or with
    /// <paramref name="edgeList"/> (implied by a p-value or distance filter in the program) the edge list <c>name TAB name TAB distance TAB p-value TAB shared/total</c> of the pairs that pass the filters.
    /// </summary>
    /// <param name="inputs">The input files.</param>
    /// <param name="cancellationToken">Polled per input, per sketched window block and per compared pair.</param>
    /// <exception cref="ArgumentException">As <see cref="Dist"/>.</exception>
    public static string Triangle(IReadOnlyList<MashInput> inputs, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, bool edgeList = false, double maxPValue = 1.0, double maxDistance = 1.0, CancellationToken cancellationToken = default, MashReadsOptions? reads = null, bool showComments = false, MashAlphabet? alphabet = null)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        RequireFraction(maxPValue, nameof(maxPValue), "-v");
        RequireFraction(maxDistance, nameof(maxDistance), "-d");
        InheritFromSketchFile(inputs, ref k, ref sketchSize, ref seed, ref canonical, ref alphabet);
        canonical &= !(alphabet?.ImpliesNoncanonical ?? false);
        RequireSketchParameters(k, sketchSize);
        edgeList |= maxPValue < 1.0 || maxDistance < 1.0;
        var sketches = Sketches(inputs, nameof(inputs), k, sketchSize, seed, canonical, individual, cancellationToken, reads, false, alphabet);
        string Label(Named n) => showComments ? n.Comment : n.Name;                               // -C: the comment instead of the id
        if (sketches.Count == 0 && !edgeList)                                                      // mash triangle -i of sequences all shorter than k: the matrix writer aborts (std::out_of_range)
            throw new ArgumentException($"No sequence is at least k = {k} symbols long (mash triangle -i aborts on an empty matrix).", nameof(inputs));
        var sb = new StringBuilder();
        if (!edgeList)
            sb.Append('\t').Append(sketches.Count).Append('\n');
        for (int index = 0; index < sketches.Count; index++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!edgeList)
                sb.Append(Label(sketches[index]));
            for (int i = 0; i < index; i++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var c = KmerAnalyzer.CompareMinHashSketches(sketches[i].Sketch, sketches[index].Sketch);
                if (edgeList)
                {
                    if (Passes(c, maxPValue, maxDistance))
                        sb.Append(Label(sketches[index])).Append('\t').Append(Label(sketches[i])).Append('\t').Append(Num(c.Distance)).Append('\t').Append(Num(c.PValue)).Append('\t').Append(c.SharedHashes).Append('/').Append(c.Denominator).Append('\n');
                }
                else
                {
                    sb.Append('\t').Append(Num(c.Distance));
                }
            }

            if (!edgeList)
                sb.Append('\n');
        }

        return sb.ToString();
    }

    private static void RequireFraction(double value, string name, string option)
    {
        if (!(value >= 0 && value <= 1))
            throw new ArgumentOutOfRangeException(name, value, $"{option} must be in 0..1 (mash).");
    }

    private static double ScreenIdentity(ulong common, ulong denominator, int k)
    {
        // CommandScreen.cpp estimateIdentity: 1 when every hash is shared (avoids -0), 0 for none (avoids inf), else jaccard^(1/k)
        if (common == denominator)
            return 1.0;
        if (common == 0)
            return 0.0;
        return Math.Pow((double)common / denominator, 1.0 / k);
    }

    private static string NoUsableRecordMessage(MashInput input, int k)
    {
        if (input.Records.Count > 0)
            return $"WARNING: All fasta records in input files were shorter than the k-mer size ({k}).";
        return "ERROR: Did not find fasta records in \"input files\".";
    }

    /// <summary>What the pass of <c>mash screen</c> / <c>taxscreen</c> over the mixture leaves: the query sketches, the hash table of the queries (hash -> queries that hold it), the occurrences of those hashes in the mixture and the estimated number of distinct k-mers.</summary>
    internal sealed record ScreenScan(int K, List<Named> Refs, Dictionary<ulong, List<int>> Table, Dictionary<ulong, uint> Counts, ulong SetSize, double KmerSpace)
    {
        /// <summary>The query references as (name, comment).</summary>
        public IReadOnlyList<(string Name, string Comment)> References => [.. Refs.Select(r => (r.Name, r.Comment))];

        /// <summary>The distinct hashes in the order the program's <c>robin_hood::unordered_map</c> (<c>hashCounts</c>) lists them: inserted reference by reference, each in ascending order.</summary>
        public IEnumerable<ulong> HashesInProgramOrder
        {
            get
            {
                var map = new RobinHoodMap();
                var seen = new HashSet<ulong>();
                foreach (var r in Refs)
                {
                    foreach (ulong h in r.Sketch.Hashes)
                    {
                        if (seen.Add(h))
                            map.Insert(h, 0);
                    }
                }

                return map.KeysInIterationOrder();
            }
        }

        /// <summary>The references that hold the hash in the iteration order of the program's <c>robin_hood::unordered_set</c>.</summary>
        public IEnumerable<int> OwnersInProgramOrder(ulong hash)
        {
            var set = new RobinHoodMap();
            foreach (int owner in Table[hash])
                set.Insert((ulong)owner, 0);
            return set.KeysInIterationOrder().Select(x => (int)x);
        }
    }

    internal static ScreenScan ScreenPass(IReadOnlyList<MashInput> queries, IReadOnlyList<MashInput> mixtures, int k, int sketchSize, uint seed, bool canonical, bool individual, MashAlphabet? alphabet, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentNullException.ThrowIfNull(mixtures);
        InheritFromSketchFile(queries, ref k, ref sketchSize, ref seed, ref canonical, ref alphabet);
        alphabet ??= MashAlphabet.Nucleotide();
        canonical &= !alphabet.ImpliesNoncanonical;
        RequireSketchParameters(k, sketchSize);
        var refs = Sketches(queries, nameof(queries), k, sketchSize, seed, canonical, individual, cancellationToken, null, true, alphabet);
        if (refs.Count == 0)
            throw new ArgumentException("No query sequence is at least k symbols long.", nameof(queries));
        if (mixtures.Count == 0)
            throw new ArgumentException("At least one mixture input is required.", nameof(mixtures));
        foreach (var m in mixtures)
        {
            if (m is null || m.Records is null)
                throw new ArgumentException("A mixture input is null or has no Records list.", nameof(mixtures));
        }

        var chunks = MashScreenScan.Chunks(mixtures, k, out long recordCount);
        if (recordCount == 0)
            throw new ArgumentException("ERROR: Did not find sequence records in inputs", nameof(mixtures));

        // the hash table of the queries (hash -> the queries that hold it, ascending) and the occurrence counts of those hashes in the mixture
        var table = new Dictionary<ulong, List<int>>();
        for (int i = 0; i < refs.Count; i++)
        {
            foreach (ulong h in refs[i].Sketch.Hashes)
            {
                if (!table.TryGetValue(h, out var owners))
                    table[h] = owners = [];
                owners.Add(i);
            }
        }

        var counts = new Dictionary<ulong, uint>(table.Count);
        foreach (ulong h in table.Keys)
            counts[h] = 0;
        var bottom = MashScreenScan.Scan(chunks, k, sketchSize, canonical, seed, alphabet.ScreenView(), alphabet.IsProtein, h =>
        {
            if (counts.TryGetValue(h, out uint c))
                counts[h] = c + 1;
        }, cancellationToken);
        double hashBits = alphabet.Uses64BitHash(k) ? 64.0 : 32.0;
        ulong setSize = 0;                                      // MinHashHeap::estimateSetSize: 2^bits * hashes / largest hash, truncated
        if (bottom.Length > 0)
            setSize = (ulong)(Math.Pow(2.0, hashBits) * bottom.Length / (double)bottom[^1]);
        double kmerSpace = alphabet.KmerSpace(k);
        return new ScreenScan(k, refs, table, counts, setSize, kmerSpace);
    }

    /// <summary>
    /// <c>mash screen [-w] [-i I] [-v P] queries.msh mixture...</c> for nucleotide queries (<c>CommandScreen.cpp</c>, marbl/Mash v2.3; audit A5-10b, B09 F185): the hashes of the query sketches are counted in every window of the mixture
    /// (reads or contigs, all records of all mixture files), and for every query one line <c>identity TAB shared/hashes TAB median-multiplicity TAB p-value TAB query-ID TAB query-comment</c> is written, where
    /// identity = (shared / hashes)^(1/k), the multiplicity is the median number of occurrences of the shared hashes in the mixture and the p-value is the binomial tail of the shared count with
    /// r = (estimated distinct k-mers of the mixture) / 4^k. With <paramref name="winnerTakesAll"/> a hash found in several queries is kept only by the query with the best identity (ties: the longer query).
    /// </summary>
    /// <param name="queries">The query files (what <c>mash sketch</c> writes into the .msh: one sketch per file, or per record with <paramref name="individual"/>); their names are the query IDs.</param>
    /// <param name="mixtures">The mixture files.</param>
    /// <param name="k">K-mer size of the query sketch (<c>mash sketch -k</c>).</param>
    /// <param name="sketchSize">Sketch size of the query sketch (<c>-s</c>).</param>
    /// <param name="seed">Hash seed (<c>-S</c>).</param>
    /// <param name="canonical">False when the query sketch was made with <c>-n</c>.</param>
    /// <param name="individual">The query sketch was made with <c>-i</c> (one reference per sequence).</param>
    /// <param name="winnerTakesAll"><c>-w</c>.</param>
    /// <param name="minIdentity"><c>-i</c>: minimum identity (-1 - 1; -1 writes every query, 0 only those with a shared hash).</param>
    /// <param name="maxPValue"><c>-v</c>: maximum p-value (0 - 1).</param>
    /// <param name="cancellationToken">Polled per file and per window block.</param>
    /// <exception cref="ArgumentException">A query without a record of k symbols, no mixture record (the program's exit status 1 with "ERROR: Did not find sequence records in inputs"), a parameter outside the program's range.</exception>
    public static string Screen(IReadOnlyList<MashInput> queries, IReadOnlyList<MashInput> mixtures, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, bool winnerTakesAll = false, double minIdentity = 0.0, double maxPValue = 1.0, CancellationToken cancellationToken = default, MashAlphabet? alphabet = null)
    {
        if (double.IsNaN(minIdentity) || minIdentity < -1.0 || minIdentity > 1.0)
            throw new ArgumentOutOfRangeException(nameof(minIdentity), minIdentity, "-i must be in -1..1 (mash screen).");
        RequireFraction(maxPValue, nameof(maxPValue), "-v");
        var pass = ScreenPass(queries, mixtures, k, sketchSize, seed, canonical, individual, alphabet, cancellationToken);
        k = pass.K;
        var refs = pass.Refs;
        var table = pass.Table;
        var counts = pass.Counts;
        ulong setSize = pass.SetSize;
        double kmerSpace = pass.KmerSpace;

        var shared = new ulong[refs.Count];
        var depths = new List<uint>[refs.Count];
        for (int i = 0; i < refs.Count; i++)
            depths[i] = [];
        foreach (var (hash, count) in counts)
        {
            if (count < 1)
                continue;
            foreach (int owner in table[hash])
            {
                shared[owner]++;
                depths[owner].Add(count);
            }
        }

        if (winnerTakesAll)
        {
            var scores = new double[refs.Count];
            for (int i = 0; i < refs.Count; i++)
                scores[i] = ScreenIdentity(shared[i], (ulong)refs[i].Sketch.Hashes.Count, k);
            Array.Clear(shared);
            foreach (var d in depths)
                d.Clear();
            foreach (var (hash, owners) in table)
            {
                if (counts[hash] < 1)
                    continue;
                double maxScore = 0;
                long maxLength = 0;
                int maxIndex = -1;
                foreach (int owner in owners)
                {
                    if (scores[owner] > maxScore)
                    {
                        maxScore = scores[owner];
                        maxIndex = owner;
                        maxLength = refs[owner].Sketch.Length;
                    }
                    else if (scores[owner] == maxScore && refs[owner].Sketch.Length > maxLength)
                    {
                        maxIndex = owner;
                        maxLength = refs[owner].Sketch.Length;
                    }
                }

                if (maxIndex < 0)
                    maxIndex = owners[0];                 // the program leaves the index uninitialised when no query has a positive length; the first owner is what a zero-initialised stack gives
                shared[maxIndex]++;
                depths[maxIndex].Add(counts[hash]);
            }
        }

        var sb = new StringBuilder();
        for (int i = 0; i < refs.Count; i++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            depths[i].Sort();
            if (shared[i] == 0 && minIdentity >= 0.0)
                continue;
            ulong denominator = (ulong)refs[i].Sketch.Hashes.Count;
            double identity = ScreenIdentity(shared[i], denominator, k);
            if (identity < minIdentity)
                continue;
            double pValue = shared[i] == 0 ? 1.0 : StatisticsHelper.BinomialUpperTail((long)shared[i], (long)denominator, setSize / kmerSpace);
            if (pValue > maxPValue)
                continue;
            sb.Append(Num(identity)).Append('\t').Append(shared[i]).Append('/').Append(denominator).Append('\t').Append(shared[i] > 0 ? depths[i][(int)(shared[i] / 2)] : 0).Append('\t').Append(Num(pValue)).Append('\t')
                .Append(refs[i].Name).Append('\t').Append(refs[i].Comment).Append('\n');
        }

        return sb.ToString();
    }
}
