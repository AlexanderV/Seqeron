using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Analysis;

/// <summary>The options of <see cref="GenomicAnalyzer.Nucmer"/>; defaults are those of <c>nucmer</c> 3.23.</summary>
public sealed record NucmerOptions
{
    /// <summary>Anchor matches: <see cref="CommonRegionMode.MaximalUnique"/> (<c>--mum</c>), <see cref="CommonRegionMode.MaximalUniqueReference"/> (<c>--mumreference</c> / <c>--mumcand</c>, the default) or <see cref="CommonRegionMode.MaximalExactMatches"/> (<c>--maxmatch</c>).</summary>
    public CommonRegionMode MatchMode { get; init; } = CommonRegionMode.MaximalUniqueReference;

    /// <summary>Minimum length of a single match (<c>-l|--minmatch</c>, default 20).</summary>
    public int MinMatch { get; init; } = 20;

    /// <summary>Maximum gap between two adjacent matches of a cluster (<c>-g|--maxgap</c>, default 90).</summary>
    public long MaxGap { get; init; } = 90;

    /// <summary>Minimum length of a cluster of matches (<c>-c|--mincluster</c>, default 65).</summary>
    public long MinCluster { get; init; } = 65;

    /// <summary>Maximum diagonal difference of two adjacent anchors (<c>-D|--diagdiff</c>, default 5).</summary>
    public int DiagonalDifference { get; init; } = 5;

    /// <summary>Diagonal difference as a fraction of the gap length (<c>-d|--diagfactor</c>, default 0.12).</summary>
    public double DiagonalFactor { get; init; } = 0.12;

    /// <summary>Distance an extension tries to extend poor scoring regions before giving up (<c>-b|--breaklen</c>, default 200).</summary>
    public int BreakLength { get; init; } = 200;

    /// <summary>Enforce absolute banding of the dynamic programming matrix with the diagonal difference (<c>--banded</c>; default off).</summary>
    public bool Banded { get; init; }

    /// <summary>Use only the forward strand of the query (<c>-f|--forward</c>).</summary>
    public bool ForwardOnly { get; init; }

    /// <summary>Use only the reverse complement of the query (<c>-r|--reverse</c>); with <see cref="ForwardOnly"/> both strands (as the script).</summary>
    public bool ReverseOnly { get; init; }

    /// <summary>Extend the clusters (<c>--extend</c>, default; false = <c>--noextend</c>).</summary>
    public bool Extend { get; init; } = true;

    /// <summary>Create the delta (<c>--delta</c>, default; false = <c>--nodelta</c>, the cluster file is produced instead).</summary>
    public bool Delta { get; init; } = true;

    /// <summary>Alignment score optimisation (<c>--optimize</c>, default; false = <c>--nooptimize</c> = <c>postnuc -t</c>).</summary>
    public bool Optimize { get; init; } = true;

    /// <summary>Remove shadowed alignments (<c>--simplify</c>, default; false = <c>--nosimplify</c> = <c>postnuc -s</c>).</summary>
    public bool Simplify { get; init; } = true;

    /// <summary>Also produce the NUCmer 1.1 coords file (<c>-o|--coords</c>: <c>show-coords -r out.delta</c>); needs the delta (not with <c>--nodelta</c>).</summary>
    public bool Coords { get; init; }

    /// <summary>
    /// The output prefix (<c>-p|--prefix</c>, default <c>out</c>): the script writes <c>&lt;prefix&gt;.delta</c> (or <c>&lt;prefix&gt;.cluster</c> with <c>--nodelta</c>) and, with <c>-o</c>, <c>&lt;prefix&gt;.coords</c>; its temporary files
    /// <c>&lt;prefix&gt;.ntref</c> (the joined reference, read by <c>mummer</c>) and <c>&lt;prefix&gt;.mgaps</c> are deleted at the end. The prefix names the keys of <see cref="NucmerResult.Files"/> and the files in <see cref="NucmerResult.Log"/>.
    /// </summary>
    public string Prefix { get; init; } = "out";
}

/// <summary>Output of <see cref="GenomicAnalyzer.Nucmer"/>.</summary>
/// <param name="Delta">Text of <c>out.delta</c> (null with <c>--nodelta</c>).</param>
/// <param name="Cluster">Text of <c>out.cluster</c> (only with <c>--nodelta</c>).</param>
/// <param name="Messages">Warnings the programs write to standard error.</param>
/// <param name="Coords">Text of <c>out.coords</c> (only with <see cref="NucmerOptions.Coords"/>).</param>
/// <param name="Log">Everything the script and its programs write to standard error, in order (<c>1: PREPARING DATA</c>, the warnings of <c>prenuc</c>, <c>2,3: RUNNING mummer AND CREATING CLUSTERS</c>, the reading lines of <c>mummer</c> - the lines with its timings are not reproduced -, <c>4: FINISHING DATA</c>, the warnings of <c>postnuc</c>, <c>5: GENERATING COORDS FILE</c>).</param>
/// <param name="Prefix">The output prefix (<see cref="NucmerOptions.Prefix"/>).</param>
public sealed record NucmerResult(string? Delta, string? Cluster, string Messages, string? Coords = null, string Log = "", string Prefix = "out")
{
    /// <summary>The files the script leaves behind, by name: <c>&lt;prefix&gt;.delta</c> or <c>&lt;prefix&gt;.cluster</c> (with <c>--nodelta</c>) and <c>&lt;prefix&gt;.coords</c> (with <c>-o</c>).</summary>
    public IReadOnlyDictionary<string, string> Files
    {
        get
        {
            var files = new Dictionary<string, string>();
            if (Delta != null)
                files[Prefix + ".delta"] = Delta;
            if (Cluster != null)
                files[Prefix + ".cluster"] = Cluster;
            if (Coords != null)
                files[Prefix + ".coords"] = Coords;
            return files;
        }
    }
}

public static partial class GenomicAnalyzer
{
    #region nucmer (mummer | mgaps | postnuc)

    /// <summary>
    /// MUMmer 3.23 <c>nucmer</c>: <c>prenuc</c> (the reference records joined with an <c>x</c> after each), <c>mummer</c> (<see cref="NucmerOptions.MatchMode"/>, both strands unless restricted, only A/C/G/T match)
    /// piped into <c>mgaps</c> (<see cref="NucmerClusters"/>), then <c>postnuc</c> (<see cref="NucmerPostnuc"/>: banded Smith-Waterman extension of the clusters). Returns the text of the <c>.delta</c> file
    /// (the alignments between every pair of reference record and query record, as the distances between insertions and deletions) - byte for byte what the installed program writes.
    /// </summary>
    /// <param name="referenceFasta">The reference FASTA text (multi-record).</param>
    /// <param name="queryFasta">The query FASTA text (multi-record).</param>
    /// <param name="referenceName">The reference file name printed in the first line of the output (<c>nucmer</c> prints the absolute path).</param>
    /// <param name="queryName">The query file name printed in the first line of the output.</param>
    /// <param name="options">The options (null = defaults).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="FormatException">A program aborts on the input (no reference sequence, unparsable FASTA header ...).</exception>
    public static NucmerResult Nucmer(
        string referenceFasta, string queryFasta, string referenceName = "reference.fa", string queryName = "query.fa",
        NucmerOptions? options = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(referenceName);
        ArgumentNullException.ThrowIfNull(queryName);
        var o = options ?? new NucmerOptions();
        var refs = NucmerPostnuc.ReadFasta(referenceFasta, out string messages);
        var log = new StringBuilder("1: PREPARING DATA\n").Append(messages);
        if (refs.Count == 0)
            throw new FormatException("mummer: empty sequence in multiple fasta file");   // prenuc's file has no sequence
        var prenuc = new StringBuilder();
        foreach (var (_, seq) in refs)
            prenuc.Append(seq).Append('x');

        var strands = StrandSelection.Both;
        if (o.ForwardOnly != o.ReverseOnly)
            strands = o.ForwardOnly ? StrandSelection.Forward : StrandSelection.ReverseComplement;
        var clusterOptions = new NucmerClusterOptions { MaxSeparation = o.MaxGap, MinClusterLength = o.MinCluster, DiagonalDifference = o.DiagonalDifference, DiagonalFraction = o.DiagonalFactor };
        string concatenated = prenuc.ToString();
        string pfx = o.Prefix ?? "out";

        // mummer | mgaps: a header line per query record and strand; mgaps prints it before the first cluster of the block and a '#' line before the others
        var mgaps = new StringBuilder();
        var queries = NucmerPostnuc.ReadFasta(queryFasta, out string queryMessages);
        CheckMummerInput(queryFasta, queries, "qry");
        long residues = MummerMultiSeq.Scan(queryName, queryFasta, MummerMultiSeq.QueryReplacement).Text.Length;
        log.Append("2,3: RUNNING mummer AND CREATING CLUSTERS\n# reading input file \"").Append(pfx).Append(".ntref\" of length ").Append(concatenated.Length.ToString(CultureInfo.InvariantCulture)).Append('\n')
            .Append("# construct suffix tree for sequence of length ").Append(concatenated.Length.ToString(CultureInfo.InvariantCulture)).Append("\n# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n")
            .Append("# reading input file \"").Append(queryName).Append("\" of length ").Append(residues.ToString(CultureInfo.InvariantCulture)).Append('\n')
            .Append("# matching query-file \"").Append(queryName).Append("\"\n# against subject-file \"").Append(pfx).Append(".ntref\"\n");
        foreach (var (id, seq) in queries)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (concatenated.Length == 0 || seq.Length == 0)
                continue;
            var clusters = FindNucmerClusters(concatenated, seq, o.MinMatch, o.MatchMode, clusterOptions, strands, cancellationToken);
            foreach (bool reverse in new[] { false, true })
                MummerMgapsText.AppendBlock(mgaps, "> " + id + (reverse ? " Reverse" : string.Empty), clusters.Where(x => x.IsReverseComplement == reverse).Select(x => x.Cluster));
        }

        var post = NucmerPostnuc.Run(referenceName, queryName, referenceFasta, queryFasta, mgaps.ToString(), new NucmerPostnucOptions
        {
            BreakLength = o.BreakLength,
            Banding = o.Banded ? o.DiagonalDifference : 0,
            Delta = o.Delta,
            Extend = o.Extend,
            ToSequenceEnd = !o.Optimize,
            KeepShadowed = !o.Simplify,
        });
        string? coords = null;
        if (o.Coords)
        {
            if (post.Delta == null)
                throw new FormatException("ERROR: show-coords returned non-zero (there is no delta file with --nodelta)");
            coords = MummerShowCoords.Run(post.Delta, new MummerShowCoordsOptions { SortByReference = true }).Text;
        }

        log.Append("4: FINISHING DATA\n").Append(post.Messages);
        if (o.Coords)
            log.Append("5: GENERATING COORDS FILE\n");
        return new NucmerResult(post.Delta, post.Cluster, messages + queryMessages + post.Messages, coords, log.ToString(), pfx);

    }

    // The errors of mummer's FASTA reader: an empty file, a file without a header line and a file whose only record is empty (several records may be empty)
    private static void CheckMummerInput(string fasta, IReadOnlyList<(string Id, string Sequence)> records, string name)
    {
        if (fasta.Length == 0)
            throw new FormatException("mummer: cannot open file \"" + name + ".fa\" or file \"" + name + ".fa\" is empty");
        if (records.Count == 0)
            throw new FormatException("mummer: no sequences in multiple fasta file");
        if (records.Count == 1 && records[0].Sequence.Length == 0)
            throw new FormatException("mummer: empty sequence in multiple fasta file");
    }

    #endregion
}
