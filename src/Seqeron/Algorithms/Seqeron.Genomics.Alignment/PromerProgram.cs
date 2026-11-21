using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The options of <see cref="PromerProgram.Run"/>; defaults are those of <c>promer</c> 3.07 (MUMmer 3.23); all lengths are in amino acids.</summary>
public sealed record PromerOptions
{
    /// <summary>Anchor matches: <see cref="MummerMatchMode.Mum"/> (<c>--mum</c>), <see cref="MummerMatchMode.MumReference"/> (<c>--mumreference</c> / <c>--mumcand</c>, the default) or <see cref="MummerMatchMode.MaxMatch"/> (<c>--maxmatch</c>).</summary>
    public MummerMatchMode MatchMode { get; init; } = MummerMatchMode.MumReference;

    /// <summary>Minimum length of a single match (<c>-l|--minmatch</c>, default 6).</summary>
    public int MinMatch { get; init; } = 6;

    /// <summary>Maximum gap between two adjacent matches of a cluster (<c>-g|--maxgap</c>, default 30).</summary>
    public long MaxGap { get; init; } = 30;

    /// <summary>Minimum length of a cluster of matches (<c>-c|--mincluster</c>, default 20).</summary>
    public long MinCluster { get; init; } = 20;

    /// <summary>Diagonal difference as a fraction of the gap length (<c>-d|--diagfactor</c>, default 0.11).</summary>
    public double DiagonalFactor { get; init; } = 0.11;

    /// <summary>Distance an extension tries to extend poor scoring regions before giving up (<c>-b|--breaklen</c>, default 60).</summary>
    public int BreakLength { get; init; } = 60;

    /// <summary>The alignment matrix (<c>-x|--matrix</c>): 1 = BLOSUM 45, 2 = BLOSUM 62 (default), 3 = BLOSUM 80.</summary>
    public int Matrix { get; init; } = 2;

    /// <summary>The maximum bookend masking length (<c>-m|--masklen</c>, default 8).</summary>
    public int MaskLength { get; init; } = 8;

    /// <summary>Extend the clusters (<c>--extend</c>, default; false = <c>--noextend</c>).</summary>
    public bool Extend { get; init; } = true;

    /// <summary>Create the delta (<c>--delta</c>, default; false = <c>--nodelta</c>, the cluster file is produced instead).</summary>
    public bool Delta { get; init; } = true;

    /// <summary>Alignment score optimisation (<c>--optimize</c>, default; false = <c>--nooptimize</c> = <c>postpro -t</c>).</summary>
    public bool Optimize { get; init; } = true;

    /// <summary>Also produce the PROmer 1.1 coords file (<c>-o|--coords</c>: <c>show-coords -r out.delta</c>); needs the delta.</summary>
    public bool Coords { get; init; }

    /// <summary>The output prefix (<c>-p|--prefix</c>, default <c>out</c>): it names the temporary files <c>.aaref .aaqry .mgaps</c> shown in the messages and the keys of <see cref="PromerResult.Files"/>.</summary>
    public string Prefix { get; init; } = "out";
}

/// <summary>Output of <see cref="PromerProgram.Run"/>.</summary>
/// <param name="Delta">Text of <c>&lt;prefix&gt;.delta</c> (null with <c>--nodelta</c>).</param>
/// <param name="Cluster">Text of <c>&lt;prefix&gt;.cluster</c> (only with <c>--nodelta</c>).</param>
/// <param name="Messages">What the programs write to standard error: the progress lines of the script, the warnings and the deterministic lines of <c>mummer</c>.</param>
/// <param name="Coords">Text of <c>&lt;prefix&gt;.coords</c> (only with <see cref="PromerOptions.Coords"/>).</param>
/// <param name="AminoReference">The temporary <c>&lt;prefix&gt;.aaref</c> (the reference in six frames), deleted by the script.</param>
/// <param name="AminoQuery">The temporary <c>&lt;prefix&gt;.aaqry</c> (the query in six frames), deleted by the script.</param>
/// <param name="Mgaps">The temporary <c>&lt;prefix&gt;.mgaps</c> (the clusters of the matches), deleted by the script.</param>
public sealed record PromerResult(string? Delta, string? Cluster, string Messages, string? Coords, string AminoReference, string AminoQuery, string Mgaps)
{
    /// <summary>The files the script leaves behind, by name (<c>out.delta</c> or <c>out.cluster</c>, and <c>out.coords</c> with <see cref="PromerOptions.Coords"/>); <see cref="Prefix"/> names them.</summary>
    public IReadOnlyDictionary<string, string> Files
    {
        get
        {
            var files = new Dictionary<string, string>();
            if (Delta != null)
                files[Prefix + ".delta"] = Delta;
            if (Coords != null)
                files[Prefix + ".coords"] = Coords;
            return files;
        }
    }

    /// <summary>The prefix the files are named with (<see cref="PromerOptions.Prefix"/>).</summary>
    public string Prefix { get; init; } = "out";
}

/// <summary>
/// <c>promer</c> of MUMmer 3.23 (the perl script, version 3.07): amino-acid alignments between two multi-FASTA DNA files. The script runs <c>prepro -r</c> / <c>prepro -q</c> (<see cref="PromerPrepro"/>: the six-frame
/// translation with the masking of short stop-to-stop stretches), <c>mummer -mumreference -l 6</c> over the translations (<see cref="MummerProgram"/>), <c>mgaps -l 20 -s 30 -f 0.11</c> (<see cref="NucmerClusters"/>)
/// and <c>postpro</c> (<see cref="PromerPostpro"/>: the extension of the clusters with a BLOSUM matrix), and with <c>-o</c> <c>show-coords -r</c> (<see cref="MummerShowCoords"/>). The result is the <c>.delta</c> file
/// (<c>PROMER</c>): the coordinates are those of the original DNA, one delta value is three nucleotides.
/// </summary>
public static class PromerProgram
{
    /// <summary>Runs <c>promer</c>.</summary>
    /// <param name="referenceFasta">The reference FASTA text (multi-record DNA).</param>
    /// <param name="queryFasta">The query FASTA text (multi-record DNA).</param>
    /// <param name="referenceName">The reference file name printed in the first line of the delta file (<c>promer</c> prints the absolute path).</param>
    /// <param name="queryName">The query file name printed in the first line of the delta file.</param>
    /// <param name="options">The options (null = defaults).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="FormatException">A program of the pipeline aborts on the input (no sequence in a file, a non-positive match length ...).</exception>
    public static PromerResult Run(string referenceFasta, string queryFasta, string referenceName = "reference.fa", string queryName = "query.fa", PromerOptions? options = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(referenceName);
        ArgumentNullException.ThrowIfNull(queryName);
        var o = options ?? new PromerOptions();
        string pfx = o.Prefix ?? "out";
        var msg = new StringBuilder();

        msg.Append("1: PREPARING DATA\n");
        string aaRef = PromerPrepro.Run(referenceFasta, true, o.MaskLength, referenceName, out string m1);
        msg.Append(m1);
        string aaQry = PromerPrepro.Run(queryFasta, false, o.MaskLength, queryName, out string m2);
        msg.Append(m2);

        msg.Append("2,3: RUNNING mummer AND CREATING CLUSTERS\n");
        var mgaps = new StringBuilder();
        try
        {
            var mm = MummerMultiSeq.Scan(pfx + ".aaref", aaRef, 0);
            msg.Append("# reading input file \"").Append(pfx).Append(".aaref\" of length ").Append(mm.Text.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
            msg.Append("# construct suffix tree for sequence of length ").Append(mm.Text.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
            msg.Append("# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n");
            if (o.MinMatch <= 0)
                throw new ArgumentException("mummer: argument for option -l is not a positive integer");
            var queries = MummerMultiSeq.Scan(pfx + ".aaqry", aaQry, 0);
            msg.Append("# reading input file \"").Append(pfx).Append(".aaqry\" of length ").Append(queries.Text.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
            msg.Append("# matching query-file \"").Append(pfx).Append(".aaqry\"\n# against subject-file \"").Append(pfx).Append(".aaref\"\n");
            var index = new MummerIndex(mm.Text);
            var clusterOptions = new NucmerClusterOptions { MaxSeparation = o.MaxGap, MinClusterLength = o.MinCluster, DiagonalDifference = 5, DiagonalFraction = o.DiagonalFactor };
            var mo = new MummerOptions { Mode = o.MatchMode, MinLength = o.MinMatch };
            foreach (var block in MummerProgram.Enumerate(index, queries, mo, o.MinMatch, forward: true, reverse: false, cancellationToken))
            {
                var matches = block.Matches.Select(m => ((long)m.ReferenceStart + 1, (long)m.QueryStart + 1, (long)m.Length)).ToList();
                MummerMgapsText.AppendBlock(mgaps, "> " + block.QueryName, NucmerClusters.Cluster(matches, clusterOptions, cancellationToken));
            }
        }
        catch (Exception e) when (e is FormatException or ArgumentException)
        {
            throw new FormatException("ERROR: mummer and/or mgaps returned non-zero\n" + e.Message);
        }

        msg.Append("4: FINISHING DATA\n");
        var post = PromerPostpro.Run(referenceName, queryName, referenceFasta, queryFasta, mgaps.ToString(), new PromerPostproOptions
        {
            BreakLength = o.BreakLength, Delta = o.Delta, Extend = o.Extend, ToSequenceEnd = !o.Optimize, Matrix = o.Matrix,
        });
        msg.Append(post.Messages);

        string? coords = null;
        if (o.Coords)
        {
            msg.Append("5: GENERATING COORDS FILE\n");
            if (post.Delta == null)
                throw new FormatException("ERROR: show-coords returned non-zero (there is no delta file with --nodelta)");
            coords = MummerShowCoords.Run(post.Delta, new MummerShowCoordsOptions { SortByReference = true }).Text;
        }

        return new PromerResult(post.Delta, post.Cluster, msg.ToString(), coords, aaRef, aaQry, mgaps.ToString()) { Prefix = pfx };
    }
}
