using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The files <c>run-mummer1</c> leaves behind.</summary>
/// <param name="Out">The matches (<c>&lt;prefix&gt;.out</c>).</param>
/// <param name="Gaps">The consistent matches and the other matches (<c>&lt;prefix&gt;.gaps</c>, <see cref="MummerGaps"/>).</param>
/// <param name="Align">The gaps with their alignments (<c>&lt;prefix&gt;.align</c>, <see cref="MummerAnnotate"/>).</param>
/// <param name="ErrorsGaps">The gaps file with the number of errors of each gap (<c>&lt;prefix&gt;.errorsgaps</c>, the <c>witherrors.gaps</c> of the annotation).</param>
/// <param name="Prefix">The prefix the files are named with.</param>
public sealed record RunMummer1Result(string Out, string Gaps, string Align, string ErrorsGaps, string Prefix)
{
    /// <summary>The four files by name (<c>&lt;prefix&gt;.out</c>, <c>.gaps</c>, <c>.align</c>, <c>.errorsgaps</c>).</summary>
    public IReadOnlyDictionary<string, string> Files => new Dictionary<string, string>
    {
        [Prefix + ".out"] = Out, [Prefix + ".gaps"] = Gaps, [Prefix + ".align"] = Align, [Prefix + ".errorsgaps"] = ErrorsGaps,
    };
}

/// <summary>The files and messages <c>run-mummer3</c> leaves behind.</summary>
/// <param name="Out">The maximal unique matches of both strands (<c>&lt;prefix&gt;.out</c>, <see cref="MummerProgram"/>).</param>
/// <param name="Gaps">The clusters of the matches (<c>&lt;prefix&gt;.gaps</c>, <see cref="MummerMgaps"/>).</param>
/// <param name="Align">The clusters extended and combined (<c>&lt;prefix&gt;.align</c>, standard output of <see cref="MummerCombineMums"/>).</param>
/// <param name="ErrorsGaps">The gaps file with the number of errors of each gap (<c>&lt;prefix&gt;.errorsgaps</c>).</param>
/// <param name="Prefix">The prefix the files are named with.</param>
/// <param name="Stdout">The standard output of the script: the three progress lines.</param>
/// <param name="Log">The standard error of the script: the lines of <c>mummer</c> and the summary of <c>combineMUMs</c>.</param>
public sealed record RunMummer3Result(string Out, string Gaps, string Align, string ErrorsGaps, string Prefix, string Stdout, string Log)
{
    /// <summary>The four files by name (<c>&lt;prefix&gt;.out</c>, <c>.gaps</c>, <c>.align</c>, <c>.errorsgaps</c>).</summary>
    public IReadOnlyDictionary<string, string> Files => new Dictionary<string, string>
    {
        [Prefix + ".out"] = Out, [Prefix + ".gaps"] = Gaps, [Prefix + ".align"] = Align, [Prefix + ".errorsgaps"] = ErrorsGaps,
    };
}

/// <summary>
/// The shell scripts of MUMmer 3.23 that run its programs one after the other: <c>run-mummer1</c> (the "severely antiquated" mummer 1 suite: <c>mummer -mum -l 20 [-r]</c> without its first line,
/// <c>gaps</c>, <c>mummer-annotate</c>) and <c>run-mummer3</c> (the basic mummer 3 suite: <c>mummer -mumreference -b -l 20</c>, <c>mgaps -l 100 -f .12 -s 600</c>, <c>combineMUMs -x -e .10 -W</c>).
/// </summary>
public static class MummerRunPrograms
{
    /// <summary><c>run-mummer1 &lt;fasta reference&gt; &lt;fasta query&gt; &lt;prefix&gt; [-r]</c>.</summary>
    /// <param name="referenceFasta">The reference (a single sequence).</param>
    /// <param name="queryFasta">The query (a single sequence).</param>
    /// <param name="referenceName">The file name of the reference as the script receives it (printed in the gaps file and used to find the reference again).</param>
    /// <param name="prefix">The prefix of the output files.</param>
    /// <param name="reverse">The fourth argument <c>-r</c>: only reverse-complement matches (the coordinates are relative to the reversed query). Without it the script passes an empty argument to <c>mummer</c>, which fails.</param>
    /// <exception cref="FormatException">A program of the script fails (<c>mummer</c> without <c>-r</c>, an empty match list: <c>gaps</c> crashes, a second query record: <c>gaps</c> never terminates).</exception>
    public static RunMummer1Result RunMummer1(string referenceFasta, string queryFasta, string referenceName = "reference.fa", string prefix = "out", bool reverse = true)
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(referenceName);
        ArgumentNullException.ThrowIfNull(prefix);
        if (prefix.Length == 0)
            throw new FormatException("Usage: run-mummer1 <fasta reference> <fasta query> <prefix> [-r]");
        if (!reverse)
            throw new FormatException("/usr/bin/mummer: cannot open file \"\" or file \"\" is empty (without -r the script passes an empty argument to mummer; gaps then crashes on the empty list)");
        string mummer = MummerProgram.Run(referenceFasta, [queryFasta], new MummerOptions { Mode = MummerMatchMode.Mum, MinLength = 20, OnlyReverse = true, ReferenceName = referenceName, QueryNames = ["query.fa"] }).Text;
        int nl = mummer.IndexOf('\n', StringComparison.Ordinal);
        string outText = nl < 0 ? string.Empty : mummer[(nl + 1)..];   // tail +2
        string gaps = MummerGaps.Run(referenceName, true, outText);
        var annotated = MummerAnnotate.Run(gaps, queryFasta, new Dictionary<string, string> { [referenceName] = referenceFasta });
        return new RunMummer1Result(outText, gaps, annotated.Text, annotated.WithErrors, prefix);
    }

    /// <summary><c>run-mummer3 &lt;fasta reference&gt; &lt;multi-fasta query&gt; &lt;prefix&gt;</c>.</summary>
    /// <param name="referenceFasta">The reference (only its first record is aligned by <c>combineMUMs</c>).</param>
    /// <param name="queryFasta">The multi-FASTA query file.</param>
    /// <param name="referenceName">The file name of the reference as the script receives it (printed in the messages of <c>mummer</c>).</param>
    /// <param name="queryName">The file name of the query as the script receives it.</param>
    /// <param name="prefix">The prefix of the output files.</param>
    /// <exception cref="FormatException">A program of the script fails (an empty input file, a query without a sequence, an assertion of <c>combineMUMs</c>); the script stops there (<c>sh -e</c>).</exception>
    public static RunMummer3Result RunMummer3(string referenceFasta, string queryFasta, string referenceName = "reference.fa", string queryName = "query.fa", string prefix = "out")
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(referenceName);
        ArgumentNullException.ThrowIfNull(queryName);
        ArgumentNullException.ThrowIfNull(prefix);
        if (prefix.Length == 0)
            throw new FormatException("Usage: run-mummer3 <fasta reference> <multi-fasta query> <prefix>");
        var mummer = MummerProgram.Run(referenceFasta, [queryFasta], new MummerOptions { Mode = MummerMatchMode.MumReference, Both = true, MinLength = 20, ReferenceName = referenceName, QueryNames = [queryName] });
        string gaps = MummerMgaps.Run(mummer.Text, new NucmerClusterOptions { MinClusterLength = 100, DiagonalFraction = .12, MaxSeparation = 600 });
        var combined = MummerCombineMums.Run(referenceFasta, queryFasta, gaps, new MummerCombineMumsOptions
        {
            NoCoverFiles = true, ErrorRate = .10, ErrorsFileName = prefix + ".errorsgaps", ReferenceName = referenceName, QueryName = queryName,
        });
        return new RunMummer3Result(mummer.Text, gaps, combined.Text, combined.ErrorsFile, prefix, "Find MUMs\nDetermine gaps\nAlign gaps\n", mummer.Messages + combined.Log);
    }
}
