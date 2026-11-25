using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The input contract shared by the multiple-alignment entry points (audit A2-38 / A2-39):
/// <list type="bullet">
/// <item>a <b>null</b> collection or a null element is an <see cref="ArgumentNullException"/> on EVERY overload (the Seqeron aligners and the
/// ClustalW / MUSCLE / T-Coffee program ports alike; it used to be a <see cref="NullReferenceException"/> on the ClustalW and MUSCLE DNA overloads);</item>
/// <item>malformed text (a user tree, a user matrix) is an <see cref="ArgumentException"/> ("User matrix": <see cref="FormatException"/> by design of
/// <c>ClustalWUserMatrix.Parse</c>, wrapped into an <see cref="ArgumentException"/> where the text is an option of an aligner);</item>
/// <item>an <b>empty</b> member: the Seqeron aligners (<c>MultipleAlign</c>, <c>MultipleAlignProgressive</c>, <c>MultipleAlignIterative</c>,
/// <c>MultipleAlignConsistency</c>, <c>MultipleAlignExact</c>, <c>MultipleAlignAnchored</c>) return an all-gap row for it (like the pairwise
/// aligners, F39), the three PROGRAM ports (<c>MultipleAlignClustalW</c>, <c>MultipleAlignMuscle</c>, <c>MultipleAlignTCoffeeProgram</c>) refuse it with an
/// <see cref="ArgumentException"/> because the programs do (ClustalW cannot read an empty sequence, MUSCLE ignores empty FASTA records, t_coffee stops):
/// all three share one message (<see cref="EmptyMemberMessage"/>).</item>
/// </list>
/// </summary>
public static partial class SequenceAligner
{
    /// <summary>The one message of the program ports for an empty member (the programs refuse empty input; the Seqeron aligners do not).</summary>
    public static string EmptyMemberMessage(string program, int index) =>
        $"{program} cannot align an empty sequence (sequence {index} is empty): the program refuses empty input. " +
        "MultipleAlign, MultipleAlignProgressive, MultipleAlignIterative, MultipleAlignConsistency, MultipleAlignExact and MultipleAlignAnchored accept empty members and return an all-gap row for them.";

    /// <summary>
    /// The one message of the program ports for an EMPTY LIST of sequences (audit A3-4). Installed programs on an empty FASTA file: <c>clustalw</c> "No sequences in file. No alignment!", <c>muscle3</c>
    /// "No sequences in input file", <c>t_coffee</c> exit 1 ("The format of the file is not supported"), <c>kalign</c> exit 1 ("No alignment"), <c>poa</c> exit 1 ("Error reading sequence file"),
    /// <c>probcons</c> exit 1 ("No sequences read."): the ports refuse it with an <see cref="ArgumentException"/>. <c>mafft</c> alone exits 0 with an empty output
    /// ("Warning: Only 0 sequence found."), so <c>MultipleAlignMafft</c> returns an empty result, as do the Seqeron aligners.
    /// </summary>
    internal static string EmptyInputMessage(string program, string programsOwnText) =>
        $"{program} cannot align an empty list of sequences: the program refuses an empty input ({programsOwnText}). " +
        "MultipleAlignMafft and the Seqeron aligners (MultipleAlign, MultipleAlignProgressive, MultipleAlignIterative, MultipleAlignConsistency, MultipleAlignExact, MultipleAlignAnchored) return an empty result.";

    /// <summary>Throws the shared <see cref="EmptyInputMessage"/> <see cref="ArgumentException"/> for an empty list.</summary>
    public static void RequireProgramInput(int count, string program, string programsOwnText, string paramName)
    {
        if (count == 0)
            throw new ArgumentException(EmptyInputMessage(program, programsOwnText), paramName);
    }

    /// <summary>Materialises <paramref name="sequences"/>; a null element is an <see cref="ArgumentNullException"/> naming its index.</summary>
    public static List<T> ToListNoNull<T>(IEnumerable<T> sequences, string paramName) where T : class
    {
        if (sequences is null)
            throw new ArgumentNullException(paramName);
        var list = new List<T>();
        foreach (var s in sequences)
        {
            if (s is null)
                throw new ArgumentNullException(paramName, $"Sequence {list.Count} is null.");
            list.Add(s);
        }

        return list;
    }

    /// <summary>The residue strings of <paramref name="sequences"/> (null element -> <see cref="ArgumentNullException"/>).</summary>
    public static List<string> DnaStringsNoNull(IEnumerable<DnaSequence> sequences, string paramName) =>
        ToListNoNull(sequences, paramName).Select(s => s.Sequence).ToList();

    /// <summary>A null list or a null element is an <see cref="ArgumentNullException"/> naming the parameter (and the index of the element): the one check of the entry points that take rows / arguments as given (audit A5-3).</summary>
    public static void RequireNoNullElements<T>(IReadOnlyList<T?>? items, string paramName, string what = "Row") where T : class
    {
        if (items is null)
            throw new ArgumentNullException(paramName);
        for (int i = 0; i < items.Count; i++)
            if (items[i] is null)
                throw new ArgumentNullException(paramName, $"{what} {i} is null.");
    }

    /// <summary>The argument list and the text "files" of a program runner (<c>KalignProgram</c>, <c>PoaProgram</c>, <c>ProbConsProgram</c>): a null list, null argument or null file text is an <see cref="ArgumentNullException"/> (an argv entry of a real program is never null).</summary>
    public static void RequireProgramCall(IReadOnlyList<string?>? arguments, IReadOnlyDictionary<string, string>? files)
    {
        RequireNoNullElements(arguments, "arguments", "Argument");
        if (files is null)
            return;
        foreach (var kv in files)
            if (kv.Value is null)
                throw new ArgumentNullException(nameof(files), $"The text of file \"{kv.Key}\" is null.");
    }

    /// <summary>Null element -> <see cref="ArgumentNullException"/>, empty element -> the shared <see cref="EmptyMemberMessage"/> <see cref="ArgumentException"/>.</summary>
    public static void RequireProgramMembers(IReadOnlyList<string?> sequences, string program, string paramName)
    {
        for (int i = 0; i < sequences.Count; i++)
        {
            if (sequences[i] is null)
                throw new ArgumentNullException(paramName, $"Sequence {i} is null.");
            if (sequences[i]!.Length == 0)
                throw new ArgumentException(EmptyMemberMessage(program, i), paramName);
        }
    }
}
