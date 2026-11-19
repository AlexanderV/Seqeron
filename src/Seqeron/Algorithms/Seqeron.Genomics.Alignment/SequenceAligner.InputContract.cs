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
    internal static string EmptyMemberMessage(string program, int index) =>
        $"{program} cannot align an empty sequence (sequence {index} is empty): the program refuses empty input. " +
        "MultipleAlign, MultipleAlignProgressive, MultipleAlignIterative, MultipleAlignConsistency, MultipleAlignExact and MultipleAlignAnchored accept empty members and return an all-gap row for them.";

    /// <summary>Materialises <paramref name="sequences"/>; a null element is an <see cref="ArgumentNullException"/> naming its index.</summary>
    internal static List<T> ToListNoNull<T>(IEnumerable<T> sequences, string paramName) where T : class
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
    internal static List<string> DnaStringsNoNull(IEnumerable<DnaSequence> sequences, string paramName) =>
        ToListNoNull(sequences, paramName).Select(s => s.Sequence).ToList();

    /// <summary>Null element -> <see cref="ArgumentNullException"/>, empty element -> the shared <see cref="EmptyMemberMessage"/> <see cref="ArgumentException"/>.</summary>
    internal static void RequireProgramMembers(IReadOnlyList<string?> sequences, string program, string paramName)
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
