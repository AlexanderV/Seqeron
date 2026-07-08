using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Predictor-agnostic well-formedness checks shared by the H-type
/// (<see cref="RnaPseudoknotPredictFuzzTests"/>) and recursive
/// (<see cref="RnaPseudoknotRecursiveFuzzTests"/>) pseudoknot fuzz suites. Both derive the SAME
/// structural invariants from the pseudoknot algorithm docs; the pieces that do not depend on WHICH
/// predictor produced the <see cref="PseudoknotStructure"/> live here once. Each suite keeps its own
/// <c>AssertWellFormed</c> (which calls its specific predictor) and delegates these shared checks.
/// </summary>
internal static class PseudoknotFuzzAssertions
{
    /// <summary>True iff the pair set contains a crossing pair (i,j),(k,l) with i &lt; k &lt; j &lt; l (§2.1).</summary>
    public static bool HasCrossingPair(IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        for (int a = 0; a < pairs.Count; a++)
        for (int b = a + 1; b < pairs.Count; b++)
        {
            int i = Math.Min(pairs[a].Position1, pairs[a].Position2);
            int j = Math.Max(pairs[a].Position1, pairs[a].Position2);
            int k = Math.Min(pairs[b].Position1, pairs[b].Position2);
            int l = Math.Max(pairs[b].Position1, pairs[b].Position2);
            if (k < i) (i, j, k, l) = (k, l, i, j);
            if (i < k && k < j && j < l) return true;
        }
        return false;
    }

    /// <summary>
    /// Asserts the two-layer dot-bracket is well formed and EXACTLY encodes the base-pair set:
    /// length n; only {(,),[,],.}; both families balanced; every column is '.' iff unpaired and a
    /// matching open/close bracket iff it is the 5'/3' end of a reported pair (§3.2).
    /// </summary>
    public static void AssertDotBracketConsistent(PseudoknotStructure pk, int n)
    {
        string db = pk.DotBracket;
        db.Length.Should().Be(n, "dot-bracket length must equal the sequence length");
        db.Should().MatchRegex(@"^[\(\)\[\]\.]*$", "dot-bracket uses only ( ) [ ] .");

        Balanced(db, '(', ')').Should().BeTrue("the () family must be balanced");
        Balanced(db, '[', ']').Should().BeTrue("the [] family must be balanced");

        var open = new bool[n];
        var close = new bool[n];
        foreach (var (a, b) in pk.BasePairs)
        {
            int i = Math.Min(a, b), j = Math.Max(a, b);
            open[i] = true;
            close[j] = true;
        }
        for (int p = 0; p < n; p++)
        {
            bool isOpen = db[p] is '(' or '[';
            bool isClose = db[p] is ')' or ']';
            bool isDot = db[p] == '.';
            if (open[p])
                isOpen.Should().BeTrue($"column {p} is a 5' end of a pair → must be an opening bracket");
            else if (close[p])
                isClose.Should().BeTrue($"column {p} is a 3' end of a pair → must be a closing bracket");
            else
                isDot.Should().BeTrue($"column {p} is unpaired → must be '.'");
        }
    }

    /// <summary>Asserts the documented empty pseudoknot-free result: no pairs, all dots, ΔG 0, no knot.</summary>
    public static void AssertEmptyStructure(PseudoknotStructure pk)
    {
        pk.BasePairs.Should().BeEmpty("a degenerate / too-short input yields no base pairs (§6.1)");
        pk.HasPseudoknot.Should().BeFalse("no pseudoknot is possible (§6.1)");
        pk.FreeEnergy.Should().Be(0.0, "the empty structure has ΔG = 0 (§6.1)");
        pk.DotBracket.Should().Be(new string('.', pk.Sequence.Length), "all positions are unpaired dots");
    }

    private static bool Balanced(string s, char open, char close)
    {
        int depth = 0;
        foreach (char c in s)
        {
            if (c == open) depth++;
            else if (c == close) depth--;
            if (depth < 0) return false;
        }
        return depth == 0;
    }
}
