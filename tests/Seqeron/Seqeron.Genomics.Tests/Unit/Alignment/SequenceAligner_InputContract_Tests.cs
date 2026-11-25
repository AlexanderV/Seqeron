using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-38, A2-39, A2-41, A2-42, A2-43; B09 F81, F82, F83, F84, F85): the input contract of the multiple
/// aligners and the overload-resolution contract of the scoring argument.
///
/// A2-38: null collection / null element = <see cref="ArgumentNullException"/> on EVERY overload (the ClustalW and MUSCLE DNA overloads used to throw
///   <see cref="NullReferenceException"/>); malformed text of a user tree / a T-Coffee method parameter = <see cref="ArgumentException"/>, never <see cref="FormatException"/>.
/// A2-39: an empty member is refused by the three program ports with ONE message (the programs refuse it), the Seqeron aligners return an all-gap row.
/// A2-41: the T-Coffee proba pair-HMM and the MUSCLE profile DP refuse too large inputs with an <see cref="ArgumentException"/> (no int overflow / OOM).
/// A2-42: an exception inside the parallel pair stage reaches the caller with its original type for every k, never wrapped in <see cref="AggregateException"/>.
/// A2-43: <c>null</c> for the scoring argument compiles (no CS0121) and selects the default scoring.
/// </summary>
[TestFixture]
public class SequenceAligner_InputContract_Tests
{
    // ------------------------------------------------------------------ A2-38 null elements

    private static readonly string[] WithNullString = { "ACGTACGT", null! };
    private static readonly DnaSequence[] WithNullDna = { new("ACGTACGT"), null! };
    private static readonly ProteinSequence[] WithNullProtein = { new("MKVLA"), null! };

    [Test]
    public void NullElement_IsArgumentNullException_OnEveryMuscleOverload()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle(WithNullString));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle(WithNullProtein));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle((IEnumerable<DnaSequence>)null!));
    }

    // ------------------------------------------------------------------ A2-38 malformed text

    [Test]
    public void MalformedUserTreeBranchLength_IsArgumentException_NotFormatException()
    {
        var seqs = new[] { "ACGTACGTAC", "ACGTACGAAC", "ACGGACGTAC" };
        var ex = Assert.Throws<ArgumentException>(() =>
            SequenceAligner.MultipleAlignMuscle(seqs, new MuscleOptions { UserTree = "(0:abc,1:1,2:1);" }));
        Assert.That(ex!.Message, Does.Contain("abc"));
    }

    [Test]
    public void MuscleProfileDp_RefusesAnOversizedProfilePair_WithAClearArgumentException()
    {
        Assert.That(SequenceAligner.MuscleProfileMatrixMaxCells, Is.EqualTo(1_000_000_000L));
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.MuscleEngine.RequireMuscleProfileMatrixFits(40_000, 40_000));
        Assert.That(ex!.Message, Does.Contain("MuscleProfileMatrixMaxCells"));
        Assert.DoesNotThrow(() => SequenceAligner.MuscleEngine.RequireMuscleProfileMatrixFits(30_000, 30_000));
    }

    // ------------------------------------------------------------------ A2-42 AggregateException

    [TestCase(2)]
    [TestCase(3)]
    [TestCase(5)]
    public void ExceptionInTheParallelPairStage_KeepsItsType_ForEveryK(int k)
    {
        // match 1e6 x (m + n + 2) = 4.02e8 >= 2^28: the range guard fires INSIDE the pair dynamic program of a worker thread
        string s = new('A', 200);
        var seqs = Enumerable.Range(0, k).Select(_ => new DnaSequence(s)).ToArray();
        var scoring = new ScoringMatrix(1_000_000, -1, -2, -1);
        foreach (var gap in new[] { MsaGapModel.Linear, MsaGapModel.Affine })
        {
            var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlign(seqs, scoring, new MsaOptions { GapModel = gap }), $"k={k} {gap}");
            Assert.That(ex!.Message, Does.StartWith("Score range too large"));
        }
    }

    [Test]
    public void MaxDegreeOfParallelism_ChangesNothing_AndIsValidated()
    {
        var seqs = new[] { "ACGTACGTAC", "ACGTACGAAC", "ACGGACGTAC", "ACGTTCGTAC", "TCGTACGTAC" }.Select(x => new DnaSequence(x)).ToArray();
        var reference = SequenceAligner.MultipleAlign(seqs);
        foreach (int dop in new[] { 1, 2, 8, -1 })
        {
            var r = SequenceAligner.MultipleAlign(seqs, null, new MsaOptions { MaxDegreeOfParallelism = dop });
            Assert.That(r.AlignedSequences, Is.EqualTo(reference.AlignedSequences), $"dop {dop}");
            Assert.That(r.TotalScore, Is.EqualTo(reference.TotalScore));
        }

        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlign(seqs, null, new MsaOptions { MaxDegreeOfParallelism = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.MultipleAlign(seqs, null, new MsaOptions { MaxDegreeOfParallelism = -2 }));
    }

    [Test]
    public void NullScoring_TheTypedMatrixOverloadsStillBind()
    {
        var matrix = new MatrixScoring(SubstitutionMatrix.Blosum62, -11, -1);
        // Biopython PairwiseAligner (BLOSUM62, open = extend = -1): global 17 (MKV-A: 5 + 5 + 4 - 1 + 4), local 17; Waterman-Eggert scores 17, 2 (L / M), 1 (independent Python WE of the same matrix)
        Assert.That(SequenceAligner.GlobalAlign("MKVLA", "MKVA", matrix).Score, Is.EqualTo(17));
        var top = SequenceAligner.LocalAlignTopN("MKVLA", "MKVA", matrix, 3, 1);
        Assert.That(top.Select(r => r.ExactScore), Is.EqualTo(new[] { 17.0, 2.0, 1.0 }));
        Assert.That((top[0].AlignedSequence1, top[0].AlignedSequence2), Is.EqualTo(("MKVLA", "MKV-A")));
    }

    [Test]
    public void EveryScoringMatrixOverloadWithAMatrixScoringTwin_CarriesTheResolutionPriority()
    {
        var methods = typeof(SequenceAligner).GetMethods(BindingFlags.Public | BindingFlags.Static);
        var violations = new List<string>();
        foreach (var m in methods)
        {
            var ps = m.GetParameters();
            for (int i = 0; i < ps.Length; i++)
            {
                if (ps[i].ParameterType != typeof(ScoringMatrix))
                    continue;
                bool hasTwin = methods.Any(t => t.Name == m.Name && t != m && t.GetParameters().Length == ps.Length
                    && t.GetParameters()[i].ParameterType == typeof(MatrixScoring));
                bool hasPriority = m.GetCustomAttribute<System.Runtime.CompilerServices.OverloadResolutionPriorityAttribute>() is { Priority: >= 1 };
                if (hasTwin && !hasPriority)
                    violations.Add($"{m.Name}({string.Join(", ", ps.Select(p => p.ParameterType.Name))})");
            }
        }

        Assert.That(violations, Is.Empty, "ScoringMatrix? overloads with a MatrixScoring twin need [OverloadResolutionPriority(1)] so that `null` binds to the default scoring");
    }
}
