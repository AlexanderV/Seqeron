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
    public void NullElement_IsArgumentNullException_OnEveryClustalWOverload()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW(WithNullString));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW(WithNullProtein));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW(WithNullString, new ClustalWProteinOptions()));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW((IEnumerable<string>)null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW((IEnumerable<DnaSequence>)null!));
    }

    [Test]
    public void NullElement_IsArgumentNullException_OnEveryMuscleOverload()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle(WithNullString));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle(WithNullProtein));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMuscle((IEnumerable<DnaSequence>)null!));
    }

    [Test]
    public void NullElement_IsArgumentNullException_OnTCoffeeProgram_AndTheSeqeronAligners()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(WithNullString));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignProgressive(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignIterative(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignConsistency(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignExact(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignAnchored(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign(WithNullProtein, new MatrixScoring(SubstitutionMatrix.Blosum62, -2, -1)));
    }

    [Test]
    public void NullElement_Message_NamesTheIndex()
    {
        var ex = Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignClustalW(WithNullString));
        Assert.That(ex!.Message, Does.Contain("Sequence 1 is null"));
        Assert.That(ex.ParamName, Is.EqualTo("sequences"));
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
    public void MalformedTCoffeeMethodAndWeightValues_AreArgumentException_NotFormatException()
    {
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-method slow_pair@EP@GOP@abc"));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-method slow_pair@EP@GEP@1.5"));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-seq_weight @1,x,3"));
    }

    [Test]
    public void ClustalWUserMatrixParse_StaysFormatException_ByDocumentedDesign()
    {
        Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("not a matrix", "text"));
    }

    // ------------------------------------------------------------------ A2-39 empty members

    [Test]
    public void EmptyMember_ThreeProgramPorts_ShareOneMessage()
    {
        var withEmpty = new[] { "ACGTACGT", "", "ACGTACCT" };
        var c = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignClustalW(withEmpty));
        var m = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMuscle(withEmpty));
        var t = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(withEmpty));
        Assert.That(c!.Message, Does.StartWith("ClustalW cannot align an empty sequence (sequence 1 is empty)"));
        Assert.That(m!.Message, Does.StartWith("MUSCLE cannot align an empty sequence (sequence 1 is empty)"));
        Assert.That(t!.Message, Does.StartWith("t_coffee cannot align an empty sequence (sequence 1 is empty)"));
        foreach (var ex in new[] { c, m, t })
            Assert.That(ex.Message, Does.Contain("MultipleAlign, MultipleAlignProgressive").And.Contain("return an all-gap row"));
    }

    [Test]
    public void EmptyMember_SeqeronAligners_ReturnAnAllGapRow()
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence(""), new DnaSequence("ACGTACCT") };
        foreach (var r in new[]
        {
            SequenceAligner.MultipleAlign(seqs), SequenceAligner.MultipleAlignProgressive(seqs),
            SequenceAligner.MultipleAlignIterative(seqs), SequenceAligner.MultipleAlignConsistency(seqs),
        })
        {
            Assert.That(r.AlignedSequences[1], Is.EqualTo(new string('-', r.AlignedSequences[0].Length)));
        }
    }

    // ------------------------------------------------------------------ A2-41 size guards

    [Test]
    public void TCoffeeProba_RefusesAnOversizedPair_BeforeAllocating()
    {
        // 9 000 x 9 000 protein: (9 001)^2 x 5 = 405 M floats > 2^28; nothing is allocated (the guard runs first)
        string big = new('A', 9000);
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.TcProbaPairEntries(big, big, TcSeqKind.Protein));
        Assert.That(ex!.Message, Does.Contain("proba_pair").And.Contain("slow_pair"));
        Assert.That(SequenceAligner.TCoffeeProbaMaxMatrixFloats, Is.EqualTo(1L << 28));
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

    // ------------------------------------------------------------------ A2-43 null scoring overloads

    [Test]
    public void NullScoring_Compiles_AndSelectsTheDefault_OnEveryPairwiseFamily()
    {
        const string a = "ACGTACGTTAGC", b = "ACGTCGTTAGC";
        // each line below was CS0121 (ambiguous between ScoringMatrix? and MatrixScoring) before OverloadResolutionPriority
        Assert.That(SequenceAligner.GlobalAlign(a, b, null).Score, Is.EqualTo(SequenceAligner.GlobalAlign(a, b).Score));
        Assert.That(SequenceAligner.GlobalAlignAffine(a, b, null).Score, Is.EqualTo(SequenceAligner.GlobalAlignAffine(a, b).Score));
        Assert.That(SequenceAligner.LocalAlign(a, b, null).Score, Is.EqualTo(SequenceAligner.LocalAlign(a, b).Score));
        Assert.That(SequenceAligner.LocalAlignAffine(a, b, null).Score, Is.EqualTo(SequenceAligner.LocalAlignAffine(a, b).Score));
        Assert.That(SequenceAligner.SemiGlobalAlign(a, b, null).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b).Score));
        Assert.That(SequenceAligner.SemiGlobalAlignAffine(a, b, null).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b).Score));
        Assert.That(SequenceAligner.SemiGlobalAlign(a, b, null, EndGaps.All).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, SequenceAligner.SimpleDna, EndGaps.All).Score));
        Assert.That(SequenceAligner.LocalAlignTopN(a, b, null, 3, 1), Has.Count.GreaterThan(0));
        Assert.That(SequenceAligner.LocalAlignAffineTopN(a, b, null, 3, 1), Has.Count.GreaterThan(0));
        Assert.That(SequenceAligner.LocalAlignCoOptimal(a, b, null, 5), Has.Count.GreaterThan(0));
        Assert.That(SequenceAligner.GlobalAlignCoOptimal(a, b, null, 5), Has.Count.GreaterThan(0));
        Assert.That(SequenceAligner.GlobalAlignLinearSpace(a, b, null).Score, Is.EqualTo(SequenceAligner.GlobalAlign(a, b).Score));
        Assert.That(SequenceAligner.LocalAlignTopNLinearSpace(a, b, null, 3, 1), Has.Count.GreaterThan(0));
        var dna = new[] { new DnaSequence(a), new DnaSequence(b) };
        Assert.That(SequenceAligner.MultipleAlign(dna, null, null).AlignedSequences, Has.Length.EqualTo(2));
        Assert.That(SequenceAligner.MultipleAlignProgressive(dna, null, null).AlignedSequences, Has.Length.EqualTo(2));
        Assert.That(SequenceAligner.MultipleAlignClustalW(new[] { a, b }, null, null).AlignedSequences, Has.Length.EqualTo(2));
        Assert.That(SequenceAligner.MultipleAlignExact(dna, null, null).AlignedSequences, Has.Length.EqualTo(2));
    }

    [Test]
    public void NullScoring_TheTypedMatrixOverloadsStillBind()
    {
        var matrix = new MatrixScoring(SubstitutionMatrix.Blosum62, -11, -1);
        Assert.That(SequenceAligner.GlobalAlign("MKVLA", "MKVA", matrix).Score, Is.Not.EqualTo(int.MinValue));
        Assert.That(SequenceAligner.LocalAlignTopN("MKVLA", "MKVA", matrix, 3, 1), Has.Count.GreaterThan(0));
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
