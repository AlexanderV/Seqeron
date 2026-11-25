using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / ALIGN-STATS-001 (audit A5-3, B09 F169): the SECONDARY public entry points of the alignment classes refuse a null list / null element / out-of-range value with a typed
/// <see cref="ArgumentException"/> family exception naming the parameter (the library's contract; the main entry points were fixed in F81 / F82). Before: <see cref="NullReferenceException"/> /
/// <see cref="IndexOutOfRangeException"/> / <see cref="KeyNotFoundException"/> found by a reflection fuzzer over every public static method.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_NullContract_Tests
{
    private static readonly string[] Rows = ["ACGTACGTAC", "ACGTCGTAAC", "ACGTACGTTC"];
    private static readonly string[] NullRow = ["ACGTACGTAC", null!, "ACGTACGTTC"];

    private static void ThrowsNull(TestDelegate call, string paramName) =>
        Assert.That(Assert.Throws<ArgumentNullException>(call)!.ParamName, Is.EqualTo(paramName));

    // ------------------------------------------------------------------ (1) the seven MUSCLE stage methods: a null list / null row

    [Test]
    public void MuscleStageMethods_NullListOrRow_IsAnArgumentNullException()
    {
        Assert.Multiple(() =>
        {
            ThrowsNull(() => SequenceAligner.MuscleGuideTree(null!), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleGuideTree(NullRow), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleRefine(null!), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleRefine(NullRow), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleRefineWindows(null!), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleRefineWindows(NullRow), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleScoreAlignment(null!), "alignment");
            ThrowsNull(() => SequenceAligner.MuscleScoreAlignment(NullRow), "alignment");
            foreach (var call in new Func<IReadOnlyList<string>, IReadOnlyList<string>, object>[]
            {
                (a, b) => SequenceAligner.MuscleProfileAlign(a, b), (a, b) => SequenceAligner.MuscleProfileScore(a, b), (a, b) => SequenceAligner.MuscleLocalAlign(a, b),
            })
            {
                ThrowsNull(() => call(null!, Rows), "alignment1");
                ThrowsNull(() => call(NullRow, Rows), "alignment1");
                ThrowsNull(() => call(Rows, null!), "alignment2");
                ThrowsNull(() => call(Rows, NullRow), "alignment2");
            }
        });
    }

    [Test]
    public void MuscleStageMethods_ACharacterAboveTheByteRange_IsAnArgumentException_NotAnIndexOutOfRange()
    {
        // muscle reads one byte per residue of an alignment file: U+4E2D has no byte
        string[] wide = ["ACGTACGTAC", "ACGT中CGTAC"];
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.MuscleScoreAlignment(wide))!;
        Assert.That(ex.Message, Does.Contain("U+4E2D").And.Contain("Row 1, column 4"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MuscleRefine(wide));
        // a Latin-1 letter is a byte: the wildcard, as for any letter that is not a residue
        Assert.DoesNotThrow(() => SequenceAligner.MuscleScoreAlignment(["ACGTACGTAC", "ACGTéCGTAC"]));
        ThrowsNull(() => SequenceAligner.MuscleRun(new[] { "ACGT", null! }), "sequences");
    }

    // ------------------------------------------------------------------ (3) minimap2 Map: a null target

    [Test]
    public void Minimap2Map_NullTargetElement_IsAnArgumentNullException()
    {
        ThrowsNull(() => Minimap2Chaining.Map(["ACGTACGTACGTACGTACGTAC", null!], "ACGTACGTACGTACGTACGTAC"), "targets");
        ThrowsNull(() => Minimap2Chaining.Map(null!, "ACGT"), "targets");
        ThrowsNull(() => Minimap2Chaining.Map(["ACGT"], null!), "query");
    }

    // ------------------------------------------------------------------ (6) ksw2 SimpleMatrix

    [Test]
    public void Ksw2SimpleMatrix_OutOfRangeAlphabetOrScores_IsAnArgumentOutOfRangeException()
    {
        Assert.Multiple(() =>
        {
            foreach (int m in new[] { 0, -1, int.MinValue, 257, int.MaxValue })
                Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.SimpleMatrix(m, 2, -4, 1), $"m = {m}");
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.SimpleMatrix(5, 200, -4, 1))!.ParamName, Is.EqualTo("match"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.SimpleMatrix(5, 2, int.MinValue, 1))!.ParamName, Is.EqualTo("mismatch"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.SimpleMatrix(5, 2, -4, int.MinValue))!.ParamName, Is.EqualTo("ambiguous"));
            // the exact matrix of ksw_gen_simple_mat(5, mat, 2, 4, 1) of minimap2 is untouched
            Assert.That(Ksw2Aligner.SimpleMatrix(5, 2, -4, 1), Is.EqualTo(new sbyte[] { 2, -4, -4, -4, -1, -4, 2, -4, -4, -1, -4, -4, 2, -4, -1, -4, -4, -4, 2, -1, -1, -1, -1, -1, -1 }));
            Assert.That(Ksw2Aligner.SimpleMatrix(1, 2, -4, 1), Is.EqualTo(new sbyte[] { -1 }));
            Assert.That(Ksw2Aligner.SimpleMatrix(256, 127, -128, -127), Has.Length.EqualTo(65536));
        });
    }

    // ------------------------------------------------------------------ (7) AlignmentFormats: a2m without the state annotation

    [Test]
    public void AlignmentFormats_A2mWithoutTheStateAnnotation_IsAnArgumentException()
    {
        var alignment = new BioAlignment
        {
            Records = [new BioAlignmentRecord { Id = "a", Row = "AC-T" }, new BioAlignmentRecord { Id = "b", Row = "ACGT" }],
        };
        var ex = Assert.Throws<ArgumentException>(() => AlignmentFormats.Format(alignment, AlignmentFileFormat.A2m))!;
        Assert.That(ex.ParamName, Is.EqualTo("alignment"));
        Assert.That(ex.Message, Does.Contain("state").And.Contain("KeyError"));
        Assert.Throws<ArgumentException>(() => AlignmentFormats.Write([alignment], AlignmentFileFormat.A2m));
        // a state that is not a string is the same refusal, a state string is written as before
        var number = new BioAlignment { Records = alignment.Records, ColumnAnnotations = new JsonObject { ["state"] = 5 } };
        Assert.Throws<ArgumentException>(() => AlignmentFormats.Format(number, AlignmentFileFormat.A2m));
        var ok = new BioAlignment { Records = alignment.Records, ColumnAnnotations = new JsonObject { ["state"] = "DDDD" } };
        Assert.That(AlignmentFormats.Format(ok, AlignmentFileFormat.A2m), Is.EqualTo(">a\nAC-T\n>b\nACGT\n"));
    }
}
