// Mutation killers strengthening the MFE STRUCTURE TRACEBACK of CalculateMfeStructure (RNA-STRUCT-001).
//
// The canonical RnaSecondaryStructure_MfeStructure_Tests pin the exact dot-bracket of simple hairpins and
// one nested stem, but the Zuker traceback's INTERNAL-LOOP / BULGE cases (V(i,j) = vij2/vij3/vij4) and the
// EXTERNAL/multibranch BIFURCATION split are not exercised by those structures, leaving the traceback's
// matrix-index arithmetic (i*n+j, i-1, j+1, k+1) and per-case energy-match (|target - e| < EPS) mutants
// alive. These designed sequences fold to MFE structures that DO traverse those cases, so a corrupted
// index or case-match yields a different (or invalid) reconstructed structure.
//
// Each expected structure was obtained from the validated Zuker-Stiegler folder and is a valid,
// thermodynamically sensible minimum-free-energy secondary structure; every test also asserts the
// faithful structural invariants (valid balanced dot-bracket; every reconstructed pair is a legal
// Watson-Crick/wobble pair; bracket count == 2 x base-pair count) so the pin is a correctness anchor,
// not a bare snapshot. Theory: Zuker & Stiegler (1981) NAR 9:133-148; NNDB Turner 2004.

using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class RnaSecondaryStructure_MfeTraceback_MutationKillers_Tests
{
    private static void AssertValidStructure(MfeStructure s)
    {
        Assert.That(ValidateDotBracket(s.DotBracket), Is.True, "dot-bracket must be balanced/valid");
        Assert.That(s.DotBracket.Count(c => c == '('), Is.EqualTo(s.BasePairs.Count),
            "one '(' per reconstructed base pair");
        foreach (var (a, b) in s.BasePairs)
            Assert.That(CanPair(s.Sequence[a], s.Sequence[b]), Is.True,
                $"reconstructed pair ({a},{b}) = {s.Sequence[a]}-{s.Sequence[b]} must be a legal pair");
    }

    // Internal loop / bulge: the MFE structure is an outer 4-bp stem enclosing a 1x0 bulge and an inner
    // 2-bp stem, i.e. "((((.((....)).))))". This drives the V(i,j) internal-loop traceback cases that a
    // plain hairpin never reaches.
    [Test]
    public void CalculateMfeStructure_InternalLoopHairpin_TracesExactBulgedStem()
    {
        var s = CalculateMfeStructure("GGGGAACAAAAGUACCCC");

        Assert.That(s.DotBracket, Is.EqualTo("((((.((....)).))))"),
            "outer 4-bp stem + bulge + inner 2-bp stem (internal-loop traceback)");
        AssertValidStructure(s);
    }

    // Two independent hairpins separated by an unpaired linker: the MFE structure is
    // "((((....))))...((((....))))", which forces the EXTERNAL-loop bifurcation split (the k-split that a
    // single hairpin never exercises).
    [Test]
    public void CalculateMfeStructure_TwoHairpins_TracesExternalBifurcation()
    {
        var s = CalculateMfeStructure("GGGGAAAACCCCAAAGGGGAAAACCCC");

        Assert.That(s.DotBracket, Is.EqualTo("((((....))))...((((....))))"),
            "two separate 4-bp hairpins joined by an unpaired linker (external bifurcation)");
        Assert.That(s.BasePairs, Has.Count.EqualTo(8), "4 bp per hairpin x 2 hairpins");
        AssertValidStructure(s);
    }

    // A longer uninterrupted 6-bp stem: pins the stacking-extension traceback over more steps than the
    // canonical 3-/4-bp stems.
    [Test]
    public void CalculateMfeStructure_SixBasePairStem_TracesFullStack()
    {
        var s = CalculateMfeStructure("GGGGGCAAAAGCCCCC");

        Assert.That(s.DotBracket, Is.EqualTo("((((((....))))))"), "6-bp stem closing a 4-nt loop");
        Assert.That(s.BasePairs, Has.Count.EqualTo(6));
        AssertValidStructure(s);
    }
}
