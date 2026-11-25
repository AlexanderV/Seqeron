// MIRNA-PRECURSOR-001 / MIRNA-TARGET-001 / MIRNA-CLEAVAGE-001 — closed-form differential killers for
// the pure structural/cleavage kernels of MiRnaAnalyzer: the single-hairpin and largest-enclosed-loop
// dot-bracket parsers, the Turner terminal AU/GU pair predicate, the Drosha/Dicer fold-based cut geometry,
// and canonical seed-site finding. These are observable only through fold/score pipelines that
// re-converge under an internal index mutation, so the parsers/geometry mutants survive. They are
// asserted against hand-derived expected values (independent of the implementation).

using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
[Category("MIRNA-PRECURSOR-001")]
public class MiRnaAnalyzer_StructureCleavageKernels_Tests
{
    // ── IsTerminalPenaltyPair (canonical Turner2004Parameters; MiRnaAnalyzer's internal forwarder was removed
    //    as dead duplication): A-U / U-A / G-U / U-G are penalised; G-C / C-G / A-A are not ──
    [Test]
    public void IsTerminalPenaltyPair_OnlyAuAndGu()
    {
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('A', 'U'), Is.True);
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('U', 'A'), Is.True);
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('G', 'U'), Is.True);
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('U', 'G'), Is.True);
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('G', 'C'), Is.False);
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('C', 'G'), Is.False);
        Assert.That(Seqeron.Genomics.Core.Turner2004Parameters.IsTerminalPenaltyPair('A', 'A'), Is.False);
    }

    // ── TryDescribeSingleHairpin ──
    [Test]
    public void TryDescribeSingleHairpin_SimpleHairpin()
    {
        bool ok = MiRnaAnalyzer.TryDescribeSingleHairpin("(((...)))", out int bp, out int loopStart, out int loopSize);
        Assert.That(ok, Is.True);
        Assert.That(bp, Is.EqualTo(3));
        Assert.That(loopStart, Is.EqualTo(3));
        Assert.That(loopSize, Is.EqualTo(3));
    }

    [Test]
    public void TryDescribeSingleHairpin_BulgedStem_LoopBetweenLastOpenAndFirstClose()
    {
        // ((..((...))..))  opens at 0,1,4,5; first ')' at 10; last '(' at 5; loop = idx 6,7,8 = 3 dots.
        bool ok = MiRnaAnalyzer.TryDescribeSingleHairpin("((..((...))..))", out int bp, out int loopStart, out int loopSize);
        Assert.That(ok, Is.True);
        Assert.That(bp, Is.EqualTo(4), "4 base pairs");
        Assert.That(loopStart, Is.EqualTo(6));
        Assert.That(loopSize, Is.EqualTo(3));
    }

    [Test]
    public void TryDescribeSingleHairpin_BranchedOrUnbalanced_False()
    {
        Assert.That(MiRnaAnalyzer.TryDescribeSingleHairpin("(((...)))(((...)))", out _, out _, out _),
            Is.False, "a '(' after a ')' ⇒ branched");
        Assert.That(MiRnaAnalyzer.TryDescribeSingleHairpin("(((...))", out _, out _, out _),
            Is.False, "opens (3) != closes (2) ⇒ not balanced");
        Assert.That(MiRnaAnalyzer.TryDescribeSingleHairpin("(((...))))", out _, out _, out _),
            Is.False, "opens (3) != closes (4) ⇒ not balanced");
        Assert.That(MiRnaAnalyzer.TryDescribeSingleHairpin("...", out _, out _, out _),
            Is.False, "no stem (opens == 0) ⇒ false");
    }

    // ── LargestEnclosedLoop ──
    [Test]
    public void LargestEnclosedLoop_LongestRunBetweenFirstOpenAndLastClose()
    {
        Assert.That(MiRnaAnalyzer.LargestEnclosedLoop("(((...)))"), Is.EqualTo(3));
        Assert.That(MiRnaAnalyzer.LargestEnclosedLoop("(..(...)..)"), Is.EqualTo(3), "the inner 3-dot run, not the outer 2-dot ones");
        Assert.That(MiRnaAnalyzer.LargestEnclosedLoop("....."), Is.EqualTo(0), "no enclosing pair → 0");
        Assert.That(MiRnaAnalyzer.LargestEnclosedLoop("()"), Is.EqualTo(0), "no enclosed unpaired base");
    }

    // ── PredictDroshaDicerCleavage: fold-based cut geometry (Han 2006 / Park 2011 / miRDeep2 find_star) ──
    // Perfect hairpin: 5-nt flank, 35-bp stem, 8-nt loop, 35-bp stem, 25-nt flank (108 nt); pair k = (5+k, 82−k).
    private static readonly string KernelPri =
        new string('A', 5) + new string('G', 35) + new string('A', 8) + new string('C', 35) + new string('A', 25);
    private static readonly string KernelFold =
        new string('.', 5) + new string('(', 35) + new string('.', 8) + new string(')', 35) + new string('.', 25);

    [Test]
    public void PredictDroshaDicerCleavage_ExactFoldGeometry()
    {
        // junction 5: 12th pair = (16, 71) ⇒ Drosha 5' cut 16, 3' cut 71+2 = 73; 5p [16..37] (22 nt);
        // 3p start = partner(37−2 = 35) = 82−30 = 52 ⇒ 3p [52..73].
        var c = MiRnaAnalyzer.PredictDroshaDicerCleavage(KernelPri, 5, KernelFold);
        Assert.That(c, Is.Not.Null);
        Assert.That(c!.Value.DroshaCut5Prime, Is.EqualTo(16));
        Assert.That(c.Value.MatureStart, Is.EqualTo(16));
        Assert.That(c.Value.MatureEnd, Is.EqualTo(37));
        Assert.That(c.Value.StarStart, Is.EqualTo(52));
        Assert.That(c.Value.StarEnd, Is.EqualTo(73));
        Assert.That(c.Value.DroshaCut3Prime, Is.EqualTo(73));
    }

    [Test]
    public void PredictDroshaDicerCleavage_OutOfRangeBasalJunction_Null()
    {
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(KernelPri, -1, KernelFold), Is.Null);
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(KernelPri, KernelPri.Length, KernelFold), Is.Null);
    }

    [Test]
    public void PredictDroshaDicerCleavage_CnncMotifDetection()
    {
        // CNNC window starts at DroshaCut3Prime + {16,17,18} = {89,90,91} (Auyeung 2013, downstream of the 3p cut).
        var with = KernelPri.ToCharArray(); with[89] = 'C'; with[92] = 'C';
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(new string(with), 5, KernelFold)!.Value.HasCnncMotif, Is.True);
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(KernelPri, 5, KernelFold)!.Value.HasCnncMotif, Is.False);
    }

    // ── FindTargetSites: canonical 8mer seed match at the scan boundary ──
    [Test]
    public void FindTargetSites_Finds8merBuiltFromSeedReverseComplement()
    {
        // miRNA seed = positions 2-8; an 8mer target is [pos8 RC][6mer core = RC of 2-7][A opposite 1].
        // Build the mRNA so the only 8mer sits with its 6mer core at a known offset.
        var miRna = CreateMiRna("test-mir", "AACGUACGUAGCUAGCUAGCUA");
        string seedRC = GetReverseComplement(miRna.SeedSequence); // 7 nt: [RC pos8]...[RC pos2]
        // mRNA = GG + seedRC + A + GG → core at index 3, pos8 match at index 2, A1 at index 9.
        string mrna = "GG" + seedRC + "A" + "GG";

        var sites = FindTargetSites(mrna, miRna).ToList();
        Assert.That(sites.Any(s => s.Type == TargetSiteType.Seed8mer && s.Start == 2),
            $"expected an 8mer at start 2; got [{string.Join(", ", sites.Select(s => $"{s.Type}@{s.Start}"))}]");
    }

    [Test]
    public void FindTargetSites_7merM8_AtScanBoundary()
    {
        // mRNA = "AAAA" + seedRC: the 6mer core sits at index 5 = mrna.Length-6 (the LAST scan index),
        // preceded by seedRC[0] (pos-8 match) with no room for A1 → a 7mer-m8 at start 4. This pins the
        // inclusive scan bound `i <= mrna.Length - 6`.
        var miRna = CreateMiRna("test-mir", "AACGUACGUAGCUAGCUAGCUA");
        string seedRC = GetReverseComplement(miRna.SeedSequence);
        string mrna = "AAAA" + seedRC; // length 11; core at i = 5 = 11-6
        var sites = FindTargetSites(mrna, miRna, minScore: 0.0).ToList(); // boundary test, not score filter
        Assert.That(sites.Any(s => s.Type == TargetSiteType.Seed7merM8 && s.Start == 4),
            $"expected a boundary 7mer-m8 at start 4; got [{string.Join(", ", sites.Select(s => $"{s.Type}@{s.Start}"))}]");
    }

    [Test]
    public void FindTargetSites_NoSeed_NoSites()
    {
        var miRna = CreateMiRna("test-mir", "AACGUACGUAGCUAGCUAGCUA");
        // An mRNA with no seed complement → no canonical sites.
        Assert.That(FindTargetSites(new string('A', 50), miRna).Any(), Is.False);
    }
}
