using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// MOTIF-REGULATORY-001 — bacterial σ70 promoters.
/// <see cref="MotifFinder.FindSigma70Promoters"/>: Harley &amp; Reynolds (1987) consensus TTGACA / TATAAT, spacer 15–21
/// (17 ± 1 in 92 %); oracle = independent Python brute force (300 random cases, 358 candidates, identical).
/// <see cref="MotifFinder.PredictSigma70Promoters"/>: La Fleur, Hossain &amp; Salis (2022) Promoter Calculator v1.0;
/// oracle = the authors' Python (hsalis/SalisLabCode Promoter_Calculator_v1_0.py + util.py + .npy coefficients,
/// print statements ported to Python 3, sklearn OneHotEncoder(sparse_output=False)) — 60 random sequences, 4634 per-TSS
/// predictions on both strands, every field bit-identical.
/// </summary>
[TestFixture]
public class MotifFinder_Sigma70Promoters_Tests
{
    // E. coli lac operon control region with the lacUV5 −10 (TATAAT): −35 TTTACA at 68, 18-bp spacer, −10 at 92.
    private const string LacUv5 =
        "TCAGCATTCGAGCTTACGGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCT";

    #region Consensus pairing

    [Test]
    public void FindSigma70_LacUv5_PairsTttacaTataatWithSpacer18()
    {
        var hits = MotifFinder.FindSigma70Promoters(new DnaSequence(LacUv5), maxMismatches35: 1, maxMismatches10: 0).ToList();
        Assert.That(hits, Has.Count.EqualTo(1));
        var h = hits[0];
        Assert.Multiple(() =>
        {
            Assert.That(h.Strand, Is.EqualTo('+'));
            Assert.That(h.Minus35Start, Is.EqualTo(68));
            Assert.That(h.Minus35, Is.EqualTo("TTTACA"));
            Assert.That(h.Minus10Start, Is.EqualTo(92));
            Assert.That(h.Minus10, Is.EqualTo("TATAAT"));
            Assert.That(h.Spacer, Is.EqualTo(18));
            Assert.That(h.Mismatches35, Is.EqualTo(1));
            Assert.That(h.Mismatches10, Is.EqualTo(0));
            Assert.That(h.TotalMismatches, Is.EqualTo(1));
            Assert.That(h.SpacerDeviation, Is.EqualTo(1));
        });
    }

    [Test]
    public void FindSigma70_ConsensusPromoter_SpacerBounds()
    {
        string spacer17 = new('C', 17);
        var seq = new DnaSequence("GG" + "TTGACA" + spacer17 + "TATAAT" + "GG");
        var hit = MotifFinder.FindSigma70Promoters(seq, 0, 0).Single();
        Assert.That((hit.Minus35Start, hit.Minus10Start, hit.Spacer, hit.SpacerDeviation), Is.EqualTo((2, 25, 17, 0)));

        Assert.That(MotifFinder.FindSigma70Promoters(seq, 0, 0, minSpacer: 15, maxSpacer: 16), Is.Empty);
        Assert.That(MotifFinder.FindSigma70Promoters(seq, 0, 0, minSpacer: 17, maxSpacer: 17).Count(), Is.EqualTo(1));

        var far = new DnaSequence("TTGACA" + new string('C', 22) + "TATAAT");
        Assert.That(MotifFinder.FindSigma70Promoters(far, 0, 0), Is.Empty, "spacer 22 > 21");
        Assert.That(MotifFinder.FindSigma70Promoters(far, 0, 0, maxSpacer: 22).Single().Spacer, Is.EqualTo(22));
    }

    [Test]
    public void FindSigma70_MinusStrand_ForwardCoordinatesAndOrder()
    {
        string plus = "TTGACA" + new string('C', 17) + "TATAAT";
        string rc = DnaSequence.GetReverseComplementString(plus); // ATTATA GGG… TGTCAA
        var seq = new DnaSequence("AA" + rc + "AA");
        Assert.That(MotifFinder.FindSigma70Promoters(seq, 0, 0), Is.Empty, "given strand only by default");

        var hit = MotifFinder.FindSigma70Promoters(seq, 0, 0, bothStrands: true).Single();
        Assert.Multiple(() =>
        {
            Assert.That(hit.Strand, Is.EqualTo('-'));
            Assert.That(hit.Minus35, Is.EqualTo("TTGACA"));
            Assert.That(hit.Minus10, Is.EqualTo("TATAAT"));
            Assert.That(hit.Minus10Start, Is.EqualTo(2), "−10 left of −35 on the minus strand");
            Assert.That(hit.Minus35Start, Is.EqualTo(2 + 6 + 17));
            Assert.That(DnaSequence.GetReverseComplementString(seq.Sequence.Substring(hit.Minus35Start, 6)), Is.EqualTo("TTGACA"));
        });
    }

    [Test]
    public void FindSigma70_MismatchTolerance_EqualsBruteForceEnumeration()
    {
        // Every (−35, spacer) pair within the mismatch limits is reported, in −35 / spacer order (brute-force enumeration).
        var seq = new DnaSequence("TTGACATTTTTTTTTTTTTTTTTTATAATTTGACG");
        var hits = MotifFinder.FindSigma70Promoters(seq, 1, 1, 15, 21).ToList();
        var expected = new List<(int, int, int, int)>();
        string s = seq.Sequence;
        for (int i = 0; i + 6 <= s.Length; i++)
        {
            int d35 = Ham(s.Substring(i, 6), "TTGACA");
            if (d35 > 1) continue;
            for (int sp = 15; sp <= 21 && i + 12 + sp <= s.Length; sp++)
            {
                int d10 = Ham(s.Substring(i + 6 + sp, 6), "TATAAT");
                if (d10 <= 1) expected.Add((i, sp, d35, d10));
            }
        }

        Assert.That(hits.Select(h => (h.Minus35Start, h.Spacer, h.Mismatches35, h.Mismatches10)), Is.EqualTo(expected));
        Assert.That(expected, Is.Not.Empty);
    }

    private static int Ham(string a, string b) => a.Zip(b).Count(p => p.First != p.Second);

    [Test]
    public void FindSigma70_Guards()
    {
        var seq = new DnaSequence("ACGT");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.FindSigma70Promoters(null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindSigma70Promoters(seq, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindSigma70Promoters(seq, 0, 7));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindSigma70Promoters(seq, minSpacer: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindSigma70Promoters(seq, minSpacer: 18, maxSpacer: 17));
        });
        Assert.That(MotifFinder.FindSigma70Promoters(new DnaSequence("")), Is.Empty);
    }

    #endregion

    #region Promoter Calculator v1.0

    [Test]
    public void PromoterCalculator_LacUv5_BestForwardStateEqualsReference()
    {
        var preds = MotifFinder.PredictSigma70Promoters(new DnaSequence(LacUv5));
        var fwd = preds.Where(p => p.Strand == '+').ToList();
        var rev = preds.Where(p => p.Strand == '-').ToList();
        Assert.That(fwd, Has.Count.EqualTo(65));
        Assert.That(rev, Has.Count.EqualTo(65));
        Assert.That((fwd[0].Tss, fwd[^1].Tss, rev[0].Tss, rev[^1].Tss), Is.EqualTo((58, 122, 20, 84)));

        var best = fwd.MinBy(p => p.DeltaGTotal)!;
        Assert.Multiple(() =>
        {
            Assert.That(best.Tss, Is.EqualTo(108));
            Assert.That(best.Up, Is.EqualTo("AGCTCACTCATTAGGCACCCCAGG"));
            Assert.That(best.Minus35, Is.EqualTo("TTTACA"));
            Assert.That(best.Spacer, Is.EqualTo("CTTTATGCTTCCGGCTCG"));
            Assert.That(best.Minus10, Is.EqualTo("TATAAT"));
            Assert.That(best.Discriminator, Is.EqualTo("GTGTGGAATT"));
            Assert.That(best.Itr, Is.EqualTo("GTGAGCGGATAACAATTTCA"));
            Assert.That(best.Minus35Start, Is.EqualTo(68));
            Assert.That(best.Minus10Start, Is.EqualTo(92));
            Assert.That(best.DeltaGTotal, Is.EqualTo(-3.044998584365458));
            Assert.That(best.DeltaG10, Is.EqualTo(-1.781524450740899));
            Assert.That(best.DeltaG35, Is.EqualTo(-1.095566862007599));
            Assert.That(best.DeltaGDiscriminator, Is.EqualTo(0.023900602351315935));
            Assert.That(best.DeltaGItr, Is.EqualTo(-0.08406993181999087));
            Assert.That(best.DeltaGExtended10, Is.EqualTo(0.19135161731618874));
            Assert.That(best.DeltaGSpacer, Is.EqualTo(0.116800000000012));
            Assert.That(best.DeltaGUp, Is.EqualTo(0.5424688942938479));
            Assert.That(best.DeltaGBind, Is.EqualTo(-2.026470801138449));
            Assert.That(best.TranscriptionRate, Is.EqualTo(6123.861140216731).Within(1e-9));
        });

        var tss103 = fwd.Single(p => p.Tss == 103);
        Assert.That((tss103.Minus35, tss103.Spacer, tss103.Minus10), Is.EqualTo(("TTTACA", "CTTTATGCTTCCGGC", "TCGTAT")));
        Assert.That(tss103.DeltaGTotal, Is.EqualTo(-1.6789174903686717));
    }

    [Test]
    public void PromoterCalculator_LacUv5_ReverseStrandAndInVitroEqualReference()
    {
        var rev = MotifFinder.PredictSigma70Promoters(new DnaSequence(LacUv5)).Where(p => p.Strand == '-').MinBy(p => p.DeltaGTotal)!;
        Assert.Multiple(() =>
        {
            Assert.That(rev.Tss, Is.EqualTo(63), "reference reverse key n − t");
            Assert.That(rev.Minus35, Is.EqualTo("CACACA"));
            Assert.That(rev.Minus10, Is.EqualTo("TAAAGT"));
            Assert.That(rev.DeltaGTotal, Is.EqualTo(-1.9519008907397968));
            Assert.That(rev.TranscriptionRate, Is.EqualTo(1023.9295941643023).Within(1e-9));
            string fw = LacUv5;
            Assert.That(DnaSequence.GetReverseComplementString(fw.Substring(rev.Minus35Start, 6)), Is.EqualTo(rev.Minus35));
            Assert.That(DnaSequence.GetReverseComplementString(fw.Substring(rev.Tss - 20, 20)), Is.EqualTo(rev.Itr));
        });

        var vitro = MotifFinder.PredictSigma70Promoters(new DnaSequence(LacUv5), bothStrands: false, inVitro: true);
        Assert.That(vitro.All(p => p.Strand == '+'), Is.True);
        Assert.That(vitro.Single(p => p.Tss == 108).TranscriptionRate, Is.EqualTo(504.40616830647446).Within(1e-9));
    }

    [Test]
    public void PromoterCalculator_TooShortSequence_NoPredictions()
    {
        // Shortest complete configuration: 24 + 1 + 6 + 15 + 6 + 6 + 20 = 78 nt.
        string s77 = new string('A', 77);
        Assert.That(MotifFinder.PredictSigma70Promoters(new DnaSequence(s77)), Is.Empty);
        var p78 = MotifFinder.PredictSigma70Promoters(new DnaSequence(s77 + "A"));
        Assert.That(p78.Select(p => (p.Strand, p.Tss)), Is.EqualTo(new[] { ('+', 58), ('-', 20) }));
        Assert.Throws<ArgumentNullException>(() => MotifFinder.PredictSigma70Promoters(null!));
    }

    #endregion
}
