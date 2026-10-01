// MOTIF-SHARED-001 (audit round 2, group G4) — RSAT oligo-analysis -return mseq,proba with -oneN / -onedeg.
// Evidence: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md; TestSpec: tests/TestSpecs/MOTIF-DISCOVER-001.md
// Source: rsa-tools/rsat-code master 10043f2 — perl-scripts/oligo-analysis v1.169: CountOligos mseq block (lines 635-657:
//         first occurrence per sequence, W and W' both credited with -2str), Degenerate (lines 1070-1105, incl.
//         $deg_mseq{$deg} += $patterns{$pattern_seq}->{mseq}), CalcExpected (exp_ms = S·(1 - (1 - p)^(nb_pos/S))),
//         CalcProba (ms_P = binomial_boe(exp_ms/S, S, mseq), ms_E = NPO·ms_P), lib/RSA.disco.lib NbPossibleOligos.
// Reference values: oracle copy oligo-analysis-g4 = the F31 oligo-analysis-fixdeg copy (Degenerate repaired: RSAT 1.169
//         -onedeg/-oneN return no rows) plus the matching-sequence union: CountOligos records the set of matching sequences
//         of every word, Degenerate gives a degenerate word the union of its words' sets. The unmodified RSAT sum (line 1090)
//         counts a sequence once per matching word: on the 3 sequences below, -onedeg -l 3, AYG gets mseq 4 > 3 and RSAT stops
//         ("&RSAT::stats::binomial() Successes (4) cannot be higher than trials (3)").
//         C# vs oracle: 300 random runs (2-7 sequences, k 1-5, -oneN/-onedeg, -1str/-2str, -pseudo), 68,997 patterns:
//         mseq exact, exp_freq/exp_ms/ms_P/ms_E <= 1.3e-13 relative, ms_sig <= 5.7e-14.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_SharedMotifsDegenerate_Tests
{
    private const double Rel = 1e-10; // percent units

    private static readonly DnaSequence[] D3 =
    {
        new("ACGTACGGATCC"), new("ATGCATGAAC"), new("ACGATGTT"),
    };

    private static SignificantSharedMotif Get(SharedMotifAnalysisResult r, string word) => r.Motifs.Single(m => m.Sequence == word);

    [Test]
    public void OneDegenerate_SingleStrand_MatchesRsatOracle()
    {
        // oracle: oligo-analysis-g4 -l 3 -1str -onedeg -return mseq,proba -lth mseq 1; NPO = 3·11·4^2 = 528.
        var r = MotifFinder.FindSharedMotifs(D3, 3, 1, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Single, 0.0,
            OligoDegeneracy.OneDegenerate);
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(528));
            Assert.That(r.Degeneracy, Is.EqualTo(OligoDegeneracy.OneDegenerate));
            Assert.That(r.PossiblePositions, Is.EqualTo(24));
            var ayg = Get(r, "AYG"); // ACG in s1, s3; ATG in s2, s3 — union 3 (RSAT's sum: 4 > 3 sequences)
            Assert.That(ayg.SequenceIndices, Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(ayg.ExpectedFrequency, Is.EqualTo(0.032666666666666663).Within(Rel).Percent);
            Assert.That(ayg.ExpectedMatchingSequences, Is.EqualTo(0.69998598816389546).Within(Rel).Percent);
            Assert.That(ayg.MatchingSequenceProbability, Is.EqualTo(0.01270294085234148).Within(Rel).Percent);
            Assert.That(ayg.MatchingSequenceEValue, Is.EqualTo(6.7071527700363012).Within(Rel).Percent);
            Assert.That(ayg.MatchingSequenceSignificance, Is.EqualTo(-0.82653819845756793).Within(Rel).Percent);
            var acn = Get(r, "ACN");
            Assert.That(acn.SequenceIndices, Has.Count.EqualTo(2));
            Assert.That(acn.ExpectedMatchingSequences, Is.EqualTo(1.3212545710048791).Within(Rel).Percent);
            Assert.That(acn.MatchingSequenceProbability, Is.EqualTo(0.41105031540287185).Within(Rel).Percent);
            var atn = Get(r, "ATN");
            Assert.That(atn.SequenceIndices, Has.Count.EqualTo(3));
            Assert.That(atn.MatchingSequenceProbability, Is.EqualTo(0.085427115865445352).Within(Rel).Percent);
            Assert.That(r.Motifs.All(m => m.SequenceIndices.Count <= D3.Length), Is.True);
        });
    }

    [Test]
    public void OneN_BothStrands_MatchesRsatOracle()
    {
        // oracle: -l 3 -2str -oneN -return mseq,proba -lth mseq 2; NPO = 3·1·16 − (48 − 0)/2 = 24.
        var r = MotifFinder.FindSharedMotifs(D3, 3, 2, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Both, 0.0,
            OligoDegeneracy.OneN);
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(24));
            var ang = Get(r, "ANG");
            Assert.That(ang.ReverseComplement, Is.EqualTo("CNT"));
            Assert.That(ang.SequenceIndices, Has.Count.EqualTo(3));
            Assert.That(ang.ExpectedFrequency, Is.EqualTo(0.12444444444444444).Within(Rel).Percent);
            Assert.That(ang.ExpectedMatchingSequences, Is.EqualTo(1.9639256569777699).Within(Rel).Percent);
            Assert.That(ang.MatchingSequenceProbability, Is.EqualTo(0.28055070694211726).Within(Rel).Percent);
            Assert.That(ang.MatchingSequenceEValue, Is.EqualTo(6.7332169666108141).Within(Rel).Percent);
            var acn = Get(r, "ACN");
            Assert.That(acn.ReverseComplement, Is.EqualTo("NGT"));
            Assert.That(acn.SequenceIndices, Has.Count.EqualTo(2));
            Assert.That(acn.MatchingSequenceProbability, Is.EqualTo(0.72456658149428699).Within(Rel).Percent);
            var atn = Get(r, "ATN");
            Assert.That(atn.ExpectedFrequency, Is.EqualTo(0.14222222222222222).Within(Rel).Percent);
            Assert.That(atn.MatchingSequenceProbability, Is.EqualTo(0.35326100605480787).Within(Rel).Percent);
            Assert.That(atn.MatchingSequenceSignificance, Is.EqualTo(-0.92830694314348305).Within(Rel).Percent);
        });
    }

    [Test]
    public void MatchingSequences_AreTheUnionOverMatchingWords()
    {
        // Brute force: a sequence matches D when one of its windows matches D.
        var r = MotifFinder.FindSharedMotifs(D3, 3, 1, OligoBackgroundModel.Equiprobable, OligoStrandMode.Single, 0.0,
            OligoDegeneracy.OneDegenerate);
        foreach (var m in r.Motifs)
        {
            var expected = Enumerable.Range(0, D3.Length).Where(i =>
                Enumerable.Range(0, D3[i].Length - 2).Any(p =>
                    Enumerable.Range(0, 3).All(j => IupacHelper.MatchesIupac(D3[i].Sequence[p + j], m.Sequence[j])))).ToArray();
            Assert.That(m.SequenceIndices, Is.EqualTo(expected), m.Sequence);
        }
    }

    [Test]
    public void NoDegeneracy_EqualsPlainOverload()
    {
        var a = MotifFinder.FindSharedMotifs(D3, 3, 2, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Both, 0.1);
        var b = MotifFinder.FindSharedMotifs(D3, 3, 2, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Both, 0.1, OligoDegeneracy.None);
        Assert.That(b.Motifs.Select(m => (m.Sequence, m.MatchingSequenceProbability)),
            Is.EqualTo(a.Motifs.Select(m => (m.Sequence, m.MatchingSequenceProbability))));
    }

    [Test]
    public void PseudoFrequency_UsesDegenerateNpo()
    {
        var plain = MotifFinder.FindSharedMotifs(D3, 3, 1, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Single, 0.0, OligoDegeneracy.OneN);
        var pseudo = MotifFinder.FindSharedMotifs(D3, 3, 1, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Single, 0.2, OligoDegeneracy.OneN);
        Assert.That(Get(pseudo, "ACN").ExpectedFrequency,
            Is.EqualTo(0.8 * Get(plain, "ACN").ExpectedFrequency + 0.2 / 48).Within(Rel).Percent);
    }

    [Test]
    public void InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => MotifFinder.FindSharedMotifs(D3, 3, 1, OligoBackgroundModel.Equiprobable, OligoStrandMode.Single, 0.0, (OligoDegeneracy)9),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.FindSharedMotifs(D3, 3, 1, OligoBackgroundModel.Equiprobable, OligoStrandMode.Single, 1.5, OligoDegeneracy.OneN),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.FindSharedMotifs(D3, 3, 0, OligoBackgroundModel.Equiprobable, OligoStrandMode.Single, 0.0, OligoDegeneracy.OneN),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.FindSharedMotifs(new DnaSequence[] { null! }, 3, 1, OligoBackgroundModel.Equiprobable, OligoStrandMode.Single, 0.0, OligoDegeneracy.OneN),
                NUnit.Framework.Throws.ArgumentException);
        });
    }
}
