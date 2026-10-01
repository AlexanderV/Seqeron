// MOTIF-DISCOVER-001 (audit round 2, group G2) — RSAT oligo-analysis -seqtype dna|prot|other.
// Evidence: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md; TestSpec: tests/TestSpecs/MOTIF-DISCOVER-001.md
// Source: rsa-tools/rsat-code master 10043f2 — perl-scripts/oligo-analysis v1.169 (ReadArguments -seqtype: prot/other set
//         $sum_rc = 0; alphabet / %accepted_residue / %accepted_oligo block; CountOligos: residue counts of sequences of
//         length >= k, deletion of patterns with residues outside the alphabet and $nb_possible_pos -= discarded occurrences;
//         sub alphabet: protein 20 letters, other = sort keys %residue_occ; CalcAlphabet; CalcExpected equiprobable
//         1/alphabet_size**k and -pseudo), lib/RSA.disco.lib NbPossibleOligos (reverse-complement / degenerate corrections
//         for DNA only), lib/RSA.seq.lib FoldSequence (s/\s+//g) and OverlapCoeff.
// Reference values: RSAT oligo-analysis itself run from that clone (perl 5.38) with a %.17g dump of every pattern field
//         (oracle copy oligo-analysis-g2). Three RSAT defects are bypassed in the oracle (env-guarded, stated per test):
//         (1) OverlapCoeff's fallback for an empty %residue_proba (Markov models) is "equiprobable nucleotides"
//             a,c,g,t = 1/4 and 0 for every other residue — the oracle uses 1/alphabet_size for -seqtype prot/other;
//         (2) OverlapCoeff interpolates the word unquoted into a regex (. ? * + act as operators for -seqtype other) —
//             the oracle quotes it (\Q..\E);
//         (3) -markov dies with "Illegal division by zero" when a sub-word before a discarded residue was never counted —
//             the oracle keeps exp_freq = 0 (already 0 from the uncounted (m+1)-mer).
//         Comparison over 450 random runs (150 protein, 150 other, 150 DNA with N/R/Y/-; equi / input / Markov / lexicon,
//         -ovlp/-noov, -pseudo, DNA -2str and -oneN/-onedeg): 15,066 patterns, every column <= 1.7e-13.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_OligoSequenceType_Tests
{
    private const double Rel = 1e-10; // percent units (1e-12 relative)

    private static OligoStatistics Get(OligoAnalysisReport r, string pattern) => r.Patterns.Single(p => p.Pattern == pattern);

    private static OligoAnalysisOptions Protein(OligoBackgroundModel? background = null) => new()
    {
        SequenceType = OligoSequenceType.Protein,
        Background = background ?? OligoBackgroundModel.BernoulliFromInput,
    };

    [Test]
    public void Protein_InputBernoulli_MatchesRsat()
    {
        // RSAT: oligo-analysis -seqtype prot -l 2 -return occ,proba,zscore,freq,ratio on
        // >p1 MKLLVAAGLLKLMKXLLA* / >p2 mkllvqqKLLAA: 29 windows - 3 discarded (KX, XL, A*) = 26; NPO 20^2 = 400; 14 tested.
        var r = MotifFinder.AnalyzeOligoStrings(new[] { "MKLLVAAGLLKLMKXLLA*", "mkllvqqKLLAA" }, 2, Protein());
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceType, Is.EqualTo(OligoSequenceType.Protein));
            Assert.That(r.AlphabetSize, Is.EqualTo(20));
            Assert.That(r.PossibleOligos, Is.EqualTo(400));
            Assert.That(r.TotalOccurrences, Is.EqualTo(26));
            Assert.That(r.PossiblePositions, Is.EqualTo(26));
            Assert.That(r.TestedPatterns, Is.EqualTo(14));
            Assert.That(r.Strands, Is.EqualTo(OligoStrandMode.Single));
            Assert.That(r.Patterns.Any(p => p.Pattern.Contains('X') || p.Pattern.Contains('*')), Is.False);

            var ll = Get(r, "LL");
            Assert.That(ll.Occurrences, Is.EqualTo(5));
            Assert.That(ll.ExpectedFrequency, Is.EqualTo(0.14387633769322233).Within(Rel).Percent);
            Assert.That(ll.ExpectedOccurrences, Is.EqualTo(3.7407847800237808).Within(Rel).Percent);
            Assert.That(ll.OverlapCoefficient, Is.EqualTo(1.3793103448275863).Within(Rel).Percent);
            Assert.That(ll.ExpectedVariance, Is.EqualTo(3.8875694384551549).Within(Rel).Percent);
            Assert.That(ll.ZScore, Is.EqualTo(0.63864701205606378).Within(Rel).Percent);
            Assert.That(ll.Ratio, Is.EqualTo(1.3366179275270187).Within(Rel).Percent);
            Assert.That(ll.OccurrenceProbability, Is.EqualTo(0.31619335865366638).Within(Rel).Percent);
            Assert.That(ll.OccurrenceEValue, Is.EqualTo(4.4267070211513291).Within(Rel).Percent);
            Assert.That(ll.OccurrenceSignificance, Is.EqualTo(-0.64608077942004516).Within(Rel).Percent);

            var kl = Get(r, "KL");
            Assert.That(kl.Occurrences, Is.EqualTo(4));
            Assert.That(kl.ExpectedFrequency, Is.EqualTo(0.065398335315101072).Within(Rel).Percent);
            Assert.That(kl.OccurrenceProbability, Is.EqualTo(0.0865394727065604).Within(Rel).Percent);
            Assert.That(kl.ZScore, Is.EqualTo(2.1497117755830679).Within(Rel).Percent);
        });
    }

    [Test]
    public void Other_Equiprobable_ObservedAlphabet_MatchesRsat()
    {
        // RSAT: -seqtype other -bg equi -l 2 on "Hello, World! hello world. AbC abc": white space removed, lower-cased;
        // alphabet {! , . a b c d e h l o r w} = 13, NPO 169, 28 windows, 18 tested.
        var r = MotifFinder.AnalyzeOligoStrings(new[] { "Hello, World! hello world. AbC abc" }, 2,
            new OligoAnalysisOptions { SequenceType = OligoSequenceType.Other, Background = OligoBackgroundModel.Equiprobable });
        Assert.Multiple(() =>
        {
            Assert.That(r.AlphabetSize, Is.EqualTo(13));
            Assert.That(r.PossibleOligos, Is.EqualTo(169));
            Assert.That(r.TotalOccurrences, Is.EqualTo(28));
            Assert.That(r.TestedPatterns, Is.EqualTo(18));
            var ll = Get(r, "ll");
            Assert.That(ll.Occurrences, Is.EqualTo(2));
            Assert.That(ll.ExpectedFrequency, Is.EqualTo(0.0059171597633136093).Within(Rel).Percent);
            Assert.That(ll.OverlapCoefficient, Is.EqualTo(1.0769230769230769).Within(Rel).Percent);
            Assert.That(ll.ExpectedVariance, Is.EqualTo(0.18626798781555265).Within(Rel).Percent);
            Assert.That(ll.ZScore, Is.EqualTo(4.250165852579884).Within(Rel).Percent);
            Assert.That(ll.OccurrenceProbability, Is.EqualTo(0.01194994022319412).Within(Rel).Percent);
            Assert.That(ll.OccurrenceSignificance, Is.EqualTo(0.66736176206457676).Within(Rel).Percent);
            var ab = Get(r, "ab");
            Assert.That(ab.Occurrences, Is.EqualTo(2), "AbC and abc fold to the same word");
            Assert.That(ab.ZScore, Is.EqualTo(4.5746803547824042).Within(Rel).Percent);
            Assert.That(r.Patterns.Any(p => p.Pattern.Any(char.IsUpper)), Is.False);
        });
    }

    [Test]
    public void Protein_Markov_OverlapUsesEquiprobableResidues()
    {
        // RSAT: -seqtype prot -markov 1 -l 3 on LLLLAAAAMKLLLAA (13 windows, NPO 8000, 8 tested). Oracle fix (1): RSAT's
        // OverlapCoeff fallback gives AAA 1 + 1/4 + 1/16 = 1.3125 and LLL 1 (L is not a nucleotide); the equiprobable
        // residues of the protein alphabet give 1 + 1/20 + 1/400 = 1.0525 for both.
        var r = MotifFinder.AnalyzeOligoStrings(new[] { "LLLLAAAAMKLLLAA" }, 3, Protein(OligoBackgroundModel.MarkovFromInput(1)));
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(8000));
            Assert.That(r.TestedPatterns, Is.EqualTo(8));
            var aaa = Get(r, "AAA");
            Assert.That(aaa.ExpectedFrequency, Is.EqualTo(0.2040816326530612).Within(Rel).Percent);
            Assert.That(aaa.OverlapCoefficient, Is.EqualTo(1.0525).Within(Rel).Percent);
            Assert.That(aaa.ExpectedVariance, Is.EqualTo(-0.85845481049562633).Within(Rel).Percent);
            Assert.That(aaa.ZScore, Is.NaN, "exp_var <= 0: RSAT z-score NA");
            Assert.That(aaa.OccurrenceProbability, Is.EqualTo(0.77709852180454864).Within(Rel).Percent);
            var lll = Get(r, "LLL");
            Assert.That(lll.ExpectedFrequency, Is.EqualTo(0.27332361516034986).Within(Rel).Percent);
            Assert.That(lll.OverlapCoefficient, Is.EqualTo(1.0525).Within(Rel).Percent);
            Assert.That(lll.OccurrenceProbability, Is.EqualTo(0.73334471367664145).Within(Rel).Percent);
        });
    }

    [Test]
    public void Dna_UndefinedResidues_WindowsDiscarded_BothStrands_MatchesRsat()
    {
        // RSAT: -seqtype dna -2str -l 4 on ACGTNACGTACGRTTACGTAC / ttacgtnnACGTAC: windows with N or R are discarded,
        // nb_possible_pos = 16; NPO (4^4 + 4^2)/2 = 136; 4 tested.
        var r = MotifFinder.AnalyzeOligoStrings(new[] { "ACGTNACGTACGRTTACGTAC", "ttacgtnnACGTAC" }, 4,
            new OligoAnalysisOptions { Strands = OligoStrandMode.Both });
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceType, Is.EqualTo(OligoSequenceType.Dna));
            Assert.That(r.AlphabetSize, Is.EqualTo(4));
            Assert.That(r.TotalOccurrences, Is.EqualTo(16));
            Assert.That(r.PossibleOligos, Is.EqualTo(136));
            Assert.That(r.TestedPatterns, Is.EqualTo(4));
            var acgt = Get(r, "ACGT");
            Assert.That(acgt.Occurrences, Is.EqualTo(5));
            Assert.That(acgt.ExpectedFrequency, Is.EqualTo(0.0038334266356693561).Within(Rel).Percent);
            Assert.That(acgt.ZScore, Is.EqualTo(20.294579804259762).Within(Rel).Percent);
            Assert.That(acgt.OccurrenceProbability, Is.EqualTo(3.4909208986462239e-09).Within(Rel).Percent);
            Assert.That(acgt.OccurrenceSignificance, Is.EqualTo(7.8550000004968332).Within(Rel).Percent);
            var gtac = Get(r, "GTAC");
            Assert.That(gtac.Occurrences, Is.EqualTo(3));
            Assert.That(gtac.OccurrenceProbability, Is.EqualTo(3.0388746348848393e-05).Within(Rel).Percent);
        });
    }

    [Test]
    public void Dna_Markov_UncountedSubWordBeforeDiscardedResidue_IsNotTested()
    {
        // RSAT -markov 2 -l 5 dies on this input ("Illegal division by zero", oracle fix (3)): the 3-mer and 2-mer of a
        // valid word that precede the Y / N are never counted, so its exp_freq is 0 and it is reported untested.
        var r = MotifFinder.AnalyzeOligoStrings(new[] { "tGCtYAACTaCCAACCaATcTGNaCa" }, 5,
            new OligoAnalysisOptions { Background = OligoBackgroundModel.MarkovFromInput(2) });
        var untested = r.Patterns.Where(p => p.FittedDistribution == OligoFittedDistribution.None).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(untested, Is.Not.Empty);
            Assert.That(untested.All(p => p.ExpectedFrequency == 0 && double.IsNaN(p.OccurrenceProbability)), Is.True);
            Assert.That(r.TestedPatterns, Is.EqualTo(r.Patterns.Count - untested.Count));
        });
    }

    [Test]
    public void Dna_Strings_WithoutUndefinedResidues_AreBitIdenticalToDnaSequencePath()
    {
        string[] seqs = { "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGC", "acgtacgtTTAGGC", "GG" };
        var dna = seqs.Select(s => new DnaSequence(s.ToUpperInvariant())).ToList();
        var variants = new[]
        {
            new OligoAnalysisOptions(),
            new OligoAnalysisOptions { Strands = OligoStrandMode.Both, CountOverlapping = false, PseudoFrequency = 0.1 },
            new OligoAnalysisOptions { Background = OligoBackgroundModel.MarkovFromInput(1) },
            new OligoAnalysisOptions { Background = OligoBackgroundModel.Lexicon },
            new OligoAnalysisOptions { Degeneracy = OligoDegeneracy.OneDegenerate, Strands = OligoStrandMode.Both },
        };
        foreach (var options in variants)
        {
            var a = MotifFinder.AnalyzeOligos(dna, 3, options);
            var b = MotifFinder.AnalyzeOligoStrings(seqs.Select(s => s.Insert(1, " \n")), 3, options);
            Assert.That(b.Patterns, Has.Count.EqualTo(a.Patterns.Count));
            for (int i = 0; i < a.Patterns.Count; i++)
            {
                var (x, y) = (a.Patterns[i], b.Patterns[i]);
                Assert.That(y.Pattern, Is.EqualTo(x.Pattern));
                Assert.That(y.Positions, Is.EqualTo(x.Positions));
                Assert.That(BitConverter.DoubleToInt64Bits(y.ExpectedFrequency), Is.EqualTo(BitConverter.DoubleToInt64Bits(x.ExpectedFrequency)));
                Assert.That(BitConverter.DoubleToInt64Bits(y.ZScore), Is.EqualTo(BitConverter.DoubleToInt64Bits(x.ZScore)));
                Assert.That(BitConverter.DoubleToInt64Bits(y.OccurrenceProbability), Is.EqualTo(BitConverter.DoubleToInt64Bits(x.OccurrenceProbability)));
            }
            Assert.That(b.TestedPatterns, Is.EqualTo(a.TestedPatterns));
            Assert.That(b.PossibleOligos, Is.EqualTo(a.PossibleOligos));
        }
    }

    [Test]
    public void DnaSequenceOverload_WithProteinType_EqualsStringOverload()
    {
        var a = MotifFinder.AnalyzeOligos(new[] { new DnaSequence("ACGTTGCAACGT") }, 2, Protein());
        var b = MotifFinder.AnalyzeOligoStrings(new[] { "ACGTTGCAACGT" }, 2, Protein());
        Assert.Multiple(() =>
        {
            Assert.That(a.PossibleOligos, Is.EqualTo(400));
            Assert.That(a.Patterns.Select(p => (p.Pattern, p.OccurrenceProbability)), Is.EqualTo(b.Patterns.Select(p => (p.Pattern, p.OccurrenceProbability))));
        });
    }

    [Test]
    public void Protein_CaseAndWhitespaceFolded_DiscardedWindowsShrinkN()
    {
        var a = MotifFinder.AnalyzeOligoStrings(new[] { "mk ll\tv\r\nAAXAA" }, 2, Protein());
        var b = MotifFinder.AnalyzeOligoStrings(new[] { "MKLLVAAXAA" }, 2, Protein());
        Assert.Multiple(() =>
        {
            Assert.That(a.TotalOccurrences, Is.EqualTo(7), "9 windows minus AX and XA");
            Assert.That(a.Patterns.Select(p => (p.Pattern, p.Occurrences, p.OccurrenceProbability)),
                Is.EqualTo(b.Patterns.Select(p => (p.Pattern, p.Occurrences, p.OccurrenceProbability))));
            Assert.That(Get(a, "AA").Positions, Is.EqualTo(new[] { new OligoOccurrence(0, 5), new OligoOccurrence(0, 8) }));
        });
    }

    [Test]
    public void Protein_PseudoFrequency_UsesAlphabetNpo()
    {
        // exp_freq = (1 - psi)·p + psi / 20^k (RSAT CalcExpected with nb_possible_oligos = alphabet_size**k).
        var plain = MotifFinder.AnalyzeOligoStrings(new[] { "MKLLVAAGLLKLMK" }, 2, Protein(OligoBackgroundModel.Equiprobable));
        var pseudo = MotifFinder.AnalyzeOligoStrings(new[] { "MKLLVAAGLLKLMK" }, 2,
            Protein(OligoBackgroundModel.Equiprobable) with { PseudoFrequency = 0.5 });
        Assert.That(Get(pseudo, "LL").ExpectedFrequency,
            Is.EqualTo(0.5 * Get(plain, "LL").ExpectedFrequency + 0.5 / 400).Within(Rel).Percent);
        Assert.That(Get(plain, "LL").ExpectedFrequency, Is.EqualTo(1.0 / 400).Within(Rel).Percent);
    }

    [Test]
    public void ResidueAlphabet_UnsupportedOptions_Throw()
    {
        string[] seqs = { "MKLLV" };
        Assert.Multiple(() =>
        {
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(seqs, 2, Protein() with { Strands = OligoStrandMode.Both }),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(seqs, 2, Protein() with { Degeneracy = OligoDegeneracy.OneN }),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(seqs, 2,
                Protein() with { Calibration = OligoCalibration.FromEntries(new[] { KeyValuePair.Create("AC", (1.0, 1.0)) }) }),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(seqs, 2,
                new OligoAnalysisOptions { SequenceType = OligoSequenceType.Other, Background = OligoBackgroundModel.Bernoulli(new[] { 0.25, 0.25, 0.25, 0.25 }) }),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(seqs, 2, new OligoAnalysisOptions { SequenceType = (OligoSequenceType)7 }),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(new[] { "MK", null! }, 2, Protein()), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => MotifFinder.AnalyzeOligoStrings(null!, 2, Protein()), NUnit.Framework.Throws.ArgumentNullException);
        });
    }
}
