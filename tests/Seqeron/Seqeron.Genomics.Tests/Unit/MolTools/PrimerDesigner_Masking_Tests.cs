// PRIMER-DESIGN-001 — Primer3 template masking (audit round 3, A3-5 part 2b): PRIMER_MASK_TEMPLATE,
// PRIMER_MASK_KMERLIST_PATH / _PREFIX (caller-supplied 11-mer / 16-mer genome counts), PRIMER_MASK_FAILURE_RATE,
// PRIMER_MASK_5P_DIRECTION / _3P_DIRECTION, PRIMER_WT_MASK_FAILURE_RATE.
// Source: primer3-py 2.3.1 masker.c / masker.h (Kõressaar et al. 2018, Bioinformatics 34:1937): create_default_formula_parameters
//         (lists <prefix>_11.list / <prefix>_16.list, coefficients 0.1772 / 0.239, intercept −4.336),
//         get_frequency_of_canonical_oligo (k-mer, else reverse complement, else 1), calculate_scores (logistic, 0 when the
//         sum is 0), mask_oligo_region / read_and_mask_sequence (both strands separately, soft masking, 5000-char ring
//         buffer); libprimer3.c calc_and_check_oligo_features (is_lowercase_masked at the 3' end on trimmed_masked_seq /
//         _r; failure_rate from the last 16 nt), p_obj_fn (weights.failure_rate after the size terms); thermoanalysis.pyx
//         (window_size 16; "masking template chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified").
// Expected values: the masked copies and failure rates from masker.c compiled unchanged with a 40-line driver that sets
//         libprimer3's masker parameters (scratch h9/mk.c); the designs from primer3-py 2.3.1 design_primers with the same
//         lists written as GenomeTester4 files.
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_Masking_Tests
{
    // random.Random(5), 300 nt (same template as PrimerDesigner_BoundAndPosition_Tests).
    private const string T = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG";
    private const int TargetStart = 140, TargetEnd = 160; // SEQUENCE_TARGET [140, 20]

    // masker.c output for PRIMER_MASK_FAILURE_RATE 0.1, 5P 1, 3P 0 (forward copy, reverse copy).
    private const string MaskedFwdDefault = "GGATCACAGTCTACACTGCTCAcTCCAACcCCGGCCCCTGAGTcCGAGGAGAGGGTGCTTCAGAgTATGTATACCACTgGGTAGGATACGGCGGAGGGCaCGTCAAtACGGTTcAATGCCCTACTGCaTGCTCTtGTGGTTCATCTGCaTGGAGAGGGTGGGCATGGGTgGGGGTGCTGGCCCgTGATCTGGACCTCCCATCCAcAGCTCAtTGTACCgAGTGTAGAGAGGGgCTTGTCcTTCCAGATAGCGTtTCTGTTTCGGTGTAGGTGCTaATCGACTATGCTAcTGCGGTTAACG";
    private const string MaskedRevDefault = "GGATCACAGTCTACaCTGCTCaCTCCAAcCCCGGCCCCTGAGtCCGAGGaGAGGGTGCTTCAGaGTATGTATACCACTGGGTAGgATACGGCGGAGGGcACGTCAATACGGTTCAATGCcCTACTGcATGCTCtTGTGGTTCATCTGcATGGAGaGGGTGGGCATGGGtGGGGGTGCTGGCCCGTGATCtGGACCTCCCATCCaCAGCTCATTGTACCGAGTGTaGAGAGGgGCTTGTcCTTCCAGATAGCGtTTCTGTtTCGGTGTAGGTGCtAATCGACTATGCTACTGCGGTTAACG";

    // PRIMER_MASK_FAILURE_RATE 0.05, PRIMER_MASK_5P_DIRECTION 2, PRIMER_MASK_3P_DIRECTION 1.
    private const string MaskedFwd2 = "GGATCACAGTCTACACTGCTCactCCAAcccCGGCCCCTGAGtccGAGGAGAGGGTgctTCAGagtATGTATACCACtggGTAGgatACGGcggAGGGcacGTCAataCGGTtcaATGCCCTACTGcatGCTCttgTGGTTCATCTGcatGGAGAGGGTGGgcaTGGGtggGGGTGCTGGCCcgtGATCtggACCTcccATCCacaGCTCattGTACcgaGTGTAGAGAGGggcTTGTcctTCCAGATAGCGtttCTGTTTCGGTGtagGTGCtaaTCGACTATGCTactGCGGttaACG";
    private const string MaskedRev2 = "GGATCAcagTCTAcacTGCTcacTCCAaccCCGGCCCCTGAgtcCGAGgagAGGGTGCTTCAgagTATGTATACCActgGGTAggaTACGGCGGAGGgcaCGTCaatACGGttcAATGcccTACTgcaTGCTcttGTGGTTCATCTgcaTGGAgagGGTGGGCATGGgtgGGGGTGCTGGCccgTGATctgGACCTCCCATCcacAGCTcatTGTAccgAGTGtagAGAGgggCTTGtccTTCCAGATAGCgttTCTGtttCGGTGTAGGTGctaATCGACTATGCTACTGCGGTTAACG";

    // The lists: every third 11-mer of T (reverse-complemented at odd positions) with count 2 + (37·i) mod 1000, every
    // seventh 16-mer with count 10^(1 + i mod 5) (97 / 41 entries; later positions overwrite equal keys).
    private static readonly Dictionary<string, int> K11 = BuildK11();
    private static readonly Dictionary<string, int> K16 = BuildK16();
    private static readonly PrimerMaskingKmerLists Lists = new(K11, K16);

    private static Dictionary<string, int> BuildK11()
    {
        var d = new Dictionary<string, int>();
        for (int i = 0; i <= T.Length - 11; i++)
            if (i % 3 == 0)
            {
                string k = T.Substring(i, 11);
                d[i % 2 == 1 ? DnaSequence.GetReverseComplementString(k) : k] = 2 + i * 37 % 1000;
            }
        return d;
    }

    private static Dictionary<string, int> BuildK16()
    {
        var d = new Dictionary<string, int>();
        for (int i = 0; i <= T.Length - 16; i++)
            if (i % 7 == 0)
                d[T.Substring(i, 16)] = (int)Math.Pow(10, 1 + i % 5);
        return d;
    }

    private static PrimerPairOptions Masked(double failureRate = 0.1, int m5 = 1, int m3 = 0) =>
        PrimerPairOptions.Primer3Defaults with
        {
            MaskTemplate = true,
            MaskKmerLists = Lists,
            MaskFailureRate = failureRate,
            MaskFivePrimeDirection = m5,
            MaskThreePrimeDirection = m3,
        };

    private static PrimerParameters WithMaskWeight(double w) =>
        PrimerDesigner.Primer3DefaultParameters with
        {
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { MaskFailureRate = w },
        };

    // (left start, left len, right 3'-end coordinate, right len, PAIR_PENALTY, LEFT_PENALTY, RIGHT_PENALTY).
    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int LLen, int R, int RLen, double Pen, double LPen, double RPen)[] expected)
    {
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                var p = pairs[k];
                Assert.That(p.Forward!.Position, Is.EqualTo(e.L), $"rank {k} left");
                Assert.That(p.Forward.Length, Is.EqualTo(e.LLen), $"rank {k} left length");
                Assert.That(p.Reverse!.Position + p.Reverse.Length - 1, Is.EqualTo(e.R), $"rank {k} right");
                Assert.That(p.Reverse.Length, Is.EqualTo(e.RLen), $"rank {k} right length");
                Assert.That(p.PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(p.Forward.Penalty, Is.EqualTo(e.LPen).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(p.Reverse.Penalty, Is.EqualTo(e.RPen).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
            }
        });
    }

    #region masker.c

    [Test]
    public void MaskTemplatePrimer3_Defaults_MatchesMaskerC()
    {
        var (fwd, rev) = PrimerDesigner.MaskTemplatePrimer3(T, Lists);
        Assert.Multiple(() =>
        {
            Assert.That(fwd, Is.EqualTo(MaskedFwdDefault), "forward copy (trimmed_masked_seq)");
            Assert.That(rev, Is.EqualTo(MaskedRevDefault), "reverse copy (trimmed_masked_seq_r)");
        });
    }

    [Test]
    public void MaskTemplatePrimer3_DirectionsAndRate_MatchesMaskerC()
    {
        var (fwd, rev) = PrimerDesigner.MaskTemplatePrimer3(T, Lists, 0.05, 2, 1);
        Assert.Multiple(() =>
        {
            Assert.That(fwd, Is.EqualTo(MaskedFwd2), "forward copy");
            Assert.That(rev, Is.EqualTo(MaskedRev2), "reverse copy");
            // PRIMER_MASK_FAILURE_RATE 0: nothing is masked.
            Assert.That(PrimerDesigner.MaskTemplatePrimer3(T, Lists, 0.0).Forward, Is.EqualTo(T));
        });
    }

    [TestCase("TGCTTCAGAGTATGTATACC", 0.0)]                 // neither 3'-terminal k-mer listed: s = 0 → 0
    [TestCase("GGATCACAGTCTACACTGCT", 0.035374931720555836)] // T[0..20)
    [TestCase("AGTCTACACTGCTCACTCCAACC", 0.17017170998451242)] // T[7..30)
    [TestCase("ACGT", 0.0)]                                    // shorter than the 16-nt window
    public void CalculateMaskFailureRatePrimer3_MatchesMaskerC(string primer, double expected)
    {
        Assert.That(PrimerDesigner.CalculateMaskFailureRatePrimer3(primer, Lists), Is.EqualTo(expected).Within(1e-15));
    }

    [Test]
    public void CalculateMaskFailureRatePrimer3_ReverseComplementLookup()
    {
        // rc(T[100..122)): its 3'-terminal 16-mer / 11-mer are reverse complements of unlisted k-mers → 0 (masker.c output).
        string p = DnaSequence.GetReverseComplementString(T.Substring(100, 22));
        Assert.That(PrimerDesigner.CalculateMaskFailureRatePrimer3(p, Lists), Is.EqualTo(0.0));
        // A primer ending with a listed 16-mer, given as its reverse complement in the list, scores like the k-mer itself.
        string kmer = T.Substring(14, 16);
        var rcLists = new PrimerMaskingKmerLists(K11, new Dictionary<string, int> { [DnaSequence.GetReverseComplementString(kmer)] = 1000 });
        var fwdLists = new PrimerMaskingKmerLists(K11, new Dictionary<string, int> { [kmer] = 1000 });
        Assert.That(PrimerDesigner.CalculateMaskFailureRatePrimer3("AAAA" + kmer, rcLists),
            Is.EqualTo(PrimerDesigner.CalculateMaskFailureRatePrimer3("AAAA" + kmer, fwdLists)));
    }

    [Test]
    public void FromGenomeTester4Files_ReadsGlistmakerFormat()
    {
        string dir = Path.Combine(Path.GetTempPath(), "seqeron-mask-" + Guid.NewGuid().ToString("N")) + Path.DirectorySeparatorChar;
        Directory.CreateDirectory(dir);
        try
        {
            WriteList(dir + "test_11.list", 11, K11);
            WriteList(dir + "test_16.list", 16, K16);
            var lists = PrimerMaskingKmerLists.FromGenomeTester4Files(dir, "test");
            Assert.Multiple(() =>
            {
                Assert.That(lists.Count11, Is.EqualTo(K11.Count));
                Assert.That(lists.Count16, Is.EqualTo(K16.Count));
                Assert.That(PrimerDesigner.MaskTemplatePrimer3(T, lists).Forward, Is.EqualTo(MaskedFwdDefault));
                Assert.That(() => PrimerMaskingKmerLists.FromGenomeTester4Files(dir, "missing"),
                    NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Cannot find list file"));
            });
            File.WriteAllBytes(dir + "bad_11.list", new byte[48]);
            File.Copy(dir + "test_16.list", dir + "bad_16.list");
            Assert.That(() => PrimerMaskingKmerLists.FromGenomeTester4Files(dir, "bad"),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Given file is not a list file"));
            WriteList(dir + "empty_11.list", 11, new Dictionary<string, int>());
            File.Copy(dir + "test_16.list", dir + "empty_16.list");
            Assert.That(() => PrimerMaskingKmerLists.FromGenomeTester4Files(dir, "empty"),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("List file contains no kmers"));
        }
        finally
        {
            Directory.Delete(dir, true);
        }

        // GenomeTester4 list: magic "GT4C", k at byte 12, number of k-mers at 16, header size at 32, sorted (word, count) records.
        static void WriteList(string path, int k, Dictionary<string, int> kmers)
        {
            var records = kmers.Select(kv => (Word: kv.Key.Aggregate(0UL, (w, c) => (w << 2) | (uint)"ACGT".IndexOf(c)), Count: (uint)kv.Value))
                .OrderBy(r => r.Word).ToList();
            using var bw = new BinaryWriter(File.Create(path));
            bw.Write((uint)(('G' << 24) | ('T' << 16) | ('4' << 8) | 'C'));
            bw.Write(0u); bw.Write(0u);
            bw.Write((uint)k);
            bw.Write((uint)records.Count);
            bw.Write(0u); bw.Write(0u); bw.Write(0u);
            bw.Write(48UL);
            bw.Write(0UL);
            foreach (var (word, count) in records)
            {
                bw.Write(word);
                bw.Write(count);
            }
        }
    }

    #endregion

    #region Pair designs

    [Test]
    public void DesignPrimerPairs_MaskTemplateWithFailureRateWeight_MatchPrimer3()
    {
        // PRIMER_MASK_TEMPLATE 1, PRIMER_WT_MASK_FAILURE_RATE 1 (PRIMER_LEFT_EXPLAIN: lowercase masking of 3' end 96).
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), TargetStart, TargetEnd, WithMaskWeight(1.0), Masked());
        AssertPairs(pairs,
        [
            (46, 20, 233, 20, 0.2229611791929642, 0.1161425661219937, 0.10681861307097051),
            (116, 20, 233, 20, 0.5002207880304468, 0.39340217495947627, 0.10681861307097051),
            (46, 20, 240, 20, 0.5299746133399026, 0.1161425661219937, 0.4138320472179089),
            (117, 20, 233, 20, 0.5401089974063126, 0.43329038433534206, 0.10681861307097051),
            (46, 20, 236, 20, 0.6378549819560939, 0.1161425661219937, 0.5217124158341002),
        ]);
        Assert.Multiple(() =>
        {
            // Rank 3's left primer [117,20] carries its failure rate in the penalty (0.43329… − 0.39353… without masking).
            Assert.That(pairs[3].Forward!.MaskFailureRate!.Value, Is.EqualTo(0.43329038433534206 - 0.39353735265308387).Within(1e-12));
            foreach (var p in pairs)
            {
                Assert.That(char.IsUpper(MaskedFwdDefault[p.Forward!.Position + p.Forward.Length - 1]), "left 3' base unmasked");
                Assert.That(char.IsUpper(MaskedRevDefault[p.Reverse!.Position]), "right 3' base unmasked");
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_MaskDirectionsAndRate_MatchPrimer3()
    {
        // PRIMER_MASK_FAILURE_RATE 0.05, PRIMER_MASK_5P_DIRECTION 2, PRIMER_MASK_3P_DIRECTION 1, no failure-rate weight.
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), TargetStart, TargetEnd,
            PrimerDesigner.Primer3DefaultParameters, Masked(0.05, 2, 1));
        AssertPairs(pairs,
        [
            (117, 20, 233, 20, 0.5003559657240544, 0.39353735265308387, 0.10681861307097051),
            (117, 20, 240, 20, 0.8073693998709928, 0.39353735265308387, 0.4138320472179089),
            (117, 20, 239, 20, 0.9311901037085022, 0.39353735265308387, 0.5376527510554183),
            (117, 20, 247, 20, 0.9446679478723468, 0.39353735265308387, 0.5511305952192629),
            (100, 21, 233, 20, 1.170287183560788, 1.0634685704898175, 0.10681861307097051),
        ]);
        var single = PrimerDesigner.DesignPrimers(new DnaSequence(T), TargetStart, TargetEnd, PrimerDesigner.Primer3DefaultParameters, Masked(0.05, 2, 1));
        Assert.That(single.PairPenalty!.Value, Is.EqualTo(0.5003559657240544).Within(1e-9), "DesignPrimers rank 0");
    }

    [Test]
    public void DesignPrimerPairs_WithoutMaskTemplate_WeightIsInert()
    {
        // failure_rate is 0 unless PRIMER_MASK_TEMPLATE: the default ranking, no failure rate reported.
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), TargetStart, TargetEnd, WithMaskWeight(5.0),
            PrimerPairOptions.Primer3Defaults with { MaskKmerLists = Lists });
        Assert.Multiple(() =>
        {
            Assert.That(pairs[0].PairPenalty!.Value, Is.EqualTo(0.2229611791929642).Within(1e-9));
            Assert.That(pairs[0].Forward!.MaskFailureRate, Is.Null);
            Assert.That(PrimerDesigner.EvaluatePrimer("TGCTTCAGAGTATGTATACC", 56, true, WithMaskWeight(5.0)).MaskFailureRate, Is.Null);
        });
    }

    #endregion

    #region Penalty term and validation

    [Test]
    public void CalculatePrimer3Penalty_MaskFailureRateTerm()
    {
        var w = new Primer3PenaltyWeights(0, 0, 0, 0, 0, 0, 0, 0, 0) { MaskFailureRate = 2.0, PositionPenalty = 0 };
        Assert.That(PrimerDesigner.CalculatePrimer3Penalty(new Primer3PenaltyInputs(60, 20, 50) { MaskFailureRate = 0.25 }, w),
            Is.EqualTo(0.5).Within(1e-15));
    }

    [Test]
    public void Validation_MatchesPrimer3()
    {
        var dna = new DnaSequence(T);
        Assert.Multiple(() =>
        {
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, PrimerDesigner.Primer3DefaultParameters,
                    PrimerPairOptions.Primer3Defaults with { MaskTemplate = true }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("masking template chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified"));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.MaskTemplatePrimer3(T, Lists, 0.1, -1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.MaskTemplatePrimer3(T, Lists, 0.1, 1, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd,
                PrimerDesigner.Primer3DefaultParameters, Masked(double.NaN)));
            Assert.Throws<ArgumentException>(() => new PrimerMaskingKmerLists(new Dictionary<string, int>(), K16));
            Assert.Throws<ArgumentException>(() => new PrimerMaskingKmerLists(new Dictionary<string, int> { ["ACGT"] = 3 }, K16));
            Assert.Throws<ArgumentException>(() => new PrimerMaskingKmerLists(new Dictionary<string, int> { ["ACGTNACGTAC"] = 3 }, K16));
            Assert.Throws<ArgumentException>(() => new PrimerMaskingKmerLists(new Dictionary<string, int> { ["ACGTAACGTAC"] = -1 }, K16));
        });
    }

    #endregion
}
