// PROBE-DESIGN-001 — DesignProbes ranking option (audit round 3, A3-11, B07.md F51): the default additive score is a
// library heuristic; ProbeRanking.Primer3Penalty ranks by Primer3's internal-oligo objective p_obj_fn.
// Reference: primer3-py 2.3.1 design_primers(PRIMER_TASK = pick_hyb_probe_only, PRIMER_PICK_INTERNAL_OLIGO = 1)
//            PRIMER_INTERNAL_n / _n_PENALTY / _n_TM; Primer3 libprimer3.cc p_obj_fn (OT_INTL), primer_rec_comp,
//            _pr_data_control.
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class ProbeDesigner_Primer3Ranking_Tests
{
    // Primer3 settings of the cross-check: PRIMER_INTERNAL_MIN/MAX_SIZE 18/27, _MIN/MAX_TM 57/63, _MIN/MAX_GC 20/80,
    // _MAX_POLY_X 5, default conditions (50 nM, 50 mM, no Mg2+/dNTP), PRIMER_NUM_RETURN 12. The same window here.
    private static ProbeDesigner.ProbeParameters Primer3Window(double optTm, int optSize) =>
        new(MinLength: 18, MaxLength: 27, MinTm: 57, MaxTm: 63, MinGc: 0.20, MaxGc: 0.80,
            MaxHomopolymer: 5, AvoidSecondaryStructure: true, MaxSelfComplementarity: 0.3)
        {
            Ranking = ProbeDesigner.ProbeRanking.Primer3Penalty,
            OptTm = optTm,
            OptLength = optSize,
        };

    // random.seed(20261008); 70-nt ACGT templates with ≥ 8 Primer3 probes; (start, length, PRIMER_INTERNAL_n_PENALTY,
    // PRIMER_INTERNAL_n_TM) of primer3-py 2.3.1 pick_hyb_probe_only, best first (primer_rec_comp order).
    private static IEnumerable<TestCaseData> Primer3Cases()
    {
        yield return new TestCaseData(
            "GATCCGACGCTATATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAATAAATAATCCGTAGGGAG", 60.0, 20,
            new (int, int, double, double)[]
            {
                (2, 23, 5.437062235706321, 57.56293776429368), (2, 24, 6.098516903481425, 57.901483096518575),
                (1, 24, 6.281156733776584, 57.718843266223416), (0, 25, 6.603139914120391, 58.39686008587961),
                (3, 24, 6.762497692186912, 57.23750230781309), (2, 25, 6.786483773719397, 58.2135162262806),
                (0, 24, 6.79487996687817, 57.20512003312183), (1, 25, 6.96100380817154, 58.03899619182846),
                (0, 26, 7.323833899904571, 58.67616610009543), (3, 25, 7.4283400359979055, 57.571659964002095),
                (2, 26, 7.497965856587655, 58.502034143412345), (1, 26, 7.665048443870148, 58.33495155612985),
            }).SetName("Primer3Ranking_Template1_Opt60_20");
        yield return new TestCaseData(
            "ACCTGGCACAGATGCGAGGAGAAAGAATATTCGCGCCTTCCGACCCTGCGCAGGCAATACCGTTACCCCT", 61.5, 22,
            new (int, int, double, double)[]
            {
                (27, 22, 0.007361159225695246, 61.492638840774305), (38, 22, 0.0610278155618289, 61.43897218443817),
                (42, 22, 0.17268743365735872, 61.32731256634264), (41, 22, 0.19073118454213045, 61.69073118454213),
                (36, 22, 0.44529850216031264, 61.94529850216031), (37, 22, 0.6833419718861933, 60.81665802811381),
                (39, 22, 0.8036672902066471, 62.30366729020665), (26, 23, 1.0005314976000932, 61.50053149760009),
                (47, 22, 1.0724717885976247, 60.427528211402375), (40, 22, 1.0834851014256515, 62.58348510142565),
                (47, 23, 1.107861025097236, 61.392138974902764), (36, 23, 1.112473848469108, 61.38752615153089),
            }).SetName("Primer3Ranking_Template2_Opt61.5_22");
        yield return new TestCaseData(
            "GTGATTATTCTTCGGCCCCTGGCAGTACGCACTTGGCTCGTAACTTTTGAGGTGGGGTTTAGTTGTAACG", 59.0, 19,
            new (int, int, double, double)[]
            {
                (14, 19, 0.007730506722964492, 59.007730506722964), (14, 20, 1.3425637296118111, 59.34256372961181),
                (20, 20, 1.641291179152688, 58.35870882084731), (17, 20, 1.735156030377425, 58.264843969622575),
                (14, 18, 2.245732868651146, 57.754267131348854), (21, 20, 2.338367814458195, 57.661632185541805),
                (16, 20, 2.374406691313709, 57.62559330868629), (15, 20, 2.374406691313709, 57.62559330868629),
                (17, 21, 2.3826176519736464, 59.382617651973646), (15, 21, 2.715498229491402, 59.7154982294914),
                (19, 20, 2.7202132522774036, 57.279786747722596), (20, 21, 2.721818999737593, 59.72181899973759),
            }).SetName("Primer3Ranking_Template3_Opt59_19_TieByStartDescending");
        yield return new TestCaseData(
            "GGTCACCTCCACGGGAGATTTATACCGTGAAAGGCCTAGGACACATCTAACGTCGATGAAATGGTTGGAA", 60.0, 24,
            new (int, int, double, double)[]
            {
                (23, 24, 1.5366325208229341, 58.463367479177066), (23, 25, 1.5734342926938893, 59.42656570730611),
                (24, 24, 1.777558718832836, 58.222441281167164), (32, 24, 2.1613212931652583, 57.83867870683474),
                (31, 24, 2.166261617113207, 57.83373838288679), (31, 25, 2.1916137494297914, 58.80838625057021),
                (22, 24, 2.6378890608341408, 57.36211093916586), (21, 24, 2.6378890608341408, 57.36211093916586),
                (30, 25, 2.850354959211188, 58.14964504078881), (30, 26, 2.9214924146024828, 59.07850758539752),
                (22, 25, 2.9265142017897574, 58.07348579821024), (23, 26, 2.9892959399252845, 59.010704060074715),
            }).SetName("Primer3Ranking_Template4_Opt60_24_TieByStartDescending");
    }

    [TestCaseSource(nameof(Primer3Cases))]
    public void DesignProbes_Primer3Ranking_PenaltiesAndOrderMatchPrimer3PickHybProbeOnly(
        string template, double optTm, int optSize, (int Start, int Length, double Penalty, double Tm)[] primer3)
    {
        var ours = ProbeDesigner.DesignProbes(template, Primer3Window(optTm, optSize), maxProbes: int.MaxValue).ToList();
        var byKey = ours.ToDictionary(p => (p.Start, p.Sequence.Length));

        // Every Primer3 probe of this window is also a DesignProbes candidate (score > 0).
        var common = primer3.Where(e => byKey.ContainsKey((e.Start, e.Length))).ToList();
        Assert.That(common, Has.Count.EqualTo(primer3.Length));

        Assert.Multiple(() =>
        {
            foreach (var e in common)
            {
                var p = byKey[(e.Start, e.Length)];
                Assert.That(p.Primer3Penalty, Is.EqualTo(e.Penalty).Within(1e-9), $"penalty {e.Start}/{e.Length}");
                Assert.That(p.Tm, Is.EqualTo(e.Tm).Within(1e-9), $"Tm {e.Start}/{e.Length}");
            }

            // Order: our ranking restricted to Primer3's probes = Primer3's order (primer_rec_comp, ties by start descending).
            var keys = common.Select(e => (e.Start, e.Length)).ToHashSet();
            var ourOrder = ours.Select(p => (p.Start, p.Sequence.Length)).Where(keys.Contains).ToList();
            Assert.That(ourOrder, Is.EqualTo(common.Select(e => (e.Start, e.Length)).ToList()));

            // The ranking is penalty ascending over every candidate, ties start descending then length ascending.
            for (int i = 1; i < ours.Count; i++)
            {
                var (a, b) = (ours[i - 1], ours[i]);
                Assert.That(a.Primer3Penalty!.Value, Is.LessThanOrEqualTo(b.Primer3Penalty!.Value), $"rank {i}");
                if (a.Primer3Penalty == b.Primer3Penalty)
                    Assert.That(a.Start > b.Start || (a.Start == b.Start && a.Sequence.Length < b.Sequence.Length),
                        Is.True, $"tie order at rank {i}");
            }
        });
    }

    [Test]
    public void DesignProbes_Primer3Ranking_PenaltyEqualsDesignProbesPrimer3AndCanonicalPenalty()
    {
        // Same p_obj_fn as DesignProbesPrimer3 (PRIMER_INTERNAL_n_PENALTY) and as CalculatePrimer3Penalty with
        // Primer3's default internal-oligo weights around (OptTm, OptLength).
        const string template = "ACCTGGCACAGATGCGAGGAGAAAGAATATTCGCGCCTTCCGACCCTGCGCAGGCAATACCGTTACCCCT";
        var settings = new ProbeDesigner.Primer3ProbeSettings(MinSize: 18, OptSize: 22, MaxSize: 27, OptTm: 61.5);
        var picked = ProbeDesigner.DesignProbesPrimer3(template, settings, numReturn: 12);
        var ours = ProbeDesigner.DesignProbes(template, Primer3Window(61.5, 22), int.MaxValue)
            .ToDictionary(p => (p.Start, p.Sequence.Length));

        Assert.Multiple(() =>
        {
            foreach (var q in picked)
            {
                var p = ours[(q.Start, q.Length)];
                Assert.That(p.Primer3Penalty, Is.EqualTo(q.Penalty));
                Assert.That(p.Primer3Penalty, Is.EqualTo(PrimerDesigner.CalculatePrimer3Penalty(
                    new Primer3PenaltyInputs(p.Tm, q.Length, q.GcPercent), PrimerDesigner.DefaultPrimer3Weights,
                    new Primer3Optima(61.5, 22, null))));
                Assert.That(p.Primer3Penalty, Is.EqualTo(Math.Abs(p.Tm - 61.5) + Math.Abs(q.Length - 22)).Within(1e-12));
            }
        });
    }

    [Test]
    public void DesignProbes_Primer3Ranking_ReordersTheSameCandidatesAsTheAdditiveScore()
    {
        const string template = "GTGATTATTCTTCGGCCCCTGGCAGTACGCACTTGGCTCGTAACTTTTGAGGTGGGGTTTAGTTGTAACG";
        var primer3 = Primer3Window(59, 19);
        var additive = primer3 with { Ranking = ProbeDesigner.ProbeRanking.AdditiveScore };

        var byPenalty = ProbeDesigner.DesignProbes(template, primer3, int.MaxValue).ToList();
        var byScore = ProbeDesigner.DesignProbes(template, additive, int.MaxValue).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(byPenalty.Select(p => (p.Start, p.Sequence.Length)),
                Is.EquivalentTo(byScore.Select(p => (p.Start, p.Sequence.Length))));
            Assert.That(byPenalty.Select(p => p with { Primer3Penalty = null }), Is.EquivalentTo(byScore).Using(new ProbeComparer()));
            Assert.That(byScore.All(p => p.Primer3Penalty is null), Is.True, "additive ranking reports no Primer3 penalty");
            Assert.That(byPenalty.All(p => p.Primer3Penalty is not null), Is.True);
            // Top-k = prefix of the full ranking (lazy screening does not change the order).
            Assert.That(ProbeDesigner.DesignProbes(template, primer3, 5).Select(p => (p.Start, p.Sequence.Length)),
                Is.EqualTo(byPenalty.Take(5).Select(p => (p.Start, p.Sequence.Length))));
        });
    }

    [Test]
    public void DesignProbes_DefaultRanking_IsTheAdditiveScore()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.Defaults.Microarray.Ranking, Is.EqualTo(ProbeDesigner.ProbeRanking.AdditiveScore));
            Assert.That(new ProbeDesigner.ProbeParameters(18, 27, 57, 63, 0.2, 0.8, 5, true, 0.3).OptTm, Is.EqualTo(60.0));
            Assert.That(new ProbeDesigner.ProbeParameters(18, 27, 57, 63, 0.2, 0.8, 5, true, 0.3).OptLength, Is.EqualTo(20));
        });
        const string template = "GATCCGACGCTATATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAATAAATAATCCGTAGGGAG";
        var probes = ProbeDesigner.DesignProbes(template, Primer3Window(60, 20) with { Ranking = default }, 10).ToList();
        for (int i = 1; i < probes.Count; i++)
            Assert.That(probes[i - 1].Score, Is.GreaterThanOrEqualTo(probes[i].Score));
    }

    [Test]
    public void DesignProbes_Primer3Ranking_NonComputableTmRanksLastWithoutPenalty()
    {
        // A window with N has no Tm (Primer3 rejects it: PRIMER_INTERNAL_MAX_NS_ACCEPTED = 0); it stays a candidate
        // of the additive filter and ranks after every probe with a penalty.
        const string template = "GTGATTATTCTTCGGCCCCTGGCAGTACGCACTTGGNTCGTAACTTTTGAGG";
        var probes = ProbeDesigner.DesignProbes(template, Primer3Window(59, 19), int.MaxValue).ToList();
        int firstNull = probes.FindIndex(p => p.Primer3Penalty is null);
        Assert.Multiple(() =>
        {
            Assert.That(firstNull, Is.GreaterThan(0));
            Assert.That(probes.Skip(firstNull).All(p => p.Primer3Penalty is null && p.Sequence.Contains('N')), Is.True);
            Assert.That(probes.Take(firstNull).All(p => !p.Sequence.Contains('N')), Is.True);
        });
    }

    [Test]
    public void DesignProbes_Primer3Ranking_OptimaOutsideTheWindow_ThrowLikePrimer3()
    {
        // Primer3 _pr_data_control: "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE" / "< MIN_SIZE",
        // "Optimum internal oligo Tm lower than minimum or higher than maximum". Checked eagerly (before enumeration).
        var p = Primer3Window(60, 20);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes("ACGT", p with { OptLength = 28 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes("ACGT", p with { OptLength = 17 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes("ACGT", p with { OptTm = 63.01 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes("ACGT", p with { OptTm = double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes("ACGT", p with { Ranking = (ProbeDesigner.ProbeRanking)7 }));
            // The Microarray preset (82–90 °C, 50–60 nt) needs its own optima: Primer3's defaults 60 °C / 20 nt lie outside.
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes("ACGT",
                ProbeDesigner.Defaults.Microarray with { Ranking = ProbeDesigner.ProbeRanking.Primer3Penalty }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignTilingProbes(new string('A', 100), 20, 0,
                p with { OptTm = 50 }));
            // Optima at the window edges are accepted; the additive ranking ignores them.
            Assert.DoesNotThrow(() => ProbeDesigner.DesignProbes("ACGT", p with { OptLength = 27, OptTm = 57 }));
            Assert.DoesNotThrow(() => ProbeDesigner.DesignProbes("ACGT",
                p with { Ranking = ProbeDesigner.ProbeRanking.AdditiveScore, OptLength = 99, OptTm = -5 }));
        });
    }

    [Test]
    public void DesignProbes_Primer3Ranking_MicroarrayWithOptima_RanksByPenalty()
    {
        // Microarray preset with optima inside its window (82–90 °C, 50–60 nt): penalty = |Tm − 86| + |length − 55|
        // with the preset's OligoArray Tm (1 M Na+, 1 µM, NN over the whole oligo).
        const string template =
            "AGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATGCAGTCAGGCTAGCATCGGATCCAGTTGACCAAGTCGTTGACGC";
        var param = ProbeDesigner.Defaults.Microarray with
        {
            Ranking = ProbeDesigner.ProbeRanking.Primer3Penalty, OptTm = 86, OptLength = 55,
        };
        var probes = ProbeDesigner.DesignProbes(template, param, 20).ToList();
        Assert.That(probes, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            foreach (var p in probes)
                Assert.That(p.Primer3Penalty, Is.EqualTo(Math.Abs(p.Tm - 86) + Math.Abs(p.Sequence.Length - 55)).Within(1e-12));
            for (int i = 1; i < probes.Count; i++)
                Assert.That(probes[i - 1].Primer3Penalty!.Value, Is.LessThanOrEqualTo(probes[i].Primer3Penalty!.Value));
        });
    }

    private sealed class ProbeComparer : IEqualityComparer<ProbeDesigner.Probe>
    {
        public bool Equals(ProbeDesigner.Probe a, ProbeDesigner.Probe b) =>
            a.Sequence == b.Sequence && a.Start == b.Start && a.End == b.End && a.Tm == b.Tm && a.GcContent == b.GcContent
            && a.Score == b.Score && a.Type == b.Type && a.Warnings.SequenceEqual(b.Warnings);

        public int GetHashCode(ProbeDesigner.Probe p) => HashCode.Combine(p.Sequence, p.Start);
    }
}
