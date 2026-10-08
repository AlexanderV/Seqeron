// SEQ-DINUC-001 — Dinucleotide Analysis (frequencies, O/E relative abundance, codon frequencies)
// Evidence: docs/Evidence/SEQ-DINUC-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-DINUC-001.md
// Source: Karlin S. (1998). Pervasive properties of the genomic signature. PMC126251.
//         Karlin S, Burge C (1995). Trends Genet 11(7):283-290.
//         Nakamura Y et al. Kazusa Codon Usage Database (CUTG) readme.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateDinucleotide_Tests
{
    // Expected values are exact rationals derived from the Karlin odds-ratio
    // definition and count/total frequency definitions in the Evidence, computed
    // independently of the implementation.
    private const double Tolerance = 1e-10;

    #region CalculateDinucleotideRatios

    // M1 — rho_XY = f_XY/(f_X f_Y) for ATGCGCGT (A=1,T=2,G=3,C=2; N=8; 7 dinucleotide positions).
    // Evidence: Karlin PMC126251. rho_GC=rho_CG=(2/7)/((3/8)(2/8))=64/21; rho_AT=(1/7)/((1/8)(2/8))=32/7;
    //           rho_TG=rho_GT=(1/7)/((2/8)(3/8))=(1/7)/((3/8)(2/8))=32/21.
    [Test]
    public void CalculateDinucleotideRatios_AtgcgcgtSequence_ReturnsExactOddsRatios()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGT");

        Assert.Multiple(() =>
        {
            Assert.That(ratios["GC"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
                "rho_GC = (2/7)/((3/8)(2/8)) = 64/21 per Karlin odds ratio (INV-03)");
            Assert.That(ratios["CG"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
                "rho_CG = (2/7)/((2/8)(3/8)) = 64/21");
            Assert.That(ratios["AT"], Is.EqualTo(32.0 / 7.0).Within(Tolerance),
                "rho_AT = (1/7)/((1/8)(2/8)) = 32/7");
            Assert.That(ratios["TG"], Is.EqualTo(32.0 / 21.0).Within(Tolerance),
                "rho_TG = (1/7)/((2/8)(3/8)) = 32/21");
            Assert.That(ratios["GT"], Is.EqualTo(32.0 / 21.0).Within(Tolerance),
                "rho_GT = (1/7)/((3/8)(2/8)) = 32/21");
        });
    }

    // M7 — no-bias baseline: homopolymer AAAA has f_AA=1 and f_A=1, so rho_AA = 1/(1*1) = 1.0 exactly.
    // Evidence: Karlin PMC126251 (r=1 means no bias, INV-03).
    [Test]
    public void CalculateDinucleotideRatios_Homopolymer_ReturnsOddsRatioOne()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("AAAA");

        Assert.That(ratios["AA"], Is.EqualTo(1.0).Within(Tolerance),
            "rho_AA = f_AA/(f_A*f_A) = 1/(1*1) = 1.0, the no-bias baseline (INV-03)");
    }

    // S4 — robustness when some bases are absent: a present dinucleotide always has both of its
    //      constituent bases present, so the expected product f_X*f_Y is strictly positive and the
    //      ratio is finite (no division by zero) even when other bases (here G, C) never occur.
    //      "ATAT": A=2,T=2,G=0,C=0. Dinucs: AT,TA,AT. f_AT=2/3, f_A=2/4, f_T=2/4 -> rho_AT=(2/3)/(1/4)=8/3.
    //      f_TA=1/3, rho_TA=(1/3)/((2/4)(2/4))=4/3. Independently recomputed with exact rationals.
    // Evidence: Karlin PMC126251 (odds-ratio definition); no infinity arises because expected>0 for
    //           every dinucleotide that actually occurs.
    [Test]
    public void CalculateDinucleotideRatios_AbsentBaseContext_DoesNotProduceInfinity()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATAT");

        Assert.Multiple(() =>
        {
            Assert.That(ratios["AT"], Is.EqualTo(8.0 / 3.0).Within(Tolerance),
                "rho_AT = (2/3)/((2/4)(2/4)) = 8/3");
            Assert.That(ratios["TA"], Is.EqualTo(4.0 / 3.0).Within(Tolerance),
                "rho_TA = (1/3)/((2/4)(2/4)) = 4/3");
            Assert.That(ratios.Keys, Is.EquivalentTo(new[] { "AT", "TA" }),
                "only dinucleotides over present bases appear; absent bases G and C yield no keys");
            Assert.That(ratios.Values, Has.All.Matches<double>(double.IsFinite),
                "no ratio is infinite/NaN even though G and C are absent");
        });
    }

    // S5 — non-alphabet (ambiguous) bases are excluded from ratio inputs. "ATGNCG": the N breaks
    //      the GN and NC contexts, leaving valid dinucleotides AT, TG, CG (each f = 1/3). Single-base
    //      frequencies are over the 5 unambiguous bases (A1,T1,G2,C1; N excluded), N_total = 5.
    //      rho_AT = (1/3)/((1/5)(1/5)) = 25/3; rho_TG = (1/3)/((1/5)(2/5)) = 25/6;
    //      rho_CG = (1/3)/((1/5)(2/5)) = 25/6. Independently recomputed with exact rationals.
    // Evidence: Karlin PMC126251 (normalized base/dinucleotide frequencies over the valid alphabet).
    [Test]
    public void CalculateDinucleotideRatios_AmbiguousBase_ExcludedFromCounts()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATGNCG");

        Assert.Multiple(() =>
        {
            Assert.That(ratios.Keys, Is.EquivalentTo(new[] { "AT", "TG", "CG" }),
                "GN and NC are excluded; only ACGT dinucleotides remain");
            Assert.That(ratios["AT"], Is.EqualTo(25.0 / 3.0).Within(Tolerance),
                "rho_AT = (1/3)/((1/5)(1/5)) = 25/3");
            Assert.That(ratios["TG"], Is.EqualTo(25.0 / 6.0).Within(Tolerance),
                "rho_TG = (1/3)/((1/5)(2/5)) = 25/6");
            Assert.That(ratios["CG"], Is.EqualTo(25.0 / 6.0).Within(Tolerance),
                "rho_CG = (1/3)/((1/5)(2/5)) = 25/6");
        });
    }

    // S1 — guards: null, empty, and length < 2 all return an empty dictionary.
    [Test]
    public void CalculateDinucleotideRatios_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios(null!), Is.Empty,
                "null input yields empty dictionary (input guard)");
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios(string.Empty), Is.Empty,
                "empty input yields empty dictionary (input guard)");
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios("A"), Is.Empty,
                "length < 2 yields empty dictionary (no dinucleotide positions)");
        });
    }

    // C1 — case-insensitivity: lowercase input produces the same ratios as uppercase (M1 values).
    [Test]
    public void CalculateDinucleotideRatios_LowercaseInput_MatchesUppercase()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("atgcgcgt");

        Assert.That(ratios["GC"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
            "lowercase is normalized via ToUpperInvariant; rho_GC equals the uppercase value 64/21");
    }

    #endregion

    #region CalculateDinucleotideFrequencies

    // M2 — normalized dinucleotide frequencies for ATGCGCGT (7 positions): count/(N-1).
    // Evidence: Karlin PMC126251 (f_XY normalized). GC=CG=2/7; AT=TG=GT=1/7.
    [Test]
    public void CalculateDinucleotideFrequencies_AtgcgcgtSequence_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        Assert.Multiple(() =>
        {
            Assert.That(freq["GC"], Is.EqualTo(2.0 / 7.0).Within(Tolerance),
                "f_GC = 2/7 (2 of 7 dinucleotide positions)");
            Assert.That(freq["CG"], Is.EqualTo(2.0 / 7.0).Within(Tolerance),
                "f_CG = 2/7");
            Assert.That(freq["AT"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_AT = 1/7");
            Assert.That(freq["TG"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_TG = 1/7");
            Assert.That(freq["GT"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_GT = 1/7");
        });
    }

    // M3 — INV-01: dinucleotide frequencies sum to 1.0.
    [Test]
    public void CalculateDinucleotideFrequencies_AtgcgcgtSequence_SumsToOne()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
            "normalized frequencies over all valid dinucleotide positions sum to 1.0 (INV-01)");
    }

    // S2 — guards: null, empty, length < 2 return empty.
    [Test]
    public void CalculateDinucleotideFrequencies_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(null!), Is.Empty,
                "null input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(string.Empty), Is.Empty,
                "empty input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies("A"), Is.Empty,
                "length < 2 yields empty dictionary");
        });
    }

    // S5f — non-alphabet (ambiguous) bases excluded from frequency counts. "ATGNCG": dinucleotides
    //       GN and NC contain N and are dropped, leaving AT, TG, CG over 3 valid positions, each 1/3.
    //       Independently recomputed (count/(valid positions)).
    // Evidence: Karlin PMC126251 (frequency over valid dinucleotide positions); ambiguous units excluded.
    [Test]
    public void CalculateDinucleotideFrequencies_AmbiguousBase_ExcludedFromCounts()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGNCG");

        Assert.Multiple(() =>
        {
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "AT", "TG", "CG" }),
                "GN and NC contain the ambiguous base N and are excluded");
            Assert.That(freq["AT"], Is.EqualTo(1.0 / 3.0).Within(Tolerance), "f_AT = 1/3");
            Assert.That(freq["TG"], Is.EqualTo(1.0 / 3.0).Within(Tolerance), "f_TG = 1/3");
            Assert.That(freq["CG"], Is.EqualTo(1.0 / 3.0).Within(Tolerance), "f_CG = 1/3");
            Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
                "valid-position frequencies still sum to 1.0 (INV-01)");
        });
    }

    // C2 — RNA U handling: U is a valid base; AUGCGC includes the AU dinucleotide.
    [Test]
    public void CalculateDinucleotideFrequencies_RnaInput_IncludesUracilDinucleotides()
    {
        // AUGCGC -> dinucs AU,UG,GC,CG,GC (5). AU=1/5, GC=2/5.
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("AUGCGC");

        Assert.Multiple(() =>
        {
            Assert.That(freq.ContainsKey("AU"), Is.True,
                "U is in the dinucleotide alphabet, so AU is counted (RNA support)");
            Assert.That(freq["AU"], Is.EqualTo(1.0 / 5.0).Within(Tolerance),
                "f_AU = 1/5 (1 of 5 dinucleotide positions)");
            Assert.That(freq["GC"], Is.EqualTo(2.0 / 5.0).Within(Tolerance),
                "f_GC = 2/5");
        });
    }

    #endregion

    #region Reference implementations (seqinr rho, EMBOSS compseq) and Karlin rho* (review 2026-09, F13)

    private static void AssertRatios(IReadOnlyDictionary<string, double> actual,
        IReadOnlyDictionary<string, double> expected, string source)
    {
        Assert.Multiple(() =>
        {
            Assert.That(actual.Keys, Is.EquivalentTo(expected.Keys), $"key set ({source})");
            foreach (var (k, v) in expected)
                Assert.That(actual.GetValueOrDefault(k, double.NaN), Is.EqualTo(v).Within(Tolerance), $"{k} ({source})");
        });
    }

    // EMBOSS 6.6.0 `compseq -word 2 -calcfreq` Obs/Exp (executed) = seqinr rho (cran/seqinr R/rho.R,
    // ported literally) for pure-ACGT input: 32/27 = 1.1851852, 16/9 = 1.7777778, 16/27 = 0.5925926.
    [Test]
    public void CalculateDinucleotideRatios_MatchesEmbossCompseqAndSeqinrRho()
    {
        var expected = new Dictionary<string, double>
        {
            ["AA"] = 32.0 / 27, ["AG"] = 32.0 / 27, ["AT"] = 16.0 / 9, ["CA"] = 32.0 / 27,
            ["CC"] = 16.0 / 9, ["CG"] = 16.0 / 27, ["CT"] = 16.0 / 27, ["GA"] = 16.0 / 27,
            ["GC"] = 16.0 / 9, ["GG"] = 32.0 / 27, ["TA"] = 32.0 / 27, ["TC"] = 16.0 / 27,
            ["TG"] = 16.0 / 27, ["TT"] = 16.0 / 9,
        };
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("GATCCTTAAAGGCGCATTTAGGCCCATG"),
            expected, "compseq/seqinr");
    }

    // seqinr rho/count: N, IUPAC codes and case are dropped from both numerator and denominators
    // (11 in-alphabet dinucleotides, 14 in-alphabet bases). ρ_CG = (3/11)/((3/14)(4/14)) = 49/11.
    [Test]
    public void CalculateDinucleotideRatios_AmbiguityAndLowercase_MatchesSeqinrRho()
    {
        const string seq = "ATGCGCGTNNacgtRAT";
        var expected = new Dictionary<string, double>
        {
            ["AC"] = 196.0 / 99, ["AT"] = 98.0 / 33, ["CG"] = 49.0 / 11,
            ["GC"] = 98.0 / 33, ["GT"] = 49.0 / 22, ["TG"] = 49.0 / 44,
        };
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios(seq), expected, "seqinr rho");

        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        Assert.Multiple(() =>
        {
            Assert.That(freq["CG"], Is.EqualTo(3.0 / 11).Within(Tolerance), "seqinr count(freq=TRUE): 3/11");
            Assert.That(freq["AT"], Is.EqualTo(2.0 / 11).Within(Tolerance), "seqinr count(freq=TRUE): 2/11");
        });
    }

    // Karlin rho*: f*_A=f*_T=(f_A+f_T)/2, f*_C=f*_G=(f_C+f_G)/2, f*_GT=(f_GT+f_AC)/2 (Karlin PMC126251; restated
    // in PMC1829274, WebSearch snippet). Reference = seqinr rho on seq + separator + reverse complement (executed).
    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetric_MatchesKarlinRhoStar()
    {
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGT", strandSymmetric: true),
            new Dictionary<string, double>
            {
                ["AC"] = 128.0 / 105, ["AT"] = 256.0 / 63, ["CA"] = 128.0 / 105, ["CG"] = 512.0 / 175,
                ["GC"] = 512.0 / 175, ["GT"] = 128.0 / 105, ["TG"] = 128.0 / 105,
            }, "rho* ATGCGCGT");

        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("GATCCTTAAAGGCGCATTTAGGCCCATG", strandSymmetric: true),
            new Dictionary<string, double>
            {
                ["AA"] = 40.0 / 27, ["AG"] = 8.0 / 9, ["AT"] = 16.0 / 9, ["CA"] = 8.0 / 9,
                ["CC"] = 40.0 / 27, ["CG"] = 16.0 / 27, ["CT"] = 8.0 / 9, ["GA"] = 16.0 / 27,
                ["GC"] = 16.0 / 9, ["GG"] = 40.0 / 27, ["TA"] = 32.0 / 27, ["TC"] = 16.0 / 27,
                ["TG"] = 8.0 / 9, ["TT"] = 40.0 / 27,
            }, "rho* 28-mer");

        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGTNNacgtRAT", strandSymmetric: true),
            new Dictionary<string, double>
            {
                ["AC"] = 24.0 / 11, ["AT"] = 32.0 / 11, ["CA"] = 8.0 / 11, ["CG"] = 48.0 / 11,
                ["GC"] = 32.0 / 11, ["GT"] = 24.0 / 11, ["TG"] = 8.0 / 11,
            }, "rho* with N/R/lowercase");
    }

    // rho* is strand-invariant: rho*_XY = rho*_{revcomp(XY)} and rho*(S) = rho*(revcomp(S)).
    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetric_IsReverseComplementInvariant()
    {
        const string seq = "GATCCTTAAAGGCGCATTTAGGCCCATGNNACGTTGCA";
        var forward = SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: true);
        var reverse = SequenceStatistics.CalculateDinucleotideRatios(
            DnaSequence.GetReverseComplementString(seq), strandSymmetric: true);

        AssertRatios(reverse, forward, "rho*(revcomp S)");
        Assert.Multiple(() =>
        {
            Assert.That(forward["GT"], Is.EqualTo(forward["AC"]).Within(Tolerance), "rho*_GT = rho*_AC");
            Assert.That(forward["AG"], Is.EqualTo(forward["CT"]).Within(Tolerance), "rho*_AG = rho*_CT");
        });
    }

    // rho* is a double-stranded DNA statistic: U is excluded like N. "AUGC": only GC counted
    // (symmetric count 2 of 2), bases A,G,C -> f*_C=f*_G=2/6 -> rho*_GC = 1/(1/9) = 9.
    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetric_ExcludesUracil()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("AUGC", strandSymmetric: true);
        Assert.Multiple(() =>
        {
            Assert.That(ratios.Keys, Is.EquivalentTo(new[] { "GC" }));
            Assert.That(ratios["GC"], Is.EqualTo(9.0).Within(Tolerance));
        });
    }

    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetricFalse_EqualsSingleArgumentOverload()
    {
        foreach (var seq in new[] { "ATGCGCGT", "AUGCGC", "ATGCGCGTNNacgtRAT", "A", "", "NNNN" })
            AssertRatios(SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: false),
                SequenceStatistics.CalculateDinucleotideRatios(seq), seq);
    }

    // D3: counting delegates to the canonical overlapping k-mer counter (k = 2), filtered to {A,C,G,T,U}.
    [Test]
    public void CalculateDinucleotideFrequencies_EqualsCanonicalKmerCounterFiltered()
    {
        var rng = new Random(20260928);
        const string alphabet = "ACGTUNacgtu-RY";
        for (int n = 0; n < 200; n++)
        {
            var seq = new string(Enumerable.Range(0, rng.Next(2, 60)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            var counts = KmerAnalyzer.CountKmers(seq, 2)
                .Where(kv => kv.Key.All(c => "ACGTU".Contains(c)))
                .ToDictionary(kv => kv.Key, kv => kv.Value);
            int total = counts.Values.Sum();
            var expected = counts.ToDictionary(kv => kv.Key, kv => (double)kv.Value / total);
            AssertRatios(SequenceStatistics.CalculateDinucleotideFrequencies(seq), expected, seq);
        }
    }

    #endregion
}
