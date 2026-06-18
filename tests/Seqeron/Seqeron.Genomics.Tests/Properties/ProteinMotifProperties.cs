using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for protein motif finding: common motifs, PROSITE patterns,
/// domain prediction. Uses FsCheck for invariant verification with random protein sequences.
///
/// Test Units: PROTMOTIF-FIND-001, PROTMOTIF-PROSITE-001, PROTMOTIF-DOMAIN-001, PROTMOTIF-CC-001, PROTMOTIF-COMMON-001, PROTMOTIF-LC-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class ProteinMotifProperties
{
    #region Generators

    /// <summary>
    /// Generates random protein sequences from the 20 standard amino acids.
    /// </summary>
    private static Arbitrary<string> ProteinArbitrary(int minLen = 30) =>
        Gen.Elements('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// A protein with known N-glycosylation motif N-X-[ST]-X (where X≠P).
    /// The NIS and NQT triplets in this protein match N-{P}-[ST]-{P}.
    /// </summary>
    private static readonly string TestProtein =
        "MKTLLLTLVVVTLVLSSQPVLSRELRECPRGSGKSCQACPAG" +
        "NISTYQCQSYVMSHLCSYQCNQRCFQSLENQCQTFHCRGFQF" +
        "NSTRTMPLHCRGFQFNSTRTMPLHCRG";

    /// <summary>
    /// A protein with a zinc-finger-like CxxC motif for pattern matching.
    /// </summary>
    private static readonly string ZincFingerProtein =
        "MKCPICGKSFSQSSSLERHIRTHTGEKPYVC" +
        "ELCGKRFRDQANLIRHLRSHTGERPFQCEWC" +
        "GKTFSDKSNLTRHQRTHTGEKKFAC";

    #endregion

    #region PROTMOTIF-FIND-001: R: positions valid; M: broader pattern → ≥ matches; D: deterministic

    /// <summary>
    /// INV-1: All motif positions are within protein bounds.
    /// Evidence: Start and End indices reference valid positions in the protein string.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCommonMotifs_Positions_WithinBounds()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var motifs = ProteinMotifFinder.FindCommonMotifs(seq).ToList();
            return motifs.All(m => m.Start >= 0 && m.End <= seq.Length)
                .Label("All motif positions must be within [0, seqLen]");
        });
    }

    /// <summary>
    /// INV-2: Each found motif sequence is a substring of the protein.
    /// Evidence: Motif finding extracts substrings from the input protein.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCommonMotifs_Sequence_IsSubstringOfProtein()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var motifs = ProteinMotifFinder.FindCommonMotifs(seq).ToList();
            return motifs.All(m => seq.Contains(m.Sequence))
                .Label("Motif sequence must be a substring of the protein");
        });
    }

    /// <summary>
    /// INV-3: FindCommonMotifs is deterministic.
    /// Evidence: Regex-based pattern scanning is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCommonMotifs_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var m1 = ProteinMotifFinder.FindCommonMotifs(seq).ToList();
            var m2 = ProteinMotifFinder.FindCommonMotifs(seq).ToList();
            bool same = m1.Count == m2.Count &&
                        m1.Zip(m2).All(p => p.First.Start == p.Second.Start &&
                                             p.First.End == p.Second.End &&
                                             p.First.MotifName == p.Second.MotifName);
            return same.Label("FindCommonMotifs must be deterministic");
        });
    }

    /// <summary>
    /// INV-4: CommonMotifs dictionary is non-empty (built-in patterns exist).
    /// Evidence: The library ships with standard motif definitions.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CommonMotifs_Dictionary_NotEmpty()
    {
        Assert.That(ProteinMotifFinder.CommonMotifs.Count, Is.GreaterThan(0),
            "CommonMotifs dictionary should not be empty");
    }

    /// <summary>
    /// INV-5: FindMotifByPattern with a shorter (broader) pattern finds ≥ as many matches
    /// as a longer (more specific) pattern.
    /// Evidence: A longer regex is more restrictive, matching a subset of the shorter pattern's matches.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindMotifByPattern_BroaderPattern_MoreMatches()
    {
        string broad = "C..C";       // CxxC — any 2 residues between cysteines
        string narrow = "C[PIC]..C"; // C[PIC]xxC — first residue must be P, I, or C

        var broadMatches = ProteinMotifFinder.FindMotifByPattern(ZincFingerProtein, broad).ToList();
        var narrowMatches = ProteinMotifFinder.FindMotifByPattern(ZincFingerProtein, narrow).ToList();

        Assert.That(broadMatches.Count, Is.GreaterThanOrEqualTo(narrowMatches.Count),
            $"Broader pattern got {broadMatches.Count} matches, narrower got {narrowMatches.Count}");
    }

    #endregion

    #region PROTMOTIF-PROSITE-001: R: match positions valid; P: match conforms to PROSITE pattern regex; D: deterministic

    /// <summary>
    /// INV-6: FindMotifByProsite returns matches with valid positions.
    /// Evidence: PROSITE pattern scanning via regex produces valid boundary positions.
    /// Source: PROSITE database — patterns defined using PROSITE syntax notation.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindMotifByProsite_Positions_WithinBounds()
    {
        string prositePattern = "N-{P}-[ST]-{P}";  // N-glycosylation motif
        var motifs = ProteinMotifFinder.FindMotifByProsite(TestProtein, prositePattern).ToList();

        foreach (var m in motifs)
        {
            Assert.That(m.Start, Is.GreaterThanOrEqualTo(0),
                $"Start={m.Start} must be ≥ 0");
            Assert.That(m.End, Is.LessThanOrEqualTo(TestProtein.Length),
                $"End={m.End} must be ≤ {TestProtein.Length}");
        }
    }

    /// <summary>
    /// INV-7: N-glycosylation motif match starts with N (Asparagine).
    /// Evidence: PROSITE pattern N-{P}-[ST]-{P} mandates N at position 1.
    /// Source: PROSITE PS00001 — ASN_GLYCOSYLATION.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindMotifByProsite_NGlycosylation_StartsWithN()
    {
        string prositePattern = "N-{P}-[ST]-{P}";
        var motifs = ProteinMotifFinder.FindMotifByProsite(TestProtein, prositePattern).ToList();

        foreach (var m in motifs)
            Assert.That(m.Sequence[0], Is.EqualTo('N'),
                $"N-glycosylation motif must start with N, got '{m.Sequence}'");
    }

    /// <summary>
    /// INV-8: ConvertPrositeToRegex produces a valid regex.
    /// Evidence: PROSITE-to-regex conversion follows defined translation rules.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ConvertPrositeToRegex_ProducesValidRegex()
    {
        string prositePattern = "N-{P}-[ST]-{P}";
        string regex = ProteinMotifFinder.ConvertPrositeToRegex(prositePattern);

        Assert.That(regex, Is.Not.Null.And.Not.Empty);
        Assert.DoesNotThrow(() => new System.Text.RegularExpressions.Regex(regex),
            $"Converted regex '{regex}' is not valid");
    }

    /// <summary>
    /// INV-9: FindMotifByProsite is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindMotifByProsite_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            string pattern = "N-{P}-[ST]-{P}";
            var m1 = ProteinMotifFinder.FindMotifByProsite(seq, pattern).ToList();
            var m2 = ProteinMotifFinder.FindMotifByProsite(seq, pattern).ToList();
            bool same = m1.Count == m2.Count &&
                        m1.Zip(m2).All(p => p.First.Start == p.Second.Start &&
                                             p.First.End == p.Second.End);
            return same.Label("FindMotifByProsite must be deterministic");
        });
    }

    /// <summary>
    /// INV-10: PROSITE match sequence conforms to the expected regex.
    /// Evidence: Each match should be a valid instance of the PROSITE pattern.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindMotifByProsite_MatchConformsToRegex()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            string prositePattern = "N-{P}-[ST]-{P}";
            string regex = ProteinMotifFinder.ConvertPrositeToRegex(prositePattern);
            var re = new System.Text.RegularExpressions.Regex(regex);
            var motifs = ProteinMotifFinder.FindMotifByProsite(seq, prositePattern).ToList();
            return motifs.All(m => re.IsMatch(m.Sequence))
                .Label("PROSITE match sequence must conform to the converted regex");
        });
    }

    #endregion

    #region PROTMOTIF-DOMAIN-001: R: domain start < end; D: deterministic; P: domain score above threshold

    /// <summary>
    /// INV-11: Domain start &lt; end for all found domains.
    /// Evidence: A protein domain spans at least one amino acid residue.
    /// Source: Pfam database — domains are defined as contiguous regions.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindDomains_StartLessThanEnd()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var domains = ProteinMotifFinder.FindDomains(seq).ToList();
            return domains.All(d => d.Start < d.End)
                .Label("Domain start must be < end");
        });
    }

    /// <summary>
    /// INV-12: Domain positions are within protein bounds.
    /// Evidence: Domain boundaries reference valid positions in the input protein.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindDomains_Positions_WithinBounds()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var domains = ProteinMotifFinder.FindDomains(seq).ToList();
            return domains.All(d => d.Start >= 0 && d.End <= seq.Length)
                .Label("Domain positions must be within [0, seqLen]");
        });
    }

    /// <summary>
    /// INV-13: FindDomains is deterministic.
    /// Evidence: Pattern-based domain scanning uses deterministic regex matching.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindDomains_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var d1 = ProteinMotifFinder.FindDomains(seq).ToList();
            var d2 = ProteinMotifFinder.FindDomains(seq).ToList();
            bool same = d1.Count == d2.Count &&
                        d1.Zip(d2).All(p => p.First.Start == p.Second.Start &&
                                             p.First.End == p.Second.End &&
                                             p.First.Name == p.Second.Name);
            return same.Label("FindDomains must be deterministic");
        });
    }

    /// <summary>
    /// INV-14: Domain score is finite.
    /// Evidence: Scores are computed from regex match quality, yielding finite values.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindDomains_Score_IsFinite()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var domains = ProteinMotifFinder.FindDomains(seq).ToList();
            return domains.All(d => double.IsFinite(d.Score))
                .Label("Domain scores must be finite");
        });
    }

    // Signal-peptide prediction invariants are covered by Test Unit PROTMOTIF-SP-001
    // (ProteinMotifFinder_PredictSignalPeptide_Tests.cs), including the cleavage-position
    // bound property test, so they are not duplicated here.

    /// <summary>
    /// Transmembrane helix predictions have valid positions.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictTransmembraneHelices_Positions_WithinBounds()
    {
        string protein = "AAAILLLLLLLLLLLLLLLLLLLLAAAAGGGGG" +
                          "LLLLLLLLLLLLLLLLLLLLLLLLAAAA";
        var helices = ProteinMotifFinder.PredictTransmembraneHelices(protein).ToList();

        foreach (var (start, end, score) in helices)
        {
            Assert.That(start, Is.GreaterThanOrEqualTo(0));
            Assert.That(start, Is.LessThanOrEqualTo(end));
            Assert.That(end, Is.LessThan(protein.Length)); // 0-based inclusive index (INV-02: End ≤ length-1)
            Assert.That(double.IsFinite(score), Is.True);
        }
    }

    #endregion

    #region PROTMOTIF-CC-001: R: score ∈ [0,1]; P: heptad periodicity detected; D: deterministic

    // PredictCoiledCoils scores each window by the fraction of heptad a/d core positions occupied by a
    // hydrophobic core residue (I/L/V), maximised over the 7 registers (Lupas 1991; Mason & Arndt 2004).
    // Contiguous windows ≥ threshold form a region spanning at least 3 heptads (21 residues).

    /// <summary>
    /// INV-1 (R): every reported coiled-coil region has a score in [threshold,1], valid positions,
    /// and a length of at least 3 heptads (21 residues).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CoiledCoil_Regions_AreValid()
    {
        return Prop.ForAll(ProteinArbitrary(40), seq =>
        {
            var regions = ProteinMotifFinder.PredictCoiledCoils(seq).ToList();
            bool ok = regions.All(r =>
                r.Score >= 0.5 - 1e-9 && r.Score <= 1.0 + 1e-9 &&
                r.Start >= 0 && r.Start <= r.End && r.End < seq.Length &&
                r.End - r.Start + 1 >= 21);
            return ok.Label("a coiled-coil region had an out-of-range score, position, or length");
        });
    }

    /// <summary>
    /// INV-2 (P): A sequence with hydrophobic core residues placed at every heptad a/d position is
    /// detected as a coiled coil with peak occupancy 1.0; a sequence with no core residues is not.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CoiledCoil_HeptadPeriodicity_IsDetected()
    {
        // 6 heptads (42 residues): L at positions a(0) and d(3) of each heptad, A elsewhere.
        var chars = new char[42];
        for (int i = 0; i < 42; i++) chars[i] = (i % 7 == 0 || i % 7 == 3) ? 'L' : 'A';
        string coiled = new string(chars);
        string noCore = new string('A', 42);

        var regions = ProteinMotifFinder.PredictCoiledCoils(coiled).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(regions, Is.Not.Empty, "perfect heptad pattern must be detected");
            Assert.That(regions.Max(r => r.Score), Is.EqualTo(1.0).Within(1e-9), "a/d core fully occupied → score 1");
            Assert.That(ProteinMotifFinder.PredictCoiledCoils(noCore), Is.Empty,
                "no hydrophobic core residues → no coiled coil");
        });
    }

    /// <summary>
    /// INV-3 (D): Coiled-coil prediction is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CoiledCoil_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(40), seq =>
            ProteinMotifFinder.PredictCoiledCoils(seq).SequenceEqual(ProteinMotifFinder.PredictCoiledCoils(seq))
                .Label("PredictCoiledCoils must be deterministic"));
    }

    #endregion

    #region PROTMOTIF-COMMON-001: R: positions valid; P: each match conforms to its motif pattern; M: more occurrences → ≥ support; D: deterministic

    // FindCommonMotifs scans the curated PROSITE common-motif catalogue over a protein and reports
    // each occurrence (overlapping, ScanProsite style; De Castro et al. 2006).

    private static readonly IReadOnlyDictionary<string, string> MotifNameToRegex =
        ProteinMotifFinder.CommonMotifs.Values
            .GroupBy(m => m.Name)
            .ToDictionary(g => g.Key, g => g.First().RegexPattern);

    /// <summary>
    /// INV-1 (R): every match has valid coordinates and its reported Sequence equals the substring at
    /// those coordinates.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CommonMotifs_PositionsAndSubstring_AreValid()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
        {
            var matches = ProteinMotifFinder.FindCommonMotifs(seq).ToList();
            bool ok = matches.All(m =>
                m.Start >= 0 && m.Start <= m.End && m.End < seq.Length &&
                seq.Substring(m.Start, m.End - m.Start + 1) == m.Sequence);
            return ok.Label("a common-motif match had invalid coordinates or substring");
        });
    }

    /// <summary>
    /// INV-2 (P): every reported match conforms to its motif's PROSITE-derived regex pattern.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CommonMotifs_EachMatch_ConformsToPattern()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
        {
            var matches = ProteinMotifFinder.FindCommonMotifs(seq).ToList();
            bool ok = matches.All(m =>
                MotifNameToRegex.TryGetValue(m.MotifName, out string? rx) &&
                System.Text.RegularExpressions.Regex.IsMatch(
                    m.Sequence, "^(?:" + rx + ")$",
                    System.Text.RegularExpressions.RegexOptions.IgnoreCase));
            return ok.Label("a match did not conform to its motif regex");
        });
    }

    /// <summary>
    /// INV-3 (M, support): a motif-bearing protein yields matches, and duplicating it does not reduce
    /// the number of matches (more occurrences → ≥ support).
    /// </summary>
    [Test]
    [Category("Property")]
    public void CommonMotifs_MoreOccurrences_DoNotReduceSupport()
    {
        int single = ProteinMotifFinder.FindCommonMotifs(TestProtein).Count();
        int doubled = ProteinMotifFinder.FindCommonMotifs(TestProtein + TestProtein).Count();

        Assert.Multiple(() =>
        {
            Assert.That(single, Is.GreaterThan(0), "the test protein must contain common motifs");
            Assert.That(doubled, Is.GreaterThanOrEqualTo(single), "more occurrences must not reduce support");
        });
    }

    /// <summary>
    /// INV-4 (D): Common-motif scanning is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CommonMotifs_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
        {
            var a = ProteinMotifFinder.FindCommonMotifs(seq).Select(m => (m.Start, m.End, m.MotifName)).ToList();
            var b = ProteinMotifFinder.FindCommonMotifs(seq).Select(m => (m.Start, m.End, m.MotifName)).ToList();
            return a.SequenceEqual(b).Label("FindCommonMotifs must be deterministic");
        });
    }

    #endregion

    #region PROTMOTIF-LC-001: R: region start < end; M: higher complexity threshold → ≥ coverage; D: deterministic

    // FindLowComplexityRegions implements SEG (Wootton & Federhen 1993): a window with Shannon
    // complexity ≤ K1 triggers a region, extended over adjacent windows with complexity ≤ K2. NOTE:
    // because windows are flagged when complexity ≤ threshold, RAISING the threshold flags more, so
    // coverage is monotone increasing in the threshold (the checklist's wording is in the inverse
    // sense). The reported Complexity is the minimum window entropy in the region.

    private static HashSet<int> CoveredResidues(IEnumerable<(int Start, int End, double Complexity)> regions)
    {
        var set = new HashSet<int>();
        foreach (var (s, e, _) in regions)
            for (int p = s; p <= e; p++) set.Add(p);
        return set;
    }

    /// <summary>
    /// INV-1 (R + P): every region has Start &lt; End within bounds, spans at least one window, and its
    /// reported complexity does not exceed the trigger threshold (it is genuinely low-complexity).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LowComplexity_Regions_AreValidAndLowComplexity()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
        {
            var regions = ProteinMotifFinder.FindLowComplexityRegions(seq).ToList();
            bool ok = regions.All(r =>
                r.Start >= 0 && r.Start < r.End && r.End < seq.Length &&
                r.End - r.Start + 1 >= 12 &&
                r.Complexity <= 2.2 + 1e-9);
            return ok.Label("a low-complexity region was invalid or above the trigger complexity");
        });
    }

    /// <summary>
    /// INV-2 (M): Raising the complexity thresholds never reduces the residues flagged as
    /// low-complexity — the low-threshold coverage is a subset of the high-threshold coverage.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LowComplexity_HigherThreshold_CoversMore()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
        {
            var low = CoveredResidues(ProteinMotifFinder.FindLowComplexityRegions(seq, 12, 1.0, 1.5));
            var high = CoveredResidues(ProteinMotifFinder.FindLowComplexityRegions(seq, 12, 3.0, 3.5));
            return low.IsSubsetOf(high)
                .Label($"low-threshold coverage ({low.Count}) not ⊆ high-threshold coverage ({high.Count})");
        });
    }

    /// <summary>
    /// INV-3 (P, positive control): a homopolymer run is detected as a low-complexity region with
    /// complexity 0; a maximally diverse window is not flagged at the default thresholds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LowComplexity_Homopolymer_IsDetected()
    {
        var regions = ProteinMotifFinder.FindLowComplexityRegions(new string('G', 20)).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(regions, Is.Not.Empty, "a homopolymer must be low-complexity");
            Assert.That(regions.Min(r => r.Complexity), Is.EqualTo(0.0).Within(1e-9));
            // 20 distinct amino acids in a 20-residue window → maximal complexity, not flagged.
            Assert.That(ProteinMotifFinder.FindLowComplexityRegions("ACDEFGHIKLMNPQRSTVWY"), Is.Empty);
        });
    }

    /// <summary>
    /// INV-4 (D): Low-complexity detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LowComplexity_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
            ProteinMotifFinder.FindLowComplexityRegions(seq).SequenceEqual(ProteinMotifFinder.FindLowComplexityRegions(seq))
                .Label("FindLowComplexityRegions must be deterministic"));
    }

    #endregion
}
