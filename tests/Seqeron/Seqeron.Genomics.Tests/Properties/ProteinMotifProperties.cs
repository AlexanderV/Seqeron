using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for protein motif finding: common motifs, PROSITE patterns,
/// domain prediction. Uses FsCheck for invariant verification with random protein sequences.
///
/// Test Units: PROTMOTIF-FIND-001, PROTMOTIF-PROSITE-001, PROTMOTIF-DOMAIN-001
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

    /// <summary>
    /// Signal peptide prediction: cleavage position within bounds (if found).
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictSignalPeptide_CleavagePosition_WithinBounds()
    {
        string proteinWithSignal = "MKTLLLTLVVVTLVLSSQPVLSRELRECPRGSGKSCQACPAG";
        var sp = ProteinMotifFinder.PredictSignalPeptide(proteinWithSignal);

        if (sp.HasValue)
        {
            Assert.That(sp.Value.CleavagePosition, Is.GreaterThan(0));
            Assert.That(sp.Value.CleavagePosition, Is.LessThan(proteinWithSignal.Length));
            Assert.That(sp.Value.Score, Is.InRange(0.0, 1.0));
        }
    }

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
            Assert.That(end, Is.LessThanOrEqualTo(protein.Length));
            Assert.That(double.IsFinite(score), Is.True);
        }
    }

    #endregion
}
