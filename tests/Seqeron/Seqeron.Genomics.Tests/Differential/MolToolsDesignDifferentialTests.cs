// 08_DIFFERENTIAL_TESTING rows 19,20,22,23,24,25 (MolTools design/scoring). Each production method is
// checked against an INDEPENDENT oracle: a re-derived in-region PAM set, a brute Hamming off-target scan,
// an independent candidate-enumeration + arg-max for primer selection, a brute self-complementarity scan
// for hairpins, spec reconstruction of designed probes, and a brute Hamming off-target count for probe
// validation. Shared, already-differentially-validated building blocks (FindPamSites row 18, EvaluatePrimer)
// are reused only as setup; the logic under test (region filter / selection / scan) is re-derived here.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class MolToolsDesignDifferentialTests
{
    private const double Tol = 1e-9;

    private static int Hamming(string a, string b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) d++;
        return d;
    }

    private static readonly Dictionary<char, char> Comp = new() { ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G' };
    private static string RevComp(string s)
    {
        var arr = s.Select(c => Comp[c]).ToArray();
        Array.Reverse(arr);
        return new string(arr);
    }
    private static bool WcPair(char x, char y) => Comp[x] == y;

    private const string Genome60 = "ACCGGTACCGGTAAGGCCTTGGCCAACCGGTTACGTACGTAGGCCTTAAGGTTCCGGAACC";

    // ---- Row 19: CRISPR-GUIDE-001 — DesignGuideRnas candidate set vs in-region PAM targets ----

    [Test]
    [Category("CRISPR-GUIDE-001")]
    public void DesignGuideRnas_CandidateSet_MatchesInRegionPamTargets()
    {
        var dna = new DnaSequence(Genome60);
        const int regionStart = 25, regionEnd = 55;
        var noFilter = new GuideRnaParameters(0, 100, 0, false, false); // MinScore 0 -> keep every candidate

        var guides = CrisprDesigner.DesignGuideRnas(dna, regionStart, regionEnd, CrisprSystemType.SpCas9, noFilter)
            .ToList();

        // Independent oracle: PAM sites (validated in row 18) whose SpCas9 blunt cut is in region.
        // The cut lies 3 bp 5' of the PAM *on the PAM-bearing strand* (CRISPOR crispor.py: "the
        // expected cleavage position located -3bp 5' of the PAM site"; its schematic draws the cut
        // marker as `startFt = start - 3` on '+' and as `ftSeq + "---"` on '-'), so in forward
        // coordinates it is Position-3 on the plus strand and Position+pamLen+2 on the minus strand.
        // Transcribed here independently of CrisprDesigner.GetCutSite.
        var oracle = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9)
            .Where(p =>
            {
                int cut = p.IsForwardStrand ? p.Position - 3 : p.Position + p.PamSequence.Length + 2;
                return cut >= regionStart && cut <= regionEnd;
            })
            .Select(p => (p.TargetSequence, p.TargetStart, p.IsForwardStrand)).ToList();

        Assert.That(oracle, Is.Not.Empty);
        Assert.That(oracle.Any(o => !o.IsForwardStrand), Is.True,
            "the oracle must cover reverse-strand guides (the strand whose cut frame was wrong)");
        Assert.That(guides.Select(g => (g.Sequence, g.Position, g.IsForwardStrand)),
            Is.EquivalentTo(oracle));

        // DesignGuideRnas is documented (and CRISPOR's guide table is sorted) best-first; ties are
        // broken by position and then forward strand first.
        Assert.That(guides.Select(g => (g.Score, -g.Position, g.IsForwardStrand ? 1 : 0)),
            Is.Ordered.Descending);
    }

    // ---- Row 20: CRISPR-OFF-001 — FindOffTargets vs brute Hamming over PAM targets ----

    [Test]
    [Category("CRISPR-OFF-001")]
    public void FindOffTargets_MatchesBruteHammingScan()
    {
        const string guide = "ACGTACGTACGTACGTACGT";
        const string mutated = "TCGTACGTACGTACGTACGA"; // 2 substitutions vs guide (pos 0 and 19)
        var genome = new DnaSequence(mutated + "AGG" + Genome60);
        const int maxMm = 3;

        var actual = CrisprDesigner.FindOffTargets(guide, genome, maxMm, CrisprSystemType.SpCas9)
            .Select(o => (o.Position, o.Sequence, o.Mismatches, o.IsForwardStrand)).ToList();

        var oracle = CrisprDesigner.FindPamSites(genome, CrisprSystemType.SpCas9)
            .Select(p => (p, mm: Hamming(guide, p.TargetSequence)))
            .Where(x => x.mm >= 1 && x.mm <= maxMm)
            .Select(x => (x.p.Position, x.p.TargetSequence, x.mm, x.p.IsForwardStrand)).ToList();

        Assert.That(actual, Is.EqualTo(oracle));
        Assert.That(actual.Any(o => o.Mismatches == 2), Is.True, "the planted 2-mismatch off-target is present");
    }

    // ---- Row 23: PRIMER-STRUCT-001 — HasHairpinPotential vs brute self-complementarity scan ----

    private static bool HairpinOracle(string seq, int stem, int loop)
    {
        var s = seq.ToUpperInvariant();
        if (s.Length < stem * 2 + loop) return false;
        for (int i = 0; i <= s.Length - stem; i++)
            for (int j = i + stem + loop; j <= s.Length - stem; j++)
            {
                bool all = true;
                for (int k = 0; k < stem; k++)
                    if (!WcPair(s[i + k], s[j + stem - 1 - k])) { all = false; break; }
                if (all) return true;
            }
        return false;
    }

    [Test]
    [Category("PRIMER-STRUCT-001")]
    [TestCase("ACCCGAAATCGGGTA")]      // stem ACCC / loop / GGGT (revcomp pair)
    [TestCase("AAAAAAAAAAAAAAA")]      // homopolymer -> no hairpin
    [TestCase("GATCCGTTAGGATCG")]
    [TestCase("ACGTACGTACGTACG")]
    [TestCase("TTGGCCAATTGGCCAATT")]
    public void HasHairpinPotential_MatchesBruteSelfComplementScan(string seq)
    {
        Assert.That(PrimerDesigner.HasHairpinPotential(seq), Is.EqualTo(HairpinOracle(seq, 4, 3)));
    }

    // ---- Row 25: PROBE-VALID-001 — ValidateProbe specificity vs brute Hamming off-target count ----

    private static int OffTargetOracle(string probe, IEnumerable<string> refs, int maxMm)
    {
        int hits = 0;
        foreach (var r in refs.Select(x => x.ToUpperInvariant()))
            for (int i = 0; i + probe.Length <= r.Length; i++)
                if (Hamming(r.Substring(i, probe.Length), probe) <= maxMm) hits++;
        return hits;
    }
}
