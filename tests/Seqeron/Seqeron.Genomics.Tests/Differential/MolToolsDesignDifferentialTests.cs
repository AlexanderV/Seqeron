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
            .Select(g => (g.Sequence, g.Position, g.IsForwardStrand)).ToList();

        // Independent oracle: PAM sites (validated in row 18) whose Cas9 cut (PAM.Position-3) is in region.
        var oracle = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9)
            .Where(p => { int cut = p.Position - 3; return cut >= regionStart && cut <= regionEnd; })
            .Select(p => (p.TargetSequence, p.TargetStart, p.IsForwardStrand)).ToList();

        Assert.That(guides, Is.EqualTo(oracle));
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

    // ---- Row 22: PRIMER-DESIGN-001 — DesignPrimers vs brute-force Primer3 pair search ----

    [Test]
    [Category("PRIMER-DESIGN-001")]
    public void DesignPrimers_MatchesBruteForcePrimer3PairSearch()
    {
        const string template =
            "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAGGCCAATTGGCCAATTGGCCAATTACGTACGTACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        var dna = new DnaSequence(template);
        int targetStart = 50, targetEnd = 62;
        var param = PrimerDesigner.DefaultParameters;

        var result = PrimerDesigner.DesignPrimers(dna, targetStart, targetEnd, param);

        // Independent enumeration of all valid candidates on each side.
        var fwd = new List<PrimerCandidate>();
        for (int start = Math.Max(0, targetStart - 200); start < targetStart; start++)
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= targetStart; len++)
            {
                var c = PrimerDesigner.EvaluatePrimer(template.Substring(start, len), start, true, param);
                if (c.IsValid) fwd.Add(c);
            }
        var rev = new List<PrimerCandidate>();
        int rEnd = Math.Min(template.Length, targetEnd + 200);
        for (int end = targetEnd + param.MinLength; end <= rEnd; end++)
            for (int len = param.MinLength; len <= param.MaxLength && end - len >= targetEnd; len++)
            {
                int start = end - len;
                var c = PrimerDesigner.EvaluatePrimer(RevComp(template.Substring(start, len)), start, false, param);
                if (c.IsValid) rev.Add(c);
            }

        // Exhaustive Primer3 pair search: minimise left+right penalty over all pairs meeting the pair
        // constraints; ties (1e-6) broken as libprimer3.cc compare_primer_pair.
        (PrimerCandidate F, PrimerCandidate R, double Q)? best = null;
        foreach (var f in fwd)
            foreach (var r in rev)
            {
                double dTm = Math.Abs(PrimerDesigner.CalculateMeltingTemperaturePrimer3(f.Sequence)
                                      - PrimerDesigner.CalculateMeltingTemperaturePrimer3(r.Sequence));
                if (dTm > 5.0 || PrimerDesigner.HasPrimerDimer(f.Sequence, r.Sequence)) continue;
                double q = f.Penalty + r.Penalty;
                if (best is null || Better(q, f, r, best.Value.Q, best.Value.F, best.Value.R))
                    best = (f, r, q);
            }

        Assert.That(best, Is.Not.Null, "fixture must admit a valid pair");
        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.True);
            Assert.That(result.Forward?.Sequence, Is.EqualTo(best!.Value.F.Sequence), "forward primer");
            Assert.That(result.Forward?.Position, Is.EqualTo(best.Value.F.Position), "forward position");
            Assert.That(result.Reverse?.Sequence, Is.EqualTo(best.Value.R.Sequence), "reverse primer");
            Assert.That(result.Reverse?.Position, Is.EqualTo(best.Value.R.Position), "reverse position");
            Assert.That(result.ProductSize,
                Is.EqualTo(best.Value.R.Position + best.Value.R.Length - best.Value.F.Position));
        });

        static bool Better(double q1, PrimerCandidate l1, PrimerCandidate r1, double q2, PrimerCandidate l2, PrimerCandidate r2)
        {
            if (q1 + 1e-6 < q2) return true;
            if (q1 > q2 + 1e-6) return false;
            if (l1.Position != l2.Position) return l1.Position > l2.Position;
            int e1 = r1.Position + r1.Length, e2 = r2.Position + r2.Length;
            if (e1 != e2) return e1 < e2;
            if (l1.Length != l2.Length) return l1.Length < l2.Length;
            return r1.Length < r2.Length;
        }
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

    // ---- Row 24: PROBE-DESIGN-001 — designed probes reconstructible + correctly ranked ----

    [Test]
    [Category("PROBE-DESIGN-001")]
    public void DesignProbes_ProbesReconstructibleAndRanked()
    {
        const string target =
            "ACGTGCATGCATGCGCATGCATCGATCGTAGCTAGCATGCATGCGCATGCTAGCATCGATCGTAGCTAGCATGCATGCATGCGC";
        var upper = target.ToUpperInvariant();
        var param = ProbeDesigner.Defaults.Microarray; // length 50-60, GC 0.40-0.60
        const int maxProbes = 10;

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes).ToList();
        Assert.That(probes.Count, Is.LessThanOrEqualTo(maxProbes));

        double prevScore = double.PositiveInfinity;
        foreach (var p in probes)
        {
            // Probe is an actual window of the target at its reported start.
            Assert.That(p.Sequence, Is.EqualTo(upper.Substring(p.Start, p.Sequence.Length)), "probe is a real window");
            Assert.That(p.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength), "length in range");

            // GC content (fraction) re-computed independently.
            double gc = (double)p.Sequence.Count(c => c == 'G' || c == 'C') / p.Sequence.Length;
            Assert.That(p.GcContent, Is.EqualTo(gc).Within(Tol), "GC fraction");
            Assert.That(gc, Is.InRange(param.MinGc, param.MaxGc), "GC within parameter band");

            // Returned in non-increasing score order.
            Assert.That(p.Score, Is.LessThanOrEqualTo(prevScore), "ranked by score descending");
            prevScore = p.Score;
        }
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

    [Test]
    [Category("PROBE-VALID-001")]
    public void ValidateProbe_SpecificityMatchesBruteHammingCount()
    {
        const string probe = "ACGTGCAT";
        const int maxMm = 2;

        foreach (var refs in new[]
        {
            new[] { "TTTTACGTGCATTTTT" },                 // exactly one match -> specificity 1.0
            new[] { "ACGTGCATxxxxACGTGCAT".Replace("x", "A") }, // two exact matches -> 1/2
            new[] { "TTTTTTTTTTTT" },                      // no match -> 0.0
        })
        {
            var v = ProbeDesigner.ValidateProbe(probe, refs, maxMm);
            int oracle = OffTargetOracle(probe, refs, maxMm);
            Assert.That(v.OffTargetHits, Is.EqualTo(oracle), "off-target count");

            double expectedSpec = oracle == 0 ? 0.0 : oracle == 1 ? 1.0 : 1.0 / oracle;
            Assert.That(v.SpecificityScore, Is.EqualTo(expectedSpec).Within(Tol), "specificity score");
        }
    }
}
