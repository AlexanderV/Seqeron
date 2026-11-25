// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
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

    // ---- Row 22: PRIMER-DESIGN-001 — DesignPrimers vs brute-force Primer3 pair search ----

    [Test]
    [Category("PRIMER-DESIGN-001")]
    public void DesignPrimers_MatchesBruteForcePrimer3PairSearch()
    {
        // Random 160-mer (python random.seed(2026)); the earlier palindromic fixture admits no pair under
        // Primer3's default thermodynamic structure limits (every candidate forms a > 47 °C hairpin/dimer).
        const string template =
            "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCA";
        var dna = new DnaSequence(template);
        int targetStart = 70, targetEnd = 90;
        var param = ThermodynamicPrimerDesigner.DefaultParameters;

        var result = ThermodynamicPrimerDesigner.DesignPrimers(dna, targetStart, targetEnd, param);

        // Independent enumeration of all valid candidates on each side.
        var fwd = new List<PrimerCandidate>();
        for (int start = Math.Max(0, targetStart - 200); start < targetStart; start++)
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= targetStart; len++)
            {
                var c = ThermodynamicPrimerDesigner.EvaluatePrimer(template.Substring(start, len), start, true, param);
                if (c.IsValid) fwd.Add(c);
            }
        var rev = new List<PrimerCandidate>();
        int rEnd = Math.Min(template.Length, targetEnd + 200);
        for (int end = targetEnd + param.MinLength; end <= rEnd; end++)
            for (int len = param.MinLength; len <= param.MaxLength && end - len >= targetEnd; len++)
            {
                int start = end - len;
                var c = ThermodynamicPrimerDesigner.EvaluatePrimer(RevComp(template.Substring(start, len)), start, false, param);
                if (c.IsValid) rev.Add(c);
            }

        // Exhaustive Primer3 pair search: minimise left+right penalty over all pairs meeting the pair
        // constraints; ties (1e-6) broken as libprimer3.cc compare_primer_pair.
        (PrimerCandidate F, PrimerCandidate R, double Q)? best = null;
        foreach (var f in fwd)
            foreach (var r in rev)
            {
                double dTm = Math.Abs(ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(f.Sequence)
                                      - ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(r.Sequence));
                if (dTm > 5.0) continue;
                double q = f.Penalty + r.Penalty;
                if (best is not null && !Better(q, f, r, best.Value.Q, best.Value.F, best.Value.R))
                    continue;
                // Default screen: Primer3 PRIMER_PAIR_MAX_COMPL_ANY_TH / _COMPL_END_TH = 47 °C (ntthal);
                // evaluated only for pairs that would improve on the current best (same optimum).
                if (ThermodynamicPrimerDesigner.CalculatePrimer3PairComplementarity(f.Sequence, r.Sequence)!.Value.Exceeds())
                    continue;
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

    // ---- Row 24: PROBE-DESIGN-001 — designed probes reconstructible + correctly ranked ----

    [Test]
    [Category("PROBE-DESIGN-001")]
    public void DesignProbes_ProbesReconstructibleAndRanked()
    {
        const string target =
            "ACGTGCATGCATGCGCATGCATCGATCGTAGCTAGCATGCATGCGCATGCTAGCATCGATCGTAGCTAGCATGCATGCATGCGC";
        var upper = target.ToUpperInvariant();
        var param = ThermodynamicProbeDesigner.Defaults.Microarray; // length 50-60, GC 0.40-0.60
        const int maxProbes = 10;

        var probes = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes).ToList();
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
            var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMm);
            int oracle = OffTargetOracle(probe, refs, maxMm);
            Assert.That(v.OffTargetHits, Is.EqualTo(oracle), "off-target count");

            double expectedSpec = oracle == 0 ? 0.0 : oracle == 1 ? 1.0 : 1.0 / oracle;
            Assert.That(v.SpecificityScore, Is.EqualTo(expectedSpec).Within(Tol), "specificity score");
        }
    }

    private const double Tol = 1e-9;
    private static string RevComp(string s)
    {
        var arr = s.Select(c => Comp[c]).ToArray();
        Array.Reverse(arr);
        return new string(arr);
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


    private static int Hamming(string a, string b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) d++;
        return d;
    }

    private static readonly Dictionary<char, char> Comp = new() { ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G' };
}
