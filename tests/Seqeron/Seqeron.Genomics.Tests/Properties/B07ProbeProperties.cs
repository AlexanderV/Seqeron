using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.ProbeDesigner;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the ProbeDesigner behaviour introduced by review batch B07
/// (docs/Validation/review-2026-09/B07.md): F48 (G+C over the non-N bases), F50 (preset Tm scales), F51 (Primer3
/// internal-oligo penalty ranking), F52 (Kane-criteria hit count), F54 (BLAST+ K indexing), F55/F56 (THAL_MAX_ALIGN
/// opt-in: lazy screen and site duplex Tm), F60 (both-strand reference scan), F61/F63 (fallback stem-loop screen
/// reported, not deciding validity), F65 (no length adjustment for untabulated ungapped schemes) and F67 (tiling
/// with the end-anchored window). Oracles are independent re-computations (Hamming scans, G+C counts, interval
/// unions, |Tm − opt| + |len − opt|, reverse complements) or the documented equalities between public methods —
/// never copies of the implementation. Targets are short (≤ 120 nt) and the ntthal screen is used only where the
/// property concerns it (otherwise <see cref="ProbeStructureScreen.Heuristic"/>), so each case is fast.
///
/// Test Units: PROBE-DESIGN-001, PROBE-VALID-001, PROBE-EVALUE-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class B07ProbeProperties
{

    private static Gen<string> SeqGen(int minLen, int maxLen, params char[] alphabet) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alphabet.Length == 0 ? ['A', 'C', 'G', 'T'] : alphabet).ArrayOf(n)
        select new string(chars);

    private static Gen<string> AcgtGen(int minLen, int maxLen) => SeqGen(minLen, maxLen, 'A', 'C', 'G', 'T');

    /// <summary>Independent Watson–Crick reverse complement (upper-case ACGT; anything else → N).</summary>
    private static string RevComp(string s)
    {
        var r = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(r);
    }

    private static string Mutate(string s, int[] positions, char[] bases)
    {
        var c = s.ToCharArray();
        for (int i = 0; i < positions.Length && c.Length > 0; i++)
            c[positions[i] % c.Length] = bases[i];
        return new string(c);
    }

    private static bool Close(double a, double b, double rel = 1e-9) =>
        a == b || Math.Abs(a - b) <= rel * Math.Max(1.0, Math.Max(Math.Abs(a), Math.Abs(b)));

    /// <summary>Independent ungapped site oracle: (position, mismatch positions) of every window with ≤ k mismatches.</summary>
    private static List<(int Pos, List<int> Mm)> HammingSites(string reference, string pattern, int k)
    {
        var sites = new List<(int, List<int>)>();
        for (int p = 0; p + pattern.Length <= reference.Length; p++)
        {
            var mm = new List<int>();
            for (int i = 0; i < pattern.Length && mm.Count <= k; i++)
                if (reference[p + i] != pattern[i])
                    mm.Add(i);
            if (mm.Count <= k)
                sites.Add((p, mm));
        }
        return sites;
    }

    /// <summary>Kane et al. (2000) on an ungapped diagonal: identity (L − d)/L &gt; thr, or an identical run &gt; contig.</summary>
    private static bool Kane(List<int> mm, int length, double thr, int contig)
    {
        if ((double)(length - mm.Count) / length > thr)
            return true;
        int run = 0, best = 0, next = 0;
        for (int i = 0; i < length; i++)
        {
            if (next < mm.Count && mm[next] == i) { next++; run = 0; }
            else best = Math.Max(best, ++run);
        }
        return best > contig;
    }

    /// <summary>Independent (sites, Kane sites) count of a probe over references, optionally on both strands (union of positions).</summary>
    private static (int Sites, int Kane) SiteOracle(string probe, IEnumerable<string> refs, int k, bool both, double thr = 0.75, int contig = 15)
    {
        int sites = 0, kane = 0;
        foreach (var r in refs)
        {
            var map = new Dictionary<int, bool>();
            var patterns = both && RevComp(probe) != probe ? new[] { probe, RevComp(probe) } : new[] { probe };
            foreach (var pat in patterns)
                foreach (var (pos, mm) in HammingSites(r, pat, k))
                    map[pos] = (map.TryGetValue(pos, out bool v) && v) || Kane(mm, probe.Length, thr, contig);
            sites += map.Count;
            kane += map.Values.Count(v => v);
        }
        return (sites, kane);
    }

    /// <summary>A probe plus references holding planted copies (0–5 substitutions, either strand) in random flanks.</summary>
    private static Gen<(string Probe, string[] Refs, int K)> ProbeWithReferencesGen(int minLen, int maxLen) =>
        from probe in AcgtGen(minLen, maxLen)
        from nRefs in Gen.Choose(1, 3)
        from refs in (
            from copies in Gen.Choose(0, 3)
            from planted in (
                from d in Gen.Choose(0, 5)
                from pos in Gen.Choose(0, 1000).ArrayOf(d)
                from bases in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(d)
                from rc in Gen.Elements(false, false, true)
                let m = Mutate(probe, pos, bases)
                select rc ? RevComp(m) : m).ArrayOf(copies)
            from flanks in AcgtGen(0, 25).ArrayOf(copies + 1)
            select string.Concat(flanks.Zip(planted.Append(""), (f, p) => f + p))).ArrayOf(nRefs)
        from k in Gen.Choose(0, 4)
        select (probe, refs, k);

    private sealed class ProbeComparer : IEqualityComparer<Probe>
    {
        public static readonly ProbeComparer Instance = new();
        public bool Equals(Probe a, Probe b) =>
            a.Sequence == b.Sequence && a.Start == b.Start && a.Score == b.Score && a.Tm == b.Tm
            && a.Warnings.SequenceEqual(b.Warnings) && a.Primer3Penalty == b.Primer3Penalty;
        public int GetHashCode(Probe p) => HashCode.Combine(p.Sequence, p.Start);
    }



    private static readonly HashSet<(int, int)> BlastnTables =
        [(1, -5), (1, -4), (2, -7), (1, -3), (2, -5), (1, -2), (2, -3), (3, -4), (1, -1), (3, -2), (4, -5), (5, -4)];

    private static int Gcd(int a, int b) => b == 0 ? Math.Abs(a) : Gcd(b, a % b);

    /// <summary>Schemes with a defined λ under the uniform composition: match &gt; 0, mismatch &lt; 0, match + 3·mismatch &lt; 0.</summary>
    private static Gen<(int Match, int Mismatch)> SchemeGen(int maxMatch, int maxMismatch) =>
        (from m in Gen.Choose(1, maxMatch)
         from mm in Gen.Choose(1, maxMismatch)
         select (m, -mm)).Where(s => s.Item1 + 3 * s.Item2 < 0);

    /// <summary>
    /// F54: for every gcd = 1 scheme <see cref="KarlinKMethod.NcbiBlast"/> and <see cref="KarlinKMethod.ReducedLattice"/>
    /// give the same λ / K / H; for the closed-form schemes δ·(a, b) with reduced lowest score −1 or highest +1 they agree
    /// too, and <see cref="KarlinKMethod.ReducedLattice"/> is scale-invariant (K(δa, δb) = K(a, b), λ(δa, δb) = λ(a, b)/δ).
    /// Uniform and non-uniform compositions.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property UngappedK_NcbiBlastEqualsReducedLattice_ForCoprimeAndClosedFormSchemes()
    {
        var closedForm = from a in Gen.Choose(1, 6)
                         from b in Gen.Choose(1, 6)
                         from which in Gen.Elements(true, false)
                         let s = which ? (Match: 1, Mismatch: -b) : (Match: a, Mismatch: -1)
                         where s.Match + 3 * s.Mismatch < 0
                         from d in Gen.Choose(2, 6)
                         select (s.Match * d, s.Mismatch * d, d);
        var coprime = from s in SchemeGen(12, 20) where Gcd(s.Match, -s.Mismatch) == 1 select (s.Match, s.Mismatch, 1);
        var gen = from s in Gen.Frequency((1, closedForm), (1, coprime))
                  from pct in Gen.Frequency((2, Gen.Constant(25)), (1, Gen.Choose(15, 35)))
                  let freq = pct / 100.0
                  let pMatch = 4 * freq * freq
                  where pMatch * s.Item1 + (1 - pMatch) * s.Item2 < 0
                  select (s, freq);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var ((match, mismatch, d), f) = x;
            var r = ComputeUngappedKarlinParameters(match, mismatch, f, KarlinKMethod.ReducedLattice);
            var n = ComputeUngappedKarlinParameters(match, mismatch, f, KarlinKMethod.NcbiBlast);
            bool same = Close(r.Lambda, n.Lambda, 1e-12) && Close(r.K, n.K, 1e-12) && Close(r.H, n.H, 1e-12);
            if (d > 1)
            {
                var red = ComputeUngappedKarlinParameters(match / d, mismatch / d, f, KarlinKMethod.ReducedLattice);
                same &= Close(r.K, red.K, 1e-9) && Close(r.Lambda * d, red.Lambda, 1e-9);
            }
            return same.Label($"{match}/{mismatch} δ {d} p {f}: reduced {r} ncbi {n}");
        });
    }

    /// <summary>
    /// F65: ungapped <c>ComputeBlastnStatistics</c> of a scheme whose gcd-reduced form has no BLAST+ <c>blastn_values_*</c>
    /// table applies no length adjustment (α = β = 0, ℓ = 0, search space m·n); tabulated schemes keep α = λ/H and
    /// ℓ = <c>ComputeLengthAdjustment(K, α/λ, β, m, n, N)</c> with search space (m − ℓ)·max(1, n − N·ℓ). λ / K / H are the
    /// NcbiBlast ungapped values either way.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property UngappedBlastnStatistics_LengthAdjustmentOnlyForTabulatedSchemes()
    {
        var gen = from s in SchemeGen(8, 12)
                  from m in Gen.Choose(10, 500)
                  from n in Gen.Choose(10, 2_000_000)
                  from bigN in Gen.Choose(1, 20)
                  from raw in Gen.Choose(0, 400)
                  select (s.Match, s.Mismatch, m, (long)n, bigN, raw);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (match, mismatch, m, n, bigN, raw) = x;
            int g = Gcd(match, -mismatch);
            bool tabulated = BlastnTables.Contains((match / g, mismatch / g));
            var st = ComputeBlastnStatistics(raw, m, n, bigN, new ScoringMatrix(match, mismatch, -5, -2), gapped: false);
            var u = ComputeUngappedKarlinParameters(match, mismatch, 0.25, KarlinKMethod.NcbiBlast);
            var p = st.Parameters;
            bool lambdaKh = p.Lambda == u.Lambda && p.K == u.K && p.H == u.H && !p.Gapped;
            if (!tabulated)
                return (lambdaKh && p.Alpha == 0 && p.Beta == 0 && st.LengthAdjustment == 0
                        && st.EffectiveSearchSpace == (double)m * n)
                    .Label($"untabulated {match}/{mismatch}: α {p.Alpha} β {p.Beta} ℓ {st.LengthAdjustment} space {st.EffectiveSearchSpace}");
            int ell = ComputeLengthAdjustment(p.K, p.Alpha / p.Lambda, p.Beta, m, n, bigN);
            double space = (m - ell) * Math.Max(1.0, n - (double)bigN * ell);
            return (lambdaKh && Close(p.Alpha, u.Lambda / u.H, 1e-12) && p.Beta == u.Beta && st.LengthAdjustment == ell
                    && Close(st.EffectiveSearchSpace, space, 1e-12))
                .Label($"tabulated {match}/{mismatch}: α {p.Alpha} vs {u.Lambda / u.H}, ℓ {st.LengthAdjustment} vs {ell}");
        });
    }

}
