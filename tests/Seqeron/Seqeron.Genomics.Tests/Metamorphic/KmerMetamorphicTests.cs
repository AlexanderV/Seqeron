using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the K-mer area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-COUNT-001 — k-mer counting (K-mer).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 32.
///
/// API under test (KmerAnalyzer.CountKmers):
///   CountKmers(sequence, k) → Dictionary&lt;kmer,count&gt; over the (|seq| − k + 1) sliding
///   windows of length k (case-insensitive; empty when k &gt; |seq|).
///
/// Relations (derived from the definition, NOT from output):
///   • INV (conservation): the counts sum to the number of length-k windows, |seq| − k + 1.
///   • MON  (occurrence subsumption): every occurrence of a (k+1)-mer w implies an
///          occurrence of its length-k prefix AND its length-k suffix at the same/adjacent
///          position, so count_{k+1}(w) ≤ count_k(prefix) and ≤ count_k(suffix); and the
///          TOTAL instance count strictly decreases by exactly 1 when k grows by 1.
///   • INV (reverse / canonical): reversing the sequence maps each k-mer occurrence to its
///          reverse, so count_s(w) = count_{reverse(s)}(reverse(w)); hence the CANONICAL
///          histogram keyed by min(w, reverse(w)) is identical for s and reverse(s).
///
/// ── Reconciliation with the checklist wording ──
///   The checklist row reads "MON: k+1 → ≤ distinct k-mers". The literal claim that the
///   number of DISTINCT k-mers is non-increasing in k is FALSE in general — e.g. "AACAAC"
///   has 2 distinct 1-mers {A,C} but 3 distinct 2-mers {AA,AC,CA}. The rigorously TRUE
///   monotone facts are the per-k-mer occurrence subsumption (a longer word cannot occur
///   more often than its prefix/suffix) and the strict −1 step in the TOTAL instance count;
///   those are what is asserted here, rather than rubber-stamping a false relation.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class KmerMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20260619);

    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    private static string Reverse(string s) => new(s.Reverse().ToArray());

    /// <summary>Watson–Crick complement WITHOUT reversal — a position-preserving relabel A↔T, C↔G.</summary>
    private static string Complement(string s)
    {
        var chars = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            chars[i] = s[i] switch
            {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => throw new ArgumentException($"non-ACGT base '{s[i]}'", nameof(s)),
            };
        return new string(chars);
    }

    /// <summary>Canonical histogram keyed by min(kmer, reverse(kmer)) — invariant under sequence reversal.</summary>
    private static Dictionary<string, int> CanonicalHistogram(Dictionary<string, int> counts)
    {
        var hist = new Dictionary<string, int>();
        foreach (var (kmer, count) in counts)
        {
            string canon = string.CompareOrdinal(kmer, Reverse(kmer)) <= 0 ? kmer : Reverse(kmer);
            hist[canon] = hist.GetValueOrDefault(canon) + count;
        }
        return hist;
    }

    /// <summary>Uppercase bodies (repetitive, homopolymer, and fixed-seed random).</summary>
    private static IEnumerable<string> KmerBodies()
    {
        yield return "ACGTACGTACGT";
        yield return "AAAAAAAA";
        yield return "AACAACAACAAC";
        yield return "GATTACAGATTACA";
        yield return RandomDna(60);
        yield return RandomDna(100);
    }

    private static readonly int[] KValues = { 1, 2, 3, 4, 6 };

    #endregion

    #region INV — counts sum to the number of length-k windows (|seq| − k + 1)

    [Test]
    [Description("INV: the total k-mer instance count equals |seq| − k + 1 (the number of sliding windows), and 0 when k exceeds the sequence length.")]
    public void CountKmers_TotalInstances_EqualsWindowCount()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                var counts = KmerAnalyzer.CountKmers(body, k);
                int expected = body.Length >= k ? body.Length - k + 1 : 0;

                counts.Values.Sum().Should().Be(expected,
                    because: $"every length-{k} window contributes exactly one k-mer instance, so the counts sum to |seq| − k + 1 = {expected}");
            }
        }
    }

    #endregion

    #region MON — extending k by 1 cannot increase a k-mer's count; total instances drop by 1

    [Test]
    [Description("MON: every (k+1)-mer occurs no more often than its length-k prefix or suffix (occurrence subsumption), and the total instance count strictly decreases by exactly 1 as k grows by 1.")]
    public void CountKmers_ExtendingK_SubsumesOccurrences_AndTotalDropsByOne()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                if (body.Length < k + 1)
                    continue;

                var ck = KmerAnalyzer.CountKmers(body, k);
                var ck1 = KmerAnalyzer.CountKmers(body, k + 1);

                foreach (var (w, count) in ck1)
                {
                    ck[w[..k]].Should().BeGreaterThanOrEqualTo(count,
                        because: $"every occurrence of the (k+1)-mer '{w}' contains an occurrence of its length-{k} prefix, so the prefix occurs at least as often");
                    ck[w[1..]].Should().BeGreaterThanOrEqualTo(count,
                        because: $"every occurrence of '{w}' contains an occurrence of its length-{k} suffix, so the suffix occurs at least as often");
                }

                ck1.Values.Sum().Should().Be(ck.Values.Sum() - 1,
                    because: "there is exactly one fewer length-(k+1) window than length-k window, so the total instance count drops by 1");
            }
        }
    }

    #endregion

    #region INV — reversing the sequence maps each k-mer to its reverse (canonical histogram preserved)

    [Test]
    [Description("INV: reversing the sequence maps each k-mer to its reverse with the same count, so the canonical histogram (keyed by min(kmer, reverse(kmer))) is identical.")]
    public void CountKmers_ReverseSequence_PreservesCanonicalHistogram()
    {
        foreach (var body in KmerBodies())
        {
            string reversed = Reverse(body);

            foreach (int k in KValues)
            {
                var forward = KmerAnalyzer.CountKmers(body, k);
                var rev = KmerAnalyzer.CountKmers(reversed, k);

                rev.Count.Should().Be(forward.Count,
                    because: "the bijection kmer ↔ reverse(kmer) pairs distinct k-mers, so the number of distinct k-mers is preserved under reversal");

                foreach (var (w, count) in forward)
                {
                    rev.GetValueOrDefault(Reverse(w)).Should().Be(count,
                        because: $"every occurrence of '{w}' in s corresponds to an occurrence of '{Reverse(w)}' in reverse(s), so their counts match");
                }

                CanonicalHistogram(rev).Should().BeEquivalentTo(CanonicalHistogram(forward),
                    because: "min(kmer, reverse(kmer)) is invariant under reversal, so the canonical k-mer histogram is the same for s and reverse(s)");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-FREQ-001 — normalized k-mer frequencies (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 33.
    //
    // API under test (KmerAnalyzer.GetKmerFrequencies):
    //   GetKmerFrequencies(sequence, k) = CountKmers(sequence, k) normalized by the total
    //   instance count, i.e. freq(w) = count(w) / (|seq| − k + 1). Case-insensitive; the
    //   empty map when k > |seq| (no windows ⇒ nothing to normalize).
    //
    // Relations (derived from the definition — a probability distribution over k-mers — NOT
    // from output):
    //   • INV (normalization / sum = 1): for any non-empty result the frequencies form a
    //          probability distribution, so they sum to 1 and each lies in (0, 1]; when
    //          k > |seq| the result is empty (sum 0).
    //   • INV (duplicate ⇒ same freqs): frequency is scale-invariant — multiplying every
    //          k-mer count by the same constant leaves count/total unchanged. Self-
    //          concatenation S·S is the boundary-free realization of this at k = 1 (each
    //          mononucleotide count doubles, the total doubles), and likewise any m-fold
    //          repeat S^m. For k ≥ 2 the m−1 junctions inject k−1 extra windows each, so the
    //          relation is exact only at k = 1 — that boundary effect is asserted explicitly
    //          rather than pretending the equality holds verbatim for all k.
    //   • SYM (complement ⇒ relabelled profile): the Watson–Crick complement (A↔T, C↔G,
    //          WITHOUT reversal) is a position-preserving bijection on the sequence, so
    //          window i of comp(S) is comp(window i of S). Hence freq_{comp(S)}(comp(w)) =
    //          freq_S(w) exactly for every k, and the multiset of frequency VALUES is
    //          identical between S and comp(S).
    // ───────────────────────────────────────────────────────────────────────────

    #region INV — frequencies form a probability distribution (sum to 1; empty when k > |seq|)

    [Test]
    [Description("INV: normalized k-mer frequencies sum to 1 and each lies in (0, 1] whenever there is at least one window; the result is empty when k exceeds the sequence length.")]
    public void GetKmerFrequencies_FormProbabilityDistribution()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                var freqs = KmerAnalyzer.GetKmerFrequencies(body, k);

                if (body.Length < k)
                {
                    freqs.Should().BeEmpty(
                        because: $"there are no length-{k} windows in a sequence of length {body.Length}, so there is nothing to normalize");
                    continue;
                }

                freqs.Values.Sum().Should().BeApproximately(1.0, 1e-9,
                    because: "every window contributes one k-mer instance and frequencies divide by the total instance count, so the frequencies partition probability mass 1");
                freqs.Values.Should().OnlyContain(f => f > 0.0 && f <= 1.0,
                    because: "a frequency is a strictly positive count over the total, capped at 1 when a single k-mer fills every window");
            }
        }
    }

    #endregion

    #region INV — frequency is scale-invariant: duplicating the sequence preserves it (exact at k = 1)

    [Test]
    [Description("INV: self-concatenation and m-fold repetition leave the k=1 frequency distribution unchanged, because every base count scales by the same factor and frequency is count/total.")]
    public void GetKmerFrequencies_RepeatingSequence_PreservesMononucleotideFrequencies()
    {
        foreach (var body in KmerBodies())
        {
            var single = KmerAnalyzer.GetKmerFrequencies(body, 1);

            foreach (int m in new[] { 2, 3, 5 })
            {
                string repeated = string.Concat(Enumerable.Repeat(body, m));
                var repeatedFreqs = KmerAnalyzer.GetKmerFrequencies(repeated, 1);

                repeatedFreqs.Keys.Should().BeEquivalentTo(single.Keys,
                    because: "repeating a sequence introduces no new bases, so the support of the mononucleotide distribution is unchanged");
                foreach (var (baseLetter, freq) in single)
                {
                    repeatedFreqs[baseLetter].Should().BeApproximately(freq, 1e-9,
                        because: $"repeating S {m}× multiplies the count of '{baseLetter}' and the total by {m}, so count/total — the frequency — is unchanged");
                }
            }
        }
    }

    [Test]
    [Description("MON/boundary: for k ≥ 2 self-concatenation perturbs frequencies only through the k−1 junction windows, so the total absolute frequency change is bounded by 2·(k−1)/(2|seq|−k+1) — exactly 0 at k = 1.")]
    public void GetKmerFrequencies_SelfConcatenation_DeviationBoundedByJunctionWindows()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                if (body.Length < k)
                    continue;

                var single = KmerAnalyzer.GetKmerFrequencies(body, k);
                var doubled = KmerAnalyzer.GetKmerFrequencies(body + body, k);

                double l1 = single.Keys.Union(doubled.Keys)
                    .Sum(w => Math.Abs(single.GetValueOrDefault(w) - doubled.GetValueOrDefault(w)));

                // The doubled sequence has 2n−k+1 windows; exactly k−1 of them straddle the junction
                // and can carry mass the single copy lacks (and the renormalization touches the rest).
                int doubledWindows = 2 * body.Length - k + 1;
                double bound = 2.0 * (k - 1) / doubledWindows + 1e-12;

                l1.Should().BeLessThanOrEqualTo(bound,
                    because: $"only the {k - 1} junction window(s) differ between S·S and two independent copies of S, so the frequency vectors agree up to that bounded perturbation (exactly when k = 1)");
            }
        }
    }

    #endregion

    #region SYM — complement is a position-preserving relabel: freq_{comp(S)}(comp(w)) = freq_S(w)

    [Test]
    [Description("SYM: the Watson–Crick complement (A↔T, C↔G, no reversal) relabels each k-mer to its complement at the same positions, so freq_{comp(S)}(comp(w)) = freq_S(w) and the multiset of frequency values is identical.")]
    public void GetKmerFrequencies_Complement_RelabelsProfileByComplement()
    {
        foreach (var body in KmerBodies())
        {
            string comp = Complement(body);

            foreach (int k in KValues)
            {
                if (body.Length < k)
                    continue;

                var forward = KmerAnalyzer.GetKmerFrequencies(body, k);
                var complemented = KmerAnalyzer.GetKmerFrequencies(comp, k);

                foreach (var (w, freq) in forward)
                {
                    complemented[Complement(w)].Should().BeApproximately(freq, 1e-9,
                        because: $"complementing is positional, so every occurrence of '{w}' in S is an occurrence of '{Complement(w)}' in comp(S) at the same window — the frequencies must match exactly");
                }

                complemented.Values.OrderBy(f => f).Should().Equal(
                    forward.Values.OrderBy(f => f),
                    because: "k-mer ↔ complement is a bijection, so the complement merely relabels the distribution and the multiset of frequency values is preserved");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-FIND-001 — frequent / recurrent k-mer selection (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 34.
    //
    // APIs under test:
    //   FindKmersWithMinCount(seq, k, minCount) = { (w, count(w)) : count(w) ≥ minCount },
    //     ordered by count DESCENDING (Compeau & Pevzner "Count(Text,Pattern) ≥ t").
    //   FindMostFrequentKmers(seq, k) = the k-mers whose count equals the maximum count
    //     (the top tier; ties all returned). Both empty when k > |seq|.
    //
    // Relations (derived from the definitions, NOT from output):
    //   • MON (threshold antitone ⇒ result monotone): the qualifying set
    //          {w : count(w) ≥ t} shrinks as the threshold t rises, so a LOWER minCount
    //          yields a SUPERSET and a non-decreasing result size.
    //   • SUB (most-frequent ⊆ any lower threshold): the maximal-count k-mers clear every
    //          threshold ≤ maxCount, so FindMostFrequentKmers ⊆ FindKmersWithMinCount(·,t)
    //          for all t ≤ maxCount, and equals exactly the tier at t = maxCount. The
    //          descending-ranked list also nests by prefix — top-1 ⊆ top-5 — with counts
    //          non-increasing along it (the "top-N" reading of the checklist).
    //   • INV (repeat ⇒ unit k-mer dominates): for a tandem repeat U^m with k ≤ |U| the
    //          string has period |U|, so a window's content depends only on its start
    //          residue mod |U|. By pigeonhole the maximum count is ≥ m−1 (it grows with m),
    //          and every most-frequent k-mer is a substring of U·U — i.e. a k-mer of the
    //          repeat unit, not an artefact of the flanks.
    // ───────────────────────────────────────────────────────────────────────────

    #region MON — lower minCount yields a superset (more recurrent k-mers)

    [Test]
    [Description("MON: FindKmersWithMinCount with a lower threshold returns a superset of the higher-threshold result, and the result size is non-increasing as the threshold rises.")]
    public void FindKmersWithMinCount_LowerThreshold_YieldsSuperset()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                if (body.Length < k)
                    continue;

                int maxCount = KmerAnalyzer.CountKmers(body, k).Values.Max();

                for (int t = 1; t < maxCount; t++)
                {
                    var lower = KmerAnalyzer.FindKmersWithMinCount(body, k, t).Select(x => x.Kmer).ToHashSet();
                    var higher = KmerAnalyzer.FindKmersWithMinCount(body, k, t + 1).Select(x => x.Kmer).ToHashSet();

                    higher.Should().BeSubsetOf(lower,
                        because: $"requiring count ≥ {t + 1} is strictly stronger than count ≥ {t}, so the qualifying set can only shrink");
                    lower.Count.Should().BeGreaterThanOrEqualTo(higher.Count,
                        because: "a superset has at least as many elements");
                }
            }
        }
    }

    #endregion

    #region SUB — most-frequent k-mers sit atop every lower threshold; ranked prefixes nest

    [Test]
    [Description("SUB: the most-frequent k-mers equal the tier at minCount = maxCount and are a subset of every lower-threshold result.")]
    public void FindMostFrequentKmers_AreTheTopThresholdTier()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                if (body.Length < k)
                    continue;

                var mostFrequent = KmerAnalyzer.FindMostFrequentKmers(body, k).ToHashSet();
                int maxCount = KmerAnalyzer.CountKmers(body, k).Values.Max();

                var atMax = KmerAnalyzer.FindKmersWithMinCount(body, k, maxCount).Select(x => x.Kmer).ToHashSet();
                mostFrequent.Should().BeEquivalentTo(atMax,
                    because: "the most-frequent k-mers are exactly those whose count reaches the maximum, i.e. those clearing the threshold minCount = maxCount");

                for (int t = 1; t <= maxCount; t++)
                {
                    var tier = KmerAnalyzer.FindKmersWithMinCount(body, k, t).Select(x => x.Kmer).ToHashSet();
                    mostFrequent.Should().BeSubsetOf(tier,
                        because: $"a maximal-count k-mer clears every threshold t = {t} ≤ {maxCount}, so it appears in that tier");
                }
            }
        }
    }

    [Test]
    [Description("SUB: FindKmersWithMinCount is ranked by count descending, so its prefixes nest (top-1 ⊆ top-5) and the counts are non-increasing along the list.")]
    public void FindKmersWithMinCount_RankedPrefixes_Nest()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                if (body.Length < k)
                    continue;

                var ranked = KmerAnalyzer.FindKmersWithMinCount(body, k, 1).ToList();

                for (int i = 1; i < ranked.Count; i++)
                    ranked[i].Count.Should().BeLessThanOrEqualTo(ranked[i - 1].Count,
                        because: "the API documents a count-descending order, so each rank's count is ≤ the previous rank's");

                var top1 = ranked.Take(1).Select(x => x.Kmer).ToHashSet();
                var top5 = ranked.Take(5).Select(x => x.Kmer).ToHashSet();
                top1.Should().BeSubsetOf(top5,
                    because: "the top-1 of a ranked list is a prefix of the top-5, hence a subset");
            }
        }
    }

    #endregion

    #region INV — a tandem repeat makes its unit k-mers dominate (count grows with copies)

    [Test]
    [Description("INV: for a tandem repeat U^m with k ≤ |U|, the maximum k-mer count is ≥ m−1 and grows with m, and every most-frequent k-mer is a substring of U·U (a k-mer of the repeat unit).")]
    public void FindMostFrequentKmers_TandemRepeat_UnitKmerDominates()
    {
        string[] units = { "ACGT", "ATG", "GATTACA", "AAC" };

        foreach (var unit in units)
        {
            string doubledUnit = unit + unit;

            foreach (int k in new[] { 1, 2, 3 })
            {
                if (k > unit.Length)
                    continue;

                int prevMax = 0;
                foreach (int m in new[] { 2, 3, 5, 8 })
                {
                    string repeat = string.Concat(Enumerable.Repeat(unit, m));

                    int maxCount = KmerAnalyzer.CountKmers(repeat, k).Values.Max();
                    maxCount.Should().BeGreaterThanOrEqualTo(m - 1,
                        because: $"U^{m} has period |U|={unit.Length}, so by pigeonhole over {unit.Length} start residues some length-{k} k-mer recurs ≥ m−1 times");
                    maxCount.Should().BeGreaterThanOrEqualTo(prevMax,
                        because: "appending another copy of the unit adds one window per residue class, so the maximum count cannot decrease as m grows");
                    prevMax = maxCount;

                    foreach (var w in KmerAnalyzer.FindMostFrequentKmers(repeat, k))
                        doubledUnit.Should().Contain(w,
                            because: $"every length-{k} window of U^{m} (k ≤ |U|) is a substring of U·U, so a dominating k-mer is a k-mer of the repeat unit — not a flank artefact");
                }
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-ASYNC-001 — asynchronous k-mer counting (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 156.
    //
    // API under test (KmerAnalyzer.CountKmersAsync vs CountKmers):
    //   The async entry point offloads the same deterministic counting to a task.
    //
    // Relations (derived from the counting being a pure function, NOT from output):
    //   • INV  (async = sync): the asynchronous result equals the synchronous result.
    //   • INV  (execution order independent): running many counts concurrently yields the identical
    //          result regardless of completion order — concurrency does not affect the counts.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-ASYNC-001 INV — async equals sync

    [Test]
    [Description("INV: CountKmersAsync delegates to the same deterministic counting, so its result equals CountKmers for every sequence and k.")]
    public void KmerAsync_EqualsSync()
    {
        foreach (var seq in new[] { "ACGTACGTTGGCCAATAC", "AAAAAAA", "GCGCGCGCGC" })
            foreach (int k in new[] { 1, 2, 3, 4 })
            {
                var async = KmerAnalyzer.CountKmersAsync(seq, k).GetAwaiter().GetResult();
                async.Should().BeEquivalentTo(KmerAnalyzer.CountKmers(seq, k),
                    because: $"the async count of '{seq}' (k={k}) must equal the synchronous count");
            }
    }

    #endregion

    #region KMER-ASYNC-001 INV — concurrent execution order does not affect the result

    [Test]
    [Description("INV: counting is a pure function, so launching many async counts concurrently yields the identical result regardless of the order in which they complete.")]
    public void KmerAsync_ConcurrentExecution_OrderIndependent()
    {
        const string seq = "ACGTACGTTGGCCAATACGTACGT";
        const int k = 3;
        var expected = KmerAnalyzer.CountKmers(seq, k);

        var tasks = Enumerable.Range(0, 16).Select(_ => KmerAnalyzer.CountKmersAsync(seq, k)).ToArray();
        System.Threading.Tasks.Task.WaitAll(tasks);

        foreach (var task in tasks)
            task.Result.Should().BeEquivalentTo(expected,
                because: "concurrent, out-of-order completion cannot change a pure counting result");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-BOTH-001 — both-strand k-mer counting (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 157.
    //
    // API under test (KmerAnalyzer.CountKmersBothStrands):
    //   Counts each k-mer on the forward strand plus the reverse-complement strand.
    //
    // Relations (derived from strand-symmetric counting, NOT from output):
    //   • SYM  (reverse-complement invariance): counting both strands of a sequence and of its
    //          reverse complement scans the same two strands, so the result is identical.
    //   • ADD  (additive on concatenation): with k = 1 there are no junction k-mers, so the
    //          both-strand counts of a+b equal the element-wise sum of the parts' counts.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-BOTH-001 SYM — both-strand counts are reverse-complement invariant

    [Test]
    [Description("SYM: counting both strands of a sequence and of its reverse complement scans the same forward+reverse strands, so the count dictionaries are identical.")]
    public void KmerBothStrands_ReverseComplement_Invariant()
    {
        const string seq = "ACGTACGTTGGCCAATAC";
        string rc = Seqeron.Genomics.Core.DnaSequence.GetReverseComplementString(seq);

        foreach (int k in new[] { 1, 2, 3, 4 })
            KmerAnalyzer.CountKmersBothStrands(rc, k).Should().BeEquivalentTo(KmerAnalyzer.CountKmersBothStrands(seq, k),
                because: $"both-strand counting (k={k}) already includes the reverse-complement strand");
    }

    #endregion

    #region KMER-BOTH-001 ADD — k=1 both-strand counts are additive on concatenation

    [Test]
    [Description("ADD: at k=1 no k-mer spans the junction, so the both-strand counts of a+b equal the element-wise sum of the both-strand counts of a and b.")]
    public void KmerBothStrands_Concatenation_AdditiveAtK1()
    {
        foreach (var (a, b) in new[] { ("ACGT", "GGCC"), ("AACC", "GGTTAC"), ("G", "ACGTAC") })
        {
            var ca = KmerAnalyzer.CountKmersBothStrands(a, 1);
            var cb = KmerAnalyzer.CountKmersBothStrands(b, 1);
            var cab = KmerAnalyzer.CountKmersBothStrands(a + b, 1);

            var summed = new Dictionary<string, int>(ca);
            foreach (var kvp in cb)
                summed[kvp.Key] = summed.GetValueOrDefault(kvp.Key) + kvp.Value;

            cab.Should().BeEquivalentTo(summed,
                because: $"single-base counting of '{a}'+'{b}' on both strands is the sum of the parts (no junction 1-mer)");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-DIST-001 — k-mer frequency distance (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 158.
    //
    // API under test (KmerAnalyzer.KmerDistance):
    //   Euclidean distance between the two sequences' k-mer relative-frequency vectors.
    //
    // Relations (derived from the metric definition, NOT from output):
    //   • SYM  (d(a,b)=d(b,a)): the Euclidean distance is symmetric in its arguments.
    //   • INV  (d(x,x)=0): a sequence is at zero distance from itself.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-DIST-001 SYM — the distance is symmetric

    [Test]
    [Description("SYM: the k-mer frequency distance is a Euclidean metric, so d(a,b) equals d(b,a).")]
    public void KmerDistance_Symmetric()
    {
        foreach (var (a, b) in new[] { ("ACGTACGT", "ACGTTGCA"), ("AAAACCCC", "GGGGTTTT"), ("ACGTACGTAC", "TTTTAAAA") })
            foreach (int k in new[] { 1, 2, 3 })
                KmerAnalyzer.KmerDistance(b, a, k).Should().BeApproximately(KmerAnalyzer.KmerDistance(a, b, k), 1e-12,
                    because: $"the Euclidean k-mer distance (k={k}) does not depend on argument order");
    }

    #endregion

    #region KMER-DIST-001 INV — a sequence is at zero distance from itself

    [Test]
    [Description("INV: identical frequency vectors give a Euclidean distance of 0, so d(x,x)=0.")]
    public void KmerDistance_SelfDistance_Zero()
    {
        foreach (var seq in new[] { "ACGTACGT", "AAAAAAAA", "GCGCGCGCGC" })
            foreach (int k in new[] { 1, 2, 3 })
                KmerAnalyzer.KmerDistance(seq, seq, k).Should().Be(0.0,
                    because: $"'{seq}' has identical k-mer frequencies to itself (k={k})");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-GENERATE-001 — exhaustive k-mer generation (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 159.
    //
    // API under test (KmerAnalyzer.GenerateAllKmers):
    //   Enumerates every length-k string over the alphabet (the Cartesian product alphabetᵏ).
    //
    // Relations (derived from the Cartesian-product definition, NOT from output):
    //   • INV  (order independent): the generated SET depends only on the alphabet's members and k,
    //          not on the order of the alphabet letters.
    //   • P    (set closed under all k-mers): the result is complete — exactly |alphabet|ᵏ distinct
    //          k-mers — and closed: it equals { a + s : a ∈ alphabet, s ∈ generate(k-1) }.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-GENERATE-001 INV — the generated set is independent of alphabet order

    [Test]
    [Description("INV: the set of all k-mers depends only on the alphabet's members and k, so permuting the alphabet order yields the same set.")]
    public void GenerateAllKmers_AlphabetOrder_Invariant()
    {
        foreach (int k in new[] { 1, 2, 3 })
            KmerAnalyzer.GenerateAllKmers(k, "TGCA").ToHashSet()
                .Should().BeEquivalentTo(KmerAnalyzer.GenerateAllKmers(k, "ACGT").ToHashSet(),
                    because: $"alphabet order does not change which k-mers exist (k={k})");
    }

    #endregion

    #region KMER-GENERATE-001 P — the set is complete and closed under all k-mers

    [Test]
    [Description("P: the generated set has exactly |alphabet|^k distinct k-mers and is closed — it equals { a + s : a in alphabet, s in generate(k-1) }.")]
    public void GenerateAllKmers_CompleteAndClosed()
    {
        const string alphabet = "ACGT";
        foreach (int k in new[] { 1, 2, 3, 4 })
        {
            var all = KmerAnalyzer.GenerateAllKmers(k, alphabet).ToList();

            all.Should().OnlyHaveUniqueItems(because: "every generated k-mer is distinct");
            all.Count.Should().Be((int)System.Math.Pow(alphabet.Length, k), because: $"there are |alphabet|^{k} k-mers");

            if (k >= 2)
            {
                var closure = (from a in alphabet
                               from s in KmerAnalyzer.GenerateAllKmers(k - 1, alphabet)
                               select a + s).ToHashSet();
                all.ToHashSet().Should().BeEquivalentTo(closure,
                    because: "the k-mer set is the alphabet prepended to every (k-1)-mer — closed under extension");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-POSITIONS-001 — k-mer occurrence positions (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 160.
    //
    // API under test (KmerAnalyzer.FindKmerPositions):
    //   Reports every (overlapping) start index where the k-mer occurs, scanning left to right.
    //
    // Relations (derived from positional matching, NOT from output):
    //   • SHIFT (prepend flank shifts positions): prepending an f-base flank that does not create a
    //          new occurrence shifts every reported position by f.
    //   • INV  (order/scan independent): the result is the complete, deterministic, strictly
    //          ascending set of occurrences — independent of how the scan is run.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-POSITIONS-001 SHIFT — a prepended flank shifts the positions

    [Test]
    [Description("SHIFT: prepending an f-base flank that introduces no new occurrence shifts every reported position by exactly f.")]
    public void FindKmerPositions_PrependFlank_ShiftsPositions()
    {
        const string seq = "ACGTACGTTACG";
        const string kmer = "ACG";
        var original = KmerAnalyzer.FindKmerPositions(seq, kmer).ToList();
        original.Should().NotBeEmpty();

        foreach (var flank in new[] { "TT", "GGGTTT" }) // neither creates an ACG at the junction
        {
            var shifted = KmerAnalyzer.FindKmerPositions(flank + seq, kmer).ToList();
            shifted.Should().Equal(original.Select(p => p + flank.Length),
                because: $"the {flank.Length}-base flank relocates every occurrence by {flank.Length} without adding one");
        }
    }

    #endregion

    #region KMER-POSITIONS-001 INV — the result is the complete, deterministic ascending set

    [Test]
    [Description("INV: the occurrence list is the complete set of overlapping matches in strictly ascending order, and is identical on repeated calls.")]
    public void FindKmerPositions_Deterministic_AscendingComplete()
    {
        const string seq = "AAAA"; // overlapping occurrences of "AA" at 0,1,2
        const string kmer = "AA";

        var positions = KmerAnalyzer.FindKmerPositions(seq, kmer).ToList();
        positions.Should().Equal(new[] { 0, 1, 2 });
        positions.Should().BeInAscendingOrder(because: "the scan reports positions left to right");
        KmerAnalyzer.FindKmerPositions(seq, kmer).Should().Equal(positions, because: "the search is a pure function");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-STATS-001 — k-mer composition statistics (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 161.
    //
    // API under test (KmerAnalyzer.AnalyzeKmers):
    //   Aggregates total/distinct/max/min/average multiplicity and Shannon entropy of the k-mer
    //   counts.
    //
    // Relations (derived from the count aggregation, NOT from output):
    //   • INV  (permutation changes positions not counts): at k = 1 a permutation preserves the
    //          single-base multiset, so every statistic is unchanged while the k-mers' positions move.
    //   • ADD  (counts additive on concatenation): the total window count satisfies
    //          TotalKmers(a+b) = TotalKmers(a) + TotalKmers(b) + (k−1) — the (k−1) junction k-mers.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-STATS-001 INV — k=1 statistics are permutation invariant

    [Test]
    [Description("INV: at k=1 a permutation preserves the single-base multiset, so all k-mer statistics are unchanged (only positions move).")]
    public void KmerStats_Permutation_PreservesCountsAtK1()
    {
        const string seq = "AACGTACGTTGGCCAATAC";
        string permuted = new string(seq.OrderBy(c => c).ToArray());

        KmerAnalyzer.AnalyzeKmers(permuted, 1).Should().Be(KmerAnalyzer.AnalyzeKmers(seq, 1),
            because: "the 1-mer count distribution depends only on the base multiset, which a permutation preserves");
    }

    #endregion

    #region KMER-STATS-001 ADD — total k-mer count is additive (with junction term) on concatenation

    [Test]
    [Description("ADD: concatenation creates k−1 junction k-mers, so TotalKmers(a+b) = TotalKmers(a) + TotalKmers(b) + (k−1).")]
    public void KmerStats_TotalKmers_AdditiveWithJunction()
    {
        foreach (var (a, b) in new[] { ("ACGTAC", "GGTTAC"), ("AAACCC", "GGGTTTAC") })
            foreach (int k in new[] { 1, 2, 3 })
            {
                int total = KmerAnalyzer.AnalyzeKmers(a + b, k).TotalKmers;
                int parts = KmerAnalyzer.AnalyzeKmers(a, k).TotalKmers + KmerAnalyzer.AnalyzeKmers(b, k).TotalKmers;
                total.Should().Be(parts + (k - 1),
                    because: $"concatenating '{a}' and '{b}' adds exactly {k - 1} junction-spanning {k}-mers");
            }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: KMER-UNIQUE-001 — unique (count-1) k-mers (K-mer).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 162.
    //
    // API under test (KmerAnalyzer.FindUniqueKmers):
    //   Returns the k-mers occurring exactly once in the sequence.
    //
    // Relations (derived from the count-1 definition, NOT from output):
    //   • MON  (duplicating a k-mer removes it from the unique set): a unique k-mer that gains a
    //          second occurrence no longer has count 1, so it leaves the set; the rest stay.
    //   • INV  (order independent): at k = 1 the unique-base set depends only on the multiset of
    //          bases, so a permutation yields the same unique set.
    // ───────────────────────────────────────────────────────────────────────────

    #region KMER-UNIQUE-001 MON — duplicating a unique k-mer drops it from the set

    [Test]
    [Description("MON: appending another copy of a currently-unique base raises its count to 2, so it leaves the unique set while every other unique base remains.")]
    public void UniqueKmers_DuplicateUniqueKmer_Removed()
    {
        const string seq = "AACGT"; // 1-mer counts: A=2, C=1, G=1, T=1 ⇒ unique {C,G,T}
        var before = KmerAnalyzer.FindUniqueKmers(seq, 1).ToHashSet();
        before.Should().Contain("C", because: "C occurs exactly once");

        var after = KmerAnalyzer.FindUniqueKmers(seq + "C", 1).ToHashSet(); // C now occurs twice
        after.Should().NotContain("C", because: "duplicating C makes its count 2, so it is no longer unique");
        after.Should().BeEquivalentTo(before.Where(km => km != "C"),
            because: "only the duplicated k-mer leaves the unique set; the others are unaffected at k=1");
    }

    #endregion

    #region KMER-UNIQUE-001 INV — the unique set is permutation invariant at k=1

    [Test]
    [Description("INV: at k=1 the unique-base set depends only on the base multiset, so permuting the sequence yields the same unique set.")]
    public void UniqueKmers_Permutation_Invariant()
    {
        const string seq = "AACGTACGTT";
        string permuted = new string(seq.OrderBy(c => c).ToArray());

        KmerAnalyzer.FindUniqueKmers(permuted, 1).ToHashSet()
            .Should().BeEquivalentTo(KmerAnalyzer.FindUniqueKmers(seq, 1).ToHashSet(),
                because: "the count-1 base set is a function of the base multiset, which a permutation preserves");
    }

    #endregion
}
