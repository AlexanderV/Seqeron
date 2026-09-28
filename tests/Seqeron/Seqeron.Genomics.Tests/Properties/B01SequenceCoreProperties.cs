using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the behaviours introduced by review batch B01
/// (docs/Validation/review-2026-09/B01.md, fixes F2–F13): IUPAC validation mode and
/// IndexOfInvalid*, linguistic complexity suffix-tree path, Shannon entropy over U,
/// GC-skew / GC-analysis guards and RNA handling, Phred rounding, quality-string contract,
/// IUPAC code algebra and k-mer enumeration.
///
/// Test Units: SEQ-VALID-001, SEQ-COMPLEX-001, SEQ-ENTROPY-001, SEQ-GCSKEW-001,
/// SEQ-GC-ANALYSIS-001, SEQ-REPLICATION-001, B01-SWEEP.
/// Only mathematically exact invariants are asserted.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Composition")]
public class B01SequenceCoreProperties
{
    private const string IupacDnaCodes = "ACGTRYSWKMBDHVN";
    private const string IupacRnaCodes = "ACGURYSWKMBDHVN";

    // Mixed alphabet: strict bases, IUPAC codes, U, both cases, gaps, junk and the non-ASCII 'ſ' (U+017F).
    private static readonly char[] MixedAlphabet =
        "ACGTUacgturyswkmbdhvnRYSWKMBDHVN-. XZxz1@ſKı".ToCharArray();

    private static Arbitrary<string> StringOver(char[] alphabet, int minLen = 0, int maxLen = 60) =>
        (from n in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements(alphabet).ArrayOf(n)
         select new string(chars)).ToArbitrary();

    private static Arbitrary<string> Acgt(int minLen, int maxLen = 80) => StringOver("ACGT".ToCharArray(), minLen, maxLen);

    private static string SwapCaseAscii(string s) =>
        new(s.Select(c => c is >= 'a' and <= 'z' ? (char)(c - 32) : c is >= 'A' and <= 'Z' ? (char)(c + 32) : c).ToArray());

    #region SEQ-VALID-001 (F2/F3): IUPAC validation mode and IndexOfInvalid*

    /// <summary>IndexOfInvalidDna/Rna == −1 ⇔ IsValidDna/Rna; otherwise it points at the FIRST offending char.</summary>
    [FsCheck.NUnit.Property]
    public Property IndexOfInvalid_AgreesWithIsValid_AndIsFirstOffender()
    {
        return Prop.ForAll(StringOver(MixedAlphabet), s =>
        {
            int iDna = s.AsSpan().IndexOfInvalidDna();
            int iRna = s.AsSpan().IndexOfInvalidRna();
            static bool Dna(char c) => "ACGTacgt".Contains(c);
            static bool Rna(char c) => "ACGUacgu".Contains(c);
            int refDna = Array.FindIndex(s.ToCharArray(), c => !Dna(c));
            int refRna = Array.FindIndex(s.ToCharArray(), c => !Rna(c));
            bool ok = (iDna == -1) == s.AsSpan().IsValidDna()
                      && (iRna == -1) == s.AsSpan().IsValidRna()
                      && iDna == refDna && iRna == refRna;
            return ok.Label($"'{s}': IndexOfInvalidDna={iDna} (ref {refDna}), IndexOfInvalidRna={iRna} (ref {refRna})");
        });
    }

    /// <summary>Strict ⊆ IUPAC: IsValidDna ⇒ IsValidIupacDna and IsValidRna ⇒ IsValidIupacRna.</summary>
    [FsCheck.NUnit.Property]
    public Property StrictValid_ImpliesIupacValid()
    {
        return Prop.ForAll(StringOver(MixedAlphabet), s =>
        {
            var span = s.AsSpan();
            bool ok = (!span.IsValidDna() || span.IsValidIupacDna())
                      && (!span.IsValidRna() || span.IsValidIupacRna());
            return ok.Label($"strict-valid '{s}' must be IUPAC-valid");
        });
    }

    /// <summary>IUPAC validity equals an independent set-membership reference (ASCII case-folding only).</summary>
    [FsCheck.NUnit.Property]
    public Property IupacValidity_MatchesAlphabetReference()
    {
        return Prop.ForAll(StringOver(MixedAlphabet), s =>
        {
            bool refDna = s.All(c => IupacDnaCodes.Contains(c) || (c is >= 'a' and <= 'z' && IupacDnaCodes.Contains((char)(c - 32))));
            bool refRna = s.All(c => IupacRnaCodes.Contains(c) || (c is >= 'a' and <= 'z' && IupacRnaCodes.Contains((char)(c - 32))));
            bool ok = s.AsSpan().IsValidIupacDna() == refDna && s.AsSpan().IsValidIupacRna() == refRna;
            return ok.Label($"'{s}': IUPAC DNA {s.AsSpan().IsValidIupacDna()} (ref {refDna}), RNA {s.AsSpan().IsValidIupacRna()} (ref {refRna})");
        });
    }

    /// <summary>All four validators are invariant under ASCII case swapping.</summary>
    [FsCheck.NUnit.Property]
    public Property Validators_AreAsciiCaseInsensitive()
    {
        return Prop.ForAll(StringOver(MixedAlphabet), s =>
        {
            string t = SwapCaseAscii(s);
            bool ok = s.AsSpan().IsValidDna() == t.AsSpan().IsValidDna()
                      && s.AsSpan().IsValidRna() == t.AsSpan().IsValidRna()
                      && s.AsSpan().IsValidIupacDna() == t.AsSpan().IsValidIupacDna()
                      && s.AsSpan().IsValidIupacRna() == t.AsSpan().IsValidIupacRna()
                      && s.AsSpan().IndexOfInvalidDna() == t.AsSpan().IndexOfInvalidDna()
                      && s.AsSpan().IndexOfInvalidRna() == t.AsSpan().IndexOfInvalidRna();
            return ok.Label($"case swap changed a validity verdict: '{s}' vs '{t}'");
        });
    }

    /// <summary>T↔U swap maps IUPAC-DNA validity onto IUPAC-RNA validity (the alphabets differ only in T/U).</summary>
    [FsCheck.NUnit.Property]
    public Property IupacDnaValidity_EqualsIupacRnaValidity_UnderTUSwap()
    {
        return Prop.ForAll(StringOver(MixedAlphabet), s =>
        {
            string swapped = new(s.Select(c => c switch
            {
                'T' => 'U', 'U' => 'T', 't' => 'u', 'u' => 't', _ => c
            }).ToArray());
            return (s.AsSpan().IsValidIupacDna() == swapped.AsSpan().IsValidIupacRna())
                .Label($"IsValidIupacDna('{s}') != IsValidIupacRna('{swapped}')");
        });
    }

    /// <summary>Any string containing a non-ASCII char (incl. 'ſ' U+017F, 'K' U+212A, 'ı' U+0131) is invalid in every mode.</summary>
    [FsCheck.NUnit.Property]
    public Property NonAsciiChar_AlwaysInvalid()
    {
        var gen = from prefix in StringOver("ACGTSacgts".ToCharArray(), 0, 20).Generator
                  from bad in Gen.Elements('ſ', 'K', 'ı', 'Å', 'Ａ')
                  from suffix in StringOver("ACGTSacgts".ToCharArray(), 0, 20).Generator
                  select prefix + bad + suffix;
        return Prop.ForAll(gen.ToArbitrary(), s =>
        {
            var span = s.AsSpan();
            int badIndex = s.ToList().FindIndex(c => c > 127);
            bool ok = !span.IsValidDna() && !span.IsValidRna() && !span.IsValidIupacDna() && !span.IsValidIupacRna()
                      && span.IndexOfInvalidDna() <= badIndex && span.IndexOfInvalidRna() <= badIndex
                      && !new IupacDnaSequence(s).IsValid();
            return ok.Label($"non-ASCII input '{s}' was accepted");
        });
    }

    /// <summary>IsValidIupacDna(s) ⇒ IupacDnaSequence(s).IsValid() (container alphabet ⊇ the 15 codes, both cases).</summary>
    [FsCheck.NUnit.Property]
    public Property IupacValid_ImpliesIupacDnaSequenceValid()
    {
        return Prop.ForAll(StringOver((IupacDnaCodes + IupacDnaCodes.ToLowerInvariant()).ToCharArray()), s =>
            (s.AsSpan().IsValidIupacDna() && new IupacDnaSequence(s).IsValid())
                .Label($"IUPAC-valid '{s}' rejected by IupacDnaSequence"));
    }

    #endregion

    #region SEQ-COMPLEX-001 (F4): linguistic complexity — suffix-tree path

    /// <summary>Independent reference: Σ_{i≤m} |distinct i-substrings| / Σ_{i≤m} min(4^i, N−i+1).</summary>
    private static double ReferenceLc(string s, int maxWordLength)
    {
        int n = s.Length;
        int m = Math.Min(maxWordLength, n);
        long obs = 0, pos = 0;
        for (int i = 1; i <= m; i++)
        {
            var set = new HashSet<string>();
            for (int j = 0; j + i <= n; j++) set.Add(s.Substring(j, i));
            obs += set.Count;
            pos += Math.Min(i < 31 ? 1L << (2 * i) : long.MaxValue, n - i + 1);
        }
        return (double)obs / pos;
    }

    /// <summary>For m &gt; 12 (suffix-tree path), both overloads equal the direct hash-set reference exactly.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property LinguisticComplexity_SuffixTreePath_EqualsDirectReference()
    {
        var gen = from s in Acgt(14, 90).Generator
                  from m in Gen.Choose(13, 120)
                  select (s, m);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double expected = ReferenceLc(t.s, t.m);
            double viaString = SequenceComplexity.CalculateLinguisticComplexity(t.s, t.m);
            double viaDna = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(t.s), t.m);
            return (viaString == expected && viaDna == expected)
                .Label($"LC(m={t.m}) string={viaString:R}, dna={viaDna:R}, reference={expected:R} for '{t.s}'");
        });
    }

    /// <summary>For m ≤ 12 (hash path) the result also equals the reference, so both paths share one definition.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property LinguisticComplexity_HashPath_EqualsDirectReference()
    {
        var gen = from s in Acgt(1, 60).Generator
                  from m in Gen.Choose(1, 12)
                  select (s, m);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double expected = ReferenceLc(t.s, t.m);
            double actual = SequenceComplexity.CalculateLinguisticComplexity(t.s, t.m);
            return (actual == expected).Label($"LC(m={t.m})={actual:R} vs reference {expected:R} for '{t.s}'");
        });
    }

    /// <summary>0 &lt; LC ≤ 1 on non-empty ACGT input for any m ≥ 1 (V_1 ≥ 1; V_i ≤ min(4^i, N−i+1)).</summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property LinguisticComplexity_StrictlyPositive_AtMostOne()
    {
        var gen = from s in Acgt(1, 80).Generator
                  from m in Gen.Elements(1, 2, 5, 10, 12, 13, 20, 50, int.MaxValue)
                  select (s, m);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double lc = SequenceComplexity.CalculateLinguisticComplexity(t.s, t.m);
            return (lc > 0 && lc <= 1.0).Label($"LC={lc:R} (m={t.m}) must lie in (0,1]");
        });
    }

    /// <summary>
    /// LC is invariant under string reversal (distinct substrings of rev(s) are the reverses of those of s,
    /// so every V_i is unchanged) and under reverse complement (reversal ∘ bijective ACGT relabeling).
    /// Checked on both paths (m = 8 hash, m = N suffix tree).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property LinguisticComplexity_ReverseAndReverseComplement_Invariant()
    {
        return Prop.ForAll(Acgt(1, 80), s =>
        {
            string rev = new(s.Reverse().ToArray());
            string rc = new(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToArray());
            bool ok = true;
            foreach (int m in new[] { 8, s.Length })
            {
                double lc = SequenceComplexity.CalculateLinguisticComplexity(s, m);
                ok &= SequenceComplexity.CalculateLinguisticComplexity(rev, m) == lc
                      && SequenceComplexity.CalculateLinguisticComplexity(rc, m) == lc;
            }
            return ok.Label($"LC changed under reverse / reverse-complement for '{s}'");
        });
    }

    #endregion

    #region SEQ-ENTROPY-001 (F5): Shannon entropy over {A,C,G,T/U}

    /// <summary>H(RNA) == H(DNA with U→T) exactly, and 0 ≤ H ≤ 2.</summary>
    [FsCheck.NUnit.Property]
    public Property ShannonEntropy_RnaEqualsDnaTranscript_AndBounded()
    {
        return Prop.ForAll(StringOver("ACGUacgu".ToCharArray(), 0, 80), rna =>
        {
            string dna = rna.Replace('U', 'T').Replace('u', 't');
            double hr = SequenceComplexity.CalculateShannonEntropy(rna);
            double hd = SequenceComplexity.CalculateShannonEntropy(dna);
            return (hr == hd && hr >= 0 && hr <= 2.0 + 1e-12)
                .Label($"H('{rna}')={hr:R} vs H('{dna}')={hd:R}");
        });
    }

    /// <summary>Entropy ignores non-nucleotide symbols: H(s) == H(s restricted to A/C/G/T/U) for mixed input.</summary>
    [FsCheck.NUnit.Property]
    public Property ShannonEntropy_IgnoresNonNucleotideSymbols()
    {
        return Prop.ForAll(StringOver(MixedAlphabet), s =>
        {
            string kept = new(s.Where(c => "ACGTUacgtu".Contains(c)).ToArray());
            double h = SequenceComplexity.CalculateShannonEntropy(s);
            double hk = SequenceComplexity.CalculateShannonEntropy(kept);
            return (h == hk && h >= 0 && h <= 2.0 + 1e-12).Label($"H('{s}')={h:R} vs H(kept)={hk:R}");
        });
    }

    #endregion

    #region SEQ-GCSKEW-001 / SEQ-REPLICATION-001 / SEQ-GC-ANALYSIS-001 (F6–F8)

    /// <summary>String overloads reject window/step &lt; 1 eagerly (no enumeration needed), for any input incl. null/empty.</summary>
    [FsCheck.NUnit.Property(MaxTest = 50)]
    public Property GcSkew_StringOverloads_RejectSizesBelowOne()
    {
        var gen = from s in Gen.Elements<string?>(null, "", "ACGT", "GGGCCCATAT")
                  from bad in Gen.Choose(-5, 0)
                  select (s, bad);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            bool Throws(Action a)
            {
                try { a(); return false; }
                catch (ArgumentOutOfRangeException) { return true; }
            }
            bool ok = Throws(() => GcSkewCalculator.CalculateWindowedGcSkew(t.s!, t.bad, 1))
                      && Throws(() => GcSkewCalculator.CalculateWindowedGcSkew(t.s!, 1, t.bad))
                      && Throws(() => GcSkewCalculator.CalculateCumulativeGcSkew(t.s!, t.bad))
                      && Throws(() => GcSkewCalculator.AnalyzeGcContent(t.s!, t.bad, 1))
                      && Throws(() => GcSkewCalculator.AnalyzeGcContent(t.s!, 1, t.bad));
            return ok.Label($"size {t.bad} not rejected eagerly for '{t.s ?? "null"}'");
        });
    }

    /// <summary>
    /// Cumulative skew is the running sum of the non-overlapping windowed skews: point-for-point equal
    /// GcSkew/Position to CalculateWindowedGcSkew(s, w, w), and the last cumulative value equals their sum.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CumulativeGcSkew_IsRunningSumOfTiledWindows()
    {
        var gen = from s in StringOver("ACGTacgtN".ToCharArray(), 0, 120).Generator
                  from w in Gen.Choose(1, 15)
                  select (s, w);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var windowed = GcSkewCalculator.CalculateWindowedGcSkew(t.s, t.w, t.w).ToList();
            var cumulative = GcSkewCalculator.CalculateCumulativeGcSkew(t.s, t.w).ToList();
            bool ok = windowed.Count == cumulative.Count && windowed.Count == t.s.Length / t.w;
            double running = 0;
            for (int i = 0; ok && i < windowed.Count; i++)
            {
                running += windowed[i].GcSkew;
                ok = cumulative[i].GcSkew == windowed[i].GcSkew
                     && cumulative[i].Position == windowed[i].Position
                     && cumulative[i].CumulativeGcSkew == running;
            }
            return ok.Label($"cumulative != running sum of windows (w={t.w}, len={t.s.Length})");
        });
    }

    /// <summary>
    /// Replication origin/terminus = first argmin/argmax of the per-base cumulative skew
    /// CalculateCumulativeGcSkew(s, 1), as prefix indices with Skew_0 = 0 (index 0 if never below/above 0).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReplicationOrigin_EqualsFirstExtremumOfPerBaseCumulativeSkew()
    {
        return Prop.ForAll(StringOver("ACGTacgt".ToCharArray(), 1, 150), s =>
        {
            var cum = GcSkewCalculator.CalculateCumulativeGcSkew(s, 1).Select(p => p.CumulativeGcSkew).ToList();
            double min = 0, max = 0;
            int minPos = 0, maxPos = 0;
            for (int i = 0; i < cum.Count; i++)
            {
                if (cum[i] < min) { min = cum[i]; minPos = i + 1; }
                if (cum[i] > max) { max = cum[i]; maxPos = i + 1; }
            }
            var pred = GcSkewCalculator.PredictReplicationOrigin(s);
            bool ok = pred.PredictedOrigin == minPos && pred.PredictedTerminus == maxPos
                      && pred.OriginSkew == min && pred.TerminusSkew == max;
            return ok.Label($"origin={pred.PredictedOrigin} (ref {minPos}), terminus={pred.PredictedTerminus} (ref {maxPos})");
        });
    }

    /// <summary>
    /// AnalyzeGcContent counts RNA U in the GC denominator: an RNA string and its U→T DNA version have identical
    /// overall and per-window GC content, GC skew and GC-content variance; overall GC% == 100·CalculateGcFraction.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeGcContent_RnaEqualsDnaTranscript_ForGcMetrics()
    {
        var gen = from s in StringOver("ACGUacguN".ToCharArray(), 1, 120).Generator
                  from w in Gen.Choose(1, 12)
                  from step in Gen.Choose(1, 12)
                  select (s, w, step);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string dna = t.s.Replace('U', 'T').Replace('u', 't');
            var r = GcSkewCalculator.AnalyzeGcContent(t.s, t.w, t.step);
            var d = GcSkewCalculator.AnalyzeGcContent(dna, t.w, t.step);
            bool ok = r.OverallGcContent == d.OverallGcContent
                      && r.OverallGcSkew == d.OverallGcSkew
                      && r.GcContentVariance == d.GcContentVariance
                      && r.WindowedGcContent.Select(x => x.GcContent).SequenceEqual(d.WindowedGcContent.Select(x => x.GcContent))
                      && r.WindowedGcSkew.Select(x => x.GcSkew).SequenceEqual(d.WindowedGcSkew.Select(x => x.GcSkew))
                      && Math.Abs(r.OverallGcContent - 100.0 * t.s.AsSpan().CalculateGcFraction()) < 1e-9;
            return ok.Label($"RNA '{t.s}' GC%={r.OverallGcContent:R} vs DNA {d.OverallGcContent:R}");
        });
    }

    #endregion

    #region B01-SWEEP (F9–F13): Phred, quality strings, IUPAC code algebra, k-mers

    /// <summary>ErrorProbabilityToPhred(PhredToErrorProbability(q)) == q for every Sanger score q ∈ [0, 93].</summary>
    [Test]
    public void Phred_RoundTrip_AllSangerScores()
    {
        for (int q = 0; q <= QualitySequence.MaxSangerPhred; q++)
        {
            double p = QualitySequence.PhredToErrorProbability((byte)q);
            QualitySequence.ErrorProbabilityToPhred(p).Should().Be((byte)q, $"round(−10·log10(10^(−{q}/10))) = {q}");
        }
    }

    /// <summary>Phred is monotone non-increasing in p on [0,1], bounded by [0, 93]; p outside [0,1] or NaN throws.</summary>
    [FsCheck.NUnit.Property]
    public Property Phred_MonotoneNonIncreasing_AndBounded()
    {
        var gen = from a in Gen.Choose(0, 1_000_000)
                  from b in Gen.Choose(0, 1_000_000)
                  from e in Gen.Choose(0, 12)
                  select (Math.Min(a, b) / 1_000_000.0 * Math.Pow(10, -e), Math.Max(a, b) / 1_000_000.0 * Math.Pow(10, -e));
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            byte q1 = QualitySequence.ErrorProbabilityToPhred(t.Item1);
            byte q2 = QualitySequence.ErrorProbabilityToPhred(t.Item2);
            return (q1 >= q2 && q1 <= QualitySequence.MaxSangerPhred)
                .Label($"p1={t.Item1:R} → Q{q1}, p2={t.Item2:R} → Q{q2}");
        });
    }

    [TestCase(-1e-12)]
    [TestCase(1.0000001)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    public void Phred_OutOfDomain_Throws(double p)
    {
        var act = () => QualitySequence.ErrorProbabilityToPhred(p);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    /// <summary>A valid quality string (chars in [offset, 126]) round-trips through GetQualityString for offsets 33 and 64.</summary>
    [FsCheck.NUnit.Property]
    public Property QualityString_RoundTrip()
    {
        var gen = from offset in Gen.Elements(33, 64)
                  from n in Gen.Choose(0, 50)
                  from qs in Gen.Choose(offset, 126).Select(i => (char)i).ArrayOf(n)
                  select (offset, q: new string(qs));
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var seq = new QualitySequence(new string('A', t.q.Length), t.q, t.offset);
            bool ok = seq.GetQualityString(t.offset) == t.q
                      && seq.Qualities.Select((b, i) => b == t.q[i] - t.offset).All(x => x);
            return ok.Label($"quality round-trip failed (offset {t.offset}): '{t.q}'");
        });
    }

    /// <summary>A quality string of the wrong length, or containing a char outside [offset, 126], throws ArgumentException.</summary>
    [FsCheck.NUnit.Property]
    public Property QualityString_WrongLengthOrOutOfRangeChar_Throws()
    {
        var gen = from n in Gen.Choose(1, 30)
                  from delta in Gen.Elements(-1, 1, 2, -n)
                  from badChar in Gen.Elements('\u001F', ' ', '\u007F', 'ÿ', 'ſ')
                  from pos in Gen.Choose(0, n - 1)
                  select (n, delta, badChar, pos);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string bases = new('C', t.n);
            string wrongLength = new('I', t.n + t.delta);
            char[] q = new string('I', t.n).ToCharArray();
            q[t.pos] = t.badChar;
            bool Throws(Action a)
            {
                try { a(); return false; }
                catch (ArgumentOutOfRangeException) { return false; }
                catch (ArgumentException) { return true; }
            }
            bool ok = Throws(() => _ = new QualitySequence(bases, wrongLength))
                      && Throws(() => _ = new QualitySequence(bases, new string(q)));
            return ok.Label($"n={t.n}, delta={t.delta}, bad='{t.badChar}'@{t.pos} not rejected");
        });
    }

    /// <summary>GetIupacCode(ExpandCode(c)) == c for all 15 codes (either case); U expands to T and codes as T.</summary>
    [Test]
    public void IupacCode_ExpandThenEncode_IsIdentity()
    {
        foreach (char c in IupacDnaCodes)
        {
            IupacDnaSequence.GetIupacCode(IupacDnaSequence.ExpandCode(c)).Should().Be(c);
            IupacDnaSequence.GetIupacCode(IupacDnaSequence.ExpandCode(char.ToLowerInvariant(c))).Should().Be(c);
        }
        IupacDnaSequence.GetIupacCode(IupacDnaSequence.ExpandCode('U')).Should().Be('T');
    }

    /// <summary>GetIupacCode treats U as T and is ASCII case-insensitive for any non-empty base set.</summary>
    [FsCheck.NUnit.Property]
    public Property GetIupacCode_UAsT_AndCaseInsensitive()
    {
        return Prop.ForAll(StringOver("ACGTU".ToCharArray(), 1, 8), s =>
        {
            string asT = s.Replace('U', 'T');
            char code = IupacDnaSequence.GetIupacCode(s);
            bool ok = code == IupacDnaSequence.GetIupacCode(asT)
                      && code == IupacDnaSequence.GetIupacCode(s.ToLowerInvariant())
                      && IupacDnaSequence.ExpandCode(code).OrderBy(c => c).SequenceEqual(asT.Distinct().OrderBy(c => c));
            return ok.Label($"GetIupacCode('{s}')={code}");
        });
    }

    /// <summary>
    /// CodesMatch is symmetric, ASCII case-insensitive, and equals "expansion sets intersect" over all pairs
    /// of the 15 codes plus U.
    /// </summary>
    [Test]
    public void CodesMatch_Symmetric_CaseInsensitive_EqualsSetIntersection()
    {
        string codes = IupacDnaCodes + "U";
        foreach (char a in codes)
        foreach (char b in codes)
        {
            bool expected = IupacDnaSequence.ExpandCode(a).Intersect(IupacDnaSequence.ExpandCode(b)).Any();
            bool m = IupacDnaSequence.CodesMatch(a, b);
            m.Should().Be(expected, $"CodesMatch({a},{b})");
            IupacDnaSequence.CodesMatch(b, a).Should().Be(m, $"symmetry ({a},{b})");
            IupacDnaSequence.CodesMatch(char.ToLowerInvariant(a), b).Should().Be(m, $"case ({a},{b})");
            IupacDnaSequence.CodesMatch(a, char.ToLowerInvariant(b)).Should().Be(m, $"case ({a},{b})");
            IupacDnaSequence.CodesMatch(char.ToLowerInvariant(a), char.ToLowerInvariant(b)).Should().Be(m, $"case ({a},{b})");
        }
    }

    /// <summary>CodesMatch is symmetric for arbitrary characters (incl. non-codes and non-ASCII).</summary>
    [FsCheck.NUnit.Property]
    public Property CodesMatch_Symmetric_ForArbitraryChars()
    {
        var gen = from a in Gen.Elements(MixedAlphabet)
                  from b in Gen.Elements(MixedAlphabet)
                  select (a, b);
        return Prop.ForAll(gen.ToArbitrary(), t =>
            (IupacDnaSequence.CodesMatch(t.a, t.b) == IupacDnaSequence.CodesMatch(t.b, t.a))
                .Label($"CodesMatch not symmetric on ('{t.a}','{t.b}')"));
    }

    /// <summary>
    /// EnumerateKmers yields exactly max(0, L−k+1) windows, the i-th equal to s[i..i+k), and their upper-cased
    /// multiset equals CountKmersSpan; k ≤ 0 throws eagerly.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EnumerateKmers_CountAndContent_AgreeWithCountKmersSpan()
    {
        var gen = from s in StringOver("ACGTacgtN".ToCharArray(), 0, 60).Generator
                  from k in Gen.Choose(1, 12)
                  select (s, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var observed = new Dictionary<string, int>();
            int count = 0;
            bool windowsOk = true;
            foreach (var kmer in t.s.AsSpan().EnumerateKmers(t.k))
            {
                windowsOk &= kmer.Length == t.k && kmer.SequenceEqual(t.s.AsSpan(count, t.k));
                string key = new string(kmer).ToUpperInvariant();
                observed[key] = observed.GetValueOrDefault(key) + 1;
                count++;
            }
            var reference = t.s.AsSpan().CountKmersSpan(t.k);
            bool ok = windowsOk
                      && count == Math.Max(0, t.s.Length - t.k + 1)
                      && observed.Count == reference.Count
                      && observed.All(kv => reference.TryGetValue(kv.Key, out int v) && v == kv.Value);
            return ok.Label($"EnumerateKmers('{t.s}', {t.k}) count={count}");
        });
    }

    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void EnumerateKmers_NonPositiveK_ThrowsEagerly(int k)
    {
        var act = () => { _ = "ACGT".AsSpan().EnumerateKmers(k); };
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion
}
