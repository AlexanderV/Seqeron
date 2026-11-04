using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.CodonOptimizer;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the behaviours added by the B02 finisher
/// (docs/Validation/review-2026-09/B02.md, F31–F61 and the WP24/WP25 sweep): ASCII-only case folding
/// (F31, SequenceCase), Sec/Pyl acceptance and additive molecular weight (F32), Biopython nucleotide X
/// (F33), explicit not-calculable ENC and the ENC estimator variants (F36/F50/F51/F57), CodonW GC
/// convention (F37), RSCU pseudocount / cubar sub-families (F52/F59), Fop/CBI (F61), getorf modes
/// (F44–F46), transeq conventions and post-processing (F53/F58), getorf header ranges (F56) and
/// <see cref="GeneticCode.GetByTableNumberOrThrow"/>.
///
/// Every property asserts an exact algebraic relation, a documented range, or the equality of two
/// public paths; floating-point comparisons use 1e-9 only where a sum is re-associated.
///
/// Test Units: TRANS-PROT-001, TRANS-CODON-001, TRANS-SIXFRAME-001, CODON-ENC-001, CODON-RSCU-001,
/// CODON-STATS-001, B02-SWEEP.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Translation")]
public class B02FinisherProperties
{
    private const string Acgt = "ACGT";
    private const string IupacNucleotides = "ACGTRYSWKMBDHVN";
    private static readonly char[] NonAsciiLetters = "ſıKéÅа".ToCharArray(); // ſ ı K é Å а

    private static Gen<string> StringGen(string alphabet, int minLen, int maxLen) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alphabet.ToCharArray()).ArrayOf(n)
        select new string(chars);

    private static Gen<string> CodingGen(int minCodons, int maxCodons) =>
        from n in Gen.Choose(minCodons, maxCodons)
        from chars in Gen.Elements(Acgt.ToCharArray()).ArrayOf(3 * n)
        select new string(chars);

    private static Gen<GeneticCode> CodeGen() =>
        Gen.Elements(GeneticCode.SupportedTableNumbers.ToArray()).Select(GeneticCode.GetByTableNumber);

    // Tables whose stops are not dual-coding (the getorf / ORF API rejects 27, 28, 31).
    private static Gen<GeneticCode> OrfCodeGen() =>
        Gen.Elements(GeneticCode.SupportedTableNumbers.Where(t => t is not (27 or 28 or 31)).ToArray())
            .Select(GeneticCode.GetByTableNumber);

    private static string Dna(string s) => SequenceCase.ToUpperDna(s);

    #region F31 / WP24–25 — ASCII-only case folding (SequenceCase)

    /// <summary>SequenceCase: the char, string and span forms agree; only a–z change (by −32); length is preserved.</summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property SequenceCase_SpanStringAndCharFormsAgree()
    {
        var gen = StringGen("acgtunACGTUNxX-.*1 " + new string(NonAsciiLetters) + "\u0000￿", 0, 60);
        return Prop.ForAll(gen.ToArbitrary(), s =>
        {
            string viaString = SequenceCase.ToUpperAscii(s);
            var buffer = new char[s.Length + 3];
            SequenceCase.ToUpperAscii(s.AsSpan(), buffer);
            string viaSpan = new(buffer, 0, s.Length);
            string viaChar = new(s.Select(SequenceCase.ToUpperAscii).ToArray());
            bool onlyAsciiLowerChanged = s.Zip(viaString).All(p =>
                p.First is >= 'a' and <= 'z' ? p.Second == p.First - 32 : p.Second == p.First);
            return viaString == viaSpan && viaString == viaChar && viaString.Length == s.Length
                && onlyAsciiLowerChanged
                && SequenceCase.ToUpperDna(s) == viaString.Replace('U', 'T')
                && SequenceCase.ToUpperRna(s) == viaString.Replace('T', 'U');
        });
    }

    /// <summary>F31: a protein string with a non-ASCII letter is never accepted (ProteinSequence, TryCreate, FindMotif fold).</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property ProteinSequence_NonAsciiLetter_NeverAccepted()
    {
        var gen =
            from prefix in StringGen("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy", 0, 15)
            from bad in Gen.Elements(NonAsciiLetters)
            from suffix in StringGen("ACDEFGHIKLMNPQRSTVWY", 0, 15)
            select prefix + bad + suffix;
        return Prop.ForAll(gen.ToArbitrary(), s =>
            !ProteinSequence.TryCreate(s, out _)
            && Throws<ArgumentException>(() => new ProteinSequence(s)));
    }

    /// <summary>F31: ASCII lower case still folds — the lower-case protein equals the upper-case one.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property ProteinSequence_AsciiLowerCase_FoldsToUpper()
    {
        var gen = StringGen("ACDEFGHIKLMNPQRSTVWYUOBZJX*acdefghiklmnpqrstvwyuobzjx", 1, 40);
        return Prop.ForAll(gen.ToArbitrary(), s =>
            new ProteinSequence(s).Sequence == SequenceCase.ToUpperAscii(s));
    }

    /// <summary>F31: a codon / nucleotide string holding a non-ASCII letter is rejected by GeneticCode and Translator.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Translation_NonAsciiLetter_Rejected()
    {
        var gen =
            from code in CodeGen()
            from pos in Gen.Choose(0, 2)
            from bad in Gen.Elements(NonAsciiLetters)
            from codon in StringGen("ACGTacgt", 3, 3)
            select (code, codon.Remove(pos, 1).Insert(pos, bad.ToString()));
        return Prop.ForAll(gen.ToArbitrary(), t =>
            Throws<ArgumentException>(() => t.code.Translate(t.Item2))
            && Throws<ArgumentException>(() => t.code.Translate(t.Item2, allowNucleotideX: true))
            && Throws<ArgumentException>(() => Translator.Translate("ATG" + t.Item2, t.code)));
    }

    /// <summary>F31: ProteinPhysicochemistry skips a non-ASCII letter (EMBOSS pepstats): MW/GRAVY equal the string without it.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Physicochemistry_NonAsciiLetter_Skipped()
    {
        var gen =
            from core in StringGen("ACDEFGHIKLMNPQRSTVWY", 1, 30)
            from pos in Gen.Choose(0, 30)
            from bad in Gen.Elements(NonAsciiLetters)
            select (core, core.Insert(Math.Min(pos, core.Length), bad.ToString()));
        return Prop.ForAll(gen.ToArbitrary(), t =>
            ProteinPhysicochemistry.MolecularWeight(t.Item2) == ProteinPhysicochemistry.MolecularWeight(t.core)
            && ProteinPhysicochemistry.Gravy(t.Item2).Equals(ProteinPhysicochemistry.Gravy(t.core))
            && ProteinPhysicochemistry.IsoelectricPoint(t.Item2) == ProteinPhysicochemistry.IsoelectricPoint(t.core));
    }

    #endregion

    #region F32 — Sec / Pyl

    /// <summary>F32: U/O are valid residues; MW is additive over concatenation (minus one water per new bond).</summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property MolecularWeight_WithSecPyl_IsAdditive()
    {
        var gen =
            from a in StringGen("ACDEFGHIKLMNPQRSTVWYUO", 1, 25)
            from b in StringGen("ACDEFGHIKLMNPQRSTVWYUO", 1, 25)
            select (a, b);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var pa = new ProteinSequence(t.a);
            var pb = new ProteinSequence(t.b);
            var pab = new ProteinSequence(t.a + t.b);
            double expected = pa.MolecularWeight() + pb.MolecularWeight() - ProteinPhysicochemistry.AverageWaterMass;
            int sec = t.a.Count(c => c == 'U'), pyl = t.a.Count(c => c == 'O');
            string three = pa.ToThreeLetterCode();
            return Math.Abs(pab.MolecularWeight() - expected) < 1e-9
                && CountOccurrences(three, "Sec") == sec && CountOccurrences(three, "Pyl") == pyl;
        });
    }

    /// <summary>F32: Translator never emits U or O (UGA/UAG are stops in table 1 and residues of the code otherwise).</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Translator_NeverEmitsSecOrPyl()
    {
        var gen = from code in CodeGen() from s in CodingGen(1, 40) select (code, s);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string protein = Translator.Translate(t.s, t.code).Sequence;
            return !protein.Contains('U') && !protein.Contains('O');
        });
    }

    #endregion

    #region F33 — nucleotide X

    /// <summary>
    /// F33: with allowNucleotideX an X-codon translates exactly as the same codon with X → N, or throws —
    /// and it throws only where the N codon is a possible stop (reads 'X'). Without the flag X is rejected.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 400)]
    public Property NucleotideX_EqualsN_WhenItTranslates()
    {
        var gen =
            from code in CodeGen()
            from codon in StringGen(IupacNucleotides + "Xx", 3, 3)
            where codon.Contains('X') || codon.Contains('x')
            select (code, codon);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string asN = t.codon.Replace('X', 'N').Replace('x', 'n');
            char viaN = t.code.Translate(asN);
            bool rejectedByDefault = Throws<ArgumentException>(() => t.code.Translate(t.codon));
            try
            {
                return rejectedByDefault && t.code.Translate(t.codon, allowNucleotideX: true) == viaN
                    && t.code.IsStopCodon(t.codon, true) == t.code.IsStopCodon(asN)
                    && t.code.IsStartCodon(t.codon, true) == t.code.IsStartCodon(asN);
            }
            catch (ArgumentException)
            {
                return rejectedByDefault && viaN == 'X';
            }
        });
    }

    /// <summary>F33: Translator with allowNucleotideX equals the translation of the X → N string whenever it translates.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Translator_NucleotideX_EqualsNString()
    {
        var gen =
            from code in CodeGen()
            from s in StringGen("ACGTACGTACGTXN", 3, 45)
            from frame in Gen.Choose(0, 2)
            select (code, s, frame);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string viaN = Translator.Translate(t.s.Replace('X', 'N'), t.code, t.frame, false, false).Sequence;
            try
            {
                return Translator.Translate(t.s, t.code, t.frame, false, allowNucleotideX: true).Sequence == viaN;
            }
            catch (ArgumentException)
            {
                return viaN.Contains('X');
            }
        });
    }

    #endregion

    #region F36 / F50 / F51 / F57 — ENC

    private static int SenseCodons(GeneticCode code) => code.CodonTable.Count(kv => kv.Value != '*');

    private static int AminoAcids(GeneticCode code) => code.CodonTable.Values.Where(v => v != '*').Distinct().Count();

    /// <summary>
    /// F36: TryCalculateEnc is false exactly where the old CalculateEnc returns 0; when calculable, Wright's
    /// Nc lies in [number of amino acids, sense codons] (CodonW enc_out re-adjustment).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property TryCalculateEnc_FalseIffOldApiZero_AndInBounds()
    {
        var gen = from code in CodeGen() from s in CodingGen(0, 200) select (code, s);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            bool ok = CodonUsageAnalyzer.TryCalculateEnc(t.s, t.code, out double enc);
            double old = CodonUsageAnalyzer.CalculateEnc(t.s, t.code);
            bool viaOptions = CodonUsageAnalyzer.TryCalculateEnc(t.s, t.code, EncOptions.Default, out double encOpt);
            bool bounds = !ok || (enc >= AminoAcids(t.code) && enc <= SenseCodons(t.code));
            return ok == (old != 0) && enc == old && viaOptions == ok && encOpt == enc && bounds
                && (ok || enc == 0)
                && (t.s.Length == 0 || CodonUsageAnalyzer.GetStatistics(t.s, t.code).IsEncCalculable == ok);
        });
    }

    /// <summary>F51: the CodonW cap replaces the code's sense-codon ceiling by 61 (≤ 61; same calculability).</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property CodonWCap_AtMost61_SameCalculability()
    {
        var gen = from code in CodeGen() from s in CodingGen(0, 200) select (code, s);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            bool plain = CodonUsageAnalyzer.TryCalculateEnc(t.s, t.code, out double enc);
            bool capped = CodonUsageAnalyzer.TryCalculateEnc(t.s, t.code, new EncOptions { CodonWCap = true }, out double encCap);
            if (plain != capped) return false;
            if (!plain) return encCap == 0;
            int sense = SenseCodons(t.code);
            return encCap <= 61 && (sense == 61 ? encCap == enc : sense < 61 ? encCap >= enc : encCap <= enc);
        });
    }

    /// <summary>
    /// F50/F57: Novembre 2002, Sun–Yang–Xia 2013 and ENCprime Nc/Nc′: Try ⇔ CalculateEnc ≠ 0; calculable for
    /// every input holding a sense codon; Novembre/SYX never exceed the sense-codon count; ENCprime Nc lies in
    /// [amino acids, sense codons] (F ≤ 1 and the 1/k floor); Nc′ is positive and finite.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property EncVariants_TryConsistent_AndInBounds()
    {
        var gen =
            from code in CodeGen()
            from s in CodingGen(0, 150)
            from method in Gen.Elements(EncMethod.Novembre2002, EncMethod.SunYangXia2013, EncMethod.EncPrimeNc, EncMethod.EncPrimeNcPrime)
            select (code, s, method);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var options = new EncOptions { Method = t.method };
            bool ok = CodonUsageAnalyzer.TryCalculateEnc(t.s, t.code, options, out double enc);
            double viaCalc = CodonUsageAnalyzer.CalculateEnc(t.s, t.code, options);
            bool hasSense = Translator.SplitInFrameCodons(t.s).Any(c => c is not null && t.code.Translate(c) != '*');
            if (ok != (viaCalc != 0) || enc != viaCalc) return false;
            if (!ok) return enc == 0 && !hasSense;
            int sense = SenseCodons(t.code);
            return t.method switch
            {
                EncMethod.EncPrimeNc => enc >= AminoAcids(t.code) - 1e-9 && enc <= sense + 1e-9,
                EncMethod.EncPrimeNcPrime => enc > 0 && double.IsFinite(enc),
                _ => enc > 0 && enc <= sense,
            };
        });
    }

    #endregion

    #region F37 — CodonW GC convention

    /// <summary>
    /// F37: the CodonW convention drops the code's stop codons from GC1–3 only; counts, RSCU, ENC and GC3s are
    /// unchanged; positionwise, cusp ≤ CodonW where no stop has G/C there and cusp ≥ CodonW where every stop
    /// has; a stop-free gene gives identical values; inserting stop codons leaves the CodonW values unchanged.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property CodonWGc_StopInvariance_AndOrderRelations()
    {
        var gen =
            from code in CodeGen()
            from s in CodingGen(1, 60)
            from inserts in Gen.Choose(0, 60).ArrayOf()
            select (code, s, inserts);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var stops = t.code.CodonTable.Where(kv => kv.Value == '*').Select(kv => kv.Key.Replace('U', 'T')).ToList();
            var cusp = CodonUsageAnalyzer.GetStatistics(t.s, t.code, CodonGcConvention.EmbossCusp);
            var codonW = CodonUsageAnalyzer.GetStatistics(t.s, t.code, CodonGcConvention.CodonW);
            if (cusp.TotalCodons != codonW.TotalCodons || cusp.Gc3s != codonW.Gc3s || cusp.Enc != codonW.Enc)
                return false;

            double[] c = { cusp.Gc1, cusp.Gc2, cusp.Gc3 }, w = { codonW.Gc1, codonW.Gc2, codonW.Gc3 };
            bool hasStop = Translator.SplitInFrameCodons(t.s).Any(x => x is not null && stops.Contains(x));
            bool Relation(int k)
            {
                if (!hasStop) return c[k] == w[k];
                if (stops.All(x => x[k] is not ('G' or 'C'))) return c[k] <= w[k] + 1e-12;
                if (stops.All(x => x[k] is 'G' or 'C')) return c[k] >= w[k] - 1e-12;
                return true; // mixed stop bases at position k: no order relation
            }
            bool relations = Enumerable.Range(0, 3).All(Relation);

            // Insert stop codons at codon boundaries: CodonW GC1–3 are unchanged (unless the gene was all stops).
            string withStops = t.s;
            if (stops.Count > 0)
            {
                foreach (int i in t.inserts)
                {
                    int at = 3 * (i % (withStops.Length / 3 + 1));
                    withStops = withStops.Insert(at, stops[i % stops.Count]);
                }
            }
            var codonW2 = CodonUsageAnalyzer.GetStatistics(withStops, t.code, CodonGcConvention.CodonW);
            bool invariant = codonW.Gc1 == codonW2.Gc1 && codonW.Gc2 == codonW2.Gc2 && codonW.Gc3 == codonW2.Gc3;
            return relations && invariant;
        });
    }

    #endregion

    #region F52 / F59 — RSCU pseudocount and cubar sub-families

    private static IEnumerable<IGrouping<string, string>> Groups(GeneticCode code, RscuFamilyLevel level) =>
        code.CodonTable.GroupBy(
            kv => level == RscuFamilyLevel.AminoAcid ? kv.Value.ToString() : kv.Value + "_" + kv.Key.Replace('U', 'T')[..2],
            kv => kv.Key.Replace('U', 'T'));

    /// <summary>
    /// F52/F59: with c &gt; 0 every family (amino acid or cubar sub-family) sums to its size; with c = 0 a family
    /// present in the gene sums to its size and an absent one is 0; all 64 codons are reported.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Rscu_FamilySums_EqualFamilySize()
    {
        var gen =
            from code in CodeGen()
            from s in CodingGen(0, 80)
            from c in Gen.Elements(0.0, 0.0, 0.5, 1.0, 3.7)
            from level in Gen.Elements(RscuFamilyLevel.AminoAcid, RscuFamilyLevel.Subfamily)
            select (code, s, c, level);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var counts = Translator.CountCodons(t.s);
            var rscu = CodonUsageAnalyzer.CalculateRscu(counts, t.code, t.c, t.level);
            if (rscu.Count != 64) return false;
            return Groups(t.code, t.level).All(g =>
            {
                int size = g.Count();
                long present = g.Sum(codon => counts.GetValueOrDefault(codon));
                double sum = g.Sum(codon => rscu[codon]);
                double expected = t.c > 0 || present > 0 ? size : 0;
                return Math.Abs(sum - expected) < 1e-9;
            });
        });
    }

    /// <summary>F52/F59: c → 0 is the CodonW RSCU for every present family; c = 0 equals the historical overload bit for bit.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Rscu_PseudocountLimit_IsCodonW()
    {
        var gen = from code in CodeGen() from s in CodingGen(1, 80) select (code, s);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var counts = Translator.CountCodons(t.s);
            var codonW = CodonUsageAnalyzer.CalculateRscu(counts, t.code);
            var zero = CodonUsageAnalyzer.CalculateRscu(counts, t.code, 0.0);
            var tiny = CodonUsageAnalyzer.CalculateRscu(counts, t.code, 1e-12);
            bool bitIdentical = codonW.All(kv => zero[kv.Key] == kv.Value);
            bool limit = Groups(t.code, RscuFamilyLevel.AminoAcid)
                .Where(g => g.Sum(codon => counts.GetValueOrDefault(codon)) > 0)
                .SelectMany(g => g)
                .All(codon => Math.Abs(tiny[codon] - codonW[codon]) < 1e-6);
            return bitIdentical && limit;
        });
    }

    #endregion

    #region F61 — Fop / CBI

    private static Gen<CodonWOptimalCodonSet> OptimalSetGen() => Gen.Elements(CodonWOptimalCodonSets.All.ToArray());

    /// <summary>
    /// F61: Fop ∈ [0, 1] (modified Fop ∈ [−1, 1]) and TryCalculateFop is false exactly when no codon of a
    /// family holding an optimal (or rare) codon occurs; CBI ≤ 1 and TryCalculateCbi is false exactly when
    /// CodonW's denominator Σn − Σ n·k/d vanishes. CalculateX equals the Try value (0 when not calculable).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property FopCbi_RangesAndTryConditions()
    {
        var gen =
            from code in CodeGen()
            from s in CodingGen(0, 60)
            from set in OptimalSetGen()
            from modified in Gen.Elements(false, true)
            select (code, s, set, modified);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var counts = Translator.CountCodons(t.s);
            var rare = t.modified && t.set.RareCodons.Count > 0 ? t.set.RareCodons : null;
            var families = t.code.CodonTable.Where(kv => kv.Value != '*')
                .GroupBy(kv => kv.Value, kv => kv.Key.Replace('U', 'T'))
                .Select(g => g.ToList()).Where(g => g.Count > 1).ToList();

            long fopTotal = families
                .Where(f => f.Any(c => t.set.OptimalCodons.Contains(c) || (rare?.Contains(c) ?? false)))
                .Sum(f => f.Sum(c => (long)counts.GetValueOrDefault(c)));
            double cbiDenominator = families.Sum(f =>
            {
                int k = f.Count(t.set.OptimalCodons.Contains);
                if (k == 0) return 0.0;
                long n = f.Sum(c => (long)counts.GetValueOrDefault(c));
                return n - (double)n * k / f.Count;
            });

            bool fopOk = CodonUsageAnalyzer.TryCalculateFop(t.s, out double fop, t.set.OptimalCodons, t.code, rare);
            bool cbiOk = CodonUsageAnalyzer.TryCalculateCbi(t.s, out double cbi, t.set.OptimalCodons, t.code);
            bool fopRange = !fopOk || (rare is null ? fop is >= 0 and <= 1 : fop is >= -1 and <= 1);
            return fopOk == (fopTotal > 0) && cbiOk == (Math.Abs(cbiDenominator) >= 1e-9)
                && fopRange && (!cbiOk || cbi <= 1 + 1e-12) && (fopOk || fop == 0) && (cbiOk || cbi == 0)
                && CodonUsageAnalyzer.CalculateFop(t.s, t.set.OptimalCodons, t.code, rare) == fop
                && CodonUsageAnalyzer.CalculateCbi(t.s, t.set.OptimalCodons, t.code) == cbi;
        });
    }

    #endregion

    #region F44 / F45 — getorf modes

    private static string Strand(string upper, int frame) =>
        frame > 0 ? upper : DnaSequence.GetReverseComplementString(upper);

    private static Gen<(GeneticCode Code, string S)> OrfInput(string alphabet, int maxLen) =>
        from code in OrfCodeGen() from s in StringGen(alphabet, 0, maxLen) select (code, s);

    /// <summary>
    /// F44: with -minsize 0 the getorf -find 0 records of each frame tile that frame's complete codons
    /// stop-to-stop: contiguous, from the frame offset to the last complete codon; each protein has one residue
    /// per non-stop codon, no '*', and a terminated record ends on a codon the code reads as '*'.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfMode0_PartitionsFramesStopToStop()
    {
        return Prop.ForAll(OrfInput("ACGTACGTACGTN", 90).ToArbitrary(), t =>
        {
            string upper = Dna(t.S);
            int len = upper.Length;
            var records = Translator.FindOrfs(t.S, OrfFindMode.TranslationStopToStop, t.Code, minSize: 0).ToList();
            foreach (int frame in new[] { 1, 2, 3, -1, -2, -3 })
            {
                string strand = Strand(upper, frame);
                int offset = Math.Abs(frame) - 1;
                var inFrame = records.Where(r => r.Frame == frame).ToList();
                int codons = Math.Max(0, (len - offset) / 3);
                if (codons == 0)
                {
                    if (inFrame.Count != 0) return false;
                    continue;
                }
                int expectedStart = offset;
                foreach (var r in inFrame)
                {
                    if (r.StartPosition != expectedStart || (r.EndPosition - r.StartPosition + 1) % 3 != 0) return false;
                    int n = (r.EndPosition - r.StartPosition + 1) / 3;
                    if (r.Sequence.Length != n - (r.IsTerminated ? 1 : 0) || r.Sequence.Contains('*')) return false;
                    string last = strand.Substring(r.EndPosition - 2, 3);
                    if (r.IsTerminated != (t.Code.Translate(last) == '*')) return false;
                    expectedStart = r.EndPosition + 1;
                }
                // Covered up to the last complete codon (a frame ending on a STOP leaves no trailing record).
                if (expectedStart != offset + 3 * codons) return false;
            }
            return true;
        });
    }

    /// <summary>F44: mode 2 nucleotides are the scanned strand's bases of the record without the STOP; modes 0 and 2 give the same regions.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfMode2_SequenceIsStrandSubstring()
    {
        return Prop.ForAll(OrfInput("ACGTACGTN", 90).ToArbitrary(), t =>
        {
            string upper = Dna(t.S);
            var nuc = Translator.FindOrfs(t.S, OrfFindMode.NucleotidesStopToStop, t.Code, minSize: 0).ToList();
            var prot = Translator.FindOrfs(t.S, OrfFindMode.TranslationStopToStop, t.Code, minSize: 0).ToList();
            bool sameRegions = nuc.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated))
                .SequenceEqual(prot.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated)));
            return sameRegions && nuc.All(r =>
            {
                int length = r.EndPosition - r.StartPosition + 1 - (r.IsTerminated ? 3 : 0);
                return r.Sequence == Strand(upper, r.Frame).Substring(r.StartPosition, length);
            });
        });
    }

    /// <summary>
    /// F44: flank windows (modes 4–6) lie inside the sequence, are 3 + 2·flanking long, are centred on the
    /// codon, hold the strand's bases, and the flanked codon is a START (mode 4) or a STOP (modes 5, 6).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfFlankWindows_WithinBounds()
    {
        var gen =
            from input in OrfInput("ACGTACGTACGTN", 120)
            from mode in Gen.Elements(OrfFindMode.FlankingStart, OrfFindMode.FlankingInitialStop, OrfFindMode.FlankingEndingStop)
            from flanking in Gen.Choose(0, 12)
            select (input.Code, input.S, mode, flanking);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string upper = Dna(t.S);
            return Translator.FindOrfs(t.S, t.mode, t.Code, minSize: 0, flanking: t.flanking).All(r =>
            {
                string strand = Strand(upper, r.Frame);
                string codon = strand.Substring(r.CodonPosition, 3);
                bool codonKind = t.mode == OrfFindMode.FlankingStart
                    ? t.Code.IsStartCodon(codon)
                    : t.Code.Translate(codon) == '*';
                return r.StartPosition >= 0 && r.EndPosition < upper.Length
                    && r.EndPosition - r.StartPosition + 1 == 3 + 2 * t.flanking
                    && r.CodonPosition == r.StartPosition + t.flanking
                    && r.Sequence == strand.Substring(r.StartPosition, 3 + 2 * t.flanking)
                    && codonKind;
            });
        });
    }

    /// <summary>
    /// F44/F56: getorf -find 1 equals the legacy START→STOP overload (minSize = 3·minLength), and the two
    /// ToGetOrfRange helpers agree record by record.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfMode1_EqualsLegacy_AndRangesAgree()
    {
        var gen = from input in OrfInput("ACGTACGTACGTATGATG", 120) from minLength in Gen.Choose(0, 6) select (input.Code, input.S, minLength);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            int len = t.S.Length;
            var legacy = Translator.FindOrfs(t.S, t.Code, t.minLength)
                .Select(o => (o.StartPosition, o.EndPosition, o.Frame, Seq: o.Protein.Sequence,
                    Range: len > 0 ? o.ToGetOrfRange(len) : default))
                .Order().ToList();
            var getorf = Translator.FindOrfs(t.S, OrfFindMode.TranslationStartToStop, t.Code, minSize: 3 * t.minLength)
                .Select(o => (o.StartPosition, o.EndPosition, o.Frame, Seq: o.Sequence,
                    Range: len > 0 ? o.ToGetOrfRange(len) : default))
                .Order().ToList();
            return legacy.SequenceEqual(getorf);
        });
    }

    /// <summary>
    /// F46: -circular keeps every linear terminated START→STOP ORF, and every linear terminated STOP→STOP ORF
    /// that does not start the frame (the leading region is not an ORF on a circle).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Circular_PreservesLinearOrfsNotCrossingOrigin()
    {
        var gen = from input in OrfInput("ACGTACGTATG", 90)
                  from mode in Gen.Elements(OrfFindMode.TranslationStartToStop, OrfFindMode.TranslationStopToStop,
                      OrfFindMode.NucleotidesStartToStop, OrfFindMode.NucleotidesStopToStop)
                  select (input.Code, input.S, mode);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var circular = Translator.FindOrfs(t.S, t.mode, t.Code, minSize: 0, circular: true).ToHashSet();
            bool stopToStop = t.mode is OrfFindMode.TranslationStopToStop or OrfFindMode.NucleotidesStopToStop;
            return Translator.FindOrfs(t.S, t.mode, t.Code, minSize: 0)
                .Where(r => r.IsTerminated && (!stopToStop || r.StartPosition >= 3))
                .All(circular.Contains);
        });
    }

    #endregion

    #region F53 / F58 — transeq conventions and post-processing

    private static readonly int[] FrameKeys = { 1, 2, 3, -1, -2, -3 };

    /// <summary>
    /// F53: on A/C/G/T, every frame of the transeq -alternative convention is the Biopython frame plus the
    /// dangling-base residue (none for 0 dangling bases, X for one, the b1 b2 N codon for two); forward frames
    /// are identical in both EMBOSS conventions.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property TranseqAlternative_IsBiopythonPlusDangle()
    {
        var gen = from code in CodeGen() from s in StringGen(Acgt, 1, 50) select (code, s);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var bio = Translator.TranslateSixFrames(t.s, t.code, SixFrameConvention.Biopython);
            var alt = Translator.TranslateSixFrames(t.s, t.code, SixFrameConvention.EmbossTranseqAlternative);
            var locked = Translator.TranslateSixFrames(t.s, t.code, SixFrameConvention.EmbossTranseq);
            string rc = DnaSequence.GetReverseComplementString(t.s);
            return FrameKeys.All(k =>
            {
                string strand = k > 0 ? t.s : rc;
                int offset = Math.Abs(k) - 1;
                int dangle = Math.Max(0, strand.Length - offset) % 3;
                string extra = dangle switch
                {
                    0 => "",
                    1 => "X",
                    _ => t.code.Translate(string.Concat(strand.AsSpan(strand.Length - 2), "N")).ToString(),
                };
                return alt[k].Sequence == bio[k].Sequence + extra
                    && (k < 0 || locked[k].Sequence == alt[k].Sequence);
            });
        });
    }

    /// <summary>
    /// F58: TranseqOptions.Default equals the plain convention overload; -trim removes exactly the trailing
    /// X/* run; -clean (after trim) turns every '*' into 'X'; both are idempotent on their own output.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property TranseqTrimClean_ExactAndIdempotent()
    {
        var gen =
            from code in CodeGen()
            from s in StringGen("ACGTACGTACGTN", 1, 50)
            from convention in Gen.Elements(SixFrameConvention.EmbossTranseq, SixFrameConvention.EmbossTranseqAlternative)
            select (code, s, convention);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var plain = Translator.TranslateSixFrames(t.s, t.code, t.convention);
            var dflt = Translator.TranslateSixFrames(t.s, t.code, t.convention, TranseqOptions.Default);
            var trim = Translator.TranslateSixFrames(t.s, t.code, t.convention, new TranseqOptions { Trim = true });
            var clean = Translator.TranslateSixFrames(t.s, t.code, t.convention, new TranseqOptions { Clean = true });
            var both = Translator.TranslateSixFrames(t.s, t.code, t.convention, new TranseqOptions { Trim = true, Clean = true });
            return FrameKeys.All(k =>
            {
                string p = plain[k].Sequence;
                string trimmed = p.TrimEnd('X', '*');
                return dflt[k].Sequence == p
                    && trim[k].Sequence == trimmed
                    && clean[k].Sequence == p.Replace('*', 'X')
                    && both[k].Sequence == trimmed.Replace('*', 'X')
                    && trimmed.TrimEnd('X', '*') == trimmed
                    && both[k].Sequence.Replace('*', 'X') == both[k].Sequence;
            });
        });
    }

    #endregion

    #region GetByTableNumberOrThrow

    /// <summary>GetByTableNumberOrThrow succeeds exactly for SupportedTableNumbers and returns the GetByTableNumber instance.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetByTableNumberOrThrow_IffSupported()
    {
        return Prop.ForAll(Gen.Choose(-5, 45).ToArbitrary(), n =>
        {
            bool supported = GeneticCode.SupportedTableNumbers.Contains(n);
            try
            {
                var code = GeneticCode.GetByTableNumberOrThrow(n, "genetic_code");
                return supported && code.TableNumber == n && ReferenceEquals(code, GeneticCode.GetByTableNumber(n));
            }
            catch (ArgumentException e)
            {
                return !supported && e.ParamName == "genetic_code" && e.Message.Contains($"table: {n}.");
            }
        });
    }

    #endregion

    private static bool Throws<T>(Action action) where T : Exception
    {
        try { action(); return false; }
        catch (T) { return true; }
    }

    private static int CountOccurrences(string s, string token)
    {
        int count = 0;
        for (int i = s.IndexOf(token, StringComparison.Ordinal); i >= 0; i = s.IndexOf(token, i + 1, StringComparison.Ordinal))
            count++;
        return count;
    }
}
