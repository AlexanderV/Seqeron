using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Bacterial σ70 promoters: −35 / −10 hexamer pairing over a spacer (Harley &amp; Reynolds 1987 consensus) and the
/// La Fleur, Hossain &amp; Salis (2022) Promoter Calculator v1.0 free-energy model.
/// </summary>
public static partial class MotifFinder
{
    #region σ70 promoters — consensus pairing (Harley & Reynolds 1987)

    /// <summary>Harley &amp; Reynolds (1987) σ70 −35 consensus hexamer.</summary>
    public const string Sigma70Minus35Consensus = "TTGACA";

    /// <summary>Harley &amp; Reynolds (1987) σ70 −10 (Pribnow) consensus hexamer.</summary>
    public const string Sigma70Minus10Consensus = "TATAAT";

    /// <summary>Optimal −35/−10 spacer length (Harley &amp; Reynolds 1987: 17 ± 1 bp in 92 % of promoters).</summary>
    public const int Sigma70OptimalSpacer = 17;

    /// <summary>
    /// Pairs σ70 −35 and −10 boxes: every −35 hexamer within <paramref name="maxMismatches35"/> mismatches of
    /// TTGACA followed, after a spacer of <paramref name="minSpacer"/> … <paramref name="maxSpacer"/> bp, by a −10
    /// hexamer within <paramref name="maxMismatches10"/> mismatches of TATAAT (Harley &amp; Reynolds 1987, NAR
    /// 15:2343: consensus TTGACA / TATAAT, spacer 15–21 bp, 17 ± 1 bp in 92 % of promoters).
    /// </summary>
    /// <remarks>
    /// <para>Each candidate reports the Hamming distances of both hexamers to their consensus (canonical
    /// <see cref="SequenceExtensions.HammingDistance(ReadOnlySpan{char}, ReadOnlySpan{char})"/>) and the spacer's
    /// deviation from 17 bp. Order: strand '+' before '−', then ascending −35 start, then ascending spacer.
    /// This is the consensus definition only; for a published quantitative score of a promoter (free energy and
    /// predicted transcription rate, with UP element, extended −10, discriminator and ITR terms) use
    /// <see cref="PredictSigma70Promoters"/>.</para>
    /// <para>Coordinates are 0-based on the forward strand. On the minus strand the boxes are read 5'→3' on that
    /// strand and the −35 box lies to the right of the −10 box in forward coordinates.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="maxMismatches35">Maximum mismatches of the −35 box to TTGACA (0 … 6, default 2).</param>
    /// <param name="maxMismatches10">Maximum mismatches of the −10 box to TATAAT (0 … 6, default 2).</param>
    /// <param name="minSpacer">Minimum spacer length (default 15).</param>
    /// <param name="maxSpacer">Maximum spacer length (default 21).</param>
    /// <param name="bothStrands">Also scan the reverse-complement strand (default false).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A mismatch limit outside 0 … 6, a negative spacer or
    /// <paramref name="minSpacer"/> &gt; <paramref name="maxSpacer"/>.</exception>
    public static IEnumerable<Sigma70PromoterCandidate> FindSigma70Promoters(
        DnaSequence sequence,
        int maxMismatches35 = 2,
        int maxMismatches10 = 2,
        int minSpacer = 15,
        int maxSpacer = 21,
        bool bothStrands = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (maxMismatches35 is < 0 or > 6)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches35), maxMismatches35, "Must be in 0 … 6.");
        if (maxMismatches10 is < 0 or > 6)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches10), maxMismatches10, "Must be in 0 … 6.");
        if (minSpacer < 0)
            throw new ArgumentOutOfRangeException(nameof(minSpacer), minSpacer, "Must be non-negative.");
        if (maxSpacer < minSpacer)
            throw new ArgumentOutOfRangeException(nameof(maxSpacer), maxSpacer, "Must be ≥ minSpacer.");
        return FindSigma70PromotersCore(sequence.Sequence, maxMismatches35, maxMismatches10, minSpacer, maxSpacer, bothStrands);
    }

    private static IEnumerable<Sigma70PromoterCandidate> FindSigma70PromotersCore(
        string seq, int mm35, int mm10, int minSpacer, int maxSpacer, bool bothStrands)
    {
        foreach (var c in PairSigma70Boxes(seq, mm35, mm10, minSpacer, maxSpacer, '+'))
            yield return c;
        if (!bothStrands)
            yield break;

        int n = seq.Length;
        var minus = PairSigma70Boxes(DnaSequence.GetReverseComplementString(seq), mm35, mm10, minSpacer, maxSpacer, '-')
            .Select(c => c with
            {
                Minus35Start = n - c.Minus35Start - 6,
                Minus10Start = n - c.Minus10Start - 6,
            })
            .OrderBy(c => c.Minus35Start)
            .ThenBy(c => c.Spacer);
        foreach (var c in minus)
            yield return c;
    }

    /// <summary>Pairs boxes on one strand (coordinates of <paramref name="seq"/>).</summary>
    private static IEnumerable<Sigma70PromoterCandidate> PairSigma70Boxes(
        string seq, int mm35, int mm10, int minSpacer, int maxSpacer, char strand)
    {
        const int box = 6;
        int n = seq.Length;
        for (int i = 0; i + box <= n; i++)
        {
            int d35 = seq.AsSpan(i, box).HammingDistance(Sigma70Minus35Consensus);
            if (d35 > mm35)
                continue;

            for (int s = minSpacer; s <= maxSpacer; s++)
            {
                int j = i + box + s;
                if (j + box > n)
                    break;
                int d10 = seq.AsSpan(j, box).HammingDistance(Sigma70Minus10Consensus);
                if (d10 > mm10)
                    continue;

                yield return new Sigma70PromoterCandidate(
                    strand, i, seq.Substring(i, box), j, seq.Substring(j, box), s,
                    d35, d10, d35 + d10, Math.Abs(s - Sigma70OptimalSpacer));
            }
        }
    }

    #endregion

    #region σ70 promoters — Promoter Calculator v1.0 (La Fleur, Hossain & Salis 2022)

    /// <summary>
    /// σ70 promoter prediction with the Promoter Calculator v1.0 (La Fleur, Hossain &amp; Salis 2022, Nat Commun
    /// 13:5159, "Automated model-predictive design of synthetic promoters to control transcriptional profiles in
    /// bacteria") — a port of the authors' reference implementation (<c>Promoter_Calculator_v1_0.py</c>,
    /// <c>util.py</c>, <c>free_energy_coeffs.npy</c>, <c>model_intercept.npy</c>; github.com/hsalis/SalisLabCode).
    /// </summary>
    /// <remarks>
    /// <para>For every candidate transcription start site (TSS) the model enumerates every promoter configuration
    /// UP (24 nt) · 1 nt · −35 hexamer · spacer (15 … 20 nt) · −10 hexamer · discriminator (6 … 10 nt) · ITR (20 nt,
    /// starting at the TSS) and computes the linear free-energy model
    /// ΔG_total = ΔG_−10 + ΔG_−35 + ΔG_disc + ΔG_ITR + ΔG_ext−10 + ΔG_spacer + ΔG_UP + intercept, with
    /// ΔG_−10/ΔG_−35 = trained 3-mer energies of each hexamer half, ΔG_disc of the discriminator's first 3-mer,
    /// ΔG_ext−10 of spacer[−3:−1], ΔG_spacer = 0.1463·s² − 4.9113·s + 41.119, ΔG_ITR from the DNA:DNA − RNA:DNA
    /// hybrid energy of ITR[0:14], ΔG_UP from minor-groove width of both UP halves and the bending rigidity of
    /// UP + −35 + spacer[0:14]. The configuration with the smallest ΔG_total is kept per TSS (first found on ties,
    /// discriminator length outer loop, spacer length inner loop), and the transcription rate is
    /// K·exp(−β·ΔG_total) with K = 42 and β = 1.636217004872062 (E. coli MG1655) or 0.81632623 (in vitro).</para>
    /// <para><see cref="Sigma70PromoterPrediction.Tss"/> follows the reference: the boundary index of the first
    /// transcribed base. On '+' it is that base's forward index; on '−' the reference maps the reverse-complement
    /// TSS t to n − t, so the first transcribed base is forward index Tss − 1. Segment starts are reported in
    /// forward coordinates (segment sequences read 5'→3' on their own strand). Only TSSs with a full configuration
    /// (≥ 67 nt upstream including the discriminator, 20 nt downstream) are returned.</para>
    /// <para>Order: '+' predictions by ascending TSS, then '−' predictions by ascending TSS.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence (A/C/G/T).</param>
    /// <param name="bothStrands">Also predict on the reverse-complement strand (default true, as the reference <c>run</c>).</param>
    /// <param name="inVitro">Use the reference's in-vitro β (default false: E. coli MG1655).</param>
    /// <returns>The minimum-ΔG configuration for every TSS with at least one complete configuration.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static IReadOnlyList<Sigma70PromoterPrediction> PredictSigma70Promoters(
        DnaSequence sequence, bool bothStrands = true, bool inVitro = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        string seq = sequence.Sequence;
        int n = seq.Length;
        double beta = inVitro ? PromoterCalculatorV1.BetaInVitro : PromoterCalculatorV1.BetaEColi;

        var result = new List<Sigma70PromoterPrediction>();
        result.AddRange(PromoterCalculatorV1.Predict(seq, beta, '+', p => p));
        if (bothStrands)
        {
            var minus = PromoterCalculatorV1.Predict(DnaSequence.GetReverseComplementString(seq), beta, '-',
                p => p with
                {
                    Tss = n - p.Tss,
                    UpStart = n - p.UpStart - 24,
                    Minus35Start = n - p.Minus35Start - 6,
                    SpacerStart = n - p.SpacerStart - p.Spacer.Length,
                    Minus10Start = n - p.Minus10Start - 6,
                    DiscriminatorStart = n - p.DiscriminatorStart - p.Discriminator.Length,
                });
            minus.Reverse(); // reverse-complement TSS ascending → forward TSS ascending
            result.AddRange(minus);
        }

        return result;
    }

    /// <summary>Promoter Calculator v1.0 parameters and kernel (reference: hsalis/SalisLabCode, Promoter_Calculator).</summary>
    private static class PromoterCalculatorV1
    {
        internal const double K = 42.0;
        internal const double BetaEColi = 1.636217004872062;
        internal const double BetaInVitro = 0.81632623;
        internal const double Intercept = -0.9583584537583338;

        private const int UpLength = 24;
        private const int UpHex35Gap = 1;
        private const int HexLength = 6;
        private const int MinSpacer = 15, MaxSpacerExclusive = 21;
        private const int MinDisc = 6, MaxDiscExclusive = 11;
        private const int ItrLength = 20;

        // Feature normalisers (maximum in the training set), in the reference order
        // [UP distal groove width, UP proximal groove width, ITR hybrid ΔG, rigidity].
        private const double NormDistal = 256.0, NormProximal = 255.0, NormItr = 4.300000000000002, NormRigidity = 25.780434782608694;

        internal static List<Sigma70PromoterPrediction> Predict(
            string seq, double beta, char strand, Func<Sigma70PromoterPrediction, Sigma70PromoterPrediction> map)
        {
            var list = new List<Sigma70PromoterPrediction>();
            int n = seq.Length;
            for (int tss = 0; tss < n; tss++)
            {
                Sigma70PromoterPrediction? best = null;
                for (int disc = MinDisc; disc < MaxDiscExclusive; disc++)
                {
                    if (tss - disc < 0 || tss + ItrLength > n)
                        continue;
                    for (int spacer = MinSpacer; spacer < MaxSpacerExclusive; spacer++)
                    {
                        int upStart = tss - disc - HexLength - spacer - HexLength - UpLength - UpHex35Gap;
                        if (upStart < 0)
                            continue;
                        var p = Evaluate(seq, tss, disc, spacer, upStart, beta, strand);
                        if (best is null || p.DeltaGTotal < best.DeltaGTotal)
                            best = p;
                    }
                }

                if (best is not null)
                    list.Add(map(best));
            }

            return list;
        }

        private static Sigma70PromoterPrediction Evaluate(string seq, int tss, int disc, int spacer, int upStart, double beta, char strand)
        {
            int hex10Start = tss - disc - HexLength;
            int spacerStart = hex10Start - spacer;
            int hex35Start = spacerStart - HexLength;
            int discStart = tss - disc;

            double dg10 = Minus10First[Kmer3(seq, hex10Start)] + Minus10Last[Kmer3(seq, hex10Start + 3)];
            double dg35 = Minus35First[Kmer3(seq, hex35Start)] + Minus35Last[Kmer3(seq, hex35Start + 3)];
            double dgDisc = Discriminator[Kmer3(seq, discStart)];
            double dgExt10 = Extended10[Dimer(seq, spacerStart + spacer - 3)]; // spacer[-3:-1]

            // UP halves: distal = UP[0:12], proximal = UP[12:24]; groove width over non-overlapping dimers.
            double distal = Sum(GrooveAccess, seq, upStart, UpLength / 2);
            double proximal = Sum(GrooveAccess, seq, upStart + UpLength / 2, UpLength / 2);

            // ITR: DNA:DNA − RNA:DNA hybrid energy over the dimers of ITR[0:14].
            double dna = Sum(DnaDnaHybrid, seq, tss, 14);
            double rna = Sum(RnaDnaHybrid, seq, tss, 14);
            double hybrid = dna - rna;

            // Rigidity of UP + −35 + spacer[0:14] (44 nt), mean persistence length per nucleotide.
            const int rigidLength = UpLength + HexLength + 14;
            double rigidity = RigiditySum(seq, upStart, hex35Start, spacerStart) / rigidLength;

            double numDistal = distal / NormDistal * Numerical[0];
            double numProximal = proximal / NormProximal * Numerical[1];
            double dgItr = hybrid / NormItr * Numerical[2];
            double numRigidity = rigidity / NormRigidity * Numerical[3];

            double x = spacer;
            double dgSpacer = 0.1463 * (x * x) - 4.9113 * x + 41.119;
            double dgUp = numDistal + numProximal + numRigidity;
            double total = dg10 + dg35 + dgDisc + dgItr + dgExt10 + dgSpacer + dgUp + Intercept;
            double bind = dg10 + dg35 + dgSpacer + dgExt10 + dgUp;

            return new Sigma70PromoterPrediction(
                Strand: strand,
                Tss: tss,
                PromoterSequence: seq.Substring(upStart, tss + ItrLength - upStart),
                Up: seq.Substring(upStart, UpLength),
                Minus35: seq.Substring(hex35Start, HexLength),
                Spacer: seq.Substring(spacerStart, spacer),
                Minus10: seq.Substring(hex10Start, HexLength),
                Discriminator: seq.Substring(discStart, disc),
                Itr: seq.Substring(tss, ItrLength),
                UpStart: upStart,
                Minus35Start: hex35Start,
                SpacerStart: spacerStart,
                Minus10Start: hex10Start,
                DiscriminatorStart: discStart,
                DeltaGTotal: total,
                DeltaG10: dg10,
                DeltaG35: dg35,
                DeltaGDiscriminator: dgDisc,
                DeltaGItr: dgItr,
                DeltaGExtended10: dgExt10,
                DeltaGSpacer: dgSpacer,
                DeltaGUp: dgUp,
                DeltaGBind: bind,
                TranscriptionRate: K * Math.Exp(-beta * total));
        }

        /// <summary>Σ table[dimer] over the non-overlapping dimers of seq[start, start + length) (length even), left to right.</summary>
        private static double Sum(double[] table, string seq, int start, int length)
        {
            double s = 0;
            for (int k = 0; k < length; k += 2)
                s += table[Dimer(seq, start + k)];
            return s;
        }

        /// <summary>Σ persistence over the dimers of UP (24) + −35 (6) + spacer[0:14], accumulated left to right.</summary>
        private static double RigiditySum(string seq, int upStart, int hex35Start, int spacerStart)
        {
            double s = 0;
            for (int k = 0; k < UpLength; k += 2)
                s += Persistence[Dimer(seq, upStart + k)];
            for (int k = 0; k < HexLength; k += 2)
                s += Persistence[Dimer(seq, hex35Start + k)];
            for (int k = 0; k < 14; k += 2)
                s += Persistence[Dimer(seq, spacerStart + k)];
            return s;
        }

        // Row order A, C, G, T: the reference's one-hot categories are the sorted k-mers (AAA, AAC, …, TTT).
        private static int Dimer(string s, int i) => (AcgtIndex(s[i]) << 2) | AcgtIndex(s[i + 1]);

        private static int Kmer3(string s, int i) => (AcgtIndex(s[i]) << 4) | (AcgtIndex(s[i + 1]) << 2) | AcgtIndex(s[i + 2]);

        // util.py dimer tables in order AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT.
        private static readonly double[] GrooveAccess =
        {
            5.0, 4.0, 9.0, 0.0, 42.0, 42.0, 43.0, 9.0, 22.0, 25.0, 42.0, 4.0, 14.0, 22.0, 42.0, 5.0,
        };

        // Persistence lengths (Geggier & Vologodskii), util.persistence.
        private static readonly double[] Persistence =
        {
            50.4, 55.4, 51.0, 40.9, 46.7, 41.7, 56.0, 51.0, 54.4, 44.6, 41.7, 55.4, 44.7, 54.4, 46.7, 50.4,
        };

        private static readonly double[] DnaDnaHybrid =
        {
            -1.0, -1.44, -1.28, -0.88, -1.45, -1.42, -2.17, -1.28, -1.3, -2.24, -1.42, -1.44, -0.58, -1.3, -1.45, -1.0,
        };

        private static readonly double[] RnaDnaHybrid =
        {
            -0.2, -1.6, -1.5, -0.6, -1.1, -2.9, -2.7, -1.3, -0.9, -1.7, -2.1, -0.9, -0.9, -1.8, -2.1, -1.0,
        };

        // free_energy_coeffs.npy (343 values): [0,64) −10 first 3-mer, [64,128) −10 last 3-mer, [128,192) −35 first,
        // [192,256) −35 last, [256,320) discriminator 3-mer, [320,336) extended −10 dimer, [336,339) spacer one-hot
        // (16, 17, 18: −0.015096442164621445, −0.09696411810435467, 0.11206056026872108 — not used by the reference
        // linear_free_energy_model, which applies the quadratic spacer term), [339,343) numerical.
        private static readonly double[] Minus10First =
        {
            0.09490738746347842, 0.028821899543931544, 0.1198230376783691, -0.34434029969128693,
            0.0, 0.0, 0.0, 0.1638454358405594,
            0.0, 0.0, 0.0, 0.2147940064357312,
            0.0, 0.0, 0.0, -0.02760841459847772,
            -0.18785061878221526, -0.027935088893270586, 0.07065084457433997, -0.2999626574312279,
            0.17249598290925514, 0.22167955893809602, 0.0, 0.4295831188048784,
            0.046398931797263954, 0.0, 0.0, 0.42392556132156406,
            -0.09617404850628807, 0.0, 0.0, 0.32615138291787527,
            0.10505054947551914, -0.00040256483059495306, 0.1588895386113908, -0.2596949339348936,
            0.17726911147529376, 0.0, 0.0, 0.22783938527506442,
            0.0, 0.0, 0.0, 0.014025497244297485,
            -0.189873085272735, 0.0, 0.0, 0.23433309066242222,
            -0.7063684812690003, -0.6955332165785065, -0.5791109861977634, -0.9106392228465209,
            0.2131582968259988, 0.024649344853042127, 0.2571180467299443, 0.05923658893091885,
            0.06348749096954126, 0.021820191839442202, 0.037754073983779904, 0.03301394297027404,
            -0.16745326285599052, 0.24049535076215106, 0.03556826824357014, 0.2761609646098749,
        };

        private static readonly double[] Minus10Last =
        {
            -0.14042593692886762, -0.012384111923827271, 0.07062384545829954, -0.8708852278943782,
            -0.18378838543440634, 0.053183152602326395, 0.18086946257241995, -0.770785967815676,
            0.12034744752130364, 0.1644544995394996, -0.11797676076098046, -0.6032012647867945,
            -0.2815784985072322, 0.016760580206517762, -0.011108666369123185, -0.5839670229680539,
            0.37084711173422535, 0.2296544538769969, 0.2927349456965483, -0.5976286243429203,
            0.4160724843435361, 0.2770948968883442, 0.23833369595511988, -0.28201835991540614,
            0.0, 0.44130138079556164, 0.0, -0.22957944752798204,
            0.07217807522519586, 0.13545420990605467, 0.0, -0.21516031962995938,
            0.1703433358271497, 0.23430195379743315, 0.359782319300283, -0.5236823489952451,
            0.24985810790957888, 0.23599181501756106, 0.239697508426137, -0.3118833171984491,
            0.3054009417807217, 0.3641217865297511, 0.0, -0.2693116358365752,
            0.13490947717257368, 0.2726684004794281, 0.3607159076854683, -0.2319440535211394,
            0.06384993062886919, 0.18558906240429637, 0.0009394931280153995, -0.8503870687473457,
            0.3088283133466719, 0.30371182219728826, 0.2585312898601096, -0.389182428759729,
            0.12192168180402255, 0.42950855129388393, 0.0, -0.29282726925801883,
            0.07585053289082558, -0.035366737127009884, 0.30382604520695994, -0.255185064758909,
        };

        private static readonly double[] Minus35First =
        {
            0.0, 0.16540471136174353, 0.20399532646960403, 0.0,
            0.26202919727030116, 0.0, 0.1085144161441575, 0.0,
            0.0, 0.0, 0.11715359318273247, 0.0,
            0.23463667489026757, -0.02266078879538195, -0.42057592950756695, -0.008068733965697007,
            0.20522206154317535, 0.15203060003839258, 0.045544678282602453, 0.1500781457366836,
            0.1877433271981443, 0.2876153782860975, 0.022347560306883285, 0.1503878584281531,
            0.17716741797455585, 0.32537335947877183, 0.025433429323852983, -0.0038884919611999625,
            -0.1952143346110267, 0.1029676100034498, -0.24629712834653053, -0.11390935631664142,
            0.25840804202788786, 0.08580103150005455, 0.177621328096816, 0.0,
            0.4172176386838069, 0.252003184989742, -0.011394041424139419, 0.11882049592994029,
            0.0, 0.1897743936378747, 0.2356330850982128, 0.09817020013819935,
            -0.23849417242526405, -0.10039373057177248, -0.23783376052923252, -0.09254566213228003,
            -0.10057135366438301, 0.031594925820520314, -0.39096020475985893, -0.1518800857079465,
            0.09651004698673261, 0.0414044833258263, -0.13499031380149526, -0.2204862015355318,
            -0.12436571896176705, 0.0, -0.3332502512580711, -0.0556343867600926,
            -0.14150345445471604, -0.2610965805370817, -0.8511202866148896, -0.46946923351264314,
        };

        private static readonly double[] Minus35Last =
        {
            -0.16908705233531682, 0.043641084157395524, 0.0032413186806323083, -0.28364370837131575,
            -0.6260976284949558, -0.013027902113514216, -0.15311120173832768, -0.40719887288673423,
            -0.13692983107476275, 0.177032592979125, 0.0832860860945223, 0.13446032354401016,
            -0.05466223246930263, -0.19114400490445757, 0.03203870958479575, -0.14754411345414545,
            -0.22515898726967176, -0.00590248812297179, -0.05245668616935422, -0.06933635697725887,
            -0.23235464438478887, 0.04782476812796241, -0.013635747921151463, -0.3122991962691227,
            0.04535468125660317, 0.16364105235292595, 0.02279507617908226, -0.004260528433371471,
            -0.13420128185295316, 0.15232644152662514, 0.022378460778619363, -0.07610553119535049,
            -0.11320386467567833, 0.29318904946196794, 0.030984070093653684, 0.026543672719968405,
            -0.04432364947099053, 0.2215610511190525, -0.0032395179289541245, -0.0369869360824641,
            0.14785211603113513, 0.39444970100658056, 0.22783563226934417, 0.4579923412170347,
            -0.16901178129880748, 0.2823185034692884, 0.1206666836451072, 0.1054367402601531,
            -0.16892049808735207, 0.17037335610440976, -0.006665799606645034, -0.06804398201966647,
            -0.08322017053190033, 0.29621900270844964, 0.1303731072246929, -0.09236019035237614,
            0.07908733709760897, 0.14795245452341682, 0.0, -0.11774727064261974,
            0.05492373131928146, 0.21314019569262882, 0.0, -0.11703768409025755,
        };

        private static readonly double[] Discriminator =
        {
            -0.04487847301375619, 0.03794051524303827, 0.09488766696998044, -0.004857517483925856,
            0.02431586158030238, -0.11446373347665721, -0.020389208170715694, -0.1193189216919453,
            0.12376967412070047, 0.1658800645338869, 0.012795170885314655, -0.09222196801281188,
            -0.060202924695229726, 0.1590161608033341, 0.20286290801530318, -0.028406436170158426,
            0.09495745073130975, -0.14096411166041897, 0.13280770587783802, -0.06688516612075936,
            -0.07102722808182321, -0.10642768631174408, -0.10896914894976407, -0.15457117818056698,
            -0.1379660792586186, -0.16688377160145138, 0.06601460246286905, 0.1194498180205362,
            -0.07065392365182463, 0.024942086366553075, -0.11739483005842419, 0.2026342762130632,
            0.15203157529060746, -0.018881814329846067, 0.20430148488065647, 0.06678815959103351,
            -0.17513818723646626, -0.09485187672560347, -0.13940934892487852, -0.18001627903684744,
            -0.034709986895396124, 0.012178218133649998, 0.08046241181452651, -0.009420113359516894,
            0.1829665642586891, 0.06717619959178263, 0.023900602351315935, 0.16271416108557177,
            0.0702536687063409, 0.04647874230110192, 0.0001848609146443243, -0.1441943360313018,
            -0.10618594980487601, -0.1495728239641026, -0.12232781084798844, -0.16594764875320808,
            0.10780912735135333, -0.06396369387766375, 0.18496263427023718, -0.029965515348874584,
            0.09241381856502776, 0.0403155010901358, -0.047827908405309136, 0.15168390811051485,
        };

        private static readonly double[] Extended10 =
        {
            0.0610364204005321, 0.12238080735979254, -0.056176643129648894, -0.1710940037718414,
            -0.15950251583247507, 0.02482021652708394, 0.07186020160067112, 0.18814810202386284,
            0.11099622297853319, 0.13682442558153468, 0.027096935239089858, -0.07964961793745078,
            -0.139301362120528, 0.19135161731618874, -0.2144031377586821, -0.11438766847409058,
        };

        private static readonly double[] Numerical =
        {
            0.2546979750256038, 0.6505806574328081, 0.16210794028069989, 0.009426894866260856,
        };
    }

    #endregion
}

/// <summary>A paired σ70 −35 / −10 consensus candidate (<see cref="MotifFinder.FindSigma70Promoters"/>).</summary>
/// <param name="Strand">'+' or '-'.</param>
/// <param name="Minus35Start">0-based forward-strand start of the −35 hexamer.</param>
/// <param name="Minus35">−35 hexamer read 5'→3' on its strand.</param>
/// <param name="Minus10Start">0-based forward-strand start of the −10 hexamer.</param>
/// <param name="Minus10">−10 hexamer read 5'→3' on its strand.</param>
/// <param name="Spacer">Spacer length (bp between the hexamers).</param>
/// <param name="Mismatches35">Mismatches of the −35 box to TTGACA.</param>
/// <param name="Mismatches10">Mismatches of the −10 box to TATAAT.</param>
/// <param name="TotalMismatches">Mismatches35 + Mismatches10.</param>
/// <param name="SpacerDeviation">|Spacer − 17|.</param>
public readonly record struct Sigma70PromoterCandidate(
    char Strand,
    int Minus35Start,
    string Minus35,
    int Minus10Start,
    string Minus10,
    int Spacer,
    int Mismatches35,
    int Mismatches10,
    int TotalMismatches,
    int SpacerDeviation);

/// <summary>
/// The minimum-free-energy σ70 promoter configuration of one TSS (<see cref="MotifFinder.PredictSigma70Promoters"/>,
/// Promoter Calculator v1.0). Free energies in the model's units (kcal/mol-like, as the reference).
/// </summary>
/// <param name="Strand">'+' or '-'.</param>
/// <param name="Tss">TSS boundary index (reference convention; see <see cref="MotifFinder.PredictSigma70Promoters"/>).</param>
/// <param name="PromoterSequence">UP … ITR read 5'→3' on its strand.</param>
/// <param name="Up">UP element (24 nt).</param>
/// <param name="Minus35">−35 hexamer.</param>
/// <param name="Spacer">Spacer (15 … 20 nt).</param>
/// <param name="Minus10">−10 hexamer.</param>
/// <param name="Discriminator">Discriminator (6 … 10 nt).</param>
/// <param name="Itr">Initial transcribed region (20 nt from the TSS).</param>
/// <param name="UpStart">Forward-strand start of the UP element.</param>
/// <param name="Minus35Start">Forward-strand start of the −35 hexamer.</param>
/// <param name="SpacerStart">Forward-strand start of the spacer.</param>
/// <param name="Minus10Start">Forward-strand start of the −10 hexamer.</param>
/// <param name="DiscriminatorStart">Forward-strand start of the discriminator.</param>
/// <param name="DeltaGTotal">ΔG_total (the minimised quantity).</param>
/// <param name="DeltaG10">−10 hexamer term.</param>
/// <param name="DeltaG35">−35 hexamer term.</param>
/// <param name="DeltaGDiscriminator">Discriminator term.</param>
/// <param name="DeltaGItr">ITR (R-loop) term.</param>
/// <param name="DeltaGExtended10">Extended −10 (spacer[−3:−1]) term.</param>
/// <param name="DeltaGSpacer">Spacer-length term.</param>
/// <param name="DeltaGUp">UP element term.</param>
/// <param name="DeltaGBind">ΔG_−10 + ΔG_−35 + ΔG_spacer + ΔG_ext−10 + ΔG_UP (reference <c>dG_bind</c>).</param>
/// <param name="TranscriptionRate">K·exp(−β·ΔG_total) (reference <c>Tx_rate</c>).</param>
public sealed record Sigma70PromoterPrediction(
    char Strand,
    int Tss,
    string PromoterSequence,
    string Up,
    string Minus35,
    string Spacer,
    string Minus10,
    string Discriminator,
    string Itr,
    int UpStart,
    int Minus35Start,
    int SpacerStart,
    int Minus10Start,
    int DiscriminatorStart,
    double DeltaGTotal,
    double DeltaG10,
    double DeltaG35,
    double DeltaGDiscriminator,
    double DeltaGItr,
    double DeltaGExtended10,
    double DeltaGSpacer,
    double DeltaGUp,
    double DeltaGBind,
    double TranscriptionRate);
