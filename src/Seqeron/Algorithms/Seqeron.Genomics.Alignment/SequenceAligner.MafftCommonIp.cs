namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: commonIP, the trace-back matrices that every dynamic program of a process shares (defs.c, Salignmm.c, Galign11.c, Lalign11.c, genalign11.c, partSalignmm.c)

    /// <summary>
    /// <c>commonIP</c>, <c>commonJP</c>, <c>commonAlloc1</c>, <c>commonAlloc2</c> of <c>defs.c</c>: the programs keep ONE int matrix (and for <c>genL__align11</c> a second one) for the trace-back of every dynamic program
    /// (<c>A__align</c>, <c>G__align11</c>, <c>L__align11</c>, <c>genL__align11</c>, <c>partA__align</c>); it is allocated zeroed when a program needs more rows or columns than it has (<c>orlgth1 &gt; commonAlloc1 || orlgth2 &gt; commonAlloc2</c>:
    /// <c>AllocateIntMtx( max( orlgth1, commonAlloc1 ) + 10, max( orlgth2, commonAlloc2 ) + 10 )</c>) and is NOT cleared between calls. A cell that a dynamic program does not write and the trace-back reads therefore holds what an earlier call
    /// left there: with <c>tailgp</c> 0 the cell <c>ijp[lgth1][lgth2]</c> is only written when a finite score beats the start value of the end search, so a NaN or infinite score (a NaN branch length of a guide tree) makes the alignment depend on
    /// the earlier calls of the process. One instance is shared by the programs that run in one process (<c>tbfast</c> with its in-process <c>pairlocalalign</c>).
    /// </summary>
    internal sealed class CommonIpState
    {
        internal int[][]? Ip;
        internal int[][]? Jp;
        internal int Alloc1, Alloc2;

        private static int[][] AllocateIntMtx(int ll1, int ll2)
        {
            var mtx = new int[ll1][];
            for (int i = 0; i < ll1; i++)
                mtx[i] = new int[ll2];
            return mtx;
        }

        /// <summary>The block <c>if( orlgth1 &gt; commonAlloc1 || orlgth2 &gt; commonAlloc2 ) { ... } ijp = commonIP;</c> of the dynamic programs.</summary>
        internal int[][] IpFor(int orlgth1, int orlgth2, bool withJp = false)
        {
            if (orlgth1 > Alloc1 || orlgth2 > Alloc2)
            {
                int ll1 = Math.Max(orlgth1, Alloc1);
                int ll2 = Math.Max(orlgth2, Alloc2);
                Ip = AllocateIntMtx(ll1 + 10, ll2 + 10);
                if (withJp)
                    Jp = AllocateIntMtx(ll1 + 10, ll2 + 10);
                Alloc1 = ll1;
                Alloc2 = ll2;
            }
            return Ip!;
        }

        /// <summary><c>ijpj = commonJP</c> of <c>genL__align11</c> (the program would read a matrix that does not exist, or is too small, when another program had allocated <c>commonIP</c> alone first: allocated here).</summary>
        internal int[][] JpFor(int lgth1, int lgth2)
        {
            if (Jp == null || Jp.Length < lgth1 + 1 || Jp[0].Length < lgth2 + 1)
                Jp = AllocateIntMtx(Math.Max(Alloc1, lgth1) + 10, Math.Max(Alloc2, lgth2) + 10);
            return Jp;
        }

        /// <summary>The switch to the memory-saving dynamic program (<c>FreeIntMtx( commonIP ); commonIP = NULL; commonAlloc1 = commonAlloc2 = 0;</c>).</summary>
        internal void FreeIp()
        {
            Ip = null;
            Alloc1 = 0;
            Alloc2 = 0;
        }
    }

    /// <summary>The static capacity of one dynamic program (<c>static int orlgth1, orlgth2</c>): grown to 1.3 times the longest group seen, reset by the call with NULL arguments that the programs make at the end of a stage.</summary>
    internal sealed class DpCapacity
    {
        internal int Orlgth1, Orlgth2;

        /// <summary><c>if( lgth1 &gt; orlgth1 || lgth2 &gt; orlgth2 ) { ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100; ...; orlgth1 = ll1 - 100; }</c></summary>
        internal void Fit(int lgth1, int lgth2)
        {
            if (lgth1 > Orlgth1 || lgth2 > Orlgth2)
            {
                int ll1 = Math.Max((int)(1.3 * lgth1), Orlgth1) + 100;
                int ll2 = Math.Max((int)(1.3 * lgth2), Orlgth2) + 100;
                Orlgth1 = ll1 - 100;
                Orlgth2 = ll2 - 100;
            }
        }

        internal void Reset()
        {
            Orlgth1 = 0;
            Orlgth2 = 0;
        }
    }

    internal sealed partial class MafftGlobals
    {
        /// <summary>The trace-back matrices of the process (see <see cref="CommonIpState"/>).</summary>
        internal CommonIpState Common = new();
    }

    internal sealed partial class MafftRun
    {
        /// <summary>The trace-back read a cell of <c>commonIP</c> that the dynamic program did not write (the scores are NaN or infinite) and that an earlier program left with a value that points outside the matrix: the program stops with a segmentation fault (a message with a process id) or damages its heap.</summary>
        private static NotSupportedException StaleTraceBack() => UndefinedBehaviour("the trace-back of a dynamic program with NaN or infinite scores reads a cell of commonIP that an earlier call left behind and follows it out of the matrix (segmentation fault, process id in the message)");

        private readonly DpCapacity _gCap = new(), _lCap = new(), _genCap = new();

        /// <summary><c>true</c> for the program <c>tbfast</c> (it does not call <c>G__align11( NULL ... )</c> at the end of its progressive stage).</summary>
        private bool _isTbfast;
    }

    #endregion
}
