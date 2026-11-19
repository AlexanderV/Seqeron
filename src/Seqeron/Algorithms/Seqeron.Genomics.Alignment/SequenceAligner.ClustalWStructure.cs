using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    // ClustalW 2.1 secondary-structure gap-penalty masks of the profile commands (audit A2-14c, B09 F106): the masks read from the profile files (`ClustalFileParser::getSecStructure`, FileReader::profileInput),
    // Clustal::calcGapPenaltyMask (structure -> gap-penalty digits), ProfileBase::calcGapCoeff (the digit multiplies the gap penalties and switches the variable-region / hydrophilic / residue-specific
    // penalties off), MyersMillerProfileAlign::addGGapsMask (the masks follow the alignment), AlignmentOutput (structure lines) of the clustalw 2.1+lgpl-7 Debian source.

    private sealed partial class ClustalWMsaState
    {
        /// <summary>The structure parameters (<c>-HELIXGAP</c> ... <c>-NOSECSTR2</c>).</summary>
        public ClustalWStructureOptions StructureOptions { get; init; } = new();

        /// <summary>The structure line / mask of profile 1 (null = none).</summary>
        public ClustalWProfileStructure? Structure1Input { get; init; }

        /// <summary>The structure line / mask of profile 2.</summary>
        public ClustalWProfileStructure? Structure2Input { get; init; }

        private int _structPenalties1, _structPenalties2;   // UserParameters::structPenalties1/2: NONE 0, SECST 1, GMASK 2
        private bool _useSS1 = true, _useSS2 = true;
        private List<char> _secStruct1 = [], _gapMask1 = [], _secStruct2 = [], _gapMask2 = [];
        private string _structName1 = "", _structName2 = "";

        /// <summary>Profile 1 structure after the run (null when the profile had none).</summary>
        public ClustalWAlignmentStructure? FinalStructure1 => _structPenalties1 == 0 ? null : Snapshot(_structPenalties1, _structName1, _secStruct1, _gapMask1, _useSS1);

        /// <summary>Profile 2 structure after the run.</summary>
        public ClustalWAlignmentStructure? FinalStructure2 => _structPenalties2 == 0 ? null : Snapshot(_structPenalties2, _structName2, _secStruct2, _gapMask2, _useSS2);

        private static ClustalWAlignmentStructure Snapshot(int kind, string name, List<char> secStruct, List<char> gapMask, bool use) =>
            new(name, kind == 1, kind == 1 ? new string(secStruct.ToArray()) : "", new string(gapMask.ToArray()), use);

        // FileReader::profileInput + Clustal::profileAlign / sequencesAlignToProfile: the masks of the two profile files and their gap-penalty masks (calcGapPenaltyMask)
        private void InitStructure(int profile1NumSeqs, bool sequencesMode)
        {
            _useSS1 = StructureOptions.UseStructure1;
            _useSS2 = StructureOptions.UseStructure2;
            int len1 = _seq.Take(profile1NumSeqs).Max(r => r.Length);
            int lenAll = _seq.Max(r => r.Length);   // Alignment::getMaxAlnLength() when profile 2 is read
            (_structPenalties1, _structName1, _secStruct1, _gapMask1) = ReadStructure(Structure1Input, len1);
            (_structPenalties2, _structName2, _secStruct2, _gapMask2) = ReadStructure(Structure2Input, lenAll);

            // Clustal::sequencesAlignToProfile: "ignoring secondary structure for a list of sequences"
            if (sequencesMode && _structPenalties2 != 0 && _useSS2 && _n - profile1NumSeqs > 1)
                _useSS2 = false;

            if (_structPenalties1 == 1)
                CalcGapPenaltyMask(_seq[0].Length, _secStruct1, _gapMask1);
            if (_structPenalties2 == 1)
                CalcGapPenaltyMask(_seq[profile1NumSeqs].Length, _secStruct2, _gapMask2);
        }

        /// <summary>The state of <c>userParameters-&gt;useSS2</c> after <c>sequencesAlignToProfile</c> restores it (the structure output looks at the restored flag).</summary>
        private void RestoreStructureFlags()
        {
            _useSS1 = StructureOptions.UseStructure1;
            _useSS2 = StructureOptions.UseStructure2;
        }

        private static (int Kind, string Name, List<char> Secondary, List<char> Mask) ReadStructure(ClustalWProfileStructure? input, int length)
        {
            if (input == null)
                return (0, "", [], []);
            string name = input.Name.Replace(' ', '_');
            if (input.SecondaryStructure != null)
            {
                var secondary = Enumerable.Repeat('.', length).ToList();
                int len = 0;
                foreach (char c in input.SecondaryStructure)
                {
                    if (len >= length)
                        break;
                    secondary[len++] = c;
                }

                return (1, name, secondary, Enumerable.Repeat('.', length).ToList());
            }

            if (input.GapPenaltyMask != null)
            {
                var mask = Enumerable.Repeat('1', length).ToList();
                int len = 0;
                foreach (char c in input.GapPenaltyMask)
                {
                    if (len >= length)
                        break;
                    mask[len++] = c;
                }

                return (2, name, [], mask);
            }

            return (0, "", [], []);
        }

        // Clustal::calcGapPenaltyMask: the gap-penalty digit of every column from the structure line
        private void CalcGapPenaltyMask(int prfLength, List<char> mask, List<char> gapMask)
        {
            var structMask = new char[prfLength + 1];
            int helixEndPlus = StructureOptions.HelixEndOutside, helixEndMinus = StructureOptions.HelixEndInside;
            int strandEndPlus = StructureOptions.StrandEndOutside, strandEndMinus = StructureOptions.StrandEndInside;
            char M(int idx) => idx >= 0 && idx < mask.Count ? mask[idx] : '\0';
            static char Lower(char c) => char.ToLowerInvariant(c);
            int i = 0;
            while (i < prfLength)
            {
                if (Lower(M(i)) == 'a' || M(i) == '$')
                {
                    for (int j = -helixEndPlus; j < 0; j++)
                    {
                        if (i + j < structMask.Length && i + j >= 0 && Lower(structMask[i + j]) != 'a' && Lower(structMask[i + j]) != 'b')
                            structMask[i + j] = 'a';
                    }

                    int jj;
                    for (jj = 0; jj < helixEndMinus; jj++)
                    {
                        if (i + jj >= prfLength || (Lower(M(i + jj)) != 'a' && M(i + jj) != '$'))
                            break;
                        structMask[i + jj] = 'a';
                    }

                    i += jj;
                    while (Lower(M(i)) == 'a' || M(i) == '$')
                    {
                        if (i >= prfLength)
                            break;
                        if (M(i) == '$')
                        {
                            structMask[i] = 'A';
                            i++;
                            break;
                        }

                        structMask[i] = M(i);
                        i++;
                    }

                    for (int k = 0; k < helixEndMinus; k++)
                    {
                        if (i - k - 1 >= 0 && (Lower(M(i - k - 1)) == 'a' || M(i - k - 1) == '$'))
                            structMask[i - k - 1] = 'a';
                    }

                    for (int k = 0; k < helixEndPlus; k++)
                    {
                        if (i + k >= prfLength)
                            break;
                        structMask[i + k] = 'a';
                    }
                }
                else if (Lower(M(i)) == 'b' || M(i) == '%')
                {
                    for (int j = -strandEndPlus; j < 0; j++)
                    {
                        if (i + j >= 0 && Lower(structMask[i + j]) != 'a' && Lower(structMask[i + j]) != 'b')
                            structMask[i + j] = 'b';
                    }

                    int jj;
                    for (jj = 0; jj < strandEndPlus; jj++)
                    {
                        if (i + jj >= prfLength || (Lower(M(i + jj)) != 'b' && M(i + jj) != '%'))
                            break;
                        structMask[i + jj] = 'b';
                    }

                    i += jj;
                    while (Lower(M(i)) == 'b' || M(i) == '%')
                    {
                        if (i >= prfLength)
                            break;
                        if (M(i) == '%')
                        {
                            structMask[i] = 'B';
                            i++;
                            break;
                        }

                        structMask[i] = M(i);
                        i++;
                    }

                    for (int k = 0; k < strandEndMinus; k++)
                    {
                        if (i - k - 1 >= 0 && (Lower(M(i - k - 1)) == 'b' || M(i - k - 1) == '%'))
                            structMask[i - k - 1] = 'b';
                    }

                    for (int k = 0; k < strandEndPlus; k++)
                    {
                        if (i + k >= prfLength)
                            break;
                        structMask[i + k] = 'b';
                    }
                }
                else
                {
                    i++;
                }
            }

            for (i = 0; i < prfLength; i++)
            {
                int penalty = structMask[i] switch
                {
                    'A' => StructureOptions.HelixPenalty,
                    'a' => StructureOptions.TerminalPenalty,
                    'B' => StructureOptions.StrandPenalty,
                    'b' => StructureOptions.TerminalPenalty,
                    _ => StructureOptions.LoopPenalty,
                };
                while (gapMask.Count <= i)
                    gapMask.Add('\0');
                gapMask[i] = (char)(penalty + '0');
            }
        }

        // MyersMillerProfileAlign::addGGapsMask: the mask follows its profile: a residue column keeps its mask character, a gap column gets the gap character
        private static void AddGapsMask(List<char> mask, int len, int[] path1, int[] path2, bool switchProfiles)
        {
            var path = switchProfiles ? path2 : path1;
            var ta = new char[len];
            int ix = 0;
            for (int i = 0; i < len; i++)
            {
                if (path[i] == 2)
                {
                    ta[i] = ix < mask.Count ? mask[ix] : '\0';
                    ix++;
                }
                else if (path[i] == 1)
                {
                    ta[i] = '-';
                }
            }

            // mask->resize(len + 2): the two entries behind the alignment are zero
            while (mask.Count < len + 2)
                mask.Add('\0');
            if (mask.Count > len + 2)
                mask.RemoveRange(len + 2, mask.Count - (len + 2));
            for (int i = 0; i < len; i++)
                mask[i] = ta[i];
        }
    }
}
