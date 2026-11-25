# Third-party notices

Seqeron is distributed under the MIT licence (see [LICENSE](LICENSE)). Parts of the MIT-licensed assemblies were written by reading the source code of the programs listed here (ports, byte-for-byte checked against the installed programs) or embed their data. Each of these programs allows that under the licence reproduced below; the notices are retained here as those licences require.

The ports of GPL / LGPL programs are NOT part of the MIT assemblies: they live in separate assemblies with their own LICENSE files (see the section "Licences and assemblies" of [README.md](README.md) and [docs/Validation/review-2026-09/LICENSES.md](docs/Validation/review-2026-09/LICENSES.md), which also holds the full inventory, the black-box oracles and the items that need an owner decision).

---

## MAFFT 7.505

- Used for: Port of the alignment programs of MAFFT (src/mafft: tbfast, dvtditr, pairlocalalign, ... and the `mafft` script) in Seqeron.Genomics.Alignment (SequenceAligner.Mafft*.cs, MafftProgram*.cs).
- Copyright: Copyright 2009-2017 Kazutaka Katoh (Debian: Source http://mafft.cbrc.jp/alignment/software/source.html)
- Licence: MAFFT licence (BSD-style, 3 conditions)

    Redistribution and use in source and binary forms,
    with or without modification, are permitted provided
    that the following conditions are met:

    Redistributions of source code must retain the
    above copyright notice, this list of conditions
    and the following disclaimer.  Redistributions in
    binary form must reproduce the above copyright
    notice, this list of conditions and the following
    disclaimer in the documentation and/or other
    materials provided with the distribution.

    The name of the author may not be used to endorse
    or promote products derived from this software without
    specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
    OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
    OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
    THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
    OF SUCH DAMAGE.

The conditions: source redistributions retain the notice, list of conditions and disclaimer; binary redistributions reproduce them in the documentation; the name of the author is not used for endorsement. This file and the header comments are the reproduction of the notice.

---

## MUSCLE 3.8.1551

- Used for: Port of MUSCLE 3.8.1551 (profile score, progressive / refinement scheme, command line) in Seqeron.Genomics.Alignment (SequenceAligner.Muscle*.cs).
- Copyright: Robert C. Edgar (drive5.com)
- Licence: Public domain dedication

    MUSCLE is public domain software
    The MUSCLE software, including object and source code, is hereby donated
    to the public domain.

    Disclaimer of warranty
    THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
    EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

---

## SPOA 4.0.8

- Used for: Port of the partial-order-alignment program SPOA (SpoaEngine.cs, SpoaGraph.cs, SpoaProgram.cs).
- Copyright: Copyright (c) 2016 Robert Vaser
- Licence: MIT License

    The MIT License (MIT)

    Copyright (c) 2016 Robert Vaser

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

---

## minimap2 2.26 and ksw2

- Used for: Port of minimap2 seeding / chaining / base-level alignment and of the ksw2 kernels (Minimap2*.cs, Ksw2Aligner*.cs, SequenceAligner.DualAffine.cs).
- Copyright: Copyright (c) 2018- Dana-Farber Cancer Institute; 2017-2018 Broad Institute, Inc.; 2015-2018 Heng Li (as stated in the upstream LICENSE.txt)
- Licence: MIT License

    The MIT License

    Copyright (c) 2018-     Dana-Farber Cancer Institute
                  2017-2018 Broad Institute, Inc.

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

---

## Mash 2.3

- Used for: Port of the Mash commands sketch / dist / triangle / bounds / info / paste / screen / taxscreen (Seqeron.Genomics.Analysis: Mash*.cs, KmerAnalyzer minhash).
- Copyright: Copyright (c) 2015 Battelle National Biodefense Institute (BNBI)
- Licence: BSD 3-Clause License (as in the Debian copyright file of mash 2.3+dfsg-6)

    This Software was prepared for the Department of Homeland Security
    (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
    part of contract HSHQDC-07-C-00020 to manage and operate the National
    Biodefense Analysis and Countermeasures Center (NBACC), a Federally
    Funded Research and Development Center.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Mash also distributes kseq.h (MIT, Genome Research Ltd), MurmurHash3 (public domain), Cap'n Proto and robin_hood (MIT). The Open Bloom Filter of Mash is under the Common Public License 1.0 and is listed under NEEDS OWNER DECISION in docs/Validation/review-2026-09/LICENSES.md.

---

## Biopython 1.88

- Used for: Bio.Align / PairwiseAligner behaviour, the alignment readers and writers (Bio.Align formats), Bio.motifs, Restriction data (Restriction_Dictionary REBASE tables) and the codon / GC-skew conventions (Seqeron.Genomics.Alignment, Analysis, MolTools).
- Copyright: Copyright 1999-2024 by the Biopython contributors
- Licence: Biopython License Agreement

    Biopython is currently released under the "Biopython License Agreement"
    (given in full below). Unless stated otherwise in individual file headers,
    all Biopython's files are under the "Biopython License Agreement".

    Some files are explicitly dual licensed under your choice of the
    "Biopython License Agreement" or the "BSD 3-Clause License" (both given
    in full below). This is with the intention of later offering all of
    Biopython under this dual licensing approach.


    ------------------------------------------------------------------------------


    Biopython License Agreement
    ---------------------------

    Permission to use, copy, modify, and distribute this software and its
    documentation with or without modifications and for any purpose and
    without fee is hereby granted, provided that any copyright notices
    appear in all copies and that both those copyright notices and this
    permission notice appear in supporting documentation, and that the
    names of the contributors or copyright holders not be used in
    advertising or publicity pertaining to distribution of the software
    without specific prior permission.

    THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
    WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
    CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
    OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
    OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
    OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
    OR PERFORMANCE OF THIS SOFTWARE.


    ------------------------------------------------------------------------------


    BSD 3-Clause License
    --------------------

    Copyright (c) 1999-2025, The Biopython Contributors
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors
       may be used to endorse or promote products derived from this software
       without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## parasail 1.3.4

- Used for: Score API and variants of parasail (ParasailAligner.cs; behaviour reproduced and checked against parasail 1.3.4).
- Copyright: Copyright (c) 2015, Battelle Memorial Institute
- Licence: BSD 3-Clause License (parasail)

    404: Not Found

    Pairwise Sequence Alignment Library (parasail)

    Copyright (c) 2015, Battelle Memorial Institute

    1.  Battelle Memorial Institute (hereinafter Battelle) hereby grants
        permission to any person or entity lawfully obtaining a copy of this
        software and associated documentation files (hereinafter “the
        Software”) to redistribute and use the Software in source and binary
        forms, with or without modification.  Such person or entity may use,
        copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and may permit others to do so, subject to
        the following conditions:

        - Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimers. 

        - Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in
          the documentation and/or other materials provided with the
          distribution. 

        - Other than as used herein, neither the name Battelle Memorial
          Institute or Battelle may be used in any form whatsoever without
          the express written consent of Battelle.  

        - Redistributions of the software in any form, and publications
          based on work performed using the software should include the
          following citation as a reference:

        Daily, Jeff. (2016). Parasail: SIMD C library for global,
        semi-global, and local pairwise sequence alignments. *BMC
        Bioinformatics*, 17(1), 1-11.  doi:10.1186/s12859-016-0930-z

    2.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
        "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
        LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
        FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE
        OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
        USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
        OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
        SUCH DAMAGE.

---

## ncbi-blast+ / NCBI C++ Toolkit (DUST, ORFfinder, BLAST statistics)

- Used for: Behaviour of dustmasker (SequenceComplexity.Dustmasker.cs), the ORFfinder length rules and the BLAST statistical parameters (Karlin-Altschul tables).
- Copyright: National Center for Biotechnology Information (United States Government work)
- Licence: Public domain notice

                                     PUBLIC DOMAIN NOTICE
                        National Center for Biotechnology Information

           This software/database is a "United States Government Work" under the
           terms of the United States Copyright Act.  It was written as part of
           the author's official duties as a United States Government employee and
           thus cannot be copyrighted.  This software/database is freely available
           to the public for use. The National Library of Medicine and the U.S.
           Government have not placed any restriction on its use or reproduction.

           Although all reasonable efforts have been taken to ensure the accuracy
           and reliability of the software and data, the NLM and the U.S.
           Government do not and cannot warrant the performance or results that
           may be obtained by using this software or data. The NLM and the U.S.
           Government disclaim all warranties, express or implied, including
           warranties of performance, merchantability or fitness for any particular
           purpose.

           Please cite the author in any work or product based on this material.

---

## maxentpy (MaxEntScan tables)

- Used for: The embedded maximum-entropy splice-site tables (Seqeron.Genomics.Annotation Data/maxent_score3.txt / score5 factorisation).
- Copyright: Copyright (c) 2016 kepbod (maxentpy); Yeo & Burge 2004 (model)
- Licence: MIT License

    MIT License

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

## CNVkit

- Used for: CopyNumberMath.cs: the log2-ratio to absolute copy-number conversion of CNVkit cnvlib/call.py.
- Copyright: Copyright (c) 2013-2018 Eric Talevich; Copyright (c) 2013-2018 University of California
- Licence: Apache License 2.0

    Copyright (c) 2013-2018 Eric Talevich
    Copyright (c) 2013-2018 University of California

    The full text of the Apache License 2.0 is in docs/Validation/review-2026-09/LICENSES.md (appendix) and at https://www.apache.org/licenses/LICENSE-2.0 . Changes: the code is translated to C# and extended (the file header states the origin).

---

## MHCflurry 2.0

- Used for: MhcflurryAffinityPredictor.cs and the embedded allele pseudosequence table / network weights (Seqeron.Genomics.Oncology).
- Copyright: Copyright OpenVax (see the upstream LICENSE)
- Licence: Apache License 2.0

    See the Apache License 2.0 text in docs/Validation/review-2026-09/LICENSES.md (appendix); upstream: https://github.com/openvax/mhcflurry/blob/master/LICENSE

Changes: translated to C#; weights are shipped as an embedded resource.

---

## MobiDB-lite (thresholds)

- Used for: DisorderPredictor.cs: the consensus thresholds of MobiDB-lite v3.
- Copyright: Copyright BioComputingUP (University of Padova)
- Licence: Apache License 2.0

    See the Apache License 2.0 text in docs/Validation/review-2026-09/LICENSES.md (appendix); upstream: https://github.com/BioComputingUP/MobiDB-lite/blob/master/LICENSE

Changes: thresholds re-expressed in C#.

---

## LICHeE

- Used for: OncologyAnalyzer.PhylogenyHeterogeneity.cs: constraint network construction ported from LICHeE PHYNetwork.
- Copyright: Copyright (c) 2014 Victoria Popic
- Licence: MIT License

    Program LICHeE for multi-sample cancer phylogeny reconstruction
    by Victoria Popic (viq@stanford.edu) 2014
     *
    MIT License
     *
    Copyright (c) 2014 Victoria Popic.
    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:
     *
    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.
     *
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

---

## Pfam HMM profiles

- Used for: Embedded profile HMMs used by the protein-domain and metagenome marker units.
- Copyright: Pfam (EMBL-EBI / InterPro)
- Licence: CC0 1.0 (public domain dedication)

    Pfam is made available under the Creative Commons Zero (CC0) 1.0 Universal public domain dedication. https://creativecommons.org/publicdomain/zero/1.0/

---

## ABIS-Seq signature matrix (Monaco et al. 2019)

- Used for: Resources/ABIS_sigmatrixRNAseq.tsv (immune deconvolution).
- Copyright: (c) 2019 The Authors (Monaco et al., Cell Reports 26:1627)
- Licence: Creative Commons Attribution 4.0 (CC BY 4.0)

    Attribution: Monaco G, Lee B, Xu W, et al. RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell Reports 26(6):1627-1640 (2019). https://creativecommons.org/licenses/by/4.0/ ; the matrix is reproduced unchanged except for the file format.

