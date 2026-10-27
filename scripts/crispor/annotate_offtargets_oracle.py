#!/usr/bin/env python3
"""Runs CRISPOR's own `annotateOfftargets` on a hand-built hit list (oracle for B08 H4/H5).

Usage: annotate_offtargets_oracle.py CRISPOR_PY CFD_MISMATCH_JSON CFD_PAM_JSON
  CRISPOR_PY          crispor.py (sha256 16c3a54b454a2739... for the copy used in B08)
  CFD_*_JSON          the CFD matrices (mismatch_score.pkl / pam_scores.pkl of CRISPOR's CFD_Scoring) decoded WITHOUT
                      executing the pickles (pickletools.dis stack replay); the pickles are never loaded here.
Prints guideScore / guideCfdScore (integer-rounded, as the reference) for N perfect sites, no on-target interval
("auto-ontarget"): 100/50/20/2/0 and 100/99/96/67/17 for N = 1, 2, 5, 50, 500; and the repeat rule (totalAlnCount > MAXOCC
-> 0 / 0, empty list).
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import crispor_exec

def make(crispor_py, mm_json, pam_json):
    ns = crispor_exec.load(crispor_py)
    mm = json.load(open(mm_json)); pam = json.load(open(pam_json))
    ns['get_mm_pam_scores'] = lambda: (mm, pam)
    ns['makeAlnStr'] = lambda org, g, o, p, mit, cfd, pos, dist: ('', False, False)   # HTML only
    ns['isInPar'] = lambda db, chrom, start, end: None
    def annotate(guide, hits, pam='NGG', input_pos=''):
        """hits: [(editDist, otSeq with PAM, totalAlnCount, isRep)] in the order the aligner lists them."""
        cd = {}
        for i, (d, seq, tac, rep) in enumerate(hits):
            cd.setdefault(d, []).append(('chr1', 1000 + i * 100, 1000 + i * 100 + 23, seq, '+', '', 'g', tac, rep))
        pos, ot_desc, gs, gcfd, l12, ontarget, rep_count = ns['annotateOfftargets']('hg19', cd, guide, pam, input_pos)
        return dict(guideScore=gs, guideCfd=gcfd, offTargets=len(pos), repCount=rep_count)
    return ns, annotate

if __name__ == '__main__':
    ns, annotate = make(*sys.argv[1:4])
    g = 'GGCCGACCTGTCGCTGACGCAGG'
    for n in (1, 2, 5, 50, 500):
        print(n, 'perfect sites ->', annotate(g, [(0, g, n, False)] * n))
    print('repeat (totalAlnCount = MAXOCC + 1) ->', annotate(g, [(0, g, ns['MAXOCC'] + 1, False)] * 2))
