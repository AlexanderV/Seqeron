/* Seqeron review harness (A2-29b): runs the published ksw2 / minimap2 extension functions on the cases of cases.txt and prints every field of ksw_extz_t.
 * Input line: id variant flag w zdrop end_bonus q e q2 e2 m mat[m*m]... query target   (sequences are digit strings 0..m-1, "-" = empty)
 * Output line: id variant max zdropped max_q max_t mqe mqe_t mte mte_q score n_cigar reach_end cigar
 * variants: extz2 extd2 (SSE, minimap2 ksw2_extz2_sse.c / ksw2_extd2_sse.c) and, when compiled with -DSCALAR, extz extd (lh3/ksw2 ksw2_extz.c / ksw2_extd.c). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ksw2.h"

int main(int argc, char **argv)
{
    int start = argc > 1 ? atoi(argv[1]) : 0; /* skip the cases with id < start: used to resume after a crash of the program under test */
    setvbuf(stdout, NULL, _IOLBF, 0);
    static char variant[16], qs[1000000], ts[1000000];
    int id, flag, w, zdrop, eb, q, e, q2, e2, m;
    int8_t mat[400];
    while (scanf("%d %15s %d %d %d %d %d %d %d %d %d", &id, variant, &flag, &w, &zdrop, &eb, &q, &e, &q2, &e2, &m) == 11) {
        for (int i = 0; i < m * m; ++i) { int v; if (scanf("%d", &v) != 1) return 2; mat[i] = (int8_t)v; }
        if (scanf("%999999s %999999s", qs, ts) != 2) return 2;
        if (id < start) continue;
        fprintf(stderr, "CASE %d %s\n", id, variant); fflush(stderr); /* the last CASE line of a crashed run names the case */
        int ql = strcmp(qs, "-") ? (int)strlen(qs) : 0, tl = strcmp(ts, "-") ? (int)strlen(ts) : 0;
        uint8_t *qq = malloc(ql + 1), *tt = malloc(tl + 1);
        for (int i = 0; i < ql; ++i) qq[i] = qs[i] - '0';
        for (int i = 0; i < tl; ++i) tt[i] = ts[i] - '0';
        ksw_extz_t ez; memset(&ez, 0, sizeof ez);
        int known = 1;
#ifndef SCALAR
        if (!strcmp(variant, "extz2")) ksw_extz2_sse(0, ql, qq, tl, tt, (int8_t)m, mat, (int8_t)q, (int8_t)e, w, zdrop, eb, flag, &ez);
        else if (!strcmp(variant, "extd2")) ksw_extd2_sse(0, ql, qq, tl, tt, (int8_t)m, mat, (int8_t)q, (int8_t)e, (int8_t)q2, (int8_t)e2, w, zdrop, eb, flag, &ez);
        else known = 0;
#else
        if (!strcmp(variant, "extz")) ksw_extz(0, ql, qq, tl, tt, (int8_t)m, mat, (int8_t)q, (int8_t)e, w, zdrop, flag, &ez);
        else if (!strcmp(variant, "extd")) ksw_extd(0, ql, qq, tl, tt, (int8_t)m, mat, (int8_t)q, (int8_t)e, (int8_t)q2, (int8_t)e2, w, zdrop, flag, &ez);
        else known = 0;
#endif
        if (known) {
            printf("%d %s %u %d %d %d %d %d %d %d %d %d %d ", id, variant, (unsigned)ez.max, ez.zdropped, ez.max_q, ez.max_t, ez.mqe, ez.mqe_t, ez.mte, ez.mte_q, ez.score, ez.n_cigar, ez.reach_end);
            if (ez.n_cigar == 0) printf("*");
            for (int i = 0; i < ez.n_cigar; ++i) printf("%d%c", ez.cigar[i] >> 4, "MIDN"[ez.cigar[i] & 0xf]);
            printf("\n");
        }
        fflush(stdout);
        free(qq); free(tt); if (ez.cigar) free(ez.cigar);
    }
    return 0;
}
