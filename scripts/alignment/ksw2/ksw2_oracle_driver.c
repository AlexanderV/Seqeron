/* A3-9c2 (B09 WP45): driver of the ksw2 oracle. Compile with the ksw2 sources of minimap2 v2.26 (ksw2.h, ksw2_ll_sse.c, ksw2_exts2_sse.c) with -msse4.1 (the Debian build of the program uses the SSE4.1 code path):
 *   gcc -O1 -msse4.1 -o ksw2_oracle_driver ksw2_oracle_driver.c ksw2_ll_sse.c ksw2_exts2_sse.c
 * Input lines (stdin):
 *   LL  query target gapo gape                                    -> "score qe te"                (ksw_ll_qinit size 2 + ksw_ll_i16)
 *   S2  query target q e q2 noncan zdrop junc_bonus flag junc     -> "max zdropped max_q max_t mqe mqe_t mte mte_q score reach_end n_cigar cigar" (ksw_exts2_sse; junc = a string of digits 0-9 or "-")
 *   Z2  query target q e w zdrop end_bonus flag                   -> same fields (ksw_extz2_sse)
 *   D2  query target q e q2 e2 w zdrop end_bonus flag             -> same fields (ksw_extd2_sse)
 * query / target over ACGTN; the matrix is ksw_gen_simple_mat(5, a=2, b=4, N=1) of minimap2 (mat[i][i] = 2, mismatch -4, N -1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ksw2.h"

static int enc(char c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 4; } }

int main(void)
{
	static char mode[8], qs[200000], ts[200000], js[200000];
	int8_t mat[25];
	int i, j;
	for (i = 0; i < 4; i++) { for (j = 0; j < 4; j++) mat[i * 5 + j] = i == j ? 2 : -4; mat[i * 5 + 4] = -1; }
	for (j = 0; j < 5; j++) mat[20 + j] = -1;
	while (scanf("%7s", mode) == 1) {
		uint8_t *qv, *tv, *jv = 0;
		int ql, tl;
		if (scanf("%s %s", qs, ts) != 2) break;
		ql = strlen(qs); tl = strlen(ts);
		qv = malloc(ql + 1); tv = malloc(tl + 1);
		for (i = 0; i < ql; i++) qv[i] = enc(qs[i]);
		for (i = 0; i < tl; i++) tv[i] = enc(ts[i]);
		if (strcmp(mode, "LL") == 0) {
			int gapo, gape, qe, te, sc;
			void *qp;
			if (scanf("%d %d", &gapo, &gape) != 2) break;
			qp = ksw_ll_qinit(0, 2, ql, qv, 5, mat);
			sc = ksw_ll_i16(qp, tl, tv, gapo, gape, &qe, &te);
			printf("%d %d %d\n", sc, qe, te);
			free(qp);
		} else if (strcmp(mode, "Z2") == 0 || strcmp(mode, "D2") == 0) {
			int q, e, q2 = 0, e2 = 0, w, zdrop, eb, flag;
			ksw_extz_t ez;
			if (mode[0] == 'Z') { if (scanf("%d %d %d %d %d %d", &q, &e, &w, &zdrop, &eb, &flag) != 6) break; }
			else { if (scanf("%d %d %d %d %d %d %d %d", &q, &e, &q2, &e2, &w, &zdrop, &eb, &flag) != 8) break; }
			memset(&ez, 0, sizeof(ez));
			if (mode[0] == 'Z') ksw_extz2_sse(0, ql, qv, tl, tv, 5, mat, q, e, w, zdrop, eb, flag, &ez);
			else ksw_extd2_sse(0, ql, qv, tl, tv, 5, mat, q, e, q2, e2, w, zdrop, eb, flag, &ez);
			printf("%d %d %d %d %d %d %d %d %d %d %d", (int)ez.max, (int)ez.zdropped, ez.max_q, ez.max_t, ez.mqe, ez.mqe_t, ez.mte, ez.mte_q, ez.score, ez.reach_end, ez.n_cigar);
			for (i = 0; i < ez.n_cigar; i++) printf(" %d%c", ez.cigar[i] >> 4, "MIDN"[ez.cigar[i] & 0xf]);
			printf("\n");
			free(ez.cigar);
		} else {
			int q, e, q2, noncan, zdrop, jb, flag;
			ksw_extz_t ez;
			if (scanf("%d %d %d %d %d %d %d %s", &q, &e, &q2, &noncan, &zdrop, &jb, &flag, js) != 8) break;
			if (strcmp(js, "-") != 0) { jv = malloc(tl + 1); for (i = 0; i < tl; i++) jv[i] = js[i] - '0'; }
			memset(&ez, 0, sizeof(ez));
			ksw_exts2_sse(0, ql, qv, tl, tv, 5, mat, q, e, q2, noncan, zdrop, jb, flag, jv, &ez);
			printf("%d %d %d %d %d %d %d %d %d %d %d", (int)ez.max, (int)ez.zdropped, ez.max_q, ez.max_t, ez.mqe, ez.mqe_t, ez.mte, ez.mte_q, ez.score, ez.reach_end, ez.n_cigar);
			for (i = 0; i < ez.n_cigar; i++) printf(" %d%c", ez.cigar[i] >> 4, "MIDN"[ez.cigar[i] & 0xf]);
			printf("\n");
			free(ez.cigar); free(jv);
		}
		free(qv); free(tv);
	}
	return 0;
}
