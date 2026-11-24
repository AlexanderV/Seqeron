/* the libsais SA of a text file, written as little-endian int32 (= SA[1..n] of the container) to a file or to stdout: usage sa_dump <text> <out> ; prints the build time on stderr (ADR-005 section 8, the check at 5x10^8) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "libsais.h"
int main(int argc, char **argv) {
  if (argc < 3) return 2;
  FILE *f = fopen(argv[1], "rb"); if (!f) return 3;
  fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
  uint8_t *T = malloc(n); int32_t *SA = malloc((size_t)n * 4);
  if (!T || !SA || fread(T, 1, n, f) != (size_t)n) return 4;
  fclose(f);
  struct timespec a, b; clock_gettime(CLOCK_MONOTONIC, &a);
  int32_t rc = libsais(T, SA, (int32_t)n, 0, NULL);
  clock_gettime(CLOCK_MONOTONIC, &b);
  if (rc != 0) return 5;
  fprintf(stderr, "libsais n=%ld time=%.3f s\n", n, (b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) * 1e-9);
  FILE *o = fopen(argv[2], "wb"); if (!o) return 6;
  fwrite(SA, 4, (size_t)n, o); fclose(o);
  return 0;
}
