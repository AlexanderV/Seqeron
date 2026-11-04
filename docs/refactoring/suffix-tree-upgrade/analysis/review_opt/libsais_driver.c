#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <openssl/sha.h>
#include "libsais.h"
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
static void hex(const int32_t*a,size_t n,char*out){ unsigned char h[32]; SHA256((const unsigned char*)a,n*4,h); for(int i=0;i<32;i++) sprintf(out+2*i,"%02x",h[i]); }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"rb"); fseek(f,0,SEEK_END); long n=ftell(f); fseek(f,0,SEEK_SET);
  uint8_t*T=malloc(n); if(fread(T,1,n,f)!=(size_t)n) return 1; fclose(f);
  int32_t*SA=malloc(4*(size_t)n),*PLCP=malloc(4*(size_t)n),*LCP=malloc(4*(size_t)n);
  double t0=now(); libsais(T,SA,(int32_t)n,0,NULL); double t1=now();
  libsais_plcp(T,SA,PLCP,(int32_t)n); libsais_lcp(PLCP,SA,LCP,(int32_t)n); double t2=now();
  printf("sa=%.4f plcp_lcp=%.4f total=%.4f",t1-t0,t2-t1,t2-t0);
  if(argc>2){ char a[65],b[65]; hex(SA,n,a); hex(LCP,n,b); printf(" sa_sha=%s lcp_sha=%s",a,b); }
  printf("\n"); return 0; }
