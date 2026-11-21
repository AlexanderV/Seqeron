/* A5-8c3e (B09 F241): dumps the type and tag tables of the feature vocabularies of the installed EMBOSS 6.6.0 (libajax, ajFeatVocabGetTypes / ajFeatVocabGetTags).
   Build: gcc -o vocab_dump emboss_feature_vocab_dump.c -I/usr/lib/emboss/include -I/usr/lib/emboss/include/epcre -L/usr/lib/emboss/lib -lajax -lajaxg -lnucleus -lacd -lensembl -lajaxdb -lepcre -Wl,-rpath,/usr/lib/emboss/lib
   Run: ./vocab_dump NAME   (NAME = embl emboss gff2 gff2protein gff3 gff3protein pir protein refseqp swiss) */
#include "ajax.h"
#include <stdio.h>
static void pr(const void *key, void **value, void *cl){ printf("%s\t%s\t%s\n", (const char*)cl, ajStrGetPtr((AjPStr)key), ajStrGetPtr((AjPStr)*value)); }
int main(int argc, char **argv){
  ajNamInit("emboss");
  const AjPTable t = ajFeatVocabGetTypes(argv[1]);
  const AjPTable g = ajFeatVocabGetTags(argv[1]);
  if(!t || !g){ fprintf(stderr,"null table\n"); return 1;}
  ajTableMap((AjPTable)t, pr, "types");
  ajTableMap((AjPTable)g, pr, "tags");
  return 0; }
