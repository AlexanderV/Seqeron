#!/usr/bin/env bash
# Builds mafft-7.505-with-extensions from the official source archive into /opt/mafft-ext (PREFIX). Nothing here aborts the image build; the result is judged by in_container.sh (it lists the programs that exist).
set -u
UA="Seqeron-review/1.0"
URLS="https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz https://gitlab.com/sysimm/mafft/-/archive/v7.505/mafft-v7.505.tar.gz"
mkdir -p /opt/src && cd /opt/src
for u in $URLS; do
  echo "== try $u"
  if curl -fsSL -A "$UA" -o src.tgz "$u"; then echo "downloaded: $u"; sha256sum src.tgz; break; fi
done
[ -s src.tgz ] || { echo "NO SOURCE"; exit 0; }
mkdir -p x && tar xzf src.tgz -C x && ls x
TOP=$(find x -maxdepth 2 -type d -name core | head -1 | xargs dirname)
echo "top: $TOP"; ls "$TOP"
PREFIX=/opt/mafft-ext
( cd "$TOP/core" && make clean >/dev/null 2>&1; make PREFIX=$PREFIX 2>&1 | tail -20; make PREFIX=$PREFIX install 2>&1 | tail -20 )
if [ -d "$TOP/extensions" ]; then
  ( cd "$TOP/extensions" && make clean >/dev/null 2>&1; make PREFIX=$PREFIX 2>&1 | tail -40; make PREFIX=$PREFIX install 2>&1 | tail -40 )
else
  echo "NO extensions directory in this archive"
fi
echo "== programs"; ls -l $PREFIX/libexec/mafft 2>&1
$PREFIX/bin/mafft --version 2>&1 | head -1
exit 0
