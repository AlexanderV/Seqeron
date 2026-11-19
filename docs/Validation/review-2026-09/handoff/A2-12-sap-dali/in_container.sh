#!/usr/bin/env bash
# Runs inside the image (never aborts on a failing step: every step logs and the script goes on).
cd /work
mkdir -p out
{
  echo "date: $(date -u +%FT%TZ)"; echo "uname -m: $(uname -m)"; echo "t_coffee: $(t_coffee -version 2>&1 | head -1)"
  echo "SAP build: $( [ -f /opt/sap/BUILD_FAILED ] && echo FAILED || cat /opt/sap/sap.sha256 )"
  echo "Dali commit: $(cat /opt/dali/COMMIT 2>/dev/null) build: $( [ -f /opt/dali/BUILD_FAILED ] && echo FAILED || echo ok )"
  echo "TMalign: $(dpkg -s tm-align | grep ^Version)"; echo "mustang: $(dpkg -s mustang | grep ^Version)"
} > out/ENV.txt
cp /opt/sap/build.log out/sap_build.log 2>/dev/null; cp /opt/dali/build.log out/dali_build.log 2>/dev/null
# (1) SAP: the same synthetic CA-only structures as the committed TMalign / mustang oracle, run through t_coffee with the REAL sap on the PATH of the plug-in directory.
#     SAP_REAL=1 makes the generator link /opt/sap/bin/sap, keep the *.sap_results files of t_coffee's cache (SAP's raw output) and record rows + raw output per set.
export SAP_REAL=1 PATH=/opt/sap/bin:$PATH
python3 -I /work/generate_tcoffee_external_oracle.py /work/out sap_real > out/gen_sap.log 2>&1 || echo "sap generation failed (see gen_sap.log)"
# (2) DaliLite v5: its own command line (dali.pl / import.pl) differs from the `DaliLite -pairwise a.pdb b.pdb` call of t_coffee 13.41 (DaliLite v3). The package only records what v5 does with
#     the synthetic CA-only structures (import + pairwise alignment of the first 24 pairs); if it refuses them, the refusal is the result.
python3 -I /work/dali_probe.py /work/out > out/dali_probe.log 2>&1 || echo "dali probe failed (see dali_probe.log)"
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out
