# TEMPORARY hand-off image R14-13-aligninfo-biopython179: Python 3.9 with Biopython 1.79, the last release that has Bio.Align.AlignInfo.SummaryInfo (dumb_consensus, gap_consensus, replacement_dictionary,
# pos_specific_score_matrix, information_content, ic_vector, PSSM, print_info_content; Biopython 1.88 keeps only the deprecated get_column). It is the oracle of the C# port SummaryInfo / Pssm; the
# sandbox that builds the port may not run downloaded code. Build and run with `--platform linux/amd64` (Apple Silicon: emulated). No account data, no secrets; the only download is the PyPI wheel of Biopython 1.79.
FROM --platform=linux/amd64 python:3.9-slim
RUN pip install --no-cache-dir "biopython==1.79" "numpy<2" && python -c "import Bio; print(Bio.__version__)"
WORKDIR /work
