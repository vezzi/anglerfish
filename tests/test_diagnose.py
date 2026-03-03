from anglerfish.diagnose import run_diagnose


def test_run_diagnose_detects_motifs(tmp_path):
    seq = (
        "TTTT"
        + "AATGATACGGCGACCACCGAGATCTACAC"
        + "TAGATCGC"
        + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        + "GGGG"
        + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        + "GGACTCCT"
        + "ATCTCGTATGCCGTCTTCTGCTTG"
    )
    fq = tmp_path / "diag.fastq"
    with open(fq, "w") as f:
        f.write(f"@read1\n{seq}\n+\n{'I' * len(seq)}\n")

    outdir = tmp_path / "diag_out"
    result = run_diagnose(str(fq), max_reads=100, outdir=str(outdir))

    assert result["reads_scanned"] == 1
    assert result["motif_counts"]["i7_left_core"] == 1
    assert result["motif_counts"]["i5_left_core"] == 1
    assert result["recommendations"]["anchor_i7_left"] == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
