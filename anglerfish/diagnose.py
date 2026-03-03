import gzip
import glob
import json
import logging
import os
from pathlib import Path

log = logging.getLogger("anglerfish")


KNOWN_MOTIFS = {
    "i7_left_core": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "i7_right_core": "ATCTCGTATGCCGTCTTCTGCTTG",
    "i5_left_core": "AATGATACGGCGACCACCGAGATCTACAC",
    "i5_right_core": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "illumina_i7_short": "AGATCGGAAGAGC",
    "illumina_p5_short": "AATGATACGGCG",
}


def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _iter_fastq_reads(fastq_path: str):
    for fq in sorted(glob.glob(fastq_path)):
        path = Path(fq)
        if path.suffix == ".gz":
            handle = gzip.open(path, "rt")
        else:
            handle = open(path)
        with handle as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # plus line
                f.readline()  # quality line
                yield seq


def run_diagnose(
    fastq: str,
    max_reads: int = 50000,
    outdir: str | None = None,
):
    motif_counts = {motif: 0 for motif in KNOWN_MOTIFS}
    forward_pref = 0
    reverse_pref = 0
    reads_scanned = 0

    for seq in _iter_fastq_reads(fastq):
        reads_scanned += 1
        rc_seq = _reverse_complement(seq)

        score_fwd = 0
        score_rev = 0
        for motif_name, motif in KNOWN_MOTIFS.items():
            if motif in seq:
                motif_counts[motif_name] += 1
                score_fwd += 1
            if motif in rc_seq:
                score_rev += 1
        if score_fwd >= score_rev:
            forward_pref += 1
        else:
            reverse_pref += 1

        if reads_scanned >= max_reads:
            break

    recommended_orientation = (
        "forward" if forward_pref >= reverse_pref else "reverse"
    )
    recommendations = {
        "anchor_i7_left": KNOWN_MOTIFS["i7_left_core"],
        "anchor_i7_right": KNOWN_MOTIFS["i7_right_core"],
        "anchor_i5_left": KNOWN_MOTIFS["i5_left_core"],
        "anchor_i5_right": KNOWN_MOTIFS["i5_right_core"],
        "orientation": recommended_orientation,
    }

    result = {
        "reads_scanned": reads_scanned,
        "motif_counts": motif_counts,
        "orientation_support": {"forward": forward_pref, "reverse": reverse_pref},
        "recommendations": recommendations,
    }

    log.info(f"Scanned {reads_scanned} reads from {fastq}")
    log.info(f"Motif counts: {motif_counts}")
    log.info(f"Recommended orientation: {recommended_orientation}")
    log.info(
        "Recommended anchors: i7-left=%s i7-right=%s i5-left=%s i5-right=%s",
        recommendations["anchor_i7_left"],
        recommendations["anchor_i7_right"],
        recommendations["anchor_i5_left"],
        recommendations["anchor_i5_right"],
    )

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, "anglerfish_diagnose.json")
        with open(outfile, "w") as f:
            json.dump(result, f, indent=2)
        log.info(f"Wrote diagnose summary to {outfile}")

    return result
