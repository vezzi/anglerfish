from anglerfish.demux.adaptor import Adaptor
from anglerfish.demux.anchor import (
    AnchorConfig,
    ExtractedRead,
    extract_indices_from_fastqs,
    match_extracted_reads,
)
from types import SimpleNamespace


def _mk_entry(sample_name: str, i7: str, i5: str):
    adaptor = Adaptor(
        name="truseq_dual",
        i7_sequence_token="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC<N>ATCTCGTATGCCGTCTTCTGCTTG",
        i5_sequence_token="AATGATACGGCGACCACCGAGATCTACAC<N>ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        i7_index=i7,
        i5_index=i5,
    )
    return SimpleNamespace(sample_name=sample_name, adaptor=adaptor)


def _write_fastq(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_anchor_extracts_i7_and_i5_perfect(tmp_path):
    i7 = "GGACTCCT"
    i5 = "TAGATCGC"
    seq = (
        "TTTT"
        + "AATGATACGGCGACCACCGAGATCTACAC"
        + i5
        + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        + "NNNNNNNN"
        + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        + i7
        + "ATCTCGTATGCCGTCTTCTGCTTG"
        + "AAAA"
    )

    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, [("read1", seq)])

    config = AnchorConfig(
        i7_left="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        i7_right="ATCTCGTATGCCGTCTTCTGCTTG",
        i5_left="AATGATACGGCGACCACCGAGATCTACAC",
        i5_right="ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        max_edits=0,
        i7_index_len=8,
        i5_index_len=8,
    )
    extracted, counts = extract_indices_from_fastqs([str(fastq)], config)

    assert counts["i7_detected"] == 1
    assert counts["i5_detected"] == 1
    assert counts["both_detected"] == 1
    assert extracted["read1"].i7_index == i7
    assert extracted["read1"].i5_index == i5


def test_anchor_extracts_with_anchor_mismatches(tmp_path):
    i7 = "GGACTCCT"
    i5 = "TAGATCGC"
    # one mismatch in i7 left anchor and one mismatch in i5 left anchor
    seq = (
        "AATGATACGGCGACCACCGAGATCTACAT"
        + i5
        + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        + "CCCC"
        + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAT"
        + i7
        + "ATCTCGTATGCCGTCTTCTGCTTG"
    )
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, [("read1", seq)])

    config = AnchorConfig(
        i7_left="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        i7_right="ATCTCGTATGCCGTCTTCTGCTTG",
        i5_left="AATGATACGGCGACCACCGAGATCTACAC",
        i5_right="ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        max_edits=2,
        i7_index_len=8,
        i5_index_len=8,
    )
    extracted, counts = extract_indices_from_fastqs([str(fastq)], config)

    assert counts["both_detected"] == 1
    assert extracted["read1"].i7_index == i7
    assert extracted["read1"].i5_index == i5


def test_anchor_prefers_best_pair_in_concatemer(tmp_path):
    # first i7 pair has two anchor mismatches and wrong index, second i7 pair is perfect
    bad_i7 = "AAAAAAAA"
    good_i7 = "GGACTCCT"
    i5 = "TAGATCGC"
    seq = (
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAA"
        + bad_i7
        + "ATCTCGTATGCCGTCTTCTGCTTA"
        + "TTTTTT"
        + "AATGATACGGCGACCACCGAGATCTACAC"
        + i5
        + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        + "GGGGGG"
        + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        + good_i7
        + "ATCTCGTATGCCGTCTTCTGCTTG"
    )

    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, [("read1", seq)])

    config = AnchorConfig(
        i7_left="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        i7_right="ATCTCGTATGCCGTCTTCTGCTTG",
        i5_left="AATGATACGGCGACCACCGAGATCTACAC",
        i5_right="ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        max_edits=2,
        i7_index_len=8,
        i5_index_len=8,
    )
    extracted, _ = extract_indices_from_fastqs([str(fastq)], config)

    assert extracted["read1"].i7_index == good_i7
    assert extracted["read1"].i5_index == i5


def test_match_extracted_reads_i5_only():
    entries = [
        _mk_entry("sample1", "GGACTCCT", "TAGATCGC"),
        _mk_entry("sample2", "GGACTCCT", "CCTATCCT"),
    ]
    extracted_reads = {
        "read1": ExtractedRead(
            read_name="read1",
            read_len=120,
            i7_index="GGACTCCT",
            i5_index="TAGATCGC",
            i7_detected=True,
            i5_detected=True,
            orientation="forward",
            edits=0,
        ),
    }
    unmatched, matched = match_extracted_reads(
        entries,
        extracted_reads,
        max_distance=0,
        demux_by_i5_only=True,
    )
    assert len(unmatched) == 0
    assert len(matched) == 1
    assert matched[0][3] == "sample1"


def test_match_extracted_reads_tolerates_small_index_errors():
    entries = [_mk_entry("sample1", "GGACTCCT", "TAGATCGC")]
    extracted_reads = {
        "read1": ExtractedRead(
            read_name="read1",
            read_len=120,
            i7_index="GGACTCCT",
            i5_index="TAGATCCC",  # one mismatch from TAGATCGC
            i7_detected=True,
            i5_detected=True,
            orientation="forward",
            edits=1,
        ),
    }
    unmatched, matched = match_extracted_reads(
        entries,
        extracted_reads,
        max_distance=1,
        demux_by_i5_only=False,
    )
    assert len(unmatched) == 0
    assert len(matched) == 1
    assert matched[0][3] == "sample1"
