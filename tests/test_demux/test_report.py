from anglerfish.demux.report import AlignmentStat


def test_compute_pafstats_handles_zero_total():
    stat = AlignmentStat("test_adaptor")

    stat.compute_pafstats(
        num_fq=10,
        fragments={},
        singletons={},
        concats={},
        unknowns={},
    )

    assert stat.paf_stats["input_reads"] == [10, 1.0]
    assert stat.paf_stats["reads aligning to adaptor sequences"] == [0, 0.0]
    assert stat.paf_stats["aligned reads matching both I7 and I5 adaptor"] == [0, 0.0]
    assert stat.paf_stats["aligned reads matching only I7 or I5 adaptor"] == [0, 0.0]
    assert stat.paf_stats["aligned reads matching multiple I7/I5 adaptor pairs"] == [
        0,
        0.0,
    ]
    assert stat.paf_stats["aligned reads with uncategorized alignments"] == [0, 0.0]
