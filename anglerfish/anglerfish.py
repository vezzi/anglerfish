#!/usr/bin/env python
import glob
import gzip
import logging
import multiprocessing
import os
import sys
import uuid
from collections import Counter
from importlib.metadata import version as get_version
from itertools import groupby
from pathlib import Path

import Levenshtein as lev
import numpy as np

from .demux.demux import (
    Alignment,
    categorize_matches,
    cluster_matches,
    map_reads_to_alns,
    run_minimap2,
    write_demuxedfastq,
)
from .demux.report import AlignmentStat, Report, SampleStat
from .demux.samplesheet import SampleSheet
from .demux.anchor import AnchorConfig, extract_indices_from_fastqs, match_extracted_reads

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("anglerfish")

MAX_PROCESSES = 64  # Ought to be enough for anybody

anglerfish_logo = r"""
     ___
   ( )  \ -..__
      _.|~”~~~”…_
    ^´           `>.
(+ (+ )             “<..<^(
  `´  ``´      ___       (
   \__..~      __(   _…_(
   \                /
    “--…_     _..~%´
         ```´´
"""


def _count_fastq_reads(fastq_files: list[str]) -> int:
    num_fq_lines = 0
    for fq_file in fastq_files:
        path = Path(fq_file)
        if path.suffix == ".gz":
            handle = gzip.open(path, "rb")
        else:
            handle = open(path, "rb")
        with handle as f:
            for _ in f:
                num_fq_lines += 1
    return int(num_fq_lines / 4)


def _resolve_anchor_preset(args):
    if getattr(args, "anchors", "") == "illumina_dual_len8":
        args.anchor_i7_left = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        args.anchor_i7_right = "ATCTCGTATGCCGTCTTCTGCTTG"
        args.anchor_i5_left = "AATGATACGGCGACCACCGAGATCTACAC"
        args.anchor_i5_right = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        args.anchor_index_len_i7 = 8
        args.anchor_index_len_i5 = 8
    elif getattr(args, "anchors", "") not in ["", None]:
        log.warning(
            f" Unknown anchor preset '{args.anchors}'. Using explicit anchor arguments."
        )


def run_demux(args):
    # Boilerplate
    try:
        multiprocessing.set_start_method("spawn")
    except RuntimeError:
        pass
    version = get_version("bio-anglerfish")
    run_uuid = str(uuid.uuid4())

    # Parse arguments
    if args.debug:
        log.setLevel(logging.DEBUG)

    # Backward-compatible defaults for direct function calls/tests bypassing the CLI parser.
    for key, value in {
        "demux_mode": "auto",
        "anchors": "",
        "anchor_i7_left": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "anchor_i7_right": "ATCTCGTATGCCGTCTTCTGCTTG",
        "anchor_i5_left": "AATGATACGGCGACCACCGAGATCTACAC",
        "anchor_i5_right": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "anchor_max_edits": 3,
        "anchor_index_len_i7": 8,
        "anchor_index_len_i5": 8,
        "anchor_orientation": "both",
        "demux_by_i5_only": False,
    }.items():
        if not hasattr(args, key):
            setattr(args, key, value)

    if args.threads > MAX_PROCESSES:
        log.warning(
            f" Setting threads to {MAX_PROCESSES} as the maximum number of processes is {MAX_PROCESSES}"
        )
        args.threads = MAX_PROCESSES

    if os.path.exists(args.out_fastq):
        raise FileExistsError(
            f"Output folder '{args.out_fastq}' already exists. Please remove it or specify another --run_name"
        )
    else:
        os.mkdir(args.out_fastq)

    # Print logo and info
    sys.stderr.write(anglerfish_logo)
    log.info(f" version {version}")
    log.info(f" arguments {vars(args)}")
    log.info(f" run uuid {run_uuid}")

    # Instantiate samplesheet and report
    ss = SampleSheet(args.samplesheet, args.ont_barcodes)
    report = Report(args.run_name, run_uuid, version)
    _resolve_anchor_preset(args)

    # Calculate the minimum edit distance between indices in the samplesheet
    min_distance = ss.minimum_bc_distance()

    # Determine which edit distance to use for index matching
    if args.max_distance is None:
        # Default: Set the maximum distance for barcode matching to 0, 1 or 2
        # depending on the smallest detected edit distance between indices in the samplesheet
        args.max_distance = min(min_distance - 1, 2)
        log.info(f"Using maximum edit distance of {args.max_distance}")
    if args.max_distance >= min_distance:
        log.error(
            f" The maximum allowed edit distance for barcode matching (={args.max_distance})"
            + f"is greater than the smallest detected edit distance between indices in samplesheet (={min_distance})"
            + ", which will result in ambiguous matches."
        )
        exit()
    log.debug(f"Samplesheet bc_dist == {min_distance}")

    # Get adaptor-barcode combinations from samplesheet
    adaptor_barcode_sets = ss.get_adaptor_barcode_sets()

    # Iterate across samplesheet rows conforming to the adaptor-barcode combinations
    out_fastqs = []
    for adaptor_name, ont_barcode in adaptor_barcode_sets:
        subset_rows = ss.subset_rows(adaptor_name=adaptor_name, ont_barcode=ont_barcode)

        # Grab the fastq files for the current entries
        assert all([subset_rows[0].fastq == row.fastq for row in subset_rows])
        fastq_path = subset_rows[0].fastq
        fastq_files = sorted(glob.glob(fastq_path))
        if len(fastq_files) == 0:
            log.warning(f"No FASTQ files found for pattern: {fastq_path}")
            continue

        # If there are multiple ONT barcodes, we need to add the ONT barcode to the adaptor name
        if ont_barcode:
            adaptor_bc_name = f"{adaptor_name}_{ont_barcode}"
        else:
            adaptor_bc_name = adaptor_name

        # Easy line count in input fastq files
        stats = AlignmentStat(adaptor_bc_name)
        num_fq = _count_fastq_reads(fastq_files)

        # Demux
        flipped_i7 = False
        flipped_i5 = False
        unmatched_bed: list[list] = []
        matched_bed: list[list] = []
        flips = {
            "original": {"i7_reversed": False, "i5_reversed": False},
            "i7": {"i7_reversed": True, "i5_reversed": False},
            "i5": {"i7_reversed": False, "i5_reversed": True},
            "i7+i5": {"i7_reversed": True, "i5_reversed": True},
        }
        use_anchor_mode = args.demux_mode == "anchor"
        fragments: dict[str, list[Alignment]] = {}
        singletons: dict[str, list[Alignment]] = {}
        concats: dict[str, list[Alignment]] = {}
        unknowns: dict[str, list[Alignment]] = {}

        if not use_anchor_mode:
            # Align using predefined adaptor models
            alignment_path: str = os.path.join(args.out_fastq, f"{adaptor_bc_name}.paf")
            adaptor_fasta_path: str = os.path.join(args.out_fastq, f"{adaptor_name}.fasta")
            with open(adaptor_fasta_path, "w") as f:
                f.write(ss.get_fastastring(adaptor_name))
            for fq_file in fastq_files:
                run_minimap2(fq_file, adaptor_fasta_path, alignment_path, args.threads)

            reads_to_alns: dict[str, list[Alignment]] = map_reads_to_alns(alignment_path)
            log.info(f" Searching for adaptor hits in {adaptor_bc_name}")
            fragments, singletons, concats, unknowns = categorize_matches(
                i5_name=f"{adaptor_name}_i5",
                i7_name=f"{adaptor_name}_i7",
                reads_to_alns=reads_to_alns,
            )
            total_hits = len(fragments) + len(singletons) + len(concats) + len(unknowns)
            stats.compute_pafstats(num_fq, fragments, singletons, concats, unknowns)

            if total_hits == 0:
                log.warning(
                    f" No adaptor hits found for {adaptor_bc_name}; nothing to demultiplex for this adaptor set."
                )
                if args.demux_mode == "auto":
                    log.info(
                        f" Falling back to anchor mode for {adaptor_bc_name} because adaptor hits are zero."
                    )
                    use_anchor_mode = True

        if use_anchor_mode:
            if args.demux_mode == "anchor":
                stats.compute_pafstats(num_fq, {}, {}, {}, {})
            log.info(
                f" Running anchor-based index extraction for {adaptor_bc_name} (max_edits={args.anchor_max_edits})."
            )
            anchor_config = AnchorConfig(
                i7_left=args.anchor_i7_left,
                i7_right=args.anchor_i7_right,
                i5_left=args.anchor_i5_left,
                i5_right=args.anchor_i5_right,
                max_edits=args.anchor_max_edits,
                i7_index_len=args.anchor_index_len_i7,
                i5_index_len=args.anchor_index_len_i5,
                orientation=args.anchor_orientation,
            )
            extracted_reads, anchor_counts = extract_indices_from_fastqs(
                fastq_files, anchor_config
            )

            if args.force_rc != "original":
                log.info(
                    f" Force reverse complementing {args.force_rc} index for adaptor {adaptor_name}. Lenient mode is disabled"
                )
                unmatched_bed, matched_bed = match_extracted_reads(
                    subset_rows,
                    extracted_reads,
                    args.max_distance,
                    demux_by_i5_only=args.demux_by_i5_only,
                    **flips[args.force_rc],
                )
                flipped_i7, flipped_i5 = flips[args.force_rc].values()
            elif args.lenient:
                flipped = {}
                for flip, rev in flips.items():
                    flipped[flip] = match_extracted_reads(
                        subset_rows,
                        extracted_reads,
                        args.max_distance,
                        demux_by_i5_only=args.demux_by_i5_only,
                        i7_reversed=rev["i7_reversed"],
                        i5_reversed=rev["i5_reversed"],
                    )
                best_flip = max(flipped, key=lambda k: len(flipped[k][1]))
                if (
                    sorted([len(i[1]) for i in flipped.values()])[-1]
                    == sorted([len(i[1]) for i in flipped.values()])[-2]
                ):
                    log.warning(
                        " Lenient mode: Could not find any barcode reverse complements with unambiguously more matches. Using original index orientation for all adaptors. Please study the results carefully!"
                    )
                    unmatched_bed, matched_bed = flipped["original"]
                elif (
                    best_flip != "None"
                    and len(flipped[best_flip][1])
                    > len(flipped["original"][1]) * args.lenient_factor
                ):
                    log.info(
                        f" Lenient mode: Reverse complementing {best_flip} index for adaptor {adaptor_name} found at least {args.lenient_factor} times more matches"
                    )
                    unmatched_bed, matched_bed = flipped[best_flip]
                    flipped_i7, flipped_i5 = flips[best_flip].values()
                else:
                    log.info(
                        f" Lenient mode: using original index orientation for {adaptor_name}"
                    )
                    unmatched_bed, matched_bed = flipped["original"]
            else:
                unmatched_bed, matched_bed = match_extracted_reads(
                    subset_rows,
                    extracted_reads,
                    args.max_distance,
                    demux_by_i5_only=args.demux_by_i5_only,
                )
            stats.add_anchorstats(num_fq, anchor_counts, len(matched_bed))
            log.info(
                " Anchor detection summary for %s: i7=%s i5=%s both=%s matched=%s",
                adaptor_bc_name,
                anchor_counts["i7_detected"],
                anchor_counts["i5_detected"],
                anchor_counts["both_detected"],
                len(matched_bed),
            )
        else:
            if args.force_rc != "original":
                log.info(
                    f" Force reverse complementing {args.force_rc} index for adaptor {adaptor_name}. Lenient mode is disabled"
                )
                unmatched_bed, matched_bed = cluster_matches(
                    subset_rows,
                    fragments,
                    args.max_distance,
                    **flips[args.force_rc],
                )
                flipped_i7, flipped_i5 = flips[args.force_rc].values()
            elif args.lenient:  # Try reverse complementing the i5 and/or i7 indices and choose the best match
                flipped = {}
                results = []
                pool = multiprocessing.Pool(
                    processes=4 if args.threads >= 4 else args.threads
                )
                results = []
                for flip, rev in flips.items():
                    spawn = pool.apply_async(
                        cluster_matches,
                        args=(
                            subset_rows,
                            fragments,
                            args.max_distance,
                            rev["i7_reversed"],
                            rev["i5_reversed"],
                        ),
                    )
                    results.append((spawn, flip))
                pool.close()
                pool.join()
                flipped = {result[1]: result[0].get() for result in results}

                best_flip = max(flipped, key=lambda k: len(flipped[k][1]))

                if (
                    sorted([len(i[1]) for i in flipped.values()])[-1]
                    == sorted([len(i[1]) for i in flipped.values()])[-2]
                ):
                    log.warning(
                        " Lenient mode: Could not find any barcode reverse complements with unambiguously more matches. Using original index orientation for all adaptors. Please study the results carefully!"
                    )
                    unmatched_bed, matched_bed = flipped["original"]
                elif (
                    best_flip != "None"
                    and len(flipped[best_flip][1])
                    > len(flipped["original"][1]) * args.lenient_factor
                ):
                    log.info(
                        f" Lenient mode: Reverse complementing {best_flip} index for adaptor {adaptor_name} found at least {args.lenient_factor} times more matches"
                    )
                    unmatched_bed, matched_bed = flipped[best_flip]
                    flipped_i7, flipped_i5 = flips[best_flip].values()
                else:
                    log.info(
                        f" Lenient mode: using original index orientation for {adaptor_name}"
                    )
                    unmatched_bed, matched_bed = flipped["original"]
            else:
                unmatched_bed, matched_bed = cluster_matches(
                    subset_rows, fragments, args.max_distance
                )

        report.add_alignment_stat(stats)

        if len(matched_bed) == 0:
            log.warning(
                f" No reads matched samplesheet barcodes for {adaptor_bc_name} (max_distance={args.max_distance})."
            )

        out_pool = []
        for k, v in groupby(
            sorted(matched_bed, key=lambda x: x[3]), key=lambda y: y[3]
        ):
            # To avoid collisions in fastq filenames, we add the ONT barcode to the sample name
            fq_prefix = k
            if ont_barcode:
                fq_prefix = ont_barcode + "-" + fq_prefix
            fq_name = os.path.join(args.out_fastq, fq_prefix + ".fastq.gz")
            out_fastqs.append(fq_name)
            sample_dict = {i[0]: [i] for i in v}

            # Find read lengths
            rlens = np.array([])
            for l, w in sample_dict.items():
                for i in w:
                    rlens = np.append(rlens, i[2] - i[1])
            rmean = np.round(np.mean(rlens), 2)
            rstd = np.round(np.std(rlens), 2)

            sample_stat = SampleStat(
                k,
                len(sample_dict.keys()),
                rmean,
                rstd,
                flipped_i7,
                flipped_i5,
                ont_barcode,
            )
            report.add_sample_stat(sample_stat)
            if not args.skip_demux:
                out_pool.append((sample_dict, fastq_path, fq_name))

        # Write demuxed fastq files
        pool = multiprocessing.Pool(processes=args.threads)
        results = []
        for out in out_pool:
            log.debug(f" Writing {out[2]}")
            spawn = pool.starmap_async(write_demuxedfastq, [out])
            results.append((spawn, out[2]))
        pool.close()
        pool.join()
        for result in results:
            log.debug(
                f" PID-{result[0].get()}: wrote {result[1]}, size {os.path.getsize(result[1])} bytes"
            )

        # Top unmatched indexes
        nomatch_count = Counter([x[3] for x in unmatched_bed])
        if args.max_unknowns == 0:
            args.max_unknowns = len([sample for sample in ss]) + 10

        # We search for the closest sample in the samplesheet to the list of unknowns
        top_unknowns = []
        for i in nomatch_count.most_common(args.max_unknowns):
            if args.demux_by_i5_only:
                sample_dists = [
                    (lev.distance(i[0], (x.adaptor.i5.index_seq or "").lower()), x.sample_name)
                    for x in ss
                ]
            else:
                sample_dists = [
                    (
                        lev.distance(
                            i[0],
                            f"{x.adaptor.i7.index_seq}+{x.adaptor.i5.index_seq}".lower(),
                        ),
                        x.sample_name,
                    )
                    for x in ss
                ]
            closest_sample = min(sample_dists, key=lambda x: x[0])
            # If the distance is more than half the index length, we remove it
            if closest_sample[0] >= (len(i[0]) / 2) + 1:
                top_unknowns.append([i[0], i[1], None])
            else:
                # We might have two samples with the same distance
                all_min = [
                    x[1]
                    for x in sample_dists
                    if x[0] == closest_sample[0] and x[1] != closest_sample[1]
                ]
                # This list might be too long, so we truncate it
                if len(all_min) > 4:
                    all_min = all_min[:4]
                    all_min.append("...")
                if all_min:
                    closest_sample = (closest_sample[0], ";".join(all_min))

                top_unknowns.append(
                    [i[0], i[1], f"{closest_sample[1]} ({closest_sample[0]})"]
                )
        report.add_unmatched_stat(top_unknowns, ont_barcode, adaptor_name)

    # Check if there were samples in the samplesheet without adaptor alignments and add them to report
    for entry in ss:
        if entry.sample_name not in [
            s.sample_name for s in [stat for stat in report.sample_stats]
        ]:
            sample_stat = SampleStat(
                entry.sample_name, 0, 0, 0, False, False, entry.ont_barcode
            )
            report.add_sample_stat(sample_stat)

    report.write_report(args.out_fastq)
    report.write_json(args.out_fastq)
    report.write_dataframe(args.out_fastq, ss)
