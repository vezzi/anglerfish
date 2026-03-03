import gzip
import logging
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Literal, TYPE_CHECKING

import regex

if TYPE_CHECKING:
    from anglerfish.demux.samplesheet import SampleSheetEntry

try:
    import Levenshtein as lev
except ModuleNotFoundError:
    lev = None


def _lev_distance(a: str, b: str) -> int:
    if lev is not None:
        return lev.distance(a, b)

    # Fallback for environments missing python-Levenshtein.
    rows = len(a) + 1
    cols = len(b) + 1
    dp = [[0] * cols for _ in range(rows)]
    for i in range(rows):
        dp[i][0] = i
    for j in range(cols):
        dp[0][j] = j
    for i in range(1, rows):
        for j in range(1, cols):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + cost,
            )
    return dp[-1][-1]


@dataclass
class AnchorConfig:
    i7_left: str
    i7_right: str
    i5_left: str
    i5_right: str
    max_edits: int = 3
    i7_index_len: int = 8
    i5_index_len: int = 8
    orientation: Literal["both", "forward", "reverse"] = "both"


@dataclass
class AnchorMatch:
    start: int
    end: int
    edits: int


@dataclass
class ExtractedIndex:
    index: str | None
    detected: bool
    edits: int | None


@dataclass
class ExtractedRead:
    read_name: str
    read_len: int
    i7_index: str | None
    i5_index: str | None
    i7_detected: bool
    i5_detected: bool
    orientation: Literal["forward", "reverse"]
    edits: int


def _iter_fastq_reads(fastq_files: list[str]):
    for fq in fastq_files:
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
                yield header.strip().split()[0].lstrip("@"), seq


def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _find_anchor_matches(seq: str, anchor: str, max_edits: int) -> list[AnchorMatch]:
    if len(anchor) == 0:
        return []

    # Fast path: exact search first, even when fuzzy matching is enabled.
    # If exact matches are found they are always optimal (0 edits).
    exact_hits = []
    start = 0
    while True:
        pos = seq.find(anchor, start)
        if pos == -1:
            break
        exact_hits.append(AnchorMatch(pos, pos + len(anchor), 0))
        start = pos + 1
    if len(exact_hits) > 0:
        return exact_hits

    if max_edits == 0:
        return []

    pattern = regex.compile(
        f"({regex.escape(anchor)}){{e<={max_edits}}}", flags=regex.IGNORECASE
    )
    matches = []
    for m in pattern.finditer(seq):
        subs, ins, dels = m.fuzzy_counts
        matches.append(
            AnchorMatch(start=m.start(), end=m.end(), edits=(subs + ins + dels))
        )

    return sorted(matches, key=lambda m: (m.edits, m.start))


def _pick_best_pair(
    left_matches: list[AnchorMatch],
    right_matches: list[AnchorMatch],
    expected_index_len: int,
    read_len: int,
) -> tuple[AnchorMatch, AnchorMatch] | None:
    candidates = []
    for left, right in product(left_matches, right_matches):
        if right.start < left.end:
            continue
        index_len = right.start - left.end
        total_edits = left.edits + right.edits
        length_penalty = abs(index_len - expected_index_len)
        end_bias = min(left.start, read_len - right.end)
        candidates.append(((total_edits, length_penalty, end_bias, left.start), left, right))

    if len(candidates) == 0:
        return None
    return min(candidates, key=lambda c: c[0])[1:]


def _pick_best_left_match(matches: list[AnchorMatch], read_len: int) -> AnchorMatch | None:
    if len(matches) == 0:
        return None
    return min(matches, key=lambda m: (m.edits, min(m.start, read_len - m.end), m.start))


def _extract_index(
    seq: str, left_anchor: str, right_anchor: str, expected_index_len: int, max_edits: int
) -> ExtractedIndex:
    left_matches = _find_anchor_matches(seq, left_anchor, max_edits)
    right_matches = _find_anchor_matches(seq, right_anchor, max_edits)
    read_len = len(seq)

    best_pair = _pick_best_pair(left_matches, right_matches, expected_index_len, read_len)
    if best_pair is not None:
        left, right = best_pair
        return ExtractedIndex(
            index=seq[left.end : right.start],
            detected=True,
            edits=left.edits + right.edits,
        )

    # Loose fallback: match left anchor and take the next N bases
    left_best = _pick_best_left_match(left_matches, read_len)
    if left_best is None:
        return ExtractedIndex(index=None, detected=False, edits=None)

    i_start = left_best.end
    i_end = i_start + expected_index_len
    if i_end > read_len:
        return ExtractedIndex(index=None, detected=False, edits=None)

    return ExtractedIndex(
        index=seq[i_start:i_end],
        detected=True,
        edits=left_best.edits,
    )


def _extract_from_oriented_seq(
    read_name: str, seq: str, orientation: Literal["forward", "reverse"], config: AnchorConfig
) -> ExtractedRead:
    i7 = _extract_index(
        seq,
        config.i7_left,
        config.i7_right,
        config.i7_index_len,
        config.max_edits,
    )
    i5 = _extract_index(
        seq,
        config.i5_left,
        config.i5_right,
        config.i5_index_len,
        config.max_edits,
    )

    total_edits = 0
    if i7.edits is not None:
        total_edits += i7.edits
    if i5.edits is not None:
        total_edits += i5.edits

    return ExtractedRead(
        read_name=read_name,
        read_len=len(seq),
        i7_index=i7.index,
        i5_index=i5.index,
        i7_detected=i7.detected,
        i5_detected=i5.detected,
        orientation=orientation,
        edits=total_edits,
    )


def _pick_best_orientation(
    forward: ExtractedRead, reverse: ExtractedRead
) -> ExtractedRead:
    def score(extracted: ExtractedRead):
        both = int(extracted.i7_detected and extracted.i5_detected)
        either = int(extracted.i7_detected) + int(extracted.i5_detected)
        return (-both, -either, extracted.edits)

    return min([forward, reverse], key=score)


def extract_indices_from_fastqs(
    fastq_files: list[str],
    config: AnchorConfig,
    total_reads: int | None = None,
    progress_step_pct: int = 5,
    logger: logging.Logger | None = None,
    progress_label: str = "Anchor extraction",
) -> tuple[dict[str, ExtractedRead], dict[str, int]]:
    extracted_reads: dict[str, ExtractedRead] = {}
    counts = {
        "reads_processed": 0,
        "i7_detected": 0,
        "i5_detected": 0,
        "both_detected": 0,
        "forward_orientation": 0,
        "reverse_orientation": 0,
    }
    next_progress_pct = progress_step_pct

    for read_name, seq in _iter_fastq_reads(fastq_files):
        counts["reads_processed"] += 1
        if config.orientation == "forward":
            chosen = _extract_from_oriented_seq(read_name, seq, "forward", config)
        elif config.orientation == "reverse":
            chosen = _extract_from_oriented_seq(
                read_name, _reverse_complement(seq), "reverse", config
            )
        else:
            forward = _extract_from_oriented_seq(read_name, seq, "forward", config)
            reverse = _extract_from_oriented_seq(
                read_name, _reverse_complement(seq), "reverse", config
            )
            chosen = _pick_best_orientation(forward, reverse)

        extracted_reads[read_name] = chosen
        if chosen.i7_detected:
            counts["i7_detected"] += 1
        if chosen.i5_detected:
            counts["i5_detected"] += 1
        if chosen.i7_detected and chosen.i5_detected:
            counts["both_detected"] += 1
        if chosen.orientation == "forward":
            counts["forward_orientation"] += 1
        else:
            counts["reverse_orientation"] += 1

        if (
            logger is not None
            and total_reads is not None
            and total_reads > 0
            and next_progress_pct <= 100
        ):
            while (
                next_progress_pct <= 100
                and counts["reads_processed"] * 100 >= total_reads * next_progress_pct
            ):
                logger.info(
                    "%s progress: %s%% (%s/%s reads)",
                    progress_label,
                    next_progress_pct,
                    counts["reads_processed"],
                    total_reads,
                )
                next_progress_pct += progress_step_pct
        elif logger is not None and counts["reads_processed"] % 10000 == 0:
            logger.info(
                "%s progress: %s reads processed",
                progress_label,
                counts["reads_processed"],
            )

    if (
        logger is not None
        and total_reads is not None
        and total_reads > 0
        and next_progress_pct <= 100
    ):
        logger.info(
            "%s progress: 100%% (%s/%s reads)",
            progress_label,
            counts["reads_processed"],
            total_reads,
        )

    return extracted_reads, counts


def _flip_index(index: str | None, reverse: bool) -> str | None:
    if index is None:
        return None
    return _reverse_complement(index) if reverse else index


def match_extracted_reads(
    entries: list["SampleSheetEntry"],
    extracted_reads: dict[str, ExtractedRead],
    max_distance: int,
    demux_by_i5_only: bool = False,
    i7_reversed: bool = False,
    i5_reversed: bool = False,
) -> tuple[list[list], list[list]]:
    matched_bed: list[list] = []
    unmatched_bed: list[list] = []

    for read_name, extracted in extracted_reads.items():
        if demux_by_i5_only:
            if not extracted.i5_detected or extracted.i5_index is None:
                continue
        else:
            if (
                not extracted.i7_detected
                or extracted.i7_index is None
                or not extracted.i5_detected
                or extracted.i5_index is None
            ):
                continue

        distances = []
        for entry in entries:
            expected_i7 = _flip_index(entry.adaptor.i7.index_seq, i7_reversed)
            expected_i5 = _flip_index(entry.adaptor.i5.index_seq, i5_reversed)

            i7_dist = (
                _lev_distance((expected_i7 or "").lower(), extracted.i7_index.lower())
                if not demux_by_i5_only
                else 0
            )
            i5_dist = _lev_distance(
                (expected_i5 or "").lower(), extracted.i5_index.lower()
            )
            dist = i7_dist + i5_dist
            distances.append((dist, entry.sample_name))

        min_dist = min(distances, key=lambda d: d[0])[0]
        min_samples = [sample for dist, sample in distances if dist == min_dist]
        if len(min_samples) > 1:
            continue

        if min_dist <= max_distance:
            matched_bed.append(
                [read_name, 0, extracted.read_len, min_samples[0], "999", "."]
            )
        else:
            idx = extracted.i5_index or ""
            if not demux_by_i5_only:
                idx = f"{extracted.i7_index or ''}+{extracted.i5_index or ''}"
            unmatched_bed.append([read_name, 0, extracted.read_len, idx.lower(), "999", "."])

    return unmatched_bed, matched_bed
