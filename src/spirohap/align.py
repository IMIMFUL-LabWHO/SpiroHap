from __future__ import annotations

import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass

import edlib

from .model import Config, TrimmedRead

CIGAR_RE = re.compile(r"(\d+)([=XID])")
BASES = "ACGT"


@dataclass(frozen=True)
class AlignedRead:
    read: TrimmedRead
    bases: tuple[str, ...]
    insertions: dict[int, str]
    edit_distance: int


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if "".join(f"{length}{op}" for length, op in operations) != cigar:
        raise ValueError(f"Unsupported CIGAR returned by edlib: {cigar!r}")
    return operations


def edit_distance(query: str, target: str) -> int:
    return int(edlib.align(query, target, mode="NW", task="distance")["editDistance"])


def align_read(read: TrimmedRead, reference: str, min_base_q: int) -> AlignedRead:
    result = edlib.align(read.sequence, reference, mode="NW", task="path")
    if result["editDistance"] < 0 or not result.get("cigar"):
        raise ValueError(f"Could not align read {read.name!r}")
    bases = ["?"] * len(reference)
    insertions: dict[int, str] = {}
    query_pos = 0
    target_pos = 0
    for length, op in parse_cigar(result["cigar"]):
        if op in {"=", "X"}:
            for offset in range(length):
                qpos = query_pos + offset
                tpos = target_pos + offset
                base = read.sequence[qpos]
                bases[tpos] = base if read.qualities[qpos] >= min_base_q else "?"
            query_pos += length
            target_pos += length
        elif op == "I":
            inserted = read.sequence[query_pos : query_pos + length]
            scores = read.qualities[query_pos : query_pos + length]
            if scores and sum(scores) / len(scores) >= min_base_q:
                insertions[target_pos] = insertions.get(target_pos, "") + inserted
            query_pos += length
        elif op == "D":
            for offset in range(length):
                bases[target_pos + offset] = "-"
            target_pos += length
    if query_pos != len(read.sequence) or target_pos != len(reference):
        raise ValueError("Internal alignment-coordinate error")
    return AlignedRead(read, tuple(bases), insertions, int(result["editDistance"]))


def align_reads(
    reads: list[TrimmedRead], reference: str, min_base_q: int
) -> list[AlignedRead]:
    return [align_read(read, reference, min_base_q) for read in reads]


def choose_medoid(reads: list[TrimmedRead], cfg: Config, seed: int) -> str:
    if not reads:
        return ""
    if len(reads) == 1:
        return reads[0].sequence
    rng = random.Random(seed)
    expected_length = (
        cfg.filters.min_trimmed_length + cfg.filters.max_trimmed_length
    ) / 2
    ranked = sorted(
        reads,
        key=lambda read: (
            -read.mean_q,
            abs(len(read.sequence) - expected_length),
            read.name,
        ),
    )
    candidate_pool = ranked[
        : max(cfg.runtime.medoid_candidates * 3, cfg.runtime.medoid_candidates)
    ]
    candidates = (
        candidate_pool
        if len(candidate_pool) <= cfg.runtime.medoid_candidates
        else rng.sample(candidate_pool, cfg.runtime.medoid_candidates)
    )
    comparisons = (
        reads
        if len(reads) <= cfg.runtime.medoid_comparisons
        else rng.sample(reads, cfg.runtime.medoid_comparisons)
    )
    scored = []
    for candidate in candidates:
        total = sum(
            edit_distance(candidate.sequence, other.sequence) for other in comparisons
        )
        scored.append((total, -candidate.mean_q, candidate.name, candidate.sequence))
    return min(scored)[3]


def _insertion_consensus(strings: list[str]) -> str:
    if not strings:
        return ""
    length_counts = Counter(map(len, strings))
    modal_length = min(
        (
            length
            for length, count in length_counts.items()
            if count == max(length_counts.values())
        ),
        key=lambda value: (
            abs(value - sorted(map(len, strings))[len(strings) // 2]),
            value,
        ),
    )
    same_length = [seq for seq in strings if len(seq) == modal_length]
    if not same_length:
        return Counter(strings).most_common(1)[0][0]
    consensus = []
    for pos in range(modal_length):
        counts = Counter(seq[pos] for seq in same_length if seq[pos] in BASES)
        consensus.append(counts.most_common(1)[0][0] if counts else "N")
    return "".join(consensus)


def consensus_from_alignments(
    alignments: list[AlignedRead], reference: str, consensus_fraction: float
) -> str:
    if not alignments:
        return ""
    base_counts = [Counter() for _ in reference]
    insertion_strings: dict[int, list[str]] = defaultdict(list)
    for alignment in alignments:
        for pos, base in enumerate(alignment.bases):
            if base != "?":
                base_counts[pos][base] += 1
        for anchor, inserted in alignment.insertions.items():
            if inserted:
                insertion_strings[anchor].append(inserted)

    called_insertions: dict[int, str] = {}
    total = len(alignments)
    for anchor, strings in insertion_strings.items():
        if len(strings) / total >= consensus_fraction:
            called_insertions[anchor] = _insertion_consensus(strings)

    pieces: list[str] = []
    if 0 in called_insertions:
        pieces.append(called_insertions[0])
    for pos, reference_base in enumerate(reference):
        counts = base_counts[pos]
        observed = sum(counts.values())
        deletion_fraction = counts.get("-", 0) / observed if observed else 0.0
        if deletion_fraction < consensus_fraction:
            nucleotide_counts = {base: counts.get(base, 0) for base in BASES}
            best_count = max(nucleotide_counts.values(), default=0)
            winners = [
                base for base, count in nucleotide_counts.items() if count == best_count
            ]
            if best_count == 0 or reference_base in winners:
                called = reference_base
            else:
                called = min(winners)
            pieces.append(called)
        if pos + 1 in called_insertions:
            pieces.append(called_insertions[pos + 1])
    return "".join(pieces)


def polish_consensus(reads: list[TrimmedRead], seed_sequence: str, cfg: Config) -> str:
    consensus = seed_sequence
    for _ in range(cfg.haplotypes.polish_rounds):
        alignments = align_reads(reads, consensus, cfg.filters.min_base_q)
        updated = consensus_from_alignments(
            alignments, consensus, cfg.haplotypes.consensus_fraction
        )
        if not updated or updated == consensus:
            break
        consensus = updated
    return consensus


def mean_edit_rate(reads: list[TrimmedRead], consensus: str) -> float:
    if not reads or not consensus:
        return 0.0
    rates = [
        edit_distance(read.sequence, consensus)
        / max(len(read.sequence), len(consensus))
        for read in reads
    ]
    return sum(rates) / len(rates)


def sequence_variants(reference: str, query: str) -> list[dict[str, str | int]]:
    """Return simple 1-based SNP/indel records for query relative to reference."""
    result = edlib.align(query, reference, mode="NW", task="path")
    variants: list[dict[str, str | int]] = []
    query_pos = 0
    reference_pos = 0
    for length, op in parse_cigar(result["cigar"]):
        if op == "=":
            query_pos += length
            reference_pos += length
        elif op == "X":
            for offset in range(length):
                variants.append(
                    {
                        "type": "SNP",
                        "position": reference_pos + offset + 1,
                        "ref": reference[reference_pos + offset],
                        "alt": query[query_pos + offset],
                    }
                )
            query_pos += length
            reference_pos += length
        elif op == "I":
            variants.append(
                {
                    "type": "INS",
                    "position": reference_pos,
                    "ref": "-",
                    "alt": query[query_pos : query_pos + length],
                }
            )
            query_pos += length
        elif op == "D":
            variants.append(
                {
                    "type": "DEL",
                    "position": reference_pos + 1,
                    "ref": reference[reference_pos : reference_pos + length],
                    "alt": "-",
                }
            )
            reference_pos += length
    return variants
