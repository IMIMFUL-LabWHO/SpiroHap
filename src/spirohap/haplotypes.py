from __future__ import annotations

from collections import Counter, defaultdict

from .align import (
    AlignedRead,
    _insertion_consensus,
    align_reads,
    choose_medoid,
    edit_distance,
    mean_edit_rate,
    polish_consensus,
)
from .model import (
    Config,
    HaplotypeCall,
    ReadQC,
    Sample,
    SampleResult,
    TrimmedRead,
    VariantEvent,
)


def _candidate_events(
    alignments: list[AlignedRead], reference: str, cfg: Config
) -> tuple[list[VariantEvent], bool]:
    settings = cfg.haplotypes
    snps: list[VariantEvent] = []
    deletion_positions: list[int] = []
    deletion_support: dict[int, tuple[float, int]] = {}
    for pos, ref_base in enumerate(reference):
        counts = Counter(
            alignment.bases[pos]
            for alignment in alignments
            if alignment.bases[pos] != "?"
        )
        depth = sum(counts.values())
        if depth < settings.min_site_depth:
            continue
        base_depth = sum(counts.get(base, 0) for base in "ACGT")
        supported_alts = []
        for base in "ACGT":
            if base == ref_base or base_depth == 0:
                continue
            count = counts.get(base, 0)
            fraction = count / base_depth
            if (
                count >= settings.min_alt_reads
                and fraction >= settings.min_snp_fraction
            ):
                supported_alts.append((base, fraction, count))
        if supported_alts:
            supported_alts.sort(key=lambda item: (-item[1], item[0]))
            snps.append(
                VariantEvent(
                    "SNP",
                    pos,
                    pos + 1,
                    ref_base,
                    ",".join(item[0] for item in supported_alts),
                    supported_alts[0][1],
                    supported_alts[0][2],
                )
            )
        deletion_count = counts.get("-", 0)
        deletion_fraction = deletion_count / depth
        if (
            deletion_count >= settings.min_alt_reads
            and deletion_fraction >= settings.min_indel_fraction
        ):
            deletion_positions.append(pos)
            deletion_support[pos] = (deletion_fraction, deletion_count)

    deletions: list[VariantEvent] = []
    if deletion_positions:
        run_start = previous = deletion_positions[0]
        for pos in deletion_positions[1:] + [len(reference) + 1]:
            if pos == previous + 1:
                previous = pos
                continue
            run_end = previous + 1
            fractions = [
                deletion_support[index][0] for index in range(run_start, run_end)
            ]
            counts = [deletion_support[index][1] for index in range(run_start, run_end)]
            deletions.append(
                VariantEvent(
                    "DEL",
                    run_start,
                    run_end,
                    reference[run_start:run_end],
                    "-",
                    min(fractions),
                    min(counts),
                )
            )
            run_start = previous = pos

    insertion_strings: dict[int, list[str]] = defaultdict(list)
    for alignment in alignments:
        for anchor, inserted in alignment.insertions.items():
            if inserted:
                insertion_strings[anchor].append(inserted)
    insertions: list[VariantEvent] = []
    for anchor, strings in insertion_strings.items():
        if len(strings) < settings.min_alt_reads:
            continue
        representative = _insertion_consensus(strings)
        if not representative:
            continue
        max_distance = max(1, round(len(representative) * 0.15))
        similar = sum(
            edit_distance(value, representative) <= max_distance for value in strings
        )
        support = similar / len(alignments)
        if similar >= settings.min_alt_reads and support >= settings.min_indel_fraction:
            insertions.append(
                VariantEvent(
                    "INS", anchor, anchor, "-", representative, support, similar
                )
            )

    events = snps + deletions + insertions
    events.sort(
        key=lambda event: (event.start, {"DEL": 0, "SNP": 1, "INS": 2}[event.kind])
    )
    capped = len(events) > settings.max_variant_events
    if capped:
        events = sorted(
            events, key=lambda event: (-event.support, event.start, event.kind)
        )[: settings.max_variant_events]
        events.sort(key=lambda event: (event.start, event.kind))
    return events, capped


def _observe_event(alignment: AlignedRead, event: VariantEvent) -> str:
    if event.kind == "SNP":
        base = alignment.bases[event.start]
        if base == event.ref or base in event.alt.split(","):
            return base
        return "?"
    if event.kind == "DEL":
        values = alignment.bases[event.start : event.end]
        if values and all(value == "-" for value in values):
            return "DEL"
        if values and all(value in "ACGT" for value in values):
            return "REF"
        return "?"
    inserted = alignment.insertions.get(event.start, "")
    if not inserted:
        return "REF"
    max_distance = max(1, round(len(event.alt) * 0.15))
    return "INS" if edit_distance(inserted, event.alt) <= max_distance else "?"


def _profiles(
    alignments: list[AlignedRead], events: list[VariantEvent]
) -> list[tuple[str, ...]]:
    raw = [
        tuple(_observe_event(alignment, event) for event in events)
        for alignment in alignments
    ]
    if not events:
        return raw
    major = []
    for index in range(len(events)):
        counts = Counter(profile[index] for profile in raw if profile[index] != "?")
        major.append(counts.most_common(1)[0][0] if counts else "REF")
    return [
        tuple(
            value if value != "?" else major[index]
            for index, value in enumerate(profile)
        )
        for profile in raw
    ]


def _select_centers(
    profiles: list[tuple[str, ...]], cfg: Config
) -> list[tuple[str, ...]]:
    if not profiles or not profiles[0]:
        return [()]
    counts = Counter(profiles)
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    centers = []
    total = len(profiles)
    for profile, count in ordered:
        if not centers or (
            count >= cfg.haplotypes.min_haplotype_reads
            and count / total >= cfg.haplotypes.min_haplotype_fraction
        ):
            centers.append(profile)
        if len(centers) >= cfg.haplotypes.max_haplotypes:
            break
    return centers


def _profile_distance(left: tuple[str, ...], right: tuple[str, ...]) -> int:
    return sum(a != b for a, b in zip(left, right))


def _assign_profiles(
    profiles: list[tuple[str, ...]], centers: list[tuple[str, ...]], cfg: Config
) -> tuple[dict[int, list[int]], list[int]]:
    groups: dict[int, list[int]] = {index: [] for index in range(len(centers))}
    ambiguous: list[int] = []
    if len(centers) == 1:
        groups[0] = list(range(len(profiles)))
        return groups, ambiguous
    for read_index, profile in enumerate(profiles):
        distances = [_profile_distance(profile, center) for center in centers]
        best = min(distances)
        winners = [index for index, value in enumerate(distances) if value == best]
        if len(winners) != 1 or best > cfg.haplotypes.max_profile_mismatches:
            ambiguous.append(read_index)
        else:
            groups[winners[0]].append(read_index)
    return groups, ambiguous


def reconstruct_haplotypes(
    sample: Sample, reads: list[TrimmedRead], qc: ReadQC, cfg: Config, seed: int
) -> SampleResult:
    if len(reads) < cfg.haplotypes.min_site_depth:
        return SampleResult(
            sample=sample,
            qc=qc,
            status="FAIL_LOW_READS",
            warnings=[
                f"Only {len(reads)} reads remained; at least {cfg.haplotypes.min_site_depth} are required"
            ],
        )

    initial = choose_medoid(reads, cfg, seed)
    global_consensus = polish_consensus(reads, initial, cfg)
    alignments = align_reads(reads, global_consensus, cfg.filters.min_base_q)
    events, capped = _candidate_events(alignments, global_consensus, cfg)
    profiles = _profiles(alignments, events)
    centers = _select_centers(profiles, cfg)
    groups, ambiguous = _assign_profiles(profiles, centers, cfg)

    calls: list[HaplotypeCall] = []
    dominant_center = centers[0]
    for center_index, center in enumerate(centers):
        indices = groups.get(center_index, [])
        if len(indices) < cfg.haplotypes.min_haplotype_reads:
            continue
        fraction = len(indices) / len(reads)
        if center_index > 0 and fraction < cfg.haplotypes.min_haplotype_fraction:
            continue
        group_reads = [reads[index] for index in indices]
        hap_seed = choose_medoid(group_reads, cfg, seed + 1009 * (center_index + 1))
        sequence = polish_consensus(group_reads, hap_seed, cfg)
        forward = sum(read.orientation == "+" for read in group_reads)
        reverse = len(group_reads) - forward
        flags: list[str] = []
        if min(forward, reverse) == 0:
            flags.append("one_strand_only")
        elif (
            min(forward, reverse) / len(group_reads)
            < cfg.haplotypes.min_strand_fraction
        ):
            flags.append("strand_bias")
        if center_index > 0 and _profile_distance(center, dominant_center) == 1:
            flags.append("single_event_separation")
        residual = mean_edit_rate(group_reads, sequence)
        if residual > 0.05:
            flags.append("high_residual_error")
        calls.append(
            HaplotypeCall(
                sequence=sequence,
                assigned_reads=len(group_reads),
                fraction=fraction,
                forward_reads=forward,
                reverse_reads=reverse,
                profile=center,
                mean_edit_rate=residual,
                read_names=[read.name for read in group_reads],
                flags=flags,
            )
        )

    merged: dict[str, HaplotypeCall] = {}
    for call in calls:
        if call.sequence not in merged:
            merged[call.sequence] = call
            continue
        previous = merged[call.sequence]
        previous.assigned_reads += call.assigned_reads
        previous.fraction += call.fraction
        previous.forward_reads += call.forward_reads
        previous.reverse_reads += call.reverse_reads
        previous.read_names.extend(call.read_names)
        previous.flags = sorted(
            set(previous.flags + call.flags + ["profiles_collapsed_after_polishing"])
        )
        previous.mean_edit_rate = (
            previous.mean_edit_rate * (previous.assigned_reads - call.assigned_reads)
            + call.mean_edit_rate * call.assigned_reads
        ) / previous.assigned_reads

    warnings = []
    if capped:
        warnings.append(
            f"Candidate events exceeded max_variant_events={cfg.haplotypes.max_variant_events}; strongest retained"
        )
    if len(centers) >= cfg.haplotypes.max_haplotypes:
        warnings.append(
            "Reached max_haplotypes; inspect profile counts before interpreting rare alleles"
        )
    final_calls = sorted(
        merged.values(), key=lambda call: (-call.assigned_reads, call.sequence)
    )
    if not final_calls:
        return SampleResult(
            sample=sample,
            qc=qc,
            status="FAIL_NO_HAPLOTYPE",
            global_consensus=global_consensus,
            variant_events=events,
            ambiguous_reads=len(ambiguous),
            warnings=warnings + ["No read group passed the haplotype thresholds"],
        )
    return SampleResult(
        sample=sample,
        qc=qc,
        status="PASS",
        global_consensus=global_consensus,
        variant_events=events,
        haplotypes=final_calls,
        ambiguous_reads=len(ambiguous),
        ambiguous_read_names=[reads[index].name for index in ambiguous],
        warnings=warnings,
    )
