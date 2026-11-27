from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from . import __version__
from .align import edit_distance, sequence_variants
from .config import config_as_dict
from .haplotypes import reconstruct_haplotypes
from .io import load_sample_reads, read_fasta
from .model import Config, Sample, SampleResult
from .report import write_report


def _sample_seed(base_seed: int, sample_id: str) -> int:
    digest = hashlib.sha256(sample_id.encode("utf-8")).digest()
    return base_seed ^ int.from_bytes(digest[:8], "big")


def _process_one(payload: tuple[Sample, Config]) -> SampleResult:
    sample, cfg = payload
    seed = _sample_seed(cfg.runtime.random_seed, sample.sample_id)
    reads, qc = load_sample_reads(sample, cfg, seed)
    result = reconstruct_haplotypes(sample, reads, qc, cfg, seed)
    if sample.control_type == "negative":
        if result.status == "FAIL_LOW_READS":
            result.status = "PASS_NEGATIVE_CONTROL"
            result.warnings = []
        elif result.haplotypes:
            result.status = "WARN_NEGATIVE_CONTROL_SIGNAL"
            result.warnings.append(
                "One or more haplotypes reconstructed in a negative control"
            )
    return result


def _write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def _safe_name(value: str) -> str:
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:8]
    return f"{safe or 'sample'}_{digest}"


def _fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    row_human = a + b
    row_tick = c + d
    positives = a + c
    total = row_human + row_tick
    if total == 0:
        return 1.0
    denominator = math.comb(total, positives)

    def probability(x: int) -> float:
        return (
            math.comb(row_human, x) * math.comb(row_tick, positives - x) / denominator
        )

    observed = probability(a)
    lower = max(0, positives - row_tick)
    upper = min(row_human, positives)
    return min(
        1.0,
        sum(
            probability(x)
            for x in range(lower, upper + 1)
            if probability(x) <= observed + 1e-15
        ),
    )


def _bh_adjust(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    order = sorted(range(len(p_values)), key=lambda index: p_values[index])
    adjusted = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_index in range(total - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def _minimum_spanning_tree(
    ids: list[str], sequences: dict[str, str]
) -> list[dict[str, Any]]:
    if len(ids) < 2:
        return []
    parent = {item: item for item in ids}

    def find(item: str) -> str:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    weighted = []
    for left_index, left in enumerate(ids):
        for right in ids[left_index + 1 :]:
            weighted.append(
                (edit_distance(sequences[left], sequences[right]), left, right)
            )
    weighted.sort()
    chosen = []
    for distance, left, right in weighted:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            continue
        parent[root_right] = root_left
        chosen.append({"source": left, "target": right, "distance": distance})
        if len(chosen) == len(ids) - 1:
            break
    return chosen


def _normalize_references(path: str | Path | None) -> list[tuple[str, str]]:
    if path is None:
        return []
    records = read_fasta(path)
    if not records:
        raise ValueError("Reference FASTA contains no sequences")
    for name, sequence in records:
        if set(sequence) - set("ACGTN"):
            raise ValueError(
                f"Reference {name!r} contains non-IUPAC bases outside A/C/G/T/N"
            )
    return records


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def run_pipeline(
    samples: list[Sample],
    cfg: Config,
    output_dir: str | Path,
    jobs: int = 1,
    reference_fasta: str | Path | None = None,
    samples_path: str | Path | None = None,
    config_path: str | Path | None = None,
) -> dict[str, Any]:
    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    per_sample_dir = output / "per_sample"
    per_sample_dir.mkdir()

    payloads = [(sample, cfg) for sample in samples]
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            results = list(executor.map(_process_one, payloads))
    else:
        results = [_process_one(payload) for payload in payloads]

    sequence_occurrences: dict[str, list[tuple[SampleResult, Any]]] = {}
    for result in results:
        for call in result.haplotypes:
            sequence_occurrences.setdefault(call.sequence, []).append((result, call))
    ordered_sequences = sorted(
        sequence_occurrences,
        key=lambda sequence: (
            -len(
                {
                    result.sample.sample_id
                    for result, _ in sequence_occurrences[sequence]
                }
            ),
            -sum(call.assigned_reads for _, call in sequence_occurrences[sequence]),
            sequence,
        ),
    )
    sequence_to_id = {
        sequence: f"SIXO23S-H{index:03d}"
        for index, sequence in enumerate(ordered_sequences, 1)
    }
    for sequence, occurrences in sequence_occurrences.items():
        for _, call in occurrences:
            call.haplotype_id = sequence_to_id[sequence]

    negative_control_ids = {
        call.haplotype_id
        for result in results
        if result.sample.control_type == "negative"
        for call in result.haplotypes
    }
    replicate_support: dict[tuple[str, str], set[str]] = {}
    for result in results:
        if not result.sample.biological_id:
            continue
        for call in result.haplotypes:
            key = (result.sample.biological_id, call.haplotype_id)
            replicate_support.setdefault(key, set()).add(
                result.sample.replicate_id or result.sample.sample_id
            )

    references = _normalize_references(reference_fasta)
    hap_sequences = {
        sequence_to_id[sequence]: sequence for sequence in ordered_sequences
    }
    ids = [sequence_to_id[sequence] for sequence in ordered_sequences]

    catalog_rows: list[dict[str, Any]] = []
    for sequence in ordered_sequences:
        haplotype_id = sequence_to_id[sequence]
        occurrences = sequence_occurrences[sequence]
        sample_ids = {result.sample.sample_id for result, _ in occurrences}
        human_ids = {
            result.sample.sample_id
            for result, _ in occurrences
            if not result.sample.is_control and result.sample.source == "human"
        }
        tick_ids = {
            result.sample.sample_id
            for result, _ in occurrences
            if not result.sample.is_control and result.sample.source == "tick"
        }
        negative_ids = {
            result.sample.sample_id
            for result, _ in occurrences
            if result.sample.control_type == "negative"
        }
        neighbors = [other for other in ids if other != haplotype_id]
        if neighbors:
            nearest = min(
                neighbors,
                key=lambda other: (
                    edit_distance(sequence, hap_sequences[other]),
                    other,
                ),
            )
            nearest_distance: int | str = edit_distance(
                sequence, hap_sequences[nearest]
            )
        else:
            nearest, nearest_distance = "", ""
        closest_reference, reference_distance, reference_identity = "", "", ""
        if references:
            best_name, best_sequence = min(
                references,
                key=lambda record: (edit_distance(sequence, record[1]), record[0]),
            )
            distance = edit_distance(sequence, best_sequence)
            closest_reference = best_name
            reference_distance = distance
            reference_identity = round(
                1 - distance / max(len(sequence), len(best_sequence)), 6
            )
        catalog_rows.append(
            {
                "haplotype_id": haplotype_id,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "length": len(sequence),
                "sample_count": len(sample_ids),
                "human_samples": len(human_ids),
                "tick_samples": len(tick_ids),
                "negative_control_samples": len(negative_ids),
                "total_reads": sum(call.assigned_reads for _, call in occurrences),
                "nearest_haplotype": nearest,
                "nearest_distance": nearest_distance,
                "closest_reference": closest_reference,
                "reference_distance": reference_distance,
                "reference_identity": reference_identity,
            }
        )

    sample_haplotype_rows: list[dict[str, Any]] = []
    for result in results:
        for call in result.haplotypes:
            flags = list(call.flags)
            if call.haplotype_id in negative_control_ids:
                flags.append("detected_in_negative_control")
            if result.sample.biological_id:
                count = len(
                    replicate_support[(result.sample.biological_id, call.haplotype_id)]
                )
                replicate_status = (
                    f"replicated({count})" if count >= 2 else "single_replicate"
                )
            else:
                replicate_status = "not_available"
            event_profile = ";".join(
                f"{event.label}={allele}"
                for event, allele in zip(result.variant_events, call.profile)
            )
            sample_haplotype_rows.append(
                {
                    "sample_id": result.sample.sample_id,
                    "source": result.sample.source,
                    "biological_id": result.sample.biological_id,
                    "replicate_id": result.sample.replicate_id,
                    "is_control": str(result.sample.is_control).lower(),
                    "control_type": result.sample.control_type,
                    "haplotype_id": call.haplotype_id,
                    "assigned_reads": call.assigned_reads,
                    "fraction": round(call.fraction, 6),
                    "forward_reads": call.forward_reads,
                    "reverse_reads": call.reverse_reads,
                    "mean_read_edit_rate": round(call.mean_edit_rate, 6),
                    "replicate_status": replicate_status,
                    "event_profile": event_profile,
                    "flags": ";".join(sorted(set(flags))),
                }
            )

    qc_rows = []
    candidate_rows = []
    read_assignment_rows = []
    for result in results:
        qc = result.qc
        qc_rows.append(
            {
                "sample_id": result.sample.sample_id,
                "source": result.sample.source,
                "is_control": str(result.sample.is_control).lower(),
                "control_type": result.sample.control_type,
                "status": result.status,
                "raw_reads": qc.raw_reads,
                "passed_mean_q": qc.passed_mean_q,
                "primer_matched": qc.primer_matched,
                "passed_length": qc.passed_length,
                "downsampled_from": qc.downsampled_from,
                "retained_reads": qc.retained_reads,
                "forward_reads": qc.forward_reads,
                "reverse_reads": qc.reverse_reads,
                "candidate_events": len(result.variant_events),
                "haplotypes_called": len(result.haplotypes),
                "ambiguous_reads": result.ambiguous_reads,
                "warnings": ";".join(result.warnings),
            }
        )
        for event in result.variant_events:
            candidate_rows.append(
                {
                    "sample_id": result.sample.sample_id,
                    "event": event.label,
                    "type": event.kind,
                    "start_1based": event.start
                    if event.kind == "INS"
                    else event.start + 1,
                    "end_1based": event.start + 1 if event.kind == "INS" else event.end,
                    "ref": event.ref,
                    "alt": event.alt,
                    "support_fraction": round(event.support, 6),
                    "alt_reads": event.alt_reads,
                    "coordinate_reference": "sample_global_consensus",
                }
            )
        records = [
            (
                f"{call.haplotype_id}|sample={result.sample.sample_id}|reads={call.assigned_reads}|fraction={call.fraction:.4f}",
                call.sequence,
            )
            for call in result.haplotypes
        ]
        if records:
            _write_fasta(
                per_sample_dir
                / f"{_safe_name(result.sample.sample_id)}.haplotypes.fasta",
                records,
            )
        for call in result.haplotypes:
            read_assignment_rows.extend(
                {
                    "sample_id": result.sample.sample_id,
                    "read_id": read_name,
                    "assignment": call.haplotype_id,
                }
                for read_name in call.read_names
            )
        read_assignment_rows.extend(
            {
                "sample_id": result.sample.sample_id,
                "read_id": read_name,
                "assignment": "AMBIGUOUS",
            }
            for read_name in result.ambiguous_read_names
        )

    presence_rows: list[dict[str, Any]] = []
    for result in results:
        row: dict[str, Any] = {
            "sample_id": result.sample.sample_id,
            "source": result.sample.source,
            "status": result.status,
        }
        counts = {call.haplotype_id: call.assigned_reads for call in result.haplotypes}
        row.update({haplotype_id: counts.get(haplotype_id, 0) for haplotype_id in ids})
        presence_rows.append(row)

    association_rows: list[dict[str, Any]] = []
    human_results = [
        result
        for result in results
        if result.status == "PASS"
        and not result.sample.is_control
        and result.sample.source == "human"
    ]
    tick_results = [
        result
        for result in results
        if result.status == "PASS"
        and not result.sample.is_control
        and result.sample.source == "tick"
    ]
    p_values = []
    for haplotype_id in ids:
        human_positive = sum(
            any(call.haplotype_id == haplotype_id for call in result.haplotypes)
            for result in human_results
        )
        tick_positive = sum(
            any(call.haplotype_id == haplotype_id for call in result.haplotypes)
            for result in tick_results
        )
        a, b = human_positive, len(human_results) - human_positive
        c, d = tick_positive, len(tick_results) - tick_positive
        p_value = (
            _fisher_two_sided(a, b, c, d) if human_results and tick_results else 1.0
        )
        p_values.append(p_value)
        association_rows.append(
            {
                "haplotype_id": haplotype_id,
                "human_positive": a,
                "human_total": len(human_results),
                "tick_positive": c,
                "tick_total": len(tick_results),
                "odds_ratio_corrected": round(
                    ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)), 6
                ),
                "fisher_p": p_value,
            }
        )
    for row, q_value in zip(association_rows, _bh_adjust(p_values)):
        row["fisher_p"] = f"{row['fisher_p']:.6g}"
        row["bh_q"] = f"{q_value:.6g}"

    network_edges = _minimum_spanning_tree(ids, hap_sequences)
    if references:
        variant_reference_name, variant_reference = references[0]
    elif ids:
        variant_reference_name, variant_reference = ids[0], hap_sequences[ids[0]]
    else:
        variant_reference_name, variant_reference = "none", ""
    variant_rows = []
    if variant_reference:
        for haplotype_id in ids:
            differences = sequence_variants(
                variant_reference, hap_sequences[haplotype_id]
            )
            if not differences:
                variant_rows.append(
                    {
                        "reference": variant_reference_name,
                        "haplotype_id": haplotype_id,
                        "type": "IDENTICAL",
                        "position_1based": "",
                        "ref": "",
                        "alt": "",
                    }
                )
            for difference in differences:
                variant_rows.append(
                    {
                        "reference": variant_reference_name,
                        "haplotype_id": haplotype_id,
                        "type": difference["type"],
                        "position_1based": difference["position"],
                        "ref": difference["ref"],
                        "alt": difference["alt"],
                    }
                )

    _write_fasta(
        output / "haplotypes.fasta", [(item, hap_sequences[item]) for item in ids]
    )
    _write_tsv(
        output / "haplotype_catalog.tsv",
        catalog_rows,
        [
            "haplotype_id",
            "sequence_sha256",
            "length",
            "sample_count",
            "human_samples",
            "tick_samples",
            "negative_control_samples",
            "total_reads",
            "nearest_haplotype",
            "nearest_distance",
            "closest_reference",
            "reference_distance",
            "reference_identity",
        ],
    )
    _write_tsv(
        output / "sample_haplotypes.tsv",
        sample_haplotype_rows,
        [
            "sample_id",
            "source",
            "biological_id",
            "replicate_id",
            "is_control",
            "control_type",
            "haplotype_id",
            "assigned_reads",
            "fraction",
            "forward_reads",
            "reverse_reads",
            "mean_read_edit_rate",
            "replicate_status",
            "event_profile",
            "flags",
        ],
    )
    _write_tsv(
        output / "qc_samples.tsv",
        qc_rows,
        list(qc_rows[0]) if qc_rows else ["sample_id"],
    )
    _write_tsv(
        output / "candidate_events.tsv",
        candidate_rows,
        [
            "sample_id",
            "event",
            "type",
            "start_1based",
            "end_1based",
            "ref",
            "alt",
            "support_fraction",
            "alt_reads",
            "coordinate_reference",
        ],
    )
    _write_tsv(
        output / "haplotype_counts.tsv",
        presence_rows,
        ["sample_id", "source", "status", *ids],
    )
    _write_tsv(
        output / "read_assignments.tsv",
        read_assignment_rows,
        ["sample_id", "read_id", "assignment"],
    )
    _write_tsv(
        output / "association.tsv",
        association_rows,
        [
            "haplotype_id",
            "human_positive",
            "human_total",
            "tick_positive",
            "tick_total",
            "odds_ratio_corrected",
            "fisher_p",
            "bh_q",
        ],
    )
    _write_tsv(
        output / "haplotype_network.tsv",
        network_edges,
        ["source", "target", "distance"],
    )
    _write_tsv(
        output / "haplotype_variants.tsv",
        variant_rows,
        ["reference", "haplotype_id", "type", "position_1based", "ref", "alt"],
    )

    created_at = datetime.now(UTC).isoformat()
    summary = {
        "created_at": created_at,
        "sample_count": len(results),
        "passing_samples": sum(result.status == "PASS" for result in results),
        "haplotype_count": len(ids),
        "retained_reads": sum(result.qc.retained_reads for result in results),
        "reference_name": variant_reference_name,
    }
    write_report(
        output / "report.html",
        summary,
        qc_rows,
        catalog_rows,
        sample_haplotype_rows,
        association_rows,
        network_edges,
    )
    manifest = {
        "spirohap_version": __version__,
        "created_at": created_at,
        "python": sys.version,
        "platform": platform.platform(),
        "jobs": jobs,
        "config": config_as_dict(cfg),
        "sample_sheet": str(Path(samples_path).resolve()) if samples_path else "",
        "sample_sheet_sha256": _sha256_file(Path(samples_path).resolve())
        if samples_path
        else "",
        "config_path": str(Path(config_path).resolve()) if config_path else "defaults",
        "config_sha256": _sha256_file(Path(config_path).resolve())
        if config_path
        else "",
        "reference_fasta": str(Path(reference_fasta).resolve())
        if reference_fasta
        else "",
        "reference_fasta_sha256": _sha256_file(Path(reference_fasta).resolve())
        if reference_fasta
        else "",
        "inputs": [
            {
                "sample_id": sample.sample_id,
                "fastq": [
                    {
                        "path": str(path),
                        "bytes": os.path.getsize(path),
                        "mtime": os.path.getmtime(path),
                    }
                    for path in sample.fastq_paths
                ],
            }
            for sample in samples
        ],
        "summary": summary,
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return summary
