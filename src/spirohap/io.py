from __future__ import annotations

import csv
import glob
import gzip
import random
from collections.abc import Iterator
from pathlib import Path

import edlib

from .model import Config, ReadQC, Sample, TrimmedRead

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def read_fastq(path: Path) -> Iterator[tuple[str, str, tuple[int, ...]]]:
    with _open_text(path) as handle:
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            line_number += 4
            if not quality:
                raise ValueError(
                    f"Truncated FASTQ record ending near line {line_number} in {path}"
                )
            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n").upper()
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(
                    f"Malformed FASTQ record ending near line {line_number} in {path}"
                )
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence/quality length mismatch near line {line_number} in {path}"
                )
            yield (
                header[1:].split()[0],
                sequence,
                tuple(ord(char) - 33 for char in quality),
            )


def _parse_bool(value: str) -> bool:
    normal = value.strip().lower()
    if normal in {"", "0", "false", "no", "n"}:
        return False
    if normal in {"1", "true", "yes", "y"}:
        return True
    raise ValueError(f"Invalid Boolean value: {value!r}")


def read_sample_sheet(path: str | Path) -> list[Sample]:
    sheet = Path(path).resolve()
    with sheet.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "fastq"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Sample sheet requires columns: sample_id, fastq")
        grouped: dict[str, dict[str, object]] = {}
        order: list[str] = []
        for row_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            pattern_cell = (row.get("fastq") or "").strip()
            if not sample_id or not pattern_cell:
                raise ValueError(
                    f"Missing sample_id or fastq at sample-sheet row {row_number}"
                )
            metadata = {
                "source": (row.get("source") or "unknown").strip().lower() or "unknown",
                "biological_id": (row.get("biological_id") or "").strip(),
                "replicate_id": (row.get("replicate_id") or "").strip(),
                "is_control": _parse_bool(row.get("is_control") or ""),
                "control_type": (row.get("control_type") or "").strip().lower(),
            }
            if metadata["control_type"] not in {"", "negative", "positive"}:
                raise ValueError(
                    f"control_type must be negative, positive, or blank at row {row_number}"
                )
            if metadata["control_type"]:
                metadata["is_control"] = True
            if sample_id not in grouped:
                grouped[sample_id] = {**metadata, "paths": []}
                order.append(sample_id)
            else:
                for key, value in metadata.items():
                    if grouped[sample_id][key] != value:
                        raise ValueError(
                            f"Conflicting {key} for repeated sample_id {sample_id!r}"
                        )
            paths: list[Path] = grouped[sample_id]["paths"]  # type: ignore[assignment]
            for part in pattern_cell.split(";"):
                candidate = Path(part.strip())
                if not candidate.is_absolute():
                    candidate = sheet.parent / candidate
                matches = sorted(
                    Path(item).resolve() for item in glob.glob(str(candidate))
                )
                if not matches:
                    raise FileNotFoundError(
                        f"No FASTQ matched {part!r} for sample {sample_id!r}"
                    )
                paths.extend(matches)
        samples = []
        for sample_id in order:
            item = grouped[sample_id]
            samples.append(
                Sample(
                    sample_id=sample_id,
                    fastq_paths=tuple(dict.fromkeys(item["paths"])),  # type: ignore[arg-type]
                    source=str(item["source"]),
                    biological_id=str(item["biological_id"]),
                    replicate_id=str(item["replicate_id"]),
                    is_control=bool(item["is_control"]),
                    control_type=str(item["control_type"]),
                )
            )
    if not samples:
        raise ValueError("Sample sheet contains no samples")
    return samples


def _best_primer_match(primer: str, window: str, prefer_right: bool, max_edits: int):
    result = edlib.align(primer, window, mode="HW", task="locations", k=max_edits)
    if result["editDistance"] < 0 or not result["locations"]:
        return None
    locations = result["locations"]
    start, end = (
        max(locations, key=lambda loc: (loc[1], loc[0]))
        if prefer_right
        else min(locations)
    )
    return int(result["editDistance"]), int(start), int(end)


def _trim_forward(sequence: str, qualities: tuple[int, ...], cfg: Config):
    window_size = cfg.primers.search_window
    left_window = sequence[:window_size]
    right_offset = max(0, len(sequence) - window_size)
    right_window = sequence[right_offset:]
    left = _best_primer_match(
        cfg.primers.forward.upper(),
        left_window,
        prefer_right=False,
        max_edits=cfg.primers.max_edits,
    )
    right = _best_primer_match(
        reverse_complement(cfg.primers.reverse.upper()),
        right_window,
        prefer_right=True,
        max_edits=cfg.primers.max_edits,
    )
    if left is None or right is None:
        return None
    left_cut = left[2] + 1
    right_cut = right_offset + right[1]
    if right_cut <= left_cut:
        return None
    return (
        left[0] + right[0],
        sequence[left_cut:right_cut],
        qualities[left_cut:right_cut],
    )


def trim_and_orient(
    name: str, sequence: str, qualities: tuple[int, ...], cfg: Config
) -> TrimmedRead | None:
    forward = _trim_forward(sequence, qualities, cfg)
    rc_sequence = reverse_complement(sequence)
    rc_qualities = qualities[::-1]
    reverse = _trim_forward(rc_sequence, rc_qualities, cfg)
    choices = []
    if forward is not None:
        choices.append((forward[0], "+", forward[1], forward[2]))
    if reverse is not None:
        choices.append((reverse[0], "-", reverse[1], reverse[2]))
    if not choices:
        return None
    _, orientation, trimmed, trimmed_q = min(
        choices, key=lambda item: (item[0], item[1] == "-")
    )
    mean_q = sum(qualities) / len(qualities) if qualities else 0.0
    return TrimmedRead(name, trimmed, tuple(trimmed_q), mean_q, orientation)


def load_sample_reads(
    sample: Sample, cfg: Config, seed: int
) -> tuple[list[TrimmedRead], ReadQC]:
    qc = ReadQC()
    retained: list[TrimmedRead] = []
    rng = random.Random(seed)
    eligible = 0
    for path in sample.fastq_paths:
        for name, sequence, qualities in read_fastq(path):
            qc.raw_reads += 1
            mean_q = sum(qualities) / len(qualities) if qualities else 0.0
            if mean_q < cfg.filters.min_mean_q:
                continue
            qc.passed_mean_q += 1
            trimmed = trim_and_orient(name, sequence, qualities, cfg)
            if trimmed is None:
                continue
            qc.primer_matched += 1
            if not (
                cfg.filters.min_trimmed_length
                <= len(trimmed.sequence)
                <= cfg.filters.max_trimmed_length
            ):
                continue
            qc.passed_length += 1
            eligible += 1
            if len(retained) < cfg.filters.max_reads_per_sample:
                retained.append(trimmed)
            else:
                replacement = rng.randrange(eligible)
                if replacement < len(retained):
                    retained[replacement] = trimmed
    qc.downsampled_from = eligible
    qc.retained_reads = len(retained)
    qc.forward_reads = sum(read.orientation == "+" for read in retained)
    qc.reverse_reads = sum(read.orientation == "-" for read in retained)
    return retained, qc


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    with Path(path).open("rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0], []
            else:
                if not name:
                    raise ValueError(f"FASTA sequence before header in {path}")
                chunks.append(line)
    if name:
        records.append((name, "".join(chunks).upper()))
    return records
