from __future__ import annotations

import csv
import random
from pathlib import Path

from .io import reverse_complement

FORWARD = "GGTCAATCGGCATGAAGGATTC"
REVERSE = "AGCCGGCGAGTTACGATTA"


def _reference_sequence() -> str:
    rng = random.Random(1701)
    sequence = []
    previous = ""
    run = 0
    for _ in range(700):
        choices = list("ACGT")
        if run >= 3 and previous in choices:
            choices.remove(previous)
        base = rng.choice(choices)
        if base == previous:
            run += 1
        else:
            previous, run = base, 1
        sequence.append(base)
    return "".join(sequence)


def _mutant(reference: str) -> str:
    sequence = list(reference)
    for position in (93, 211, 418, 612):
        sequence[position] = next(base for base in "ACGT" if base != sequence[position])
    sequence[355:355] = list("GACCTTGA")
    return "".join(sequence)


def _introduce_errors(sequence: str, rng: random.Random) -> str:
    output = []
    for base in sequence:
        draw = rng.random()
        if draw < 0.002:
            continue
        if draw < 0.005:
            output.append(rng.choice("ACGT"))
        if draw < 0.011:
            output.append(
                rng.choice([candidate for candidate in "ACGT" if candidate != base])
            )
        else:
            output.append(base)
    return "".join(output)


def _fastq_record(name: str, insert: str, rng: random.Random) -> str:
    noisy_insert = _introduce_errors(insert, rng)
    sequence = FORWARD + noisy_insert + reverse_complement(REVERSE)
    if rng.random() < 0.5:
        sequence = reverse_complement(sequence)
    quality = "5" * len(sequence)  # Phred 20
    return f"@{name}\n{sequence}\n+\n{quality}\n"


def make_demo(directory: str | Path) -> dict[str, str]:
    root = Path(directory).resolve()
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Demo directory is not empty: {root}")
    fastq_dir = root / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    hap1 = _reference_sequence()
    hap2 = _mutant(hap1)
    rng = random.Random(424242)
    designs = {
        "tick_01": [hap1] * 80,
        "tick_02": [hap1] * 60 + [hap2] * 30,
        "human_01": [hap2] * 80,
        "negative_control": [],
    }
    for sample_id, inserts in designs.items():
        rng.shuffle(inserts)
        content = "".join(
            _fastq_record(f"{sample_id}_read_{index:04d}", insert, rng)
            for index, insert in enumerate(inserts, 1)
        )
        (fastq_dir / f"{sample_id}.fastq").write_text(content, encoding="utf-8")
    sample_sheet = root / "samples.csv"
    with sample_sheet.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sample_id",
                "fastq",
                "source",
                "biological_id",
                "replicate_id",
                "is_control",
                "control_type",
            ]
        )
        writer.writerow(
            ["tick_01", "fastq/tick_01.fastq", "tick", "tick_01", "pcr1", "false", ""]
        )
        writer.writerow(
            ["tick_02", "fastq/tick_02.fastq", "tick", "tick_02", "pcr1", "false", ""]
        )
        writer.writerow(
            [
                "human_01",
                "fastq/human_01.fastq",
                "human",
                "human_01",
                "pcr1",
                "false",
                "",
            ]
        )
        writer.writerow(
            [
                "negative_control",
                "fastq/negative_control.fastq",
                "control",
                "",
                "",
                "true",
                "negative",
            ]
        )
    config = root / "config.toml"
    config.write_text(
        """[primers]
forward = "GGTCAATCGGCATGAAGGATTC"
reverse = "AGCCGGCGAGTTACGATTA"
search_window = 90
max_edits = 4

[filters]
min_mean_q = 12.0
min_base_q = 8
min_trimmed_length = 620
max_trimmed_length = 820
max_reads_per_sample = 5000

[haplotypes]
min_site_depth = 20
min_alt_reads = 5
min_snp_fraction = 0.08
min_indel_fraction = 0.10
min_haplotype_reads = 6
min_haplotype_fraction = 0.08
max_haplotypes = 6
max_variant_events = 30
max_profile_mismatches = 1
consensus_fraction = 0.60
polish_rounds = 2
min_strand_fraction = 0.05

[runtime]
random_seed = 20260815
medoid_candidates = 20
medoid_comparisons = 50
""",
        encoding="utf-8",
    )
    truth = root / "truth.fasta"
    truth.write_text(
        f">demo_haplotype_1\n{hap1}\n>demo_haplotype_2\n{hap2}\n", encoding="utf-8"
    )
    return {"samples": str(sample_sheet), "config": str(config), "truth": str(truth)}
