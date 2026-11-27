from __future__ import annotations

import csv

from spirohap.config import load_config
from spirohap.demo import make_demo
from spirohap.io import read_fasta, read_sample_sheet, reverse_complement
from spirohap.pipeline import run_pipeline


def test_reverse_complement():
    assert reverse_complement("ACGTN") == "NACGT"


def test_demo_recovers_truth(tmp_path):
    demo = tmp_path / "demo"
    paths = make_demo(demo)
    samples = read_sample_sheet(paths["samples"])
    summary = run_pipeline(
        samples,
        load_config(paths["config"]),
        demo / "results",
        jobs=1,
        reference_fasta=paths["truth"],
        samples_path=paths["samples"],
        config_path=paths["config"],
    )
    truth = {sequence for _, sequence in read_fasta(paths["truth"])}
    reconstructed = {
        sequence for _, sequence in read_fasta(demo / "results" / "haplotypes.fasta")
    }
    assert reconstructed == truth
    assert summary["haplotype_count"] == 2

    with (demo / "results" / "sample_haplotypes.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mixed = [row for row in rows if row["sample_id"] == "tick_02"]
    assert len(mixed) == 2
    assert sorted(int(row["assigned_reads"]) for row in mixed) == [28, 60]

    report = (demo / "results" / "report.html").read_text(encoding="utf-8")
    assert "23S amplicon report" in report
    assert "not whole-genome strains" in report
