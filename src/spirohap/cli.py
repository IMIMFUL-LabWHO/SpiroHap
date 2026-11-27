from __future__ import annotations

import argparse
import importlib.metadata
import os
import sys
from pathlib import Path

from . import __version__
from .config import load_config
from .demo import make_demo
from .io import read_sample_sheet
from .pipeline import run_pipeline


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spirohap",
        description="Reconstruct full-length 23S amplicon haplotypes from demultiplexed Nanopore FASTQ files.",
    )
    parser.add_argument(
        "--version", action="version", version=f"spirohap {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="run the pipeline")
    run.add_argument("--samples", required=True, help="CSV sample sheet")
    run.add_argument(
        "--config", help="TOML configuration; defaults are used if omitted"
    )
    run.add_argument("--out", required=True, help="new or empty output directory")
    run.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    run.add_argument(
        "--references",
        help="optional primer-trimmed reference amplicons in FASTA; the first record sets variant coordinates",
    )

    demo = subparsers.add_parser(
        "demo", help="generate and analyze a small synthetic two-haplotype dataset"
    )
    demo.add_argument("--out", required=True, help="new or empty demo directory")
    demo.add_argument("--jobs", type=int, default=1)

    subparsers.add_parser("doctor", help="show dependency versions")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "doctor":
            print(f"spirohap\t{__version__}")
            print(f"python\t{sys.version.split()[0]}")
            print(f"edlib\t{importlib.metadata.version('edlib')}")
            return 0
        if args.command == "demo":
            paths = make_demo(args.out)
            samples = read_sample_sheet(paths["samples"])
            cfg = load_config(paths["config"])
            results_dir = Path(args.out).resolve() / "results"
            summary = run_pipeline(
                samples,
                cfg,
                results_dir,
                jobs=args.jobs,
                reference_fasta=paths["truth"],
                samples_path=paths["samples"],
                config_path=paths["config"],
            )
            print(f"Demo completed: {results_dir / 'report.html'}")
            print(f"Recovered {summary['haplotype_count']} unique haplotypes")
            return 0
        if args.jobs < 1:
            raise ValueError("--jobs must be at least 1")
        samples = read_sample_sheet(args.samples)
        cfg = load_config(args.config)
        summary = run_pipeline(
            samples,
            cfg,
            args.out,
            jobs=args.jobs,
            reference_fasta=args.references,
            samples_path=args.samples,
            config_path=args.config,
        )
        print(
            f"Completed {summary['sample_count']} samples; report: {Path(args.out).resolve() / 'report.html'}"
        )
        return 0
    except (OSError, ValueError) as error:
        print(f"spirohap: error: {error}", file=sys.stderr)
        return 2
