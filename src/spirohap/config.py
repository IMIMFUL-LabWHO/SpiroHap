from __future__ import annotations

import dataclasses
import tomllib
from pathlib import Path
from typing import Any, TypeVar

from .model import Config, FilterConfig, HaplotypeConfig, PrimerConfig, RuntimeConfig

T = TypeVar("T")


def _section(cls: type[T], raw: dict[str, Any], name: str) -> T:
    allowed = {field.name for field in dataclasses.fields(cls)}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"Unknown keys in [{name}]: {', '.join(unknown)}")
    return cls(**raw)


def _validate(cfg: Config) -> None:
    for label, seq in (
        ("forward primer", cfg.primers.forward),
        ("reverse primer", cfg.primers.reverse),
    ):
        if not seq or set(seq.upper()) - set("ACGT"):
            raise ValueError(f"{label} must contain only A, C, G, and T")
    if cfg.primers.search_window < max(
        len(cfg.primers.forward), len(cfg.primers.reverse)
    ):
        raise ValueError("primers.search_window is shorter than a primer")
    if cfg.filters.min_trimmed_length >= cfg.filters.max_trimmed_length:
        raise ValueError("min_trimmed_length must be smaller than max_trimmed_length")
    if cfg.filters.max_reads_per_sample < cfg.haplotypes.min_site_depth:
        raise ValueError("max_reads_per_sample must be >= min_site_depth")
    for name in (
        "min_snp_fraction",
        "min_indel_fraction",
        "min_haplotype_fraction",
        "consensus_fraction",
    ):
        value = getattr(cfg.haplotypes, name)
        if not 0 < value <= 1:
            raise ValueError(f"haplotypes.{name} must be in (0, 1]")
    if cfg.haplotypes.min_alt_reads < 2 or cfg.haplotypes.min_haplotype_reads < 2:
        raise ValueError("read-count thresholds must be at least 2")


def load_config(path: str | Path | None) -> Config:
    raw: dict[str, Any] = {}
    if path is not None:
        with Path(path).open("rb") as handle:
            raw = tomllib.load(handle)
    allowed = {"primers", "filters", "haplotypes", "runtime"}
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"Unknown configuration sections: {', '.join(unknown)}")
    cfg = Config(
        primers=_section(PrimerConfig, raw.get("primers", {}), "primers"),
        filters=_section(FilterConfig, raw.get("filters", {}), "filters"),
        haplotypes=_section(HaplotypeConfig, raw.get("haplotypes", {}), "haplotypes"),
        runtime=_section(RuntimeConfig, raw.get("runtime", {}), "runtime"),
    )
    _validate(cfg)
    return cfg


def config_as_dict(cfg: Config) -> dict[str, Any]:
    return dataclasses.asdict(cfg)
