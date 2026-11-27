from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PrimerConfig:
    forward: str = "GGTCAATCGGCATGAAGGATTC"
    reverse: str = "AGCCGGCGAGTTACGATTA"
    search_window: int = 90
    max_edits: int = 4


@dataclass(frozen=True)
class FilterConfig:
    min_mean_q: float = 12.0
    min_base_q: int = 8
    min_trimmed_length: int = 620
    max_trimmed_length: int = 820
    max_reads_per_sample: int = 5000


@dataclass(frozen=True)
class HaplotypeConfig:
    min_site_depth: int = 20
    min_alt_reads: int = 8
    min_snp_fraction: float = 0.05
    min_indel_fraction: float = 0.15
    min_haplotype_reads: int = 10
    min_haplotype_fraction: float = 0.05
    max_haplotypes: int = 6
    max_variant_events: int = 30
    max_profile_mismatches: int = 1
    consensus_fraction: float = 0.60
    polish_rounds: int = 2
    min_strand_fraction: float = 0.05


@dataclass(frozen=True)
class RuntimeConfig:
    random_seed: int = 20260815
    medoid_candidates: int = 30
    medoid_comparisons: int = 80


@dataclass(frozen=True)
class Config:
    primers: PrimerConfig = field(default_factory=PrimerConfig)
    filters: FilterConfig = field(default_factory=FilterConfig)
    haplotypes: HaplotypeConfig = field(default_factory=HaplotypeConfig)
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)


@dataclass(frozen=True)
class Sample:
    sample_id: str
    fastq_paths: tuple[Path, ...]
    source: str = "unknown"
    biological_id: str = ""
    replicate_id: str = ""
    is_control: bool = False
    control_type: str = ""


@dataclass(frozen=True)
class TrimmedRead:
    name: str
    sequence: str
    qualities: tuple[int, ...]
    mean_q: float
    orientation: str


@dataclass
class ReadQC:
    raw_reads: int = 0
    passed_mean_q: int = 0
    primer_matched: int = 0
    passed_length: int = 0
    retained_reads: int = 0
    downsampled_from: int = 0
    forward_reads: int = 0
    reverse_reads: int = 0


@dataclass(frozen=True)
class VariantEvent:
    kind: str
    start: int
    end: int
    ref: str
    alt: str
    support: float
    alt_reads: int

    @property
    def label(self) -> str:
        if self.kind == "SNP":
            return f"{self.start + 1}:{self.ref}>{self.alt}"
        if self.kind == "DEL":
            return f"{self.start + 1}-{self.end}:del"
        return f"{self.start}:{self.start + 1}:ins{self.alt}"


@dataclass
class HaplotypeCall:
    sequence: str
    assigned_reads: int
    fraction: float
    forward_reads: int
    reverse_reads: int
    profile: tuple[str, ...]
    mean_edit_rate: float
    read_names: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    haplotype_id: str = ""


@dataclass
class SampleResult:
    sample: Sample
    qc: ReadQC
    status: str
    global_consensus: str = ""
    variant_events: list[VariantEvent] = field(default_factory=list)
    haplotypes: list[HaplotypeCall] = field(default_factory=list)
    ambiguous_reads: int = 0
    ambiguous_read_names: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
