# SpiroHap

SpiroHap reconstructs full-length **23S amplicon haplotypes** from
demultiplexed Oxford Nanopore FASTQ files. It is configured for the published
~750-bp *Spiroplasma ixodetis* assay but keeps every threshold editable.

It is designed for the specific failure mode that matters in this dataset: a
single majority consensus can become a mosaic when one specimen contains two
real alleles. SpiroHap first finds candidate SNP/indel sites, phases those sites
on individual full-length reads, groups the read-level profiles, and then
polishes a separate consensus for each supported group.

> The output is a set of one-locus 23S alleles, not a set of whole-genome
> strains. It is a research workflow and is not a clinical diagnostic device.

## Quick start

Python 3.11 or newer is required.

```bash
unzip spirohap-pipeline.zip
cd spirohap-pipeline
python3 -m venv .venv
source .venv/bin/activate            # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install .
spirohap doctor
```

Run the included synthetic truth set first:

```bash
spirohap demo --out demo_run
```

This creates a pure tick sample, a 2-haplotype tick mixture, a human sample,
and a negative control. The automated check used during development recovered
both truth haplotypes exactly.

Then copy and edit the two templates:

```bash
cp samples.example.csv samples.csv
cp config.example.toml config.toml
spirohap run \
  --samples samples.csv \
  --config config.toml \
  --out results_23s \
  --jobs 4
```

The output folder must be new or empty; the program will not overwrite an old
analysis.

## Input

The pipeline expects demultiplexed `.fastq` or `.fastq.gz` files with the
amplicon primers still present. It accepts one or more files per sample. Repeat
the same `sample_id` on multiple rows to pool lanes/files, or separate paths in
one `fastq` cell with semicolons. Globs are accepted.

Required sample-sheet columns:

| Column | Meaning |
|---|---|
| `sample_id` | Unique analysis unit; independent PCRs must have distinct IDs. |
| `fastq` | FASTQ path relative to the sheet, absolute path, glob, or semicolon-separated paths. |

Recommended columns:

| Column | Meaning |
|---|---|
| `source` | Use exactly `human`, `tick`, or `control`; only passing human/tick samples enter the exploratory comparison. |
| `biological_id` | Original specimen ID shared by independent PCR repeats. |
| `replicate_id` | Distinct repeat, such as `pcr1` and `pcr2`. |
| `is_control` | `true`/`false`. |
| `control_type` | `negative`, `positive`, or blank. A nonblank value implies `is_control=true`. |

The default primers are:

- forward `GGTCAATCGGCATGAAGGATTC`
- reverse `AGCCGGCGAGTTACGATTA`

These are the Sixo23S sequencing primers reported for the *S. ixodetis* 23S
amplicon. The length filter applies **after both primers are removed**, which is
why the default range is 620–820 bp rather than exactly 750 bp.

If you start from POD5, see [docs/dorado.md](docs/dorado.md). Requiring barcodes
at both ends is recommended when the library design permits it, because rare
barcode leakage can mimic a minor haplotype.

## What the algorithm does

1. Filters reads by mean Q score.
2. searches both orientations for approximate forward and reverse primer
   matches, orients every accepted read identically, and removes the primers.
3. Applies the primer-trimmed length filter and deterministic reservoir
   downsampling.
4. Selects a representative read and iteratively builds a global sample
   consensus.
5. Aligns every read to that consensus and detects candidate SNPs, deletion
   blocks, and insertions above separately configurable thresholds.
6. Converts each full-length read into a phased allele profile across all
   candidate events. It does **not** create 97% OTUs, which would merge
   biologically interesting one-base variants.
7. Finds sufficiently abundant read profiles and independently polishes a
   consensus from the reads assigned to each profile.
8. Dereplicates only exact polished sequences across samples, gives them stable
   run-local IDs such as `SIXO23S-H001`, and builds a minimum-spanning network
   from pairwise edit distances.
9. Flags strand imbalance, one-event separation, negative-control detection,
   high residual read error, and replicate support.

Alignment is performed with [Edlib](https://github.com/Martinsos/edlib), using
global edit-distance alignment for full-length reads and infix alignment for
primer finding.

## Main outputs

| File | Contents |
|---|---|
| `report.html` | Self-contained overview, QC tables, sample calls, network, and exploratory comparison. |
| `haplotypes.fasta` | Exact primer-trimmed consensus sequences. |
| `sample_haplotypes.tsv` | Per-sample read counts, fractions, strand support, profiles, replicate status, and flags. |
| `haplotype_catalog.tsv` | Cohort prevalence, sequence hashes, nearest neighbor, and optional reference match. |
| `qc_samples.tsv` | Filtering counts and sample status. |
| `candidate_events.tsv` | Within-sample candidate events used for phasing; coordinates are relative to that sample's global consensus. |
| `haplotype_counts.tsv` | Sample × haplotype assigned-read matrix. |
| `read_assignments.tsv` | Read IDs assigned to each haplotype, plus ambiguous read IDs, for pileup review or FASTQ extraction. |
| `haplotype_variants.tsv` | SNP/indel differences against the first supplied reference, or against `SIXO23S-H001`. Indels are not left-normalized; an insertion position is its left reference anchor, so 0 means before base 1. |
| `haplotype_network.tsv` | Minimum-spanning-network edges and edit distances. |
| `association.tsv` | Exploratory human/tick presence comparison with Fisher p and BH q values. |
| `run_manifest.json` | Software/configuration provenance and input metadata. |

You can supply known, **primer-trimmed** amplicon sequences:

```bash
spirohap run --samples samples.csv --config config.toml \
  --references known_23s_amplicons.fasta --out results_23s
```

The first FASTA record defines the coordinate system in
`haplotype_variants.tsv`; all records are used to report the closest reference.

## Tune thresholds empirically

The defaults are conservative starting values for modern ONT reads, not
universal truth. Before biological interpretation, sequence:

- single-template controls for the dominant expected alleles;
- two-template mixtures at known fractions (for example 50%, 20%, 10%, 5%,
  2%, and 1%);
- extraction blanks and PCR no-template controls;
- independently amplified replicates.

Choose the lowest `min_haplotype_fraction` and `min_alt_reads` that recover the
mixtures without producing haplotypes in single-template controls. Indels,
especially near homopolymers, usually need a higher threshold than SNPs. If a
sample produces many candidate events or reaches `max_haplotypes`, inspect it
manually instead of merely relaxing the limits.

The pipeline never silently deletes a sequence found in a negative control; it
marks every occurrence `detected_in_negative_control` so that you can interpret
it with batch position, read fraction, and control type.

See [docs/validation.md](docs/validation.md) before claiming a human-associated
variant.

## Interpretation for the 200-tick / 3-human study

With three human samples, formal variant association is severely underpowered.
Use `association.tsv` to prioritize alleles for validation, not to label a
variant as human-adapted. A strong candidate should be exact across independent
PCRs, supported in both orientations, absent from negative controls, and
confirmed at additional independent loci or by genome-scale sequencing.

An important biological alternative is intragenomic rRNA-copy variation: one
*S. ixodetis* genome may contribute more than one 23S allele. Therefore, two
amplicon haplotypes in a tick do not by themselves prove co-infection by two
strains.

## Reproducibility and tests

```bash
python -m pip install -e '.[dev]'
pytest -q
```

Or run in Docker:

```bash
docker build -t spirohap .
docker run --rm -v "$PWD:/data" spirohap run \
  --samples /data/samples.csv --config /data/config.toml \
  --out /data/results_23s --jobs 4
```

## Scientific context

- Korva M. et al. *Emerging human pathogen: identifying Spiroplasma ixodetis
  as a frequent cause of unlocalised febrile illness.* Journal of Infection
  93, 106776 (2026). <https://doi.org/10.1016/j.jinf.2026.106776>
- Buysse M. et al. *A human-associated Spiroplasma ixodetis lineage responsible
  for infantile cataracts and adult febrile illness.* iScience (2026).
  <https://doi.org/10.1016/j.isci.2026.115233>
- Oxford Nanopore Dorado barcode documentation:
  <https://software-docs.nanoporetech.com/dorado/latest/barcoding/barcoding/>

The second study is especially relevant to follow-up design: a one-locus 23S
screen should be followed by multilocus or genome-scale work before asserting
that the observed allele represents the reported human-associated lineage.
