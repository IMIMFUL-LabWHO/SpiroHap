# Validation and interpretation checklist

## Before trusting a low-frequency allele

- Re-extract and repeat PCR from the original specimen in a separately prepared
  reaction. Give the repeats different `sample_id` values but the same
  `biological_id`, with distinct `replicate_id` values.
- Require support in both read orientations. Investigate
  `one_strand_only` and `strand_bias` flags.
- Inspect every extraction blank, PCR no-template control, and barcode control.
  A call marked `detected_in_negative_control` is not automatically removed,
  but it is not evidence of infection.
- Review the read pileup around homopolymers and every one-event separation.
  Nanopore systematic errors can recur across many reads.
- Look for recombinant profiles expected from late-cycle PCR chimeras. Reduce
  cycle number and repeat the amplification if a candidate looks like a mosaic
  of two abundant parent alleles.
- Re-run a known two-template mixture across a dilution series. Use the lowest
  frequency reliably recovered without false calls in blanks to set
  `min_haplotype_fraction`, `min_alt_reads`, and the SNP/indel thresholds.
- Confirm priority variants by an orthogonal method: allele-specific or droplet
  digital PCR, short-read sequencing, or Sanger sequencing of a clean isolate.

## Biological limits

The output sequences are **23S amplicon alleles**. They are not automatically
organismal strains or transmission clusters. Multiple rRNA operons can coexist
in one genome, and PCR can preferentially amplify one operon. A mixed call can
therefore mean co-infection, intragenomic rRNA-copy variation, contamination,
or an artefact.

For a human-associated lineage claim, add independent loci (for example 16S
plus protein-coding housekeeping loci) or, preferably, obtain isolate or
metagenomic genomes. Compare cases and ticks with temporal/geographic matching,
pre-specified QC, blinded reanalysis, and an independent validation cohort.

With only three human specimens, `association.tsv` is a descriptive triage
table. It is not powered to identify a causal human-infection variant.
