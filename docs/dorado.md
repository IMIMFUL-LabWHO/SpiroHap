# Optional Dorado basecalling and demultiplexing

SpiroHap starts from demultiplexed FASTQ. If you have POD5 files, use a current
Dorado release and substitute the exact native-barcoding kit used in your run:

```bash
dorado basecaller sup pod5/ \
  --kit-name SQK-NBD114-96 \
  --barcode-both-ends \
  > calls.bam

dorado demux calls.bam \
  --no-classify \
  --emit-fastq \
  --emit-summary \
  --output-dir demux/
```

`--barcode-both-ends` is deliberately conservative. Oxford Nanopore documents
that it sharply reduces false-positive barcode assignments, with a lower
classification rate as the trade-off. This matters here because a few
cross-assigned reads can look like a rare haplotype.

Do not copy the kit name blindly. Confirm it from the run metadata. If barcode
classification was not done during basecalling, follow Dorado's "classifying
existing datasets" route instead of using `--no-classify`.

Documentation: <https://software-docs.nanoporetech.com/dorado/latest/barcoding/barcoding/>

