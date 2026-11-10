# SpliceAI Interpret

A Python CLI that maps precomputed SpliceAI delta scores onto a selected transcript and reports candidate splice events. It supports exon skipping, full and partial intron retention, partial exon skipping, pseudoexon inclusion, and cryptic splice-site activation.

## Quick start

Use Python 3.10 or newer. The CLI has no third-party Python dependencies. Run this from the repository root:

```bash
python3 spliceai_interpret.py \
  --vcf examples/synthetic.vcf \
  --gtf examples/synthetic.gtf \
  --transcript DEMO_TX.1 \
  --fasta examples/synthetic.fa \
  --html splice-report.html > splice-events.jsonl
```

Open `splice-report.html` in a browser to see the transcript, variant marker, and gain and loss score tracks. `splice-events.jsonl` contains one JSON record per matching SpliceAI annotation, with all candidate events for that record. Use a new filename for each HTML report.

The [example files](examples/) include an invented contig, a four-exon transcript, a matching FASTA index, and three VCF records with hand-assigned scores. No downloads are needed to run the example.

| Synthetic variant | Candidate event | Length change | Modeled coding result |
| --- | --- | ---: | --- |
| `chrDemo:42 A>G` | Cryptic donor with partial intron retention | +1 nt | Frameshift, with no stop found in the modeled sequence |
| `chrDemo:60 C>T` | Exon skipping | −6 nt | In frame. The retained SNV creates a premature stop 78 nt upstream of the final junction, so NMD is predicted |
| `chrDemo:147 A>G` | Pseudoexon inclusion | +6 nt | In frame, with the reference stop retained |

The donor and pseudoexon examples also emit individual cryptic-site signals.

## Use with your data

```bash
samtools faidx GRCh38.fa
python3 spliceai_interpret.py \
  --vcf spliceai-scored.vcf.gz \
  --gtf mane.gtf.gz \
  --transcript NM_000000.1 \
  --fasta GRCh38.fa \
  --min-score 0.2 > splice-events.jsonl
```

Replace `NM_000000.1` with the exact `transcript_id` from your GTF, including its version when present.

| Input | Requirements |
| --- | --- |
| `--vcf` | Plain VCF or VCF.gz containing the standard `INFO/SpliceAI` field. The CLI uses precomputed scores. |
| `--gtf` | GTF or GTF.gz with exon records and a `gene_name` matching the SpliceAI gene symbol. Use the same assembly as the VCF. |
| `--transcript` | One exact GTF transcript ID per run. MANE Select is one way to choose a transcript. |
| `--fasta` | Optional, uncompressed reference FASTA with an adjacent `.fai`. Use the same assembly as the VCF and GTF. |
| `--min-score` | Minimum delta score for each signal used in an event. The default is 0.2. |
| `--html` | Optional standalone HTML report. The CLI creates a new file. |

## Read the output

- `sites` shows the available delta scores, their original offsets, their genomic positions, and whether they meet the threshold. Positions are 1-based and event intervals are inclusive.
- `events` lists independent splice hypotheses and their supporting signals. `reference_length_change_nt` is the length change from the modeled splice edit on the reference.
- Frame and NMD require an indexed FASTA, an SNV, a complete nuclear CDS with start and stop codon annotations, and a modeled event within the CDS. When the script finds a premature stop, it reports its position and estimates NMD using a distance of more than 50 nt from the stop's end to the final modeled exon junction.

SpliceAI's standard VCF annotation reports one peak per signal type at the gene level. This CLI projects those peaks onto the chosen transcript. Treat the resulting events as hypotheses for review. The NMD rule is a heuristic and does not account for gene- or transcript-specific exceptions.

[SpliceAI VCF format](https://github.com/Illumina/SpliceAI) · [ClinGen SVI splicing recommendations](https://doi.org/10.1016/j.ajhg.2023.06.002)
