#!/usr/bin/env python3
"""Interpret precomputed SpliceAI VCF annotations against one versioned GTF transcript.

Example:
    python spliceai_interpret.py --vcf scored.vcf.gz --gtf mane.gtf.gz \
        --transcript NM_000000.1 --fasta GRCh38.fa --html splice.html > events.jsonl

Use the same assembly for VCF, GTF, and FASTA. Run samtools faidx on the
uncompressed FASTA first. Gene matching uses the SpliceAI gene symbol and GTF
gene_name. A gene-level score is projected onto the selected transcript.

Events are independent hypotheses derived from the four reported peak scores.
Positions are 1-based and intervals inclusive. The length change describes the
reference splice edit. Coding reconstruction supports SNVs and complete nuclear
CDSs. The NMD output uses the distance from the end of the modeled premature
stop to the final modeled exon junction (>50 nt). It is a junction-based
prediction and does not account for transcript- or gene-specific exceptions.
"""

import argparse
from contextlib import ExitStack
from dataclasses import dataclass
import gzip
from html import escape
import json
import math
from pathlib import Path
import re
import sys


SITE_TYPES = ("AG", "AL", "DG", "DL")
STOPS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def contig_key(name):
    return name.removeprefix("chr")


@dataclass
class Transcript:
    id: str
    gene: str
    chrom: str
    strand: str
    # Signed, 1-based inclusive coordinates increase in transcript order.
    exons: list
    start: int | None
    stop: int | None

    @property
    def direction(self):
        return 1 if self.strand == "+" else -1


def read_transcript(path, transcript_id):
    exons = []
    starts, stops, cds = set(), set(), set()
    loci, genes = set(), set()
    unsupported_translation = False
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, low, high, _, strand, _, attrs = fields
            values = dict(re.findall(r'(\w+)\s+"([^"]+)"', attrs))
            if values.get("transcript_id") != transcript_id:
                continue
            loci.add((chrom, strand))
            if "gene_name" in values:
                genes.add(values["gene_name"])
            unsupported_translation |= (
                feature == "Selenocysteine"
                or bool(re.search(r'tag\s+"(?:cds_start_NF|cds_end_NF)"', attrs))
            )
            if feature not in {"exon", "start_codon", "stop_codon", "CDS"}:
                continue
            low, high = int(low), int(high)
            if not 1 <= low <= high:
                raise ValueError("Invalid GTF feature coordinates")
            start, end = (low, high) if strand == "+" else (-high, -low)
            if feature == "exon":
                exons.append((start, end))
            else:
                target = {"start_codon": starts, "stop_codon": stops, "CDS": cds}[feature]
                target.update(range(start, end + 1))
    if len(loci) != 1 or len(genes) != 1 or not exons:
        raise ValueError(f"Expected one locus with exons and gene_name for {transcript_id}")
    chrom, strand = loci.pop()
    exons = sorted(set(exons))
    if strand not in {"+", "-"} or any(a[1] + 1 >= b[0] for a, b in zip(exons, exons[1:])):
        raise ValueError("Expected stranded exons separated by introns")
    start, stop = None, None
    if len(starts) == len(stops) == 3 and not unsupported_translation and contig_key(chrom) not in {"M", "MT"}:
        first, last = min(starts), min(stops)
        coding = [p for lo, hi in exons for p in range(max(lo, first), min(hi, max(stops)) + 1)]
        if (len(coding) >= 6 and set(coding[:3]) == starts and set(coding[-3:]) == stops
                and set(coding[:-3]) == cds - stops):
            start, stop = first, last
    return Transcript(transcript_id, genes.pop(), chrom, strand, exons, start, stop)


class Reference:
    """Fetch subsequences from an uncompressed FASTA using its samtools .fai."""

    def __init__(self, handle, index_path, chrom):
        self.handle = handle
        if handle.read(2) == b"\x1f\x8b":
            raise ValueError("Provide an uncompressed reference FASTA")
        with open(index_path, encoding="utf-8") as index:
            matches = [row.split() for row in index if contig_key(row.split()[0]) == contig_key(chrom)]
        exact = [row for row in matches if row[0] == chrom]
        matches = exact or matches
        if len(matches) != 1:
            raise ValueError(f"Expected one FASTA index entry for {chrom}")
        self.length, self.offset, self.width, self.line_bytes = map(int, matches[0][1:5])
        if not (self.length > 0 and self.offset >= 0 and 0 < self.width <= self.line_bytes):
            raise ValueError("Invalid FASTA index")

    def fetch(self, low, high):
        if not 1 <= low <= high <= self.length:
            raise ValueError("Coordinates fall outside the reference FASTA")

        def offset(position):
            line, column = divmod(position - 1, self.width)
            return self.offset + line * self.line_bytes + column

        begin = offset(low)
        self.handle.seek(begin)
        seq = b"".join(self.handle.read(offset(high) - begin + 1).split()).decode("ascii").upper()
        if len(seq) != high - low + 1:
            raise ValueError("FASTA sequence length does not match its index")
        return seq


def parse_sites(annotation, pos, threshold):
    fields = annotation.split("|")
    if len(fields) != 10:
        raise ValueError("Expected the ten standard SpliceAI annotation fields")
    sites = {}
    for kind, raw_score, raw_offset in zip(SITE_TYPES, fields[2:6], fields[6:]):
        if raw_score == raw_offset == ".":
            continue
        score, offset = float(raw_score), int(raw_offset)
        if not math.isfinite(score) or not 0 <= score <= 1 or pos + offset < 1:
            raise ValueError(f"Invalid {kind} score or position")
        sites[kind] = {"delta_score": score, "delta_position": offset,
                       "position": pos + offset, "passes_threshold": score >= threshold}
    return sites


def predict_events(tx, sites):
    exons = tx.exons
    active = {kind: site["position"] * tx.direction for kind, site in sites.items()
              if site["passes_threshold"]}
    ag, al, dg, dl = (active.get(kind) for kind in SITE_TYPES)

    def event(name, signals, changed, low, high, delta):
        return name, signals, changed, (low, high), delta

    for i, (start, end) in enumerate(exons):
        prev_end = exons[i - 1][1] if i else None
        next_start = exons[i + 1][0] if i + 1 < len(exons) else None
        if prev_end is not None and next_start is not None and al == start and dl == end:
            yield event("exon_skipping", ("AL", "DL"), exons[:i] + exons[i + 1:],
                        start, end, -(end - start + 1))
        if next_start is not None:
            next_end = exons[i + 1][1]
            if dl == end and al == next_start:
                yield event("intron_retention", ("DL", "AL"),
                            exons[:i] + [(start, next_end)] + exons[i + 2:],
                            end + 1, next_start - 1, next_start - end - 1)
            if ag is not None and dg is not None and end < ag < dg < next_start:
                # A candidate pseudoexon must leave an intron on each side.
                if end + 1 < ag and dg + 1 < next_start:
                    yield event("pseudoexon_inclusion", ("AG", "DG"),
                                exons[:i + 1] + [(ag, dg)] + exons[i + 1:],
                                ag, dg, dg - ag + 1)
        for kind, loss, gain, is_acceptor in (("A", al, ag, True), ("D", dl, dg, False)):
            boundary = start if is_acceptor else end
            if gain is None or loss != boundary or gain == boundary:
                continue
            if is_acceptor and prev_end is not None and prev_end + 1 < gain <= end:
                changed = exons[:i] + [(gain, end)] + exons[i + 1:]
                low, high = sorted((gain, start))
                name = "partial_intron_retention" if gain < start else "partial_exon_skipping"
                yield event(name, (kind + "L", kind + "G"), changed, low, high - 1, start - gain)
            if not is_acceptor and next_start is not None and start <= gain < next_start - 1:
                changed = exons[:i] + [(start, gain)] + exons[i + 1:]
                low, high = sorted((gain, end))
                name = "partial_intron_retention" if gain > end else "partial_exon_skipping"
                yield event(name, (kind + "L", kind + "G"), changed, low + 1, high, gain - end)
    for kind, sites_to_compare, gain in (("AG", {x[0] for x in exons[1:]}, ag),
                                         ("DG", {x[1] for x in exons[:-1]}, dg)):
        if gain is not None and gain not in sites_to_compare:
            yield "cryptic_site_activation", (kind,), None, (gain, gain), None


def sequence_and_coordinates(reference, tx, exons):
    parts, coords = [], []
    for start, end in exons:
        low, high = sorted((abs(start), abs(end)))
        seq = reference.fetch(low, high)
        parts.append(seq if tx.direction == 1 else seq.translate(COMPLEMENT)[::-1])
        coords.extend(range(start, end + 1))
    return "".join(parts), coords


def stop_at(seq, start):
    return next((i for i in range(start, len(seq) - 2, 3) if seq[i:i + 3] in STOPS), None)


def coding_effect(reference, tx, changed, interval, delta, pos, ref, alt):
    unknown = {"frame": None, "nmd": "not_assessed"}
    if changed is None or tx.start is None or reference is None:
        return unknown
    if len(ref) != 1 or len(alt) != 1 or set(ref + alt) - set("ACGT"):
        return {**unknown, "reason": "Coding reconstruction supports SNVs only"}
    if not tx.start < interval[0] <= interval[1] < tx.stop:
        return {**unknown, "reason": "Event is outside the CDS or intersects a CDS boundary"}

    original, original_coords = sequence_and_coordinates(reference, tx, tx.exons)
    start = original_coords.index(tx.start)
    stop = original_coords.index(tx.stop)
    if (set(original[start:stop + 3]) - set("ACGT") or original[start:start + 3] != "ATG"
            or stop_at(original, start) != stop):
        return {**unknown, "reason": "Reference CDS failed start/stop validation"}
    if interval[0] <= original_coords[start + 2]:
        return {**unknown, "reason": "Event intersects the annotated start codon"}

    altered, coords = sequence_and_coordinates(reference, tx, changed)
    genomic_position = pos * tx.direction
    if genomic_position in coords:
        index = coords.index(genomic_position)
        base = alt if tx.direction == 1 else alt.translate(COMPLEMENT)
        altered = altered[:index] + base + altered[index + 1:]
    start, stop = coords.index(tx.start), coords.index(tx.stop)
    frame = "in_frame" if delta % 3 == 0 else "frameshift"
    if altered[start:start + 3] != "ATG":
        return {**unknown, "reason": "Altered start codon"}
    premature = stop_at(altered, start)
    translated = altered[start:premature + 3] if premature is not None else altered[start:]
    if set(translated) - set("ACGT"):
        return {"frame": frame, "nmd": "not_assessed", "reason": "Ambiguous translated sequence"}
    if premature is None or premature >= stop:
        return {"frame": frame, "nmd": "no_premature_stop_detected",
                "first_stop": ("absent" if premature is None else
                               "at_reference_stop" if premature == stop else "downstream")}
    junction = sum(end - begin + 1 for begin, end in changed[:-1])
    distance = junction - (premature + 3) if len(changed) > 1 else None
    return {"frame": frame, "ptc_cdna_position": premature + 1,
            "ptc_coding_position": premature - start + 1,
            "nt_from_stop_end_to_final_junction": distance,
            "nmd": "predicted" if distance is not None and distance > 50 else "escape_predicted"}


def results(vcf, tx, reference, threshold):
    for line_no, line in enumerate(vcf, 1):
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise ValueError(f"Malformed VCF at line {line_no}")
        chrom, position, _, ref, alt_field, _, _, info = fields[:8]
        if contig_key(chrom) != contig_key(tx.chrom):
            continue
        pos = int(position)
        if pos < 1:
            raise ValueError(f"Invalid POS at VCF line {line_no}")
        alts = set(alt_field.split(","))
        tags = [field.removeprefix("SpliceAI=") for field in info.split(";")
                if field.startswith("SpliceAI=")]
        for annotation in ",".join(tags).split(","):
            pair = annotation.split("|", 2)
            if len(pair) < 2 or pair[1] != tx.gene:
                continue
            alt = pair[0]
            if alt not in alts:
                raise ValueError(f"SpliceAI ALT is absent from VCF ALT at line {line_no}")
            if reference and reference.fetch(pos, pos + len(ref) - 1) != ref:
                raise ValueError(f"VCF REF does not match FASTA at {chrom}:{pos}")
            sites = parse_sites(annotation, pos, threshold)
            events = []
            for name, signals, changed, interval, delta in predict_events(tx, sites):
                events.append({"type": name, "signals": signals,
                               "genomic_interval": sorted(map(abs, interval)),
                               "reference_length_change_nt": delta,
                               **coding_effect(reference, tx, changed, interval, delta, pos, ref, alt)})
            yield {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "gene": tx.gene,
                   "transcript": tx.id, "strand": tx.strand, "threshold": threshold,
                   "sites": sites, "events": events}


def render_record(record, tx):
    exons = tx.exons
    points = [record["pos"] * tx.direction, exons[0][0], exons[-1][1]]
    points.extend(site["position"] * tx.direction for site in record["sites"].values()
                  if site["passes_threshold"])
    first, last = min(points), max(points)

    def x(pos):
        return 45 + 790 * (pos - first) / max(1, last - first)

    title = escape(f'{record["chrom"]}:{record["pos"]} {record["ref"]}>{record["alt"]}')
    parts = [f"<section><h2>{title}</h2>",
             '<svg viewBox="0 0 880 150" role="img" aria-label="Transcript and delta score tracks">',
             '<text x="12" y="28">5′</text><text x="851" y="28">3′</text>',
             '<line x1="45" y1="25" x2="835" y2="25" stroke="#526079"/>',
             '<line x1="45" y1="91" x2="835" y2="91" stroke="#526079"/>']
    for start, end in exons:
        parts.append(f'<rect x="{x(start):.2f}" y="19" width="{max(2, x(end)-x(start)):.2f}" '
                     f'height="12" fill="#285882"><title>{min(abs(start), abs(end))}–'
                     f'{max(abs(start), abs(end))}</title></rect>')
    variant_x = x(record["pos"] * tx.direction)
    parts.append(f'<line x1="{variant_x:.2f}" y1="12" x2="{variant_x:.2f}" y2="137" '
                 'stroke="#8b5bb6" stroke-dasharray="4 3"><title>Variant</title></line>')
    for kind, site in record["sites"].items():
        if not site["passes_threshold"]:
            continue
        pos_x = x(site["position"] * tx.direction)
        gain = kind.endswith("G")
        y = 91 + (-1 if gain else 1) * 42 * site["delta_score"]
        color = "#16734e" if gain else "#b12d3e"
        parts.append(f'<line x1="{pos_x:.2f}" x2="{pos_x:.2f}" y1="91" y2="{y:.2f}" '
                     f'stroke="{color}" stroke-width="3"><title>{kind} '
                     f'{site["delta_score"]} at {site["position"]}</title></line>')
    parts.append(f'</svg><pre>{escape(json.dumps(record, indent=2))}</pre></section>\n')
    return "".join(parts)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog="Outputs JSON Lines. Events are hypotheses. NMD is a junction-distance heuristic. "
               "See https://github.com/Illumina/SpliceAI and doi:10.1016/j.ajhg.2023.06.002.",
    )
    parser.add_argument("--vcf", required=True, help="SpliceAI INFO-tagged VCF or VCF.gz")
    parser.add_argument("--gtf", required=True, help="Matching-build GTF or GTF.gz with gene_name")
    parser.add_argument("--transcript", required=True, help="Exact GTF transcript_id, including version when present")
    parser.add_argument("--fasta", help="Matching-build uncompressed FASTA with adjacent .fai for SNV coding effects")
    parser.add_argument("--min-score", type=float, default=0.2, help="Screening threshold (default 0.2)")
    parser.add_argument("--html", type=Path, help="Optional local HTML report (new file)")
    args = parser.parse_args()
    if not 0 < args.min_score <= 1:
        parser.error("--min-score must be greater than zero and at most one")
    try:
        with ExitStack() as stack:
            tx = read_transcript(args.gtf, args.transcript)
            vcf = stack.enter_context(open_text(args.vcf))
            reference = None
            if args.fasta:
                fasta = stack.enter_context(open(args.fasta, "rb"))
                reference = Reference(fasta, args.fasta + ".fai", tx.chrom)
            report = None
            if args.html:
                report = stack.enter_context(args.html.open("x", encoding="utf-8"))
                report.write('<!doctype html><html lang="en"><meta charset="utf-8">'
                             '<title>SpliceAI interpretation</title>'
                             '<style>body{font:16px system-ui;max-width:960px;margin:2rem auto;padding:0 1rem}'
                             'pre{background:#f1f5f9;padding:1rem;overflow:auto}svg{width:100%}</style>'
                             f'<h1>{escape(tx.id)} · {escape(tx.gene)}</h1>'
                             '<p>Exons in transcript order. Purple marks the variant. '
                             'Green bars show gains and red bars show losses above the score threshold.</p>\n')
            count = 0
            for record in results(vcf, tx, reference, args.min_score):
                count += 1
                print(json.dumps(record))
                if report:
                    report.write(render_record(record, tx))
            if report:
                report.write(f"<footer>{count} matching annotations</footer></html>\n")
            if not count:
                print("No matching SpliceAI annotations for the selected gene and contig", file=sys.stderr)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
