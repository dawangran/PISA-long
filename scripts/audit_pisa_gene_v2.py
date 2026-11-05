#!/usr/bin/env python3
"""Independently audit experimental PISA v2 on a manageable BAM slice.

This is a verification utility, not a production annotator for a large BAM. It
expands only aligned reference positions into Python sets. It never uses PISA
output candidate tags to choose which genes to inspect.
"""
from __future__ import annotations

import argparse
import collections
import csv
import dataclasses
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time

import pysam


MATCH = {0, 7, 8}
ANNOTATION_TAGS = {
    "GX", "GN", "RE", "TX", "AT", "EX", "JC", "FL", "ER",
    "ZA", "ZG", "ZN", "ZB", "ZQ", "ZF", "ZC", "ZR", "ZE", "ZU", "ZX", "ZJ",
}


@dataclasses.dataclass
class Gene:
    gene_id: str
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: list[tuple[int, int]] = dataclasses.field(default_factory=list)
    junctions: set[tuple[int, int]] = dataclasses.field(default_factory=set)


def attrs(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"\s*;?', text))


def merge(intervals):
    merged = []
    for lo, hi in sorted(set(intervals)):
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    return merged


def gtf_lines(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) == 9:
                    yield fields


def load_reference(path, windows):
    """Read every GTF gene feature, then complete exons of selected genes.

    Every input alignment is contained in its chromosome window. Therefore a
    gene with actual aligned-base overlap must intersect that window. The
    second complete scan retains all exons of these genes, including exons
    outside the BAM window, so transcript-derived junctions remain complete.
    """
    genes = {}
    stats = collections.Counter()
    for fields in gtf_lines(path):
        stats["gtf_records_scanned_first_pass"] += 1
        if fields[2] != "gene":
            continue
        stats["gene_features_scanned_all_chromosomes"] += 1
        window = windows.get(fields[0])
        if window is None:
            continue
        lo, hi = int(fields[3]) - 1, int(fields[4])
        if hi <= window[0] or lo >= window[1]:
            continue
        a = attrs(fields[8])
        gene_id = a.get("gene_id")
        if not gene_id:
            raise ValueError("Overlapping GTF gene feature has no gene_id")
        gene = Gene(gene_id, a.get("gene_name", gene_id), fields[0], lo, hi, fields[6])
        if gene_id in genes:
            previous = genes[gene_id]
            if dataclasses.astuple(previous)[:6] != dataclasses.astuple(gene)[:6]:
                raise ValueError(f"Multiple incompatible GTF gene features: {gene_id}")
        genes[gene_id] = gene
    transcripts = collections.defaultdict(list)
    for fields in gtf_lines(path):
        stats["gtf_records_scanned_second_pass"] += 1
        if fields[2] != "exon" or fields[0] not in windows:
            continue
        a = attrs(fields[8])
        gene = genes.get(a.get("gene_id"))
        if gene is None:
            continue
        if fields[0] != gene.chrom:
            raise ValueError(f"Exon and gene chromosome mismatch: {gene.gene_id}")
        interval = (int(fields[3]) - 1, int(fields[4]))
        gene.exons.append(interval)
        stats["exon_records_retained"] += 1
        transcript_id = a.get("transcript_id")
        if transcript_id:
            transcripts[(gene.gene_id, transcript_id)].append(interval)
    for (gene_id, _), exons in transcripts.items():
        exons = merge(exons)
        for left, right in zip(exons, exons[1:]):
            if left[1] < right[0]:
                genes[gene_id].junctions.add((left[1], right[0]))
    for gene in genes.values():
        gene.exons = merge(gene.exons)
    stats["complete_genes_retained"] = len(genes)
    stats["transcripts_retained"] = len(transcripts)
    stats["merged_exon_intervals"] = sum(len(g.exons) for g in genes.values())
    stats["distinct_gene_junctions"] = sum(len(g.junctions) for g in genes.values())
    return genes, dict(stats)


def alignment_key(read):
    values = (
        read.query_name, read.flag, read.reference_name, read.reference_start,
        read.mapping_quality, read.cigarstring, read.next_reference_name,
        read.next_reference_start, read.template_length, read.query_sequence,
        read.qual,
    )
    return hashlib.sha256(repr(values).encode()).hexdigest()


def nonannotation_tags(read):
    return sorted((k, repr(v), t) for k, v, t in read.get_tags(with_value_type=True)
                  if k not in ANNOTATION_TAGS)


def tag(read, name, default=""):
    return read.get_tag(name) if read.has_tag(name) else default


def matched_positions(read):
    """Return M reference positions and N junctions with flanking M segments.

    I/D may neighbor N. M support must be found in each immediately adjacent
    N-delimited segment, and cannot be borrowed across another N.
    """
    ref = read.reference_start
    positions = set()
    segments = [False]
    jumps = []
    for op, length in read.cigartuples or []:
        if op in MATCH:
            positions.update(range(ref, ref + length))
            segments[-1] = True
            ref += length
        elif op == 2:
            ref += length
        elif op == 3:
            jumps.append((ref, ref + length))
            segments.append(False)
            ref += length
    supported_jumps = {jump for i, jump in enumerate(jumps)
                       if segments[i] and segments[i + 1]}
    return positions, supported_jumps


def strand_matches(read, gene, args):
    if not args.respect_strand:
        return True
    reverse = read.is_reverse ^ args.reverse_trans ^ read.is_read2
    return gene.strand == "." or gene.strand == ("-" if reverse else "+")


def union(sets):
    return set().union(*sets)


def expected(read, by_chrom, args):
    total = sum(length for op, length in read.cigartuples or [] if op in MATCH)
    quality = ("qc_fail" if read.is_qcfail else
               "unknown_mapq" if read.mapping_quality == 255 else
               "low_mapq" if read.mapping_quality < args.minimum_mapq else "")
    if read.is_unmapped:
        return {"ids": [], "kept": [], "GX": "", "GN": "", "ZA": "unmapped",
                "total": total, "body": {}, "exon": {}, "body_exclusive": {},
                "exon_exclusive": {}, "junctions": {}, "anchors": [],
                "shared_removed": [], "weak_removed": [], "coverage": 0.0}
    points, read_junctions = matched_positions(read)
    if len(points) != total:
        raise ValueError(f"Non-monotonic/reused reference positions: {read.query_name}")
    support = {}
    exonic = {}
    junctions = {}
    names = {}
    for gene in by_chrom.get(read.reference_name, []):
        if gene.end <= read.reference_start or gene.start >= read.reference_end:
            continue
        if not strand_matches(read, gene, args):
            continue
        body = {p for p in points if gene.start <= p < gene.end}
        if not body:
            continue
        support[gene.gene_id] = body
        names[gene.gene_id] = gene.name
        # Expand only intersections with actual aligned positions, not gene bodies.
        exonic[gene.gene_id] = {p for p in points
                               if any(lo <= p < hi for lo, hi in gene.exons)}
        junctions[gene.gene_id] = read_junctions & gene.junctions
    ids = sorted(support)
    body_unique = {g: support[g] - union([support[h] for h in ids if h != g]) for g in ids}
    exon_unique = {g: exonic[g] - union([exonic[h] for h in ids if h != g]) for g in ids}
    anchors = [g for g in ids
               if len(body_unique[g]) / total >= args.exclusive_frac
               and (len(support[g]) / total >= args.main_frac
                    or len(exon_unique[g]) / total >= args.exon_frac)] if total else []
    kept = set(ids)
    shared_removed = set()
    weak_removed = set()
    if len(ids) > 1 and not quality and len(anchors) == 1:
        main = anchors[0]
        shared_removed = {g for g in ids if g != main
                          and support[g] <= support[main]
                          and exonic[g] <= exonic[main]
                          and junctions[g] <= junctions[main]}
        weak = {g for g in ids if g != main and g not in shared_removed
                and not exonic[g] and not junctions[g]
                and len(support[g]) / total <= args.weak_frac}
        if len(union([support[g] for g in weak])) / total <= args.weak_frac:
            weak_removed = weak
        kept -= shared_removed | weak_removed
    kept = sorted(kept)
    gx = (kept[0] if len(kept) == 1 else
          "AMBIG:" + "|".join(kept) if kept else "")
    if len(kept) == 1:
        status = "unique" if len(support[kept[0]]) / total >= args.main_frac else "assigned_partial"
    else:
        status = "ambiguous" if kept else "unassigned"
    if quality:
        status = quality
    return {
        "ids": ids, "kept": kept, "GX": gx,
        "GN": "|".join(names[g] for g in kept), "ZA": status, "total": total,
        "body": {g: len(support[g]) for g in ids},
        "exon": {g: len(exonic[g]) for g in ids},
        "body_exclusive": {g: len(body_unique[g]) for g in ids},
        "exon_exclusive": {g: len(exon_unique[g]) for g in ids},
        "junctions": {g: len(junctions[g]) for g in ids},
        "anchors": anchors, "shared_removed": sorted(shared_removed),
        "weak_removed": sorted(weak_removed),
        "coverage": len(union(list(support.values()))) / total if total else 0.0,
    }


def record_type(read):
    return "secondary" if read.is_secondary else "supplementary" if read.is_supplementary else "primary"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--annotated", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--previous", help="Optional v1 alignment BAM for target-read comparison")
    p.add_argument("--target-read-file", type=Path,
                   help="Optional text file with one read ID per line for target_reads.tsv")
    p.add_argument("--respect-strand", action="store_true", help="Default matches PISA -ignore-strand")
    p.add_argument("--reverse-trans", action="store_true")
    p.add_argument("--minimum-mapq", type=int, default=0)
    p.add_argument("--main-frac", type=float, default=.8)
    p.add_argument("--exclusive-frac", type=float, default=.2)
    p.add_argument("--exon-frac", type=float, default=.5)
    p.add_argument("--weak-frac", type=float, default=.1)
    args = p.parse_args()
    target_reads = set()
    if args.target_read_file:
        target_reads = {line.strip() for line in args.target_read_file.read_text().splitlines()
                        if line.strip()}
        # Keep the JSON parameters serializable, as for all other path options.
        args.target_read_file = str(args.target_read_file)
    for value in (args.main_frac, args.exclusive_frac, args.exon_frac, args.weak_frac):
        if not math.isfinite(value) or not 0 <= value <= 1:
            p.error("fraction parameters must be finite and within [0,1]")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with pysam.AlignmentFile(args.input, "rb") as bam:
        originals = list(bam)
    windows = {}
    for read in originals:
        if read.is_unmapped:
            continue
        start, end = read.reference_start, read.reference_end
        old = windows.get(read.reference_name, (start, end))
        windows[read.reference_name] = (min(start, old[0]), max(end, old[1]))
    genes, reference_stats = load_reference(args.gtf, windows)
    by_chrom = collections.defaultdict(list)
    for gene in genes.values():
        by_chrom[gene.chrom].append(gene)
    annotated = collections.defaultdict(collections.deque)
    with pysam.AlignmentFile(args.annotated, "rb") as bam:
        for read in bam:
            annotated[alignment_key(read)].append(read)
    previous = collections.defaultdict(collections.deque)
    if args.previous:
        with pysam.AlignmentFile(args.previous, "rb") as bam:
            for read in bam:
                previous[alignment_key(read)].append(read)
    observed_count = sum(len(v) for v in annotated.values())
    mismatches = []
    counts = collections.defaultdict(collections.Counter)
    rows = []
    target_rows = []
    checked_tags = collections.Counter()
    for read in originals:
        key = alignment_key(read)
        if not annotated[key]:
            mismatches.append({"read_id": read.query_name, "field": "alignment_core", "expected": "present", "actual": "missing"})
            continue
        actual = annotated[key].popleft()
        if nonannotation_tags(read) != nonannotation_tags(actual):
            mismatches.append({"read_id": read.query_name, "field": "non_annotation_tags", "expected": "unchanged", "actual": "changed"})
        result = expected(read, by_chrom, args)
        ids = result["ids"]
        expected_tags = {"GX": result["GX"], "GN": result["GN"], "ZA": result["ZA"],
                         "ZG": ";".join(ids), "ZN": ";".join(genes[g].name for g in ids),
                         "ZQ": result["total"], "ZC": ";".join(result["kept"])}
        for name, field in [("ZB", "body"), ("ZE", "exon"), ("ZU", "body_exclusive"),
                            ("ZX", "exon_exclusive"), ("ZJ", "junctions")]:
            expected_tags[name] = ";".join(str(result[field][g]) for g in ids)
        for name, value in expected_tags.items():
            observed = tag(actual, name)
            checked_tags[name] += 1
            if observed != value:
                mismatches.append({"read_id": read.query_name, "flag": read.flag,
                                   "start_1based": read.reference_start + 1,
                                   "field": name, "expected": value, "actual": observed})
        if (result["shared_removed"] or result["weak_removed"] or result["ZA"] == "assigned_partial") and not tag(actual, "ZR"):
            mismatches.append({"read_id": read.query_name, "field": "ZR", "expected": "nonempty explanation", "actual": ""})
        kind = record_type(read)
        for scope in ("all", kind):
            counts[scope]["records"] += 1
            counts[scope]["expected_status_" + result["ZA"]] += 1
            counts[scope]["with_GX" if result["GX"] else "without_GX"] += 1
            counts[scope]["shared_candidates_removed"] += len(result["shared_removed"])
            counts[scope]["weak_intronic_candidates_removed"] += len(result["weak_removed"])
            if result["shared_removed"]:
                counts[scope]["records_with_shared_reduction"] += 1
            if result["weak_removed"]:
                counts[scope]["records_with_weak_intronic_reduction"] += 1
        row = {
            "read_id": read.query_name, "flag": read.flag, "alignment_type": kind,
            "chrom": read.reference_name or "*", "start_1based": read.reference_start + 1,
            "end_1based": read.reference_end or 0, "MAPQ": read.mapping_quality,
            "matched_bp": result["total"], "all_gene_coverage_fraction": f'{result["coverage"]:.9f}',
            "expected_GX": result["GX"], "actual_GX": tag(actual, "GX"),
            "actual_GN": tag(actual, "GN"), "expected_ZA": result["ZA"],
            "actual_ZA": tag(actual, "ZA"), "actual_ZR": tag(actual, "ZR"),
            "main_candidates": ";".join(result["anchors"]),
            "shared_removed": ";".join(result["shared_removed"]),
            "weak_intronic_removed": ";".join(result["weak_removed"]),
        }
        for field in ("body", "exon", "body_exclusive", "exon_exclusive", "junctions"):
            row[field] = json.dumps(result[field], sort_keys=True, separators=(",", ":"))
        if args.previous:
            prior = previous[key].popleft() if previous[key] else None
            row["previous_GX"] = tag(prior, "GX") if prior is not None else "<record_missing>"
            row["previous_GN"] = tag(prior, "GN") if prior is not None else ""
            row["previous_ZA"] = tag(prior, "ZA") if prior is not None else ""
        rows.append(row)
        if read.query_name in target_reads:
            target_rows.append(row)
    extras = sum(len(v) for v in annotated.values())
    if extras:
        mismatches.append({"field": "extra_alignment_records", "expected": 0, "actual": extras})
    fields = list(rows[0]) if rows else ["read_id"]
    for name, table in [("per_alignment_independent_audit.tsv.gz", rows), ("target_reads.tsv", target_rows)]:
        opener = gzip.open if name.endswith(".gz") else open
        with opener(out / name, "wt") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(table)
    summary = {
        "audit": "Independent actual-reference-position sets and two complete GTF scans",
        "scope": "Small BAM slice; agreement validates implementation, not biological assignment accuracy.",
        "input": str(Path(args.input).resolve()), "annotated": str(Path(args.annotated).resolve()),
        "gtf": str(Path(args.gtf).resolve()), "parameters": vars(args),
        "input_records": len(originals), "annotated_records": observed_count,
        "mapped_input_windows_0based_halfopen": windows,
        "reference_loading": reference_stats, "checked_tags": dict(checked_tags),
        "counts": {k: dict(v) for k, v in counts.items()},
        "mismatch_count": len(mismatches), "all_checks_passed": not mismatches,
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    for name, value in [("independent_validation.json", summary), ("mismatches.json", mismatches)]:
        (out / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"all_checks_passed": not mismatches, "records": len(originals),
                      "mismatch_count": len(mismatches), "elapsed_seconds": summary["elapsed_seconds"],
                      "output_dir": str(out)}, ensure_ascii=False))
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
