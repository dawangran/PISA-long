#!/usr/bin/env python3
"""Resolve PISA gene-only alignment evidence to one assignment per library/read.

This program annotates reads and prepares a primary-only UMI input. It does not
correct UMIs, remove duplicate reads, or call biologically validated fusions.
"""
from __future__ import annotations

import argparse
import collections
import csv
import dataclasses
import gzip
import itertools
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import time

import pysam

VERSION = "1.0.0"
MATCH = {0, 7, 8}
QUERY = {0, 1, 4, 7, 8}
REFERENCE = {0, 2, 3, 7, 8}
ELIGIBLE = {"unique", "ambiguous", "fusion_candidate"}
OUTPUT_STEMS = (
    "annotated", "umi_input.primary", "fusion_candidates", "ambiguous", "bypass"
)
Interval = tuple[int, int]


def merge_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if out and start <= out[-1][1]:
            out[-1] = (out[-1][0], max(end, out[-1][1]))
        else:
            out.append((start, end))
    return out


def interval_length(intervals):
    return sum(end - start for start, end in merge_intervals(intervals))


def intersection_length(left, right):
    left, right = merge_intervals(left), merge_intervals(right)
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def exclusive_length(intervals, others):
    return interval_length(intervals) - intersection_length(intervals, others)


@dataclasses.dataclass(frozen=True)
class Gene:
    gene_id: str
    name: str
    chrom: str
    start: int
    end: int
    strand: str


def read_genes(path):
    genes = collections.defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"\s*;?', fields[8]))
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            # Delimiters are part of our BAM evidence/group serialization.
            if any(c in gene_id for c in ";|\t\n"):
                raise ValueError(f"Unsupported delimiter in gene_id: {gene_id!r}")
            gene = Gene(gene_id, attrs.get("gene_name", gene_id), fields[0],
                        int(fields[3]) - 1, int(fields[4]), fields[6])
            if gene.end <= gene.start:
                raise ValueError(f"Invalid GTF gene interval: {gene_id}")
            if gene not in genes[gene_id]:
                genes[gene_id].append(gene)
    if not genes:
        raise ValueError("GTF has no gene features with gene_id; use the same full GTF as PISA")
    return dict(genes)


def get_tag(record, tag, default=None):
    return record.get_tag(tag) if record.has_tag(tag) else default


def query_blocks(record):
    """M/=/X blocks in original read coordinates, counting hard clips in length.

    Each tuple is (reference_start, reference_end, oriented_query_start,
    oriented_query_end). The last two coordinates run along the SAM alignment;
    normalize_query() converts these to the original sequencing orientation.
    """
    cigar = record.cigartuples or []
    full_length = sum(n for op, n in cigar if op in QUERY or op == 5)
    ref = record.reference_start
    query = 0
    blocks = []
    for op, n in cigar:
        if op in MATCH:
            blocks.append((ref, ref + n, query, query + n))
        if op in REFERENCE:
            ref += n
        if op in QUERY or op == 5:
            query += n
    return full_length, blocks


def normalize_query(start, end, length, reverse):
    return (length - end, length - start) if reverse else (start, end)


@dataclasses.dataclass
class Evidence:
    record: object
    length: int
    matched: list[Interval]
    genes: dict[str, list[Interval]]


def alignment_evidence(record, genes):
    if record.is_unmapped:
        return Evidence(record, 0, [], {})
    length, blocks = query_blocks(record)
    total = sum(end - start for start, end, _, _ in blocks)
    if not record.has_tag("ZQ"):
        raise ValueError(f"{record.query_name}: mapped record lacks ZQ; first run PISA anno -gene-only")
    if int(record.get_tag("ZQ")) != total:
        raise ValueError(f"{record.query_name}: ZQ does not match CIGAR M/=/X length")
    ids = str(get_tag(record, "ZG", "")).split(";") if get_tag(record, "ZG", "") else []
    encoded_bp = str(get_tag(record, "ZB", "")).split(";") if get_tag(record, "ZB", "") else []
    if len(ids) != len(encoded_bp) or len(ids) != len(set(ids)):
        raise ValueError(f"{record.query_name}: malformed ZG/ZB evidence")
    matched = [normalize_query(q0, q1, length, record.is_reverse)
               for _, _, q0, q1 in blocks]
    supports = {}
    for gene_id, bp_text in zip(ids, encoded_bp):
        if gene_id not in genes:
            raise ValueError(f"{record.query_name}: gene {gene_id} absent from supplied GTF")
        intervals = []
        for gene in genes[gene_id]:
            if gene.chrom != record.reference_name:
                continue
            for r0, r1, q0, _ in blocks:
                lo, hi = max(r0, gene.start), min(r1, gene.end)
                if lo < hi:
                    intervals.append(normalize_query(q0 + lo - r0, q0 + hi - r0,
                                                     length, record.is_reverse))
        intervals = merge_intervals(intervals)
        if interval_length(intervals) != int(bp_text):
            raise ValueError(f"{record.query_name}: {gene_id} ZB disagrees with CIGAR/GTF; use identical GTF")
        if intervals:
            supports[gene_id] = intervals
    return Evidence(record, length, merge_intervals(matched), supports)


def normalized_cigar(cigar):
    """SA may describe soft clips where the stored record uses hard clips."""
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar or "")
    normalized = []
    for n, op in tokens:
        op = "S" if op == "H" else op
        if normalized and normalized[-1][1] == op:
            normalized[-1] = (normalized[-1][0] + int(n), op)
        else:
            normalized.append((int(n), op))
    return tuple(normalized)


def missing_sa_records(main):
    observed = {(r.reference_name, r.reference_start,
                 "-" if r.is_reverse else "+", normalized_cigar(r.cigarstring))
                for r in main}
    missing = set()
    for r in main:
        sa = get_tag(r, "SA", "")
        for entry in sa.rstrip(";").split(";") if sa else []:
            fields = entry.split(",")
            if len(fields) != 6:
                raise ValueError(f"{r.query_name}: malformed SA entry {entry!r}")
            key = (fields[0], int(fields[1]) - 1, fields[2], normalized_cigar(fields[3]))
            if key not in observed:
                missing.add(key)
    return missing


def dominant_gene(support, total, args):
    """Select a dominant body while tolerating small shared-boundary hits.

    The shared exception requires support exclusive of *all* other candidates,
    not separate pairwise advantages which could cover the whole read together.
    """
    bp = {g: interval_length(v) for g, v in support.items()}
    if not bp:
        return None
    top = max(bp, key=lambda g: (bp[g], g))
    if bp[top] < args.min_gene_bp or bp[top] < args.min_gene_fraction * total:
        return None
    global_exclusive = exclusive_length(
        support[top], [i for g, spans in support.items() if g != top for i in spans])
    for gene_id, amount in bp.items():
        if gene_id == top:
            continue
        if amount > args.minor_bp:
            return None
        if amount <= args.minor_fraction * total:
            continue
        if (global_exclusive < args.min_gene_bp
                or exclusive_length(support[gene_id], support[top]) != 0):
            return None
    return top


def candidate_genes(evidence, args):
    total = interval_length(evidence.matched)
    bp = {g: interval_length(v) for g, v in evidence.genes.items()}
    top = dominant_gene(evidence.genes, total, args)
    if top is not None:
        return {top}
    substantial = {g for g, v in bp.items() if v >= args.min_gene_bp}
    # Weak alternatives do not fabricate a confident second gene.
    return substantial


def competitive_secondary(secondary, main, args):
    """Return competing alternative genes without summing alternative mappings."""
    competing = set()
    conflict = False
    for alt in secondary:
        alt_len = interval_length(alt.matched)
        alt_genes = candidate_genes(alt, args)
        if not alt_len or not alt_genes:
            continue
        for ref in main:
            ref_len = interval_length(ref.matched)
            if not ref_len:
                continue
            overlap = intersection_length(alt.matched, ref.matched)
            if (overlap < args.secondary_query_overlap * alt_len
                    or overlap < args.secondary_query_overlap * ref_len):
                continue
            if alt.record.has_tag("AS") and ref.record.has_tag("AS"):
                alt_score = float(alt.record.get_tag("AS"))
                ref_score = float(ref.record.get_tag("AS"))
                if ref_score > 0 and alt_score < args.secondary_score_ratio * ref_score:
                    continue
                if ref_score <= 0 and alt_score < ref_score:
                    continue
            ref_genes = candidate_genes(ref, args)
            if alt_genes != ref_genes:
                conflict = True
                competing.update(alt_genes)
                competing.update(ref_genes)
    return competing, conflict


@dataclasses.dataclass
class Assignment:
    status: str
    reason: str
    gx: str
    candidates: list[str]
    body_bp: dict[str, int]
    independent_bp: dict[str, int]
    total_matched: int
    full_length: int
    cb: str
    ur: str
    eligible: bool
    primary: object | None
    missing_sa_count: int = 0


def resolve_read(records, genes, args):
    primaries = [r for r in records if not (r.flag & 0x900)]
    primary = primaries[0] if len(primaries) == 1 else None
    main_records = [r for r in records if not r.is_secondary and not r.is_unmapped]
    main = [alignment_evidence(r, genes) for r in main_records]
    secondary = [alignment_evidence(r, genes) for r in records
                 if r.is_secondary and not r.is_unmapped]
    support = collections.defaultdict(list)
    reliable = collections.defaultdict(list)
    for ev in main:
        for gene_id, intervals in ev.genes.items():
            support[gene_id].extend(intervals)
            if (args.min_mapq <= ev.record.mapping_quality < 255
                    and not ev.record.is_qcfail):
                reliable[gene_id].extend(intervals)
    support = {g: merge_intervals(v) for g, v in support.items()}
    reliable = {g: merge_intervals(v) for g, v in reliable.items()}
    body = {g: interval_length(v) for g, v in support.items()}
    independent = {g: exclusive_length(reliable.get(g, []),
                                      [i for other, spans in support.items() if other != g for i in spans])
                   for g in support}
    total = interval_length([i for ev in main for i in ev.matched])
    lengths = {ev.length for ev in main if ev.length}
    full_length = max(lengths, default=0)
    candidates = sorted(body)
    cb_values = {str(get_tag(r, "CB")) for r in records if get_tag(r, "CB", "")}
    ur_values = {str(get_tag(r, "UR")) for r in records if get_tag(r, "UR", "")}
    cb = next(iter(cb_values)) if len(cb_values) == 1 else ""
    ur = next(iter(ur_values)) if len(ur_values) == 1 else ""
    missing = missing_sa_records(main_records)

    def result(status, reason, gx="", selected=None):
        chosen = candidates if selected is None else sorted(selected)
        eligible = status in ELIGIBLE
        reasons = [reason] if reason else []
        if primary is None or primary.is_unmapped:
            eligible = False
        if primary is not None and primary.is_qcfail:
            eligible = False
            reasons.append("primary_qcfail")
        if primary is not None and (not get_tag(primary, "CB", "") or not get_tag(primary, "UR", "")):
            eligible = False
            reasons.append("primary_missing_CB_or_UR")
        return Assignment(status, ";".join(reasons), gx, chosen, body, independent,
                          total, full_length, cb, ur, eligible, primary, len(missing))

    if len(primaries) != 1:
        return result("missing_primary" if not primaries else "multiple_primary",
                      "incomplete_or_nonunique_primary")
    if len(cb_values) > 1 or len(ur_values) > 1:
        return result("tag_conflict", "inconsistent_CB_or_UR")
    if primary.is_unmapped:
        return result("unmapped", "primary_unmapped")
    if len(lengths) > 1:
        return result("read_length_conflict", "inconsistent_full_query_lengths")
    if missing:
        return result("incomplete_alignment_set", "SA_partner_not_present")
    if primary.mapping_quality == 255:
        return result("unknown_mapq", "primary_MAPQ_unavailable")
    if primary.mapping_quality < args.min_mapq:
        return result("low_mapq", "primary_below_min_mapq")
    if not body or not total:
        return result("unassigned", "no_actual_gene_overlap")

    alternative, alternative_conflict = competitive_secondary(secondary, main, args)
    top = dominant_gene(support, total, args)
    significant = {g for g, bp in body.items() if bp >= args.min_gene_bp}
    covered = interval_length([i for spans in support.values() for i in spans])
    anchors = {g for g, bp in independent.items() if bp >= args.min_anchor_bp}
    if (len(anchors) >= 2 and covered >= args.min_gene_fraction * total
            and not alternative_conflict):
        # A third substantial unresolved locus prevents an over-specific partner list.
        unexplained = significant - anchors
        if not unexplained:
            fusion_ids = sorted(anchors)
            return result("fusion_candidate", "distinct_reliable_query_anchors",
                          "FUSION:" + "|".join(fusion_ids), fusion_ids)
    if alternative_conflict:
        possible = significant | alternative
        if len(possible) >= 2:
            return result("ambiguous", "competitive_secondary",
                          "AMBIG:" + "|".join(sorted(possible)), possible)
        return result("insufficient", "competitive_secondary_without_resolved_gene")
    if top is not None:
        return result("unique", "dominant_gene_body_support", top, [top])
    if len(significant) >= 2 and covered >= args.min_gene_fraction * total:
        return result("ambiguous", "overlapping_or_conflicting_gene_support",
                      "AMBIG:" + "|".join(sorted(significant)), significant)
    return result("insufficient", "gene_coverage_below_threshold")


def set_read_tags(records, assignment, genes):
    names = []
    for gene_id in assignment.candidates:
        distinct = sorted({g.name for g in genes[gene_id]})
        names.append(distinct[0] if distinct else gene_id)
    name_text = "|".join(names)
    for record in records:
        for tag in ("GX", "GN", "ZA", "ZF", "ZR", "ZC"):
            if record.has_tag(tag):
                record.set_tag(tag, None)
        if assignment.gx:
            record.set_tag("GX", assignment.gx, value_type="Z")
            record.set_tag("GN", name_text, value_type="Z")
        record.set_tag("ZA", assignment.status, value_type="Z")
        record.set_tag("ZR", assignment.reason, value_type="Z")
        if assignment.candidates:
            record.set_tag("ZC", ";".join(assignment.candidates), value_type="Z")
        if assignment.status == "fusion_candidate":
            record.set_tag("ZF", "fusion_candidate", value_type="Z")


def library_key(record, rg_libraries):
    rg = str(get_tag(record, "RG", ""))
    return rg_libraries.get(rg, f"RG:{rg}" if rg else "__input_library__")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--temp-dir", type=Path)
    parser.add_argument("--min-gene-bp", type=int, default=100)
    parser.add_argument("--min-gene-fraction", type=float, default=0.8)
    parser.add_argument("--minor-bp", type=int, default=100)
    parser.add_argument("--minor-fraction", type=float, default=0.2)
    parser.add_argument("--min-anchor-bp", type=int, default=100)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--secondary-score-ratio", type=float, default=0.95)
    parser.add_argument("--secondary-query-overlap", type=float, default=0.8)
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    for name in ("min_gene_fraction", "minor_fraction", "secondary_score_ratio", "secondary_query_overlap"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0,1]")
    if args.threads < 1 or args.min_gene_bp < 1 or args.min_anchor_bp < 1 or args.minor_bp < 0:
        parser.error("threads/min-gene-bp/min-anchor-bp must be positive; minor-bp must be nonnegative")
    if not 0 <= args.min_mapq <= 254:
        parser.error("min-mapq must be in [0,254]")
    for path in (args.input, args.gtf):
        if not path.is_file():
            parser.error(f"File not found: {path}")
    if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
        parser.error("output-dir must be absent or empty; existing results are never overwritten")
    return args


def run(args):
    started = time.monotonic()
    genes = read_genes(args.gtf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.temp_dir:
        args.temp_dir.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(args.input), "rb", check_sq=False) as inp:
        header = inp.header.to_dict()
    rg_libraries = {str(rg["ID"]): ("LB:" + str(rg["LB"]) if rg.get("LB") else "RG:" + str(rg["ID"]))
                    for rg in header.get("RG", [])}
    pg_ids = {pg["ID"] for pg in header.get("PG", [])}
    pg_id = "pisa_resolve_gene_reads"
    while pg_id in pg_ids:
        pg_id += ".1"
    header.setdefault("PG", []).append({"ID": pg_id, "PN": "pisa_resolve_gene_reads",
                                       "VN": VERSION, "CL": shlex.join(sys.argv)})
    header.setdefault("HD", {})["SO"] = "unsorted"
    rows_path = args.output_dir / "read_assignments.tsv.gz"
    columns = ["library", "read_id", "status", "reason", "GX", "candidate_gene_ids",
               "gene_body_bp", "independent_bp", "matched_query_bp", "full_read_length",
               "matched_query_fraction", "primary_flag", "primary_chrom", "primary_start_1based",
               "primary_MAPQ", "n_records", "n_primary", "n_secondary", "n_supplementary",
               "missing_SA_records", "CB", "UR", "umi_eligible"]
    status_counts = collections.Counter()
    eligible_counts = collections.Counter()
    alignment_status_counts = collections.Counter()
    file_record_counts = collections.Counter()
    counters = collections.Counter()
    library_counts = collections.Counter()
    library_eligible_counts = collections.Counter()
    with tempfile.TemporaryDirectory(prefix="pisa-read-resolve-", dir=args.temp_dir or args.output_dir) as tmp_name:
        tmp = Path(tmp_name)
        writers = {stem: pysam.AlignmentFile(str(tmp / f"{stem}.bam"), "wb", header=header)
                   for stem in OUTPUT_STEMS}
        command = [args.samtools, "collate", "-@", str(args.threads), "-u", "-O",
                   str(args.input.resolve()), str(tmp / "collate")]
        with open(args.output_dir / "collate.log", "wb") as log, gzip.open(rows_path, "wt", newline="") as table:
            writer = csv.DictWriter(table, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
            try:
                with pysam.AlignmentFile(proc.stdout, "rb", check_sq=False) as collated:
                    for read_id, group in itertools.groupby(collated, key=lambda r: r.query_name):
                        library_groups = collections.defaultdict(list)
                        for record in group:
                            library_groups[library_key(record, rg_libraries)].append(record)
                        for library, records in sorted(library_groups.items()):
                            assignment = resolve_read(records, genes, args)
                            set_read_tags(records, assignment, genes)
                            status_counts[assignment.status] += 1
                            alignment_status_counts[assignment.status] += len(records)
                            counters["read_groups"] += 1
                            library_counts[library] += 1
                            counters["alignment_records"] += len(records)
                            counters["primary_records"] += sum(not (r.flag & 0x900) for r in records)
                            counters["secondary_records"] += sum(r.is_secondary for r in records)
                            counters["supplementary_records"] += sum(r.is_supplementary for r in records)
                            for record in records:
                                writers["annotated"].write(record)
                                file_record_counts["annotated"] += 1
                                if assignment.status == "fusion_candidate":
                                    writers["fusion_candidates"].write(record)
                                    file_record_counts["fusion_candidates"] += 1
                                if assignment.status == "ambiguous":
                                    writers["ambiguous"].write(record)
                                    file_record_counts["ambiguous"] += 1
                                if not assignment.eligible:
                                    writers["bypass"].write(record)
                                    file_record_counts["bypass"] += 1
                            if assignment.eligible:
                                writers["umi_input.primary"].write(assignment.primary)
                                file_record_counts["umi_input.primary"] += 1
                                eligible_counts[assignment.status] += 1
                                library_eligible_counts[library] += 1
                            primary = assignment.primary
                            writer.writerow({
                                "library": library, "read_id": read_id,
                                "status": assignment.status, "reason": assignment.reason,
                                "GX": assignment.gx, "candidate_gene_ids": ";".join(assignment.candidates),
                                "gene_body_bp": json.dumps(assignment.body_bp, sort_keys=True, separators=(",", ":")),
                                "independent_bp": json.dumps(assignment.independent_bp, sort_keys=True, separators=(",", ":")),
                                "matched_query_bp": assignment.total_matched, "full_read_length": assignment.full_length,
                                "matched_query_fraction": f"{assignment.total_matched / assignment.full_length:.6f}" if assignment.full_length else "0",
                                "primary_flag": primary.flag if primary is not None else "",
                                "primary_chrom": primary.reference_name if primary is not None and not primary.is_unmapped else "",
                                "primary_start_1based": primary.reference_start + 1 if primary is not None and not primary.is_unmapped else "",
                                "primary_MAPQ": primary.mapping_quality if primary is not None else "",
                                "n_records": len(records), "n_primary": sum(not (r.flag & 0x900) for r in records),
                                "n_secondary": sum(r.is_secondary for r in records),
                                "n_supplementary": sum(r.is_supplementary for r in records),
                                "missing_SA_records": assignment.missing_sa_count,
                                "CB": assignment.cb, "UR": assignment.ur, "umi_eligible": int(assignment.eligible),
                            })
                if proc.wait() != 0:
                    raise RuntimeError("samtools collate failed; see collate.log")
            finally:
                if proc.poll() is None:
                    proc.terminate()
                    proc.wait()
                if proc.stdout:
                    proc.stdout.close()
                for out in writers.values():
                    out.close()
        for stem in OUTPUT_STEMS:
            destination = args.output_dir / f"{stem}.sorted.bam"
            pysam.sort("-@", str(args.threads), "-T", str(tmp / f"sort-{stem}"),
                       "-o", str(destination), str(tmp / f"{stem}.bam"))
            pysam.index("-@", str(args.threads), str(destination))
            pysam.quickcheck(str(destination))
    summary = {
        "version": VERSION, "annotation_only": True, "umi_deduplication_performed": False,
        "input": str(args.input.resolve()), "gtf": str(args.gtf.resolve()),
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "counts": dict(counters), "read_status_counts": dict(status_counts),
        "alignment_status_counts": dict(alignment_status_counts),
        "umi_eligible_read_status_counts": dict(eligible_counts),
        "library_read_counts": dict(library_counts),
        "library_umi_eligible_counts": dict(library_eligible_counts),
        "umi_input_requires_library_split": len(library_eligible_counts) > 1,
        "output_alignment_records": {stem: file_record_counts[stem] for stem in OUTPUT_STEMS},
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "notes": [
            "read_groups means library + QNAME, not independent molecules",
            "all input alignment records retained in annotated.sorted.bam",
            "fusion labels are unordered gene-combination candidates, not validated events",
            "ZG/ZN/ZB/ZQ remain alignment evidence; ZC/ZA/GX describe the read",
            "low/unknown MAPQ, incomplete SA sets, missing primary and tag conflicts bypass UMI input",
            "SA completeness can only detect missing records advertised by SA",
            "UMI-tools ordinary per-gene BAM mode does not merge across primary chromosomes",
            "if umi_input_requires_library_split is true, split UMI input by library before running UMI-tools",
        ],
    }
    with open(args.output_dir / "summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    with open(args.output_dir / "libraries.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["library", "read_groups", "umi_eligible_reads"])
        for library in sorted(library_counts):
            writer.writerow([library, library_counts[library], library_eligible_counts[library]])
    if len(library_eligible_counts) > 1:
        print("WARNING: UMI input contains multiple libraries; split by library before UMI-tools (see libraries.tsv)",
              file=sys.stderr)
    print(json.dumps({"output_dir": str(args.output_dir), "read_status_counts": dict(status_counts),
                      "umi_input_reads": sum(eligible_counts.values())}, ensure_ascii=False))


def main(argv=None):
    args = parse_args(argv)
    try:
        run(args)
    except (ValueError, RuntimeError, OSError, pysam.SamtoolsError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
