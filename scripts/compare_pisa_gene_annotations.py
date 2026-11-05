#!/usr/bin/env python3
"""Audit a manageable BAM slice before/after PISA annotation (not a full-BAM join)."""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

import pysam


ANNOTATION_TAGS = {
    "GX", "GN", "RE", "TX", "AT", "EX", "JC", "FL", "ER",
    "ZA", "ZG", "ZN", "ZB", "ZQ", "ZF", "ZC", "ZR", "ZE", "ZU", "ZX", "ZJ",
}


def tag(read, name):
    return str(read.get_tag(name)) if read.has_tag(name) else ""


def alignment_core(read):
    return (
        read.query_name, read.flag, read.reference_name, read.reference_start,
        read.mapping_quality, read.cigarstring, read.next_reference_name,
        read.next_reference_start, read.template_length, read.query_sequence,
        read.qual,
    )


def record_kind(read):
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    return "primary"


def load_bam(path):
    records = collections.defaultdict(list)
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam:
            key = hashlib.sha256(repr(alignment_core(read)).encode()).hexdigest()
            records[key].append(read)
    return records


def other_tags(read):
    # repr also handles array-valued tags without converting their numeric type.
    return sorted((k, repr(v), t) for k, v, t in read.get_tags(with_value_type=True)
                  if k not in ANNOTATION_TAGS)


def counters(records):
    counts = collections.defaultdict(collections.Counter)
    genes = collections.defaultdict(collections.Counter)
    ids = set()
    for reads in records.values():
        for read in reads:
            ids.add(read.query_name)
            gx = tag(read, "GX")
            status = tag(read, "ZA")
            for scope in ("all", record_kind(read)):
                c = counts[scope]
                c["records"] += 1
                c["with_GX" if gx else "without_GX"] += 1
                if status:
                    c["status_" + status] += 1
                if gx and not gx.startswith(("AMBIG:", "FUSION:")):
                    c["single_gene_GX" if ";" not in gx else "legacy_multigene_GX"] += 1
                    genes[scope][tag(read, "GN") or gx] += 1
                if gx.startswith("AMBIG:"):
                    c["ambiguous_set_GX"] += 1
                if gx.startswith("FUSION:"):
                    c["fusion_set_GX"] += 1
    return {"counts": dict(counts), "gene_labels": dict(genes),
            "distinct_qnames": len(ids)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--old", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--new", required=True)
    parser.add_argument("--old-categories")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=False)
    paths = {"input": args.input}
    for item in args.old:
        label, path = item.split("=", 1)
        if label in paths or label == "new":
            parser.error("duplicate/reserved label: " + label)
        paths[label] = path
    paths["new"] = args.new
    datasets = {name: load_bam(path) for name, path in paths.items()}
    baseline = datasets["input"]
    baseline_counts = {key: len(value) for key, value in baseline.items()}
    checks = {}
    for label, records in datasets.items():
        if {key: len(value) for key, value in records.items()} != baseline_counts:
            raise AssertionError(f"{label}: alignment core or multiplicity changed")
        for key, originals in baseline.items():
            if sorted(map(other_tags, originals)) != sorted(map(other_tags, records[key])):
                raise AssertionError(f"{label}: non-annotation tags changed for {originals[0].query_name}")
        checks[label] = {"all_alignment_fields_preserved": True,
                         "non_annotation_tags_preserved": True}

    categories = {}
    if args.old_categories:
        with open(args.old_categories) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                categories[row["read_id"]] = row["category"]
    transitions = collections.defaultdict(collections.Counter)
    category_counts = collections.defaultdict(collections.Counter)
    new = datasets["new"]
    old_labels = [label for label in datasets if label not in ("input", "new")]
    fields = ["read_id", "alignment_type", "flag", "chrom", "start_1based", "end_1based",
              "MAPQ", "CB", "UR", "old_unannotated_category"]
    for label in old_labels:
        fields.extend([label + "_GX", label + "_GN", label + "_RE", label + "_ZA"])
    fields.extend(["new_GX", "new_GN", "new_ZA", "new_ZF", "new_ZC", "new_ZR",
                   "alignment_candidates_ZG", "alignment_support_ZB", "alignment_matched_ZQ",
                   "alignment_exon_ZE", "alignment_exclusive_body_ZU",
                   "alignment_exclusive_exon_ZX", "alignment_known_junctions_ZJ", "CIGAR"])
    with gzip.open(out / "per_alignment_comparison.tsv.gz", "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key, reads in new.items():
            for index, read in enumerate(reads):
                kind = record_kind(read)
                gx = tag(read, "GX")
                status = tag(read, "ZA")
                category = categories.get(read.query_name, "") if kind == "primary" else ""
                row = dict(read_id=read.query_name, alignment_type=kind, flag=read.flag,
                           chrom=read.reference_name, start_1based=read.reference_start + 1,
                           end_1based=read.reference_end, MAPQ=read.mapping_quality,
                           CB=tag(read, "CB"), UR=tag(read, "UR"), old_unannotated_category=category,
                           new_GX=gx, new_GN=tag(read, "GN"), new_ZA=status,
                           new_ZF=tag(read, "ZF"), new_ZC=tag(read, "ZC"), new_ZR=tag(read, "ZR"),
                           alignment_candidates_ZG=tag(read, "ZG"),
                           alignment_support_ZB=tag(read, "ZB"),
                           alignment_matched_ZQ=tag(read, "ZQ"),
                           alignment_exon_ZE=tag(read, "ZE"),
                           alignment_exclusive_body_ZU=tag(read, "ZU"),
                           alignment_exclusive_exon_ZX=tag(read, "ZX"),
                           alignment_known_junctions_ZJ=tag(read, "ZJ"),
                           CIGAR=read.cigarstring)
                if category:
                    category_counts[category][status or "no_status"] += 1
                for label in old_labels:
                    old = datasets[label][key][index]
                    old_gx = tag(old, "GX")
                    row.update({label + "_GX": old_gx, label + "_GN": tag(old, "GN"),
                                label + "_RE": tag(old, "RE"), label + "_ZA": tag(old, "ZA")})
                    if kind == "primary":
                        if old_gx and gx == old_gx:
                            change = "same_GX"
                        elif old_gx:
                            change = "old_GX_to_" + (status or "unknown")
                        else:
                            change = "old_no_GX_to_" + (status or "unknown")
                        transitions[label][change] += 1
                writer.writerow(row)

    metrics = {label: counters(records) for label, records in datasets.items()}
    summary = {"paths": paths, "validation": checks, "metrics": metrics,
               "primary_transitions": dict(transitions),
               "old_unannotated_primary_outcomes": dict(category_counts),
               "scope": "Alignment-matched audit of supplied slice; GX presence is not accuracy."}
    with open(out / "summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    count_fields = ["dataset", "scope", "metric", "count"]
    with open(out / "counts.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(count_fields)
        for label, value in metrics.items():
            for scope, counts in value["counts"].items():
                for metric, count in sorted(counts.items()):
                    writer.writerow([label, scope, metric, count])
    print(json.dumps({"validation": checks, "primary": {
        label: value["counts"].get("primary", {}) for label, value in metrics.items()
    }}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
