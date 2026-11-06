#!/usr/bin/env python3
"""Exercise installed runtime commands with synthetic data only.

Mount this directory read-only at /tests and run:
    python3 /tests/container_smoke.py

All generated fixtures, temporary files and results live under /tmp. No source
tree, compiler, network access, writable home directory or root UID is needed.
This is a runtime smoke check; the image builder runs the unittest suite.
"""
from __future__ import annotations

import collections
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import pysam


def run(command):
    result = subprocess.run(command, text=True, capture_output=True, timeout=120)
    if result.returncode:
        raise RuntimeError(
            f"Runtime command failed ({result.returncode}): {command!r}\n"
            f"{result.stdout[-2000:]}\n{result.stderr[-5000:]}"
        )
    return result


def fixture(directory):
    # The existing fixture module derives a repository root from its filename.
    # Preserve that directory depth in /tmp, without requiring a source checkout
    # in the final image or writing into the read-only /tests bind mount.
    module_path = directory / "fixture" / "scripts" / "tests" / "test_pisa_gene_only.py"
    module_path.parent.mkdir(parents=True)
    shutil.copyfile(Path(__file__).with_name("test_pisa_gene_only.py"), module_path)
    spec = importlib.util.spec_from_file_location("pisa_runtime_fixture", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the synthetic fixture helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Use only the fixture helper; no test classes are instantiated or run.
    return module.make_fixture(directory)


def bam_records(path):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return list(bam.fetch(until_eof=True))


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    commands = ["PISA", "samtools", "pisa_resolve_gene_reads.py",
                "compare_pisa_gene_annotations.py"]
    for command in commands:
        check(shutil.which(command) is not None, f"Missing installed command: {command}")

    with tempfile.TemporaryDirectory(prefix="pisa-container-smoke-", dir="/tmp") as name:
        directory = Path(name)
        source, gtf = fixture(directory)
        original = bam_records(source)
        annotated = directory / "gene_only.bam"
        sorted_bam = directory / "gene_only.sorted.bam"
        run(["PISA", "anno", "-gene-only", "-ignore-strand", "-t", "2",
             "-gtf", str(gtf), "-o", str(annotated), "-report",
             str(directory / "annotation.csv"), str(source)])
        run(["samtools", "sort", "-@", "2", "-o", str(sorted_bam), str(annotated)])
        run(["samtools", "index", str(sorted_bam)])
        resolved = directory / "resolved"
        run(["pisa_resolve_gene_reads.py", "--input", str(sorted_bam),
             "--gtf", str(gtf), "--output-dir", str(resolved), "--threads", "2"])

        # The installed comparator checks every core alignment field and all
        # non-annotation tags, independent of output ordering.
        audit = directory / "comparison"
        run(["compare_pisa_gene_annotations.py", "--input", str(source),
             "--old", f"alignment={sorted_bam}", "--new",
             str(resolved / "annotated.sorted.bam"), "--output-dir", str(audit)])
        audit_summary = json.loads((audit / "summary.json").read_text())
        check(all(all(flags.values()) for flags in audit_summary["validation"].values()),
              "Alignment fields or non-annotation tags changed")

        summary = json.loads((resolved / "summary.json").read_text())
        check(summary["annotation_only"] is True, "Unexpected non-annotation workflow")
        check(summary["umi_deduplication_performed"] is False, "Unexpected UMI deduplication")
        check(summary["counts"]["alignment_records"] == len(original),
              "Input record count changed")
        check(summary["output_alignment_records"]["annotated"] == len(original),
              "Full annotated BAM lost records")

        stems = ("annotated", "umi_input.primary", "fusion_candidates", "ambiguous", "bypass")
        for stem in stems:
            path = resolved / f"{stem}.sorted.bam"
            check(path.is_file() and Path(str(path) + ".bai").is_file(),
                  f"Missing BAM or index: {stem}")
            run(["samtools", "quickcheck", str(path)])
            with pysam.AlignmentFile(str(path), "rb") as bam:
                check(bam.header.to_dict().get("HD", {}).get("SO") == "coordinate",
                      f"Output is not marked coordinate-sorted: {stem}")
                check(bam.has_index(), f"Output index cannot be loaded: {stem}")

        with gzip.open(resolved / "read_assignments.tsv.gz", "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        check(len(rows) == summary["counts"]["read_groups"], "Read summary row count differs")
        expected = {row["read_id"] for row in rows if row["umi_eligible"] == "1"}
        primary = bam_records(resolved / "umi_input.primary.sorted.bam")
        observed = collections.Counter(read.query_name for read in primary)
        check(expected and set(observed) == expected, "UMI-input read set differs from summary")
        check(all(count == 1 for count in observed.values()), "Repeated read in primary-only BAM")
        check(all(not (read.flag & 0x900) for read in primary), "Non-primary record in UMI input")
        check(all(read.has_tag(tag) for read in primary for tag in ("GX", "CB", "UR")),
              "UMI input is missing required tags")
        molecule_keys = collections.Counter(
            tuple(read.get_tag(tag) for tag in ("GX", "CB", "UR")) for read in primary
        )
        check(any(count > 1 for count in molecule_keys.values()),
              "Synthetic reads sharing GX/CB/UR were unexpectedly collapsed")
        unique = next(read for read in primary if read.query_name == "unique")
        check(unique.get_tag("GX") == "A" and unique.get_tag("ZA") == "unique",
              "Synthetic single-gene read was not assigned to A")
        print(json.dumps({
            "runtime_smoke_passed": True,
            "uid": os.getuid(),
            "input_alignment_records": len(original),
            "resolved_read_groups": len(rows),
            "umi_input_primary_records": len(primary),
            "umi_deduplication_performed": False,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
