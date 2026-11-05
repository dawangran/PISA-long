"""Regression tests for gene-only annotation and read-level reconciliation.

Run with:
    PISA_TEST_BINARY=/path/to/PISA python -m unittest discover -s scripts/tests -v

The fixtures deliberately distinguish an alternative mapping of the SAME query
segment from supplementary mappings of DIFFERENT query segments. No real data
or network access is required.
"""

from __future__ import annotations

import collections
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("PISA_TEST_BINARY", ROOT / "PISA"))
RESOLVER = ROOT / "scripts" / "pisa_resolve_gene_reads.py"
OWNED_TAGS = {
    "GX", "GN", "TX", "RE", "AT", "EX", "JC", "FL", "ER",
    "ZG", "ZN", "ZB", "ZQ", "ZA", "ZC", "ZR", "ZF",
}


def run_checked(argv: list[str]) -> subprocess.CompletedProcess:
    result = subprocess.run(argv, text=True, capture_output=True, timeout=120)
    if result.returncode:
        raise AssertionError(
            f"Command failed ({result.returncode}): {argv!r}\n"
            f"{result.stdout[-2000:]}\n{result.stderr[-5000:]}"
        )
    return result


def records(path: Path) -> list[pysam.AlignedSegment]:
    with pysam.AlignmentFile(str(path), "rb") as handle:
        return list(handle.fetch(until_eof=True))


def by_name(path: Path) -> dict[str, list[pysam.AlignedSegment]]:
    grouped = collections.defaultdict(list)
    for record in records(path):
        grouped[record.query_name].append(record)
    return dict(grouped)


def core(record: pysam.AlignedSegment) -> tuple[str, ...]:
    return tuple(record.to_string().split("\t")[:11])


def evidence(record: pysam.AlignedSegment) -> dict[str, int]:
    if not record.has_tag("ZG"):
        return {}
    genes = record.get_tag("ZG").split(";")
    counts = list(map(int, record.get_tag("ZB").split(";")))
    if len(genes) != len(counts):
        raise AssertionError("ZG/ZB evidence columns have different lengths")
    return dict(zip(genes, counts))


def make_fixture(directory: Path) -> tuple[Path, Path]:
    gtf = directory / "genes.gtf"
    # Coordinates here are 0-based, half-open; GTF written below is 1-based.
    genes = [
        ("chr1", "A", 100, 1000), ("chr1", "B", 2000, 3000),
        ("chr1", "C", 4000, 5000), ("chr1", "D", 4500, 5500),
        ("chr1", "HOST", 6000, 8000), ("chr1", "NEST", 6500, 7000),
        ("chr1", "G", 9000, 9100),
        ("chr1", "SMALL_HOST", 10000, 11000),
        ("chr1", "SMALL_1", 10000, 10090),
        ("chr1", "SMALL_2", 10080, 10170),
        ("chr1", "SMALL_3", 10160, 10240),
        ("chr1", "EDGE", 16383, 16483),
        ("chr2", "H", 100, 1000),
    ]
    with gtf.open("w") as handle:
        for chrom, gene, start, end in genes:
            attrs = f'gene_id "{gene}"; gene_name "Gene{gene}";'
            tx_attrs = attrs + f' transcript_id "{gene}.1";'
            for feature, left, right, attr in [
                ("gene", start, end, attrs),
                ("transcript", start, end, tx_attrs),
                ("exon", start, min(start + 100, end), tx_attrs),
                ("exon", max(start + 100, end - 100), end, tx_attrs),
            ]:
                if right > left:
                    handle.write(
                        f"{chrom}\ttest\t{feature}\t{left+1}\t{right}"
                        f"\t.\t+\t.\t{attr}\n"
                    )
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 40000}, {"SN": "chr2", "LN": 12000}],
        "RG": [{"ID": "rg1", "LB": "library1", "SM": "sample"}],
    })
    all_records = []

    def add(name, start, cigar, flag=0, chrom=0, mapq=60,
            cb="CELL1", ur="ACGTACGTAA", score=None, sa=None):
        read = pysam.AlignedSegment(header)
        read.query_name = name
        read.flag = flag
        read.reference_id = chrom if not flag & 4 else -1
        read.reference_start = start if not flag & 4 else -1
        read.mapping_quality = mapq if not flag & 4 else 0
        read.cigarstring = cigar
        query_length = read.infer_query_length() if cigar else 200
        read.query_sequence = "A" * query_length
        read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
        if cb is not None:
            read.set_tag("CB", cb)
        if ur is not None:
            read.set_tag("UR", ur)
        read.set_tag("DB", "KEEP")
        read.set_tag("RG", "rg1")
        read.set_tag("NM", 0)
        if score is not None:
            read.set_tag("AS", score)
        if sa is not None:
            read.set_tag("SA", sa)
        # Old annotation must not survive an unassigned or uncertain result.
        read.set_tag("GX", "STALE")
        read.set_tag("GN", "StaleGene")
        read.set_tag("TX", "stale_transcript")
        read.set_tag("RE", "E", value_type="A")
        read.set_tag("ZF", "fusion_candidate")
        read.set_tag("ZC", "STALE")
        read.set_tag("ZR", "stale_reason")
        all_records.append(read)

    add("unique", 200, "200M")
    add("overshoot", 850, "160M")  # 150 bp within A, 10 bp outside.
    add("novel_junction", 250, "100M100N100M")
    add("cigar_ops", 200, "40M30=30X20D100M")
    add("soft_insert", 200, "50S50M10I50M50S")
    add("deletion_only_gene", 8900, "100M100D100M")
    add("skip_only_gene", 8900, "100M100N100M")
    add("shared_region", 4600, "200M")
    add("nested_region", 6600, "200M")
    add("minor_touch", 4100, "401M")
    add("small_shared_neighbor", 4316, "250M")  # C=250, D=66, C exclusive=184.
    add("many_small_shared_neighbors", 10000, "250M")  # Global exclusive only 10.
    add("index_bin_boundary", 16383, "1M")
    add("secondary_alternative", 200, "200M", score=200)
    add("secondary_alternative", 2100, "200M", flag=256, score=198)
    add("weak_secondary", 200, "200M", score=200)
    add("weak_secondary", 2100, "200M", flag=256, score=80)
    add("same_gene_secondary", 200, "200M", score=200)
    add("same_gene_secondary", 500, "200M", flag=256, score=198)
    add("supp_fusion", 200, "200M200S", score=200)
    add("supp_fusion", 2100, "200H200M", flag=2048, score=200)
    add("sa_hard_clip_equivalence", 200, "200M200S",
        sa="chr1,2101,+,200S200M,60,0;")
    add("sa_hard_clip_equivalence", 2100, "200H200M", flag=2048,
        sa="chr1,201,+,200M200S,60,0;")
    add("secondary_supplementary", 200, "200M200S", score=200)
    add("secondary_supplementary", 2100, "200H200M", flag=2304, score=200)
    add("supp_same_query", 200, "200M200S", score=200)
    add("supp_same_query", 2100, "200M200H", flag=2048, score=200)
    add("reverse_fusion", 200, "200M200S", score=200)
    add("reverse_fusion", 2100, "200M200H", flag=2064, score=200)
    add("reverse_same_query", 200, "200M200S", score=200)
    add("reverse_same_query", 2100, "200H200M", flag=2064, score=200)
    add("same_gene_split", 200, "150M150S")
    add("same_gene_split", 500, "150H150M", flag=2048)
    add("low_quality_supp", 200, "200M200S")
    add("low_quality_supp", 2100, "200H200M", flag=2048, mapq=5)
    add("cross_chromosome", 200, "200M200S")
    add("cross_chromosome", 200, "200H200M", flag=2048, chrom=1)
    add("unresolved_third_gene", 200, "200M400S")
    add("unresolved_third_gene", 6200, "200H200M200S", flag=2048)
    add("unresolved_third_gene", 6600, "400H200M", flag=2048)
    add("missing_primary", 2100, "200H200M", flag=2048)
    add("barcode_conflict", 200, "150M150S")
    add("barcode_conflict", 500, "150H150M", flag=2048, cb="CELL2")
    add("missing_cb", 200, "200M", cb=None)
    add("missing_ur", 200, "200M", ur=None)
    add("low_mapq", 200, "200M", mapq=5)
    add("unknown_mapq", 200, "200M", mapq=255)
    add("unknown_mapq_supp", 200, "200M200S")
    add("unknown_mapq_supp", 2100, "200H200M", flag=2048, mapq=255)
    add("qc_fail", 200, "200M", flag=512)
    add("unmapped", -1, None, flag=4)
    add("multiple_primary", 200, "200M")
    add("multiple_primary", 2100, "200M")
    add("missing_sa_piece", 200, "200M200S",
        sa="chr1,2101,+,200S200M,60,0;")
    add("length_conflict", 200, "200M200S")
    add("length_conflict", 2100, "250H200M", flag=2048)
    raw = directory / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as handle:
        for read in all_records:
            handle.write(read)
    bam = directory / "input.sorted.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return bam, gtf


class PipelineFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not BINARY.is_file():
            raise RuntimeError(f"Build PISA first or set PISA_TEST_BINARY: {BINARY}")
        cls.temp = tempfile.TemporaryDirectory(prefix="pisa-gene-only-test-")
        cls.directory = Path(cls.temp.name)
        cls.bam, cls.gtf = make_fixture(cls.directory)
        cls.annotated = cls.directory / "alignment_annotated.bam"
        cls.command = [
            str(BINARY), "anno", "-gene-only", "-ignore-strand",
            "-gtf", str(cls.gtf), "-o", str(cls.annotated),
            "-report", str(cls.directory / "annotation.csv"),
            "-@", "2", "-t", "2", str(cls.bam),
        ]
        run_checked(cls.command)
        cls.original = by_name(cls.bam)
        cls.alignments = by_name(cls.annotated)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "temp"):
            cls.temp.cleanup()

    def assert_core_and_nonannotation_preserved(self, output):
        before = records(self.bam)
        after = records(output)
        self.assertEqual(collections.Counter(map(core, before)),
                         collections.Counter(map(core, after)))
        by_core = {core(read): read for read in before}
        for read in after:
            expected = {k: v for k, v in by_core[core(read)].get_tags()
                        if k not in OWNED_TAGS}
            actual = {k: v for k, v in read.get_tags() if k not in OWNED_TAGS}
            self.assertEqual(expected, actual, read.query_name)


class TestPisaGeneOnly(PipelineFixture):
    def test_exact_aligned_base_support(self):
        for name, count in [("cigar_ops", 200), ("soft_insert", 100),
                            ("novel_junction", 200), ("overshoot", 150)]:
            read = self.alignments[name][0]
            self.assertEqual(evidence(read), {"A": count}, name)
            self.assertEqual(read.get_tag("GX"), "A", name)
        self.assertEqual(self.alignments["overshoot"][0].get_tag("ZQ"), 160)
        self.assertEqual(self.alignments["cigar_ops"][0].get_tag("ZQ"), 200)
        self.assertEqual(self.alignments["soft_insert"][0].get_tag("ZQ"), 100)

    def test_skipped_reference_bases_do_not_annotate_gene(self):
        for name in ["deletion_only_gene", "skip_only_gene"]:
            read = self.alignments[name][0]
            self.assertEqual(evidence(read), {}, name)
            self.assertFalse(read.has_tag("GX"), name)
            self.assertEqual(read.get_tag("ZQ"), 200)

    def test_shared_and_nested_gene_evidence(self):
        self.assertEqual(evidence(self.alignments["shared_region"][0]),
                         {"C": 200, "D": 200})
        self.assertEqual(self.alignments["shared_region"][0].get_tag("GX"),
                         "AMBIG:C|D")
        self.assertEqual(self.alignments["nested_region"][0].get_tag("GX"),
                         "AMBIG:HOST|NEST")
        self.assertEqual(self.alignments["minor_touch"][0].get_tag("GX"), "C")

    def test_small_shared_neighbor_keeps_dominant_gene(self):
        read = self.alignments["small_shared_neighbor"][0]
        self.assertEqual(evidence(read), {"C": 250, "D": 66})
        self.assertEqual(read.get_tag("GX"), "C")
        self.assertEqual(read.get_tag("ZA"), "unique")

    def test_many_small_shared_neighbors_do_not_create_false_unique(self):
        read = self.alignments["many_small_shared_neighbors"][0]
        self.assertEqual(evidence(read), {
            "SMALL_HOST": 250, "SMALL_1": 90, "SMALL_2": 90, "SMALL_3": 80,
        })
        self.assertNotEqual(read.get_tag("ZA"), "unique")
        self.assertFalse(read.has_tag("GX"))

    def test_secondary_and_supplementary_keep_evidence(self):
        secondary = next(r for r in self.alignments["secondary_alternative"]
                         if r.is_secondary)
        supplement = next(r for r in self.alignments["supp_fusion"]
                          if r.is_supplementary)
        self.assertEqual(evidence(secondary), {"B": 200})
        self.assertEqual(evidence(supplement), {"B": 200})

    def test_all_records_and_nonannotation_content_preserved(self):
        self.assert_core_and_nonannotation_preserved(self.annotated)
        for read in records(self.annotated):
            self.assertFalse(read.has_tag("TX"), read.query_name)
            self.assertFalse(read.has_tag("RE"), read.query_name)
            self.assertFalse(read.has_tag("ZF"), read.query_name)
            self.assertFalse(read.has_tag("ZC"), read.query_name)
            self.assertFalse(read.has_tag("ZR"), read.query_name)
            if read.has_tag("GX"):
                self.assertNotEqual(read.get_tag("GX"), "STALE", read.query_name)

    def test_mapq_filter_marks_but_does_not_drop_records(self):
        output = self.directory / "mapq.bam"
        command = self.command.copy()
        command[command.index(str(self.annotated))] = str(output)
        command += ["-q", "20"]
        run_checked(command)
        self.assert_core_and_nonannotation_preserved(output)
        read = by_name(output)["low_mapq"][0]
        self.assertEqual(read.get_tag("ZA"), "low_mapq")
        self.assertFalse(read.has_tag("GX"))

    def test_gene_start_at_index_bin_boundary(self):
        output = self.directory / "boundary.bam"
        command = self.command.copy()
        command[command.index(str(self.annotated))] = str(output)
        run_checked(command + ["-gene-min-bp", "1"])
        read = by_name(output)["index_bin_boundary"][0]
        self.assertEqual(evidence(read), {"EDGE": 1})
        self.assertEqual(read.get_tag("GX"), "EDGE")

    def test_invalid_thresholds_and_destructive_mode_rejected(self):
        for option, value in [
            ("-gene-min-frac", "nan"), ("-gene-min-frac", "1.1"),
            ("-gene-minor-frac", "-0.1"), ("-gene-min-bp", "-1"),
            ("-gene-min-bp", "nonsense"), ("-gene-minor-bp", "-1"),
        ]:
            with self.subTest(option=option, value=value):
                result = subprocess.run(self.command + [option, value],
                                        capture_output=True, text=True, timeout=30)
                self.assertNotEqual(result.returncode, 0)
        result = subprocess.run(self.command + ["-anno-only"],
                                capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)

    def test_removed_v2_arguments_rejected_before_annotation(self):
        removed_options = [
            ["-gene-only-v2"], ["--gene-only-v2"],
            ["-gene-v2-main-frac", "0.8"],
            ["-gene-v2-exclusive-frac", "0.2"],
            ["-gene-v2-exon-frac", "0.5"],
            ["-gene-v2-weak-frac", "0.1"],
        ]
        # The input BAM/GTF are valid. Check both ordinary annotation and the
        # retained gene-only mode so rejection cannot be a mode conflict or a
        # missing-input error, and a removed option cannot silently fall back.
        for gene_only in (False, True):
            for index, options in enumerate(removed_options):
                with self.subTest(gene_only=gene_only, option=options[0]):
                    output = self.directory / f"removed-v2-{gene_only}-{index}.bam"
                    command = self.command.copy()
                    if not gene_only:
                        command.remove("-gene-only")
                    command[command.index(str(self.annotated))] = str(output)
                    result = subprocess.run(command + options, text=True,
                                            capture_output=True, timeout=30)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("Unknown argument", result.stderr)
                    self.assertIn(options[0], result.stderr)
                    self.assertFalse(output.exists())


class TestReadResolver(PipelineFixture):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        if not RESOLVER.is_file():
            raise RuntimeError(f"Resolver not present: {RESOLVER}")
        cls.output_dir = cls.directory / "resolved"
        cls.resolve_command = [
            sys.executable, str(RESOLVER), "--input", str(cls.annotated),
            "--gtf", str(cls.gtf), "--output-dir", str(cls.output_dir),
            "--threads", "2",
        ]
        run_checked(cls.resolve_command)
        cls.resolved_bam = cls.output_dir / "annotated.sorted.bam"
        cls.resolved = by_name(cls.resolved_bam)
        cls.umi = by_name(cls.output_dir / "umi_input.primary.sorted.bam")

    def tag(self, name, tag):
        values = {r.get_tag(tag) if r.has_tag(tag) else None
                  for r in self.resolved[name]}
        self.assertEqual(len(values), 1, f"Inconsistent read-level {tag}: {name}")
        return values.pop()

    def test_secondary_is_an_alternative_not_a_fusion(self):
        self.assertEqual(self.tag("secondary_alternative", "GX"), "AMBIG:A|B")
        self.assertEqual(self.tag("secondary_alternative", "ZA"), "ambiguous")
        self.assertIsNone(self.tag("secondary_alternative", "ZF"))
        self.assertEqual(self.tag("weak_secondary", "GX"), "A")
        self.assertEqual(self.tag("same_gene_secondary", "GX"), "A")
        self.assertEqual(self.tag("secondary_supplementary", "GX"), "A")
        self.assertIsNone(self.tag("secondary_supplementary", "ZF"))

    def test_distinct_query_segments_are_fusion_candidates(self):
        for name, genes in [("supp_fusion", "A|B"),
                            ("sa_hard_clip_equivalence", "A|B"),
                            ("reverse_fusion", "A|B"),
                            ("cross_chromosome", "A|H")]:
            with self.subTest(name=name):
                self.assertEqual(self.tag(name, "GX"), "FUSION:" + genes)
                self.assertEqual(self.tag(name, "ZA"), "fusion_candidate")
                self.assertEqual(self.tag(name, "ZF"), "fusion_candidate")

    def test_same_query_segments_cannot_supply_independent_fusion_anchors(self):
        for name in ["supp_same_query", "reverse_same_query"]:
            with self.subTest(name=name):
                self.assertNotEqual(self.tag(name, "ZA"), "fusion_candidate")
                self.assertIsNone(self.tag(name, "ZF"))
                self.assertEqual(self.tag(name, "GX"), "AMBIG:A|B")
        self.assertEqual(self.tag("same_gene_split", "GX"), "A")
        self.assertNotEqual(self.tag("low_quality_supp", "ZA"), "fusion_candidate")
        self.assertNotEqual(self.tag("unknown_mapq_supp", "ZA"), "fusion_candidate")

    def test_unresolved_third_gene_is_not_silently_dropped_from_fusion(self):
        self.assertEqual(self.tag("unresolved_third_gene", "ZA"), "ambiguous")
        self.assertEqual(self.tag("unresolved_third_gene", "GX"), "AMBIG:A|HOST|NEST")
        self.assertIsNone(self.tag("unresolved_third_gene", "ZF"))

    def test_shared_and_nested_regions_stay_ambiguous(self):
        for name, expected in [("shared_region", "AMBIG:C|D"),
                               ("nested_region", "AMBIG:HOST|NEST")]:
            self.assertEqual(self.tag(name, "GX"), expected)
            self.assertEqual(self.tag(name, "ZA"), "ambiguous")
            self.assertIsNone(self.tag(name, "ZF"))

    def test_small_shared_neighbor_keeps_dominant_gene(self):
        self.assertEqual(self.tag("small_shared_neighbor", "GX"), "C")
        self.assertEqual(self.tag("small_shared_neighbor", "ZA"), "unique")
        self.assertIsNone(self.tag("small_shared_neighbor", "ZF"))

    def test_many_small_shared_neighbors_do_not_create_false_unique(self):
        self.assertNotEqual(self.tag("many_small_shared_neighbors", "ZA"), "unique")
        self.assertIsNone(self.tag("many_small_shared_neighbors", "GX"))

    def test_incomplete_or_conflicting_reads_bypass_umi_input(self):
        statuses = {
            "missing_primary": "missing_primary",
            "multiple_primary": "multiple_primary",
            "missing_sa_piece": "incomplete_alignment_set",
            "barcode_conflict": "tag_conflict",
            "low_mapq": "low_mapq",
            "unknown_mapq": "unknown_mapq",
            "length_conflict": "read_length_conflict",
        }
        for name, status in statuses.items():
            with self.subTest(name=name):
                self.assertEqual(self.tag(name, "ZA"), status)
                self.assertIsNone(self.tag(name, "GX"))
                self.assertNotIn(name, self.umi)
        for name in ["missing_cb", "missing_ur", "qc_fail", "unmapped"]:
            self.assertNotIn(name, self.umi)

    def test_exactly_one_primary_per_umi_observation(self):
        self.assertIn("same_gene_split", self.umi)
        self.assertIn("secondary_alternative", self.umi)
        for name, aligned in self.umi.items():
            self.assertEqual(len(aligned), 1, name)
            self.assertFalse(aligned[0].flag & (4 | 256 | 2048), name)
            self.assertTrue(aligned[0].has_tag("CB"), name)
            self.assertTrue(aligned[0].has_tag("UR"), name)

    def test_all_source_records_survive_and_evidence_stays_aligned(self):
        self.assert_core_and_nonannotation_preserved(self.resolved_bam)
        original_evidence = {core(r): evidence(r) for r in records(self.annotated)}
        for read in records(self.resolved_bam):
            self.assertEqual(evidence(read), original_evidence[core(read)])
            if read.has_tag("GX"):
                self.assertNotEqual(read.get_tag("GX"), "STALE")

    def test_category_bams_include_complete_read_record_sets(self):
        paths = ["fusion_candidates.sorted.bam", "ambiguous.sorted.bam",
                 "bypass.sorted.bam"]
        for filename in paths:
            output = self.output_dir / filename
            with self.subTest(filename=filename):
                self.assertTrue(Path(str(output) + ".bai").is_file())
                with pysam.AlignmentFile(str(output), "rb") as handle:
                    self.assertTrue(handle.has_index())
                for name, selected in by_name(output).items():
                    self.assertEqual(collections.Counter(map(core, selected)),
                                     collections.Counter(map(core, self.resolved[name])))
        fusion_names = set(by_name(self.output_dir / paths[0]))
        expected = {name for name in self.resolved
                    if self.tag(name, "ZA") == "fusion_candidate"}
        self.assertEqual(fusion_names, expected)

    def test_audit_and_summary_are_present_and_valid(self):
        with gzip.open(self.output_dir / "read_assignments.tsv.gz", "rt") as handle:
            lines = [line for line in handle if line.strip()]
        self.assertEqual(len(lines) - 1, len(self.original))
        with (self.output_dir / "summary.json").open() as handle:
            self.assertIsInstance(json.load(handle), dict)

    def test_existing_output_is_not_silently_overwritten(self):
        result = subprocess.run(self.resolve_command, text=True,
                                capture_output=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)

    def test_empty_category_bams_are_valid_and_indexed(self):
        source = self.directory / "only_unique.bam"
        with pysam.AlignmentFile(str(self.annotated), "rb") as original:
            with pysam.AlignmentFile(str(source), "wb", header=original.header) as handle:
                handle.write(self.alignments["unique"][0])
        destination = self.directory / "only_unique_resolved"
        command = self.resolve_command.copy()
        command[command.index(str(self.annotated))] = str(source)
        command[command.index(str(self.output_dir))] = str(destination)
        run_checked(command)
        for stem in ["fusion_candidates", "ambiguous", "bypass"]:
            output = destination / f"{stem}.sorted.bam"
            with pysam.AlignmentFile(str(output), "rb") as handle:
                self.assertTrue(handle.has_index())
                self.assertEqual(list(handle.fetch(until_eof=True)), [])

    def test_same_read_name_in_distinct_libraries_is_kept_separate(self):
        source = self.directory / "two_libraries.bam"
        with pysam.AlignmentFile(str(self.annotated), "rb") as original:
            header = original.header.to_dict()
        header["RG"].append({"ID": "rg2", "LB": "library2", "SM": "sample2"})
        header = pysam.AlignmentHeader.from_dict(header)
        with pysam.AlignmentFile(str(source), "wb", header=header) as handle:
            for rg in ["rg1", "rg2"]:
                record = pysam.AlignedSegment.fromstring(
                    self.alignments["unique"][0].to_string(), header)
                record.set_tag("RG", rg)
                handle.write(record)
        destination = self.directory / "two_libraries_resolved"
        command = self.resolve_command.copy()
        command[command.index(str(self.annotated))] = str(source)
        command[command.index(str(self.output_dir))] = str(destination)
        run_checked(command)
        primary = records(destination / "umi_input.primary.sorted.bam")
        self.assertEqual(len(primary), 2)
        self.assertEqual({r.get_tag("RG") for r in primary}, {"rg1", "rg2"})
        self.assertTrue(all(r.get_tag("ZA") == "unique" for r in primary))
        with (destination / "summary.json").open() as handle:
            summary = json.load(handle)
        self.assertTrue(summary["umi_input_requires_library_split"])


if __name__ == "__main__":
    unittest.main()
