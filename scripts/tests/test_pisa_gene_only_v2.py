"""Integration tests for the explicitly experimental PISA gene-only v2 draft.

These tests deliberately retain the draft's documented limitations: a one-base
exon touch can prevent weak-candidate reduction, and nested bodies can prevent
selection despite gene-specific exon evidence. They do not claim biological
accuracy or silently test a revised algorithm.
"""
from __future__ import annotations

import collections
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import pysam


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("PISA_TEST_BINARY", ROOT / "PISA"))
OWNED_TAGS = {
    "GX", "GN", "TX", "RE", "AT", "EX", "JC", "FL", "ER", "ZA", "ZG",
    "ZN", "ZB", "ZQ", "ZF", "ZC", "ZR", "ZE", "ZU", "ZX", "ZJ",
}


def run_checked(argv):
    result = subprocess.run(argv, text=True, capture_output=True, timeout=120)
    if result.returncode:
        raise AssertionError(
            f"Command failed ({result.returncode}): {argv!r}\n"
            f"{result.stdout[-1500:]}\n{result.stderr[-5000:]}"
        )
    return result


def records(path):
    with pysam.AlignmentFile(str(path), "rb") as handle:
        return list(handle.fetch(until_eof=True))


def core(read):
    return tuple(read.to_string().split("\t")[:11])


def support(read, tag="ZB"):
    if not read.has_tag("ZG"):
        return {}
    ids = read.get_tag("ZG").split(";")
    counts = list(map(int, read.get_tag(tag).split(";")))
    if len(ids) != len(counts):
        raise AssertionError(f"{tag} does not have one value per ZG gene")
    return dict(zip(ids, counts))


def make_fixture(directory):
    # Genes and exons are 0-based half-open here. Every scenario has its own
    # reference so unrelated annotations cannot accidentally affect its result.
    models = []

    def gene(chrom, gid, start, end, transcripts):
        models.append((chrom, gid, start, end, transcripts))

    for chrom, prefix, b_exon_start in [
        ("partial", "PART", 1500), ("one_exonic_base", "ONE", 1402)
    ]:
        gene(chrom, prefix + "_A", 0, 1238, [[(0, 1238)]])
        gene(chrom, prefix + "_B", 1238, 1600, [[(b_exon_start, 1600)]])
    gene("threshold", "THR_A", 0, 9000, [[(0, 9000)]])
    gene("threshold", "THR_B", 9000, 9500, [[(9400, 9500)]])
    gene("shared", "SH_A", 0, 5000, [[(0, 5000)]])
    gene("shared", "SH_B", 4899, 5000, [[(4899, 5000)]])
    gene("nested", "NEST_HOST", 0, 3000, [[(0, 1000), (2000, 3000)]])
    gene("nested", "NEST_INNER", 1000, 2000, [[(1000, 2000)]])
    gene("weak_disjoint", "DIS_A", 0, 800, [[(0, 800)]])
    gene("weak_disjoint", "DIS_B", 1000, 1200, [[(1100, 1200)]])
    gene("weak_disjoint", "DIS_C", 2000, 2200, [[(2100, 2200)]])
    gene("weak_shared", "UNI_A", 0, 900, [[(0, 900)]])
    gene("weak_shared", "UNI_B", 1000, 1200, [[(1100, 1200)]])
    gene("weak_shared", "UNI_C", 1000, 1200, [[(1100, 1200)]])
    gene("intronic_main", "INT_A", 0, 1000, [[(950, 1000)]])
    gene("intronic_main", "INT_B", 1000, 1200, [[(1100, 1200)]])
    gene("junction", "JUNC_A", 0, 3000, [[(0, 3000)]])
    gene("junction", "JUNC_B", 2500, 2900, [[(2500, 2600), (2800, 2900)]])
    gene("different_junctions", "DIFF_A", 0, 5000, [[(0, 2000), (2200, 5000)]])
    gene("different_junctions", "DIFF_B", 3300, 3700, [[(3300, 3400), (3600, 3700)]])
    gene("shared_junctions", "SAME_A", 0, 5000,
         [[(0, 3400), (3600, 5000)], [(0, 2000), (2200, 3400), (3600, 5000)]])
    gene("shared_junctions", "SAME_B", 3300, 3700, [[(3300, 3400), (3600, 3700)]])
    gene("exon_union", "UNION", 0, 1000,
         [[(0, 400), (600, 1000)], [(0, 500), (600, 1000)]])
    gene("cigar", "OPS", 0, 500, [[(0, 500)]])
    gene("skipped", "SKIP", 100, 200, [[(100, 200)]])
    gene("novel", "NOVEL", 0, 1000, [[(0, 100), (900, 1000)]])

    gtf = directory / "genes.gtf"
    with gtf.open("w") as handle:
        for chrom, gid, start, end, transcripts in models:
            attrs = f'gene_id "{gid}"; gene_name "Gene{gid}";'
            handle.write(f"{chrom}\ttest\tgene\t{start+1}\t{end}\t.\t+\t.\t{attrs}\n")
            for i, exons in enumerate(transcripts):
                tx_attrs = attrs + f' transcript_id "{gid}.{i+1}";'
                handle.write(f"{chrom}\ttest\ttranscript\t{start+1}\t{end}\t.\t+\t.\t{tx_attrs}\n")
                for left, right in exons:
                    handle.write(f"{chrom}\ttest\texon\t{left+1}\t{right}\t.\t+\t.\t{tx_attrs}\n")

    chromosomes = list(dict.fromkeys(g[0] for g in models))
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": 20000} for chrom in chromosomes],
        "RG": [{"ID": "rg1", "LB": "library1", "SM": "sample"}],
    })
    all_reads = []

    def add(name, chrom="partial", start=0, cigar="1403M8597N589M", flag=0, mapq=60):
        read = pysam.AlignedSegment(header)
        read.query_name = name
        read.flag = flag
        read.reference_id = chromosomes.index(chrom) if not flag & 4 else -1
        read.reference_start = start if not flag & 4 else -1
        read.mapping_quality = mapq if not flag & 4 else 0
        read.cigarstring = cigar
        length = read.infer_query_length() if cigar else 80
        read.query_sequence = "A" * length
        read.query_qualities = pysam.qualitystring_to_array("I" * length)
        for key, value in [("CB", "CELL1"), ("UR", "ACGTACGTAA"),
                           ("DB", "KEEP"), ("RG", "rg1"), ("NM", 0)]:
            read.set_tag(key, value)
        # Includes v2 tags to make stale evidence cleanup visible.
        for key in ["GX", "GN", "ZF", "ZR", "ZC", "ZE", "ZU", "ZX", "ZJ"]:
            read.set_tag(key, "STALE")
        read.set_tag("TX", "stale_transcript")
        read.set_tag("RE", "E", value_type="A")
        all_reads.append(read)

    add("partial_main")
    add("partial_reverse", flag=16)
    add("partial_secondary", flag=256)
    add("partial_supplementary", flag=2048)
    add("partial_low_mapq", mapq=5)
    add("partial_unknown_mapq", mapq=255)
    add("partial_qcfail", flag=512)
    add("one_exonic_base", chrom="one_exonic_base")
    for count in (99, 100, 101):
        add(f"minor_{count}", "threshold", cigar=f"{9000+count}M")
    add("shared_101", "shared", cigar="5000M")
    add("all_shared", "shared", start=4900, cigar="100M")
    add("nested_exon_specific", "nested", start=1000, cigar="1000M")
    add("weak_total_too_large", "weak_disjoint", cigar="800M200N80M920N80M920N40M")
    add("weak_union_not_sum", "weak_shared", cigar="900M100N80M1920N20M")
    add("intronic_dominant_main", "intronic_main", cigar="900M100N80M1920N20M")
    add("known_junction_protects", "junction", cigar="2600M200N100M")
    add("junction_shift_is_not_exact", "junction", cigar="2600M201N99M")
    add("junction_sets_differ_despite_equal_count", "different_junctions", cigar="2000M200N1200M200N100M")
    add("shared_junction_set_is_explained", "shared_junctions", cigar="2000M200N1200M200N100M")
    add("union_unspliced", "exon_union", cigar="1000M")
    add("union_junction1", "exon_union", cigar="400M200N400M")
    add("union_junction2", "exon_union", cigar="500M100N400M")
    add("mixed_cigar", "cigar", start=100, cigar="30S40M20I30=30X20D100M40S")
    add("single_short", "cigar", start=100, cigar="89M")
    add("single_partial", "cigar", start=480, cigar="100M")
    add("single_low_mapq", "cigar", start=100, cigar="89M", mapq=5)
    add("single_qcfail", "cigar", start=100, cigar="89M", flag=512)
    add("deletion_has_no_gene_support", "skipped", cigar="100M100D100M")
    add("skip_has_no_gene_support", "skipped", cigar="100M100N100M")
    add("novel_junction", "novel", start=200, cigar="100M100N100M")
    add("unmapped", flag=4, cigar=None)
    raw = directory / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as handle:
        for read in all_reads:
            handle.write(read)
    bam = directory / "input.sorted.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return bam, gtf


class TestPisaGeneOnlyV2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not BINARY.is_file():
            raise RuntimeError(f"Build PISA first or set PISA_TEST_BINARY: {BINARY}")
        cls.temp = tempfile.TemporaryDirectory(prefix="pisa-gene-v2-test-")
        cls.directory = Path(cls.temp.name)
        cls.bam, cls.gtf = make_fixture(cls.directory)
        cls.annotated = cls.directory / "annotated.bam"
        cls.command = [
            str(BINARY), "anno", "-gene-only-v2", "-ignore-strand", "-q", "20",
            "-gtf", str(cls.gtf), "-o", str(cls.annotated),
            "-report", str(cls.directory / "annotation.csv"),
            "-@", "2", "-t", "2", str(cls.bam),
        ]
        run_checked(cls.command)
        cls.before = records(cls.bam)
        cls.after = records(cls.annotated)
        cls.by_name = {r.query_name: r for r in cls.after}

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "temp"):
            cls.temp.cleanup()

    def test_partial_gx_retained_with_all_original_evidence(self):
        read = self.by_name["partial_main"]
        self.assertEqual(read.get_tag("GX"), "PART_A")
        self.assertEqual(read.get_tag("ZA"), "assigned_partial")
        self.assertEqual(read.get_tag("ZC"), "PART_A")
        self.assertEqual(read.get_tag("ZQ"), 1992)
        self.assertEqual(support(read), {"PART_A": 1238, "PART_B": 165})
        self.assertEqual(support(read, "ZE"), {"PART_A": 1238, "PART_B": 0})
        self.assertEqual(support(read, "ZU"), {"PART_A": 1238, "PART_B": 165})
        self.assertEqual(support(read, "ZX"), {"PART_A": 1238, "PART_B": 0})
        self.assertEqual(support(read, "ZJ"), {"PART_A": 0, "PART_B": 0})
        self.assertTrue(read.get_tag("ZR"))

    def test_known_draft_limit_one_exonic_base_changes_candidate_set(self):
        read = self.by_name["one_exonic_base"]
        self.assertEqual(support(read, "ZE"), {"ONE_A": 1238, "ONE_B": 1})
        self.assertEqual(read.get_tag("GX"), "AMBIG:ONE_A|ONE_B")
        self.assertEqual(read.get_tag("ZA"), "ambiguous")

    def test_known_draft_limit_nested_gene_has_no_exclusive_body(self):
        read = self.by_name["nested_exon_specific"]
        self.assertEqual(support(read, "ZU"), {"NEST_HOST": 0, "NEST_INNER": 0})
        self.assertEqual(support(read, "ZX"), {"NEST_HOST": 0, "NEST_INNER": 1000})
        self.assertEqual(read.get_tag("GX"), "AMBIG:NEST_HOST|NEST_INNER")

    def test_fixed_100_bp_boundary_no_longer_flips_weak_intronic_candidates(self):
        for count in (99, 100, 101):
            with self.subTest(count=count):
                read = self.by_name[f"minor_{count}"]
                self.assertEqual(support(read), {"THR_A": 9000, "THR_B": count})
                self.assertEqual(read.get_tag("GX"), "THR_A")
                self.assertEqual(read.get_tag("ZA"), "unique")

    def test_shared_101_bp_does_not_veto_dominant_gene(self):
        read = self.by_name["shared_101"]
        self.assertEqual(support(read), {"SH_A": 5000, "SH_B": 101})
        self.assertEqual(support(read, "ZU"), {"SH_A": 4899, "SH_B": 0})
        self.assertEqual(read.get_tag("GX"), "SH_A")
        self.assertEqual(read.get_tag("ZA"), "unique")

    def test_all_shared_remains_ambiguous(self):
        read = self.by_name["all_shared"]
        self.assertEqual(support(read), {"SH_A": 100, "SH_B": 100})
        self.assertEqual(support(read, "ZU"), {"SH_A": 0, "SH_B": 0})
        self.assertEqual(read.get_tag("GX"), "AMBIG:SH_A|SH_B")

    def test_individually_weak_candidates_cannot_hide_large_combined_support(self):
        read = self.by_name["weak_total_too_large"]
        self.assertEqual(read.get_tag("ZQ"), 1000)
        self.assertEqual(support(read), {"DIS_A": 800, "DIS_B": 80, "DIS_C": 80})
        self.assertEqual(read.get_tag("GX"), "AMBIG:DIS_A|DIS_B|DIS_C")

    def test_weak_candidate_combination_uses_union_not_sum(self):
        read = self.by_name["weak_union_not_sum"]
        self.assertEqual(read.get_tag("ZQ"), 1000)
        self.assertEqual(support(read), {"UNI_A": 900, "UNI_B": 80, "UNI_C": 80})
        self.assertEqual(support(read, "ZU"), {"UNI_A": 900, "UNI_B": 0, "UNI_C": 0})
        self.assertEqual(read.get_tag("GX"), "UNI_A")

    def test_exact_known_junction_prevents_shared_candidate_reduction(self):
        read = self.by_name["known_junction_protects"]
        self.assertEqual(support(read, "ZJ"), {"JUNC_A": 0, "JUNC_B": 1})
        self.assertEqual(read.get_tag("GX"), "AMBIG:JUNC_A|JUNC_B")
        shifted = self.by_name["junction_shift_is_not_exact"]
        self.assertEqual(support(shifted, "ZJ"), {"JUNC_A": 0, "JUNC_B": 0})
        self.assertEqual(shifted.get_tag("GX"), "JUNC_A")

    def test_exon_union_does_not_double_count_transcripts(self):
        read = self.by_name["union_unspliced"]
        self.assertEqual(support(read), {"UNION": 1000})
        self.assertEqual(support(read, "ZE"), {"UNION": 900})
        self.assertEqual(support(read, "ZJ"), {"UNION": 0})
        for name, count in [("union_junction1", 800), ("union_junction2", 900)]:
            with self.subTest(name=name):
                read = self.by_name[name]
                self.assertEqual(support(read, "ZE"), {"UNION": count})
                self.assertEqual(support(read, "ZJ"), {"UNION": 1})

    def test_junction_subset_compares_coordinates_not_counts(self):
        different = self.by_name["junction_sets_differ_despite_equal_count"]
        self.assertEqual(support(different, "ZJ"), {"DIFF_A": 1, "DIFF_B": 1})
        self.assertEqual(different.get_tag("GX"), "AMBIG:DIFF_A|DIFF_B")
        shared = self.by_name["shared_junction_set_is_explained"]
        self.assertEqual(support(shared, "ZJ"), {"SAME_A": 2, "SAME_B": 1})
        self.assertEqual(shared.get_tag("GX"), "SAME_A")

    def test_cigar_operators_only_contribute_actual_aligned_bases(self):
        read = self.by_name["mixed_cigar"]
        self.assertEqual(read.get_tag("ZQ"), 200)
        self.assertEqual(support(read), {"OPS": 200})
        self.assertEqual(support(read, "ZE"), {"OPS": 200})
        for name in ["deletion_has_no_gene_support", "skip_has_no_gene_support"]:
            with self.subTest(name=name):
                read = self.by_name[name]
                self.assertEqual(read.get_tag("ZQ"), 200)
                self.assertFalse(read.has_tag("GX"))
                self.assertFalse(read.has_tag("ZG"))
                self.assertEqual(read.get_tag("ZA"), "unassigned")
                for key in ["ZE", "ZU", "ZX", "ZJ", "ZF", "ZC"]:
                    self.assertFalse(read.has_tag(key), key)

    def test_single_candidate_short_and_partial_keep_gx(self):
        for name, status, bp in [("single_short", "unique", 89),
                                 ("single_partial", "assigned_partial", 20)]:
            with self.subTest(name=name):
                read = self.by_name[name]
                self.assertEqual(read.get_tag("GX"), "OPS")
                self.assertEqual(read.get_tag("ZA"), status)
                self.assertEqual(support(read), {"OPS": bp})

    def test_unknown_or_low_quality_retains_gx_without_reducing_candidates(self):
        for name, status in [("partial_low_mapq", "low_mapq"),
                             ("partial_unknown_mapq", "unknown_mapq"),
                             ("partial_qcfail", "qc_fail")]:
            with self.subTest(name=name):
                read = self.by_name[name]
                self.assertEqual(read.get_tag("GX"), "AMBIG:PART_A|PART_B")
                self.assertEqual(read.get_tag("ZA"), status)
        self.assertEqual(self.by_name["single_low_mapq"].get_tag("GX"), "OPS")
        self.assertEqual(self.by_name["single_qcfail"].get_tag("GX"), "OPS")

    def test_new_junctions_do_not_require_reference_transcript_match(self):
        read = self.by_name["novel_junction"]
        self.assertEqual(read.get_tag("GX"), "NOVEL")
        self.assertEqual(read.get_tag("ZA"), "unique")
        self.assertEqual(support(read, "ZE"), {"NOVEL": 0})
        self.assertEqual(support(read, "ZJ"), {"NOVEL": 0})
        intronic = self.by_name["intronic_dominant_main"]
        self.assertEqual(support(intronic, "ZE"), {"INT_A": 0, "INT_B": 0})
        self.assertEqual(intronic.get_tag("GX"), "INT_A")

    def test_reverse_secondary_and_supplementary_are_annotated(self):
        for name in ["partial_reverse", "partial_secondary", "partial_supplementary"]:
            with self.subTest(name=name):
                read = self.by_name[name]
                self.assertEqual(read.get_tag("GX"), "PART_A")
                self.assertEqual(read.get_tag("ZA"), "assigned_partial")
                self.assertEqual(support(read), {"PART_A": 1238, "PART_B": 165})

    def test_all_records_core_fields_and_nonannotation_tags_preserved(self):
        self.assertEqual(collections.Counter(map(core, self.before)),
                         collections.Counter(map(core, self.after)))
        before = {core(read): read for read in self.before}
        for read in self.after:
            original = before[core(read)]
            self.assertEqual(sorted((k, repr(v), t) for k, v, t in original.get_tags(with_value_type=True)
                                    if k not in OWNED_TAGS),
                             sorted((k, repr(v), t) for k, v, t in read.get_tags(with_value_type=True)
                                    if k not in OWNED_TAGS), read.query_name)
            for key in OWNED_TAGS:
                if read.has_tag(key):
                    self.assertNotEqual(read.get_tag(key), "STALE", (read.query_name, key))
        unmapped = self.by_name["unmapped"]
        self.assertFalse(unmapped.has_tag("GX"))
        self.assertEqual(unmapped.get_tag("ZA"), "unmapped")

    def test_invalid_or_conflicting_parameters_rejected(self):
        for option, value in [
            ("-gene-v2-main-frac", "nan"), ("-gene-v2-exclusive-frac", "1.1"),
            ("-gene-v2-exon-frac", "-0.1"), ("-gene-v2-weak-frac", "nonsense"),
            ("-gene-min-bp", "100"), ("-gene-min-frac", "0.8"),
            ("-gene-minor-bp", "100"), ("-gene-minor-frac", "0.2"),
        ]:
            with self.subTest(option=option, value=value):
                result = subprocess.run(self.command + [option, value], capture_output=True,
                                        text=True, timeout=30)
                self.assertNotEqual(result.returncode, 0)
        result = subprocess.run(self.command + ["-anno-only"], capture_output=True,
                                text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)

    def test_v2_parameters_do_not_silently_affect_v1(self):
        command = [str(BINARY), "anno", "-gene-only", "-ignore-strand", "-gtf", str(self.gtf),
                   "-o", str(self.directory / "invalid-v1.bam"), str(self.bam),
                   "-gene-v2-weak-frac", "0.1"]
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)

    def test_configurable_weak_fraction_changes_only_supported_decision(self):
        destination = self.directory / "weak_fraction_005.bam"
        command = [str(BINARY), "anno", "-gene-only-v2", "-ignore-strand", "-q", "20",
                   "-gtf", str(self.gtf), "-o", str(destination),
                   "-gene-v2-weak-frac", "0.05", str(self.bam)]
        run_checked(command)
        reads = {read.query_name: read for read in records(destination)}
        self.assertEqual(reads["partial_main"].get_tag("GX"), "AMBIG:PART_A|PART_B")
        self.assertEqual(reads["minor_101"].get_tag("GX"), "THR_A")

    def test_configurable_main_evidence_thresholds_are_applied(self):
        for option, value, name, expected in [
            ("-gene-v2-main-frac", "0.95", "intronic_dominant_main", "AMBIG:INT_A|INT_B"),
            ("-gene-v2-exclusive-frac", "0.70", "partial_main", "AMBIG:PART_A|PART_B"),
            ("-gene-v2-exon-frac", "0.70", "partial_main", "AMBIG:PART_A|PART_B"),
        ]:
            with self.subTest(option=option):
                destination = self.directory / (option.lstrip("-") + ".bam")
                command = [str(BINARY), "anno", "-gene-only-v2", "-ignore-strand", "-q", "20",
                           "-gtf", str(self.gtf), "-o", str(destination), option, value, str(self.bam)]
                run_checked(command)
                reads = {read.query_name: read for read in records(destination)}
                self.assertEqual(reads[name].get_tag("GX"), expected)


if __name__ == "__main__":
    unittest.main()
