// annotate Gene or Peak to SAM attribution
#include "utils.h"
#include "number.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/khash_str2int.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "bam_pool.h"
#include "htslib/thread_pool.h"
#include "gtf.h"
#include "bed.h"
#include "region_index.h"
#include "read_anno.h"
#include "dict.h"
#include <zlib.h>
#include <limits.h>
#include <math.h>
#include "htslib/kseq.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)
enum gene_only_status {
    GO_UNIQUE, GO_AMBIGUOUS, GO_INSUFFICIENT, GO_UNASSIGNED,
    GO_UNMAPPED, GO_LOW_MAPQ, GO_QC_FAIL, GO_STATUS_COUNT
};
static const char *gene_only_status_names[GO_STATUS_COUNT] = {
    "unique", "ambiguous", "insufficient", "unassigned",
    "unmapped", "low_mapq", "qc_fail"
};
struct read_stat {

    uint64_t reads_in_region;
    uint64_t reads_in_region_diff_strand;
    
    // gtf
    uint64_t reads_in_intergenic;
    uint64_t reads_in_exon;
    uint64_t reads_in_exonintron;
    uint64_t reads_in_intron;
    uint64_t reads_intron_retain;
    uint64_t reads_exclude;
    uint64_t reads_antisense;
    uint64_t reads_antisenseintron;
    uint64_t reads_ambiguous;

    uint64_t reads_tss;
    uint64_t reads_anno_genes;
    uint64_t gene_only_counts[GO_STATUS_COUNT];
};

static struct args {
    const char *input_fname;
    const char *output_fname;
    const char *bed_fname;
    const char *vcf_fname;
    const char *vtag; // tag name for vcf
    const char *tag; // attribute in BAM
    const char *ctag; // tag name for TSS annotation
    
    const char *gtf_fname;
    const char *report_fname;
    const char *chr_spec_fname;

    const char *group_tag; // export summary information for each group
    
    char **chr_binding;
    const char *btag;
    
    int ignore_strand;
    int exonintron_flag;
    int intron_flag;
    int antisense;
    
    int flatten_flag;
    int exon_level;
    int reverse_trans;

    int vague_edge;
    int gene_only;
    int gene_min_bp;
    int gene_minor_bp;
    double gene_min_frac;
    double gene_minor_frac;

    // for percent spliced index (PSI)
    int psi;
    const char *psi_tags;
    
    int n_thread;
    int chunk_size;
    int tss_mode;
    int anno_only;

    int ref_alt;
    int vcf_ss;
    int phased;
    
    int input_sam;
    gzFile fp_sam;
    kstream_t *ks;
    char *preload_record;
    
    htsFile *fp;
    htsFile *out;
    bam_hdr_t *hdr;

    FILE *fp_report;

    // GTF
    struct gtf_spec *G;
    struct bed_spec *flatten; // for flatten exon 
    // bed
    struct bed_spec *B;

    // vcf
    struct bed_spec *V;
    
    uint64_t reads_input;
    uint64_t reads_pass_qc;
    int map_qual;

    struct dict *group_stat;

    int debug_mode;
} args = {
    .input_fname     = NULL,
    .output_fname    = NULL,
    .bed_fname       = NULL,
    .vcf_fname       = NULL,
    .vtag            = NULL,
    .tag             = NULL,    
    .gtf_fname       = NULL,
    .report_fname    = NULL,

    .group_tag       = NULL,
    
    .chr_binding     = NULL,
    .btag            = NULL,

    .ctag            = NULL,
    .tss_mode        = 0,
    .ignore_strand   = 0,
    .exonintron_flag = 0,
    .intron_flag     = 0,
    .antisense       = 0,
    .exon_level      = 0,
    .flatten_flag    = 0,
    .reverse_trans   = 0,
    .vague_edge      = 0,
    .gene_only       = 0,
    .gene_min_bp     = 100,
    .gene_minor_bp   = 100,
    .gene_min_frac   = 0.8,
    .gene_minor_frac = 0.2,
    .psi             = 0,
    .psi_tags        = NULL,
    
    .map_qual        = 0,
    .n_thread        = 4,
    .chunk_size      = 1000,
    .anno_only       = 0,

    .ref_alt         = 0,
    .vcf_ss          = 1,
    .phased          = 0,
    .input_sam       = 0,
    .fp_sam          = NULL,
    .ks              = NULL,
    .preload_record  = NULL,
    
    .fp              = NULL,
    .out             = NULL,
    .hdr             = NULL,
    .fp_report       = NULL,
    .G               = NULL,    
    .B               = NULL,    
    .V               = NULL,
    .reads_input     = 0,
    .reads_pass_qc   = 0,
    .group_stat      = 0,
    .debug_mode      = 0
};

struct kstring_pool {
    int n, m;
    kstring_t *str;
};
struct kstring_pool *kstring_pool_init(int size)
{
    struct kstring_pool *p = malloc(sizeof(*p));
    p->m = size;
    p->n = 0;
    p->str = malloc(p->m*sizeof(kstring_t));
    memset(p->str, 0, p->m*sizeof(kstring_t));
    return p;
}
void kstring_pool_destroy(struct kstring_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        kstring_t *s = &p->str[i];
        if (s->m) free(s->s);
    }
    free(p->str);
    free(p);
}
struct kstring_pool *kstring_pool_read(kstream_t *s, int buffer_size)
{
    struct kstring_pool *p = kstring_pool_init(buffer_size);
    
    if (args.preload_record) {
        kputs(args.preload_record, &p->str[0]);
        p->n++;
        free(args.preload_record);
        args.preload_record= NULL;
    }
    
    int ret;
    for (;;) {
        if (p->n == p->m) break;
        if (ks_getuntil(s, 2, &p->str[p->n], &ret) < 0) break;
        if (p->str[p->n].s[0] == '@') continue;
        p->n++;
    }
    
    if (p->n == 0) {
        kstring_pool_destroy(p);
        return NULL;
    }
    return p;
}

static char **chr_binding(const char *fname, bam_hdr_t *hdr)
{
    int n;
    char **list = hts_readlist(fname, 1, &n);
    if (n == 0) error("Empty file.");
    
    char **bind = malloc(n*sizeof(char*));
    struct dict *d = dict_init();
    int i;
    for (i = 0; i < n; ++i) {
        char *s = list[i];
        char *p;
        for (p = s; *p && *p != '\t'; p++);
        if (*p == '\t') *p = '\0';
        int idx = dict_push(d, s);
        if (idx != i) error("%s: Duplicated records ? %s", fname, s);
        s = p +1;
        bind[i] = strdup(s);
        free(list[i]);
    }
    free(list);

    char **ret = malloc(hdr->n_targets*sizeof(char*));
    for (i = 0; i < hdr->n_targets; ++i) {
        int idx = dict_query(d, hdr->target_name[i]);
        if (idx == -1) ret[i] = NULL;
        else {
            ret[i] = bind[idx];
            bind[idx] = NULL;
        }
    }
    
    for (i = 0; i < n; ++i)
        if (bind[i]) free(bind[i]);
    free(bind);
    return ret;
}

// default tags,
// TX for transcript name,
// GN for gene name,
// GX for gene id,
// AT for gene nane, but mapped to antisense.
// RE for region type,
// EX for exon name [chr:start-end/[+-]/Gene].
// JC for junction name
// FL for flatten exon name
// ER for exlcuded exons, usually used for PSI calculation
static char TX_tag[2] = "TX";
static char AT_tag[2] = "AT";
static char GN_tag[2] = "GN";
static char GX_tag[2] = "GX";
static char RE_tag[2] = "RE";
static char EX_tag[2] = "EX";
static char JC_tag[2] = "JC";
static char FL_tag[2] = "FL";
static char ER_tag[2] = "ER";
// default tag name for genetic variants, can change by -vtag 
static char VR_tag[2] = "VR";
// default tag name for BED, can change by -tag
static char PK_tag[2] = "PK";
// species tag
static char SP_tag[2] = "SP";
extern struct bed_spec *bed_read_vcf(const char *fn);
extern bam_hdr_t *sam_parse_header(kstream_t *s, kstring_t *line);

static int gene_only_int_arg(const char *value, const char *name, int minimum)
{
    char *end;
    errno = 0;
    long n = strtol(value, &end, 10);
    if (errno || end == value || *end || n < minimum || n > INT_MAX)
        error("%s must be an integer between %d and %d.", name, minimum, INT_MAX);
    return (int)n;
}

static double gene_only_fraction_arg(const char *value, const char *name)
{
    char *end;
    errno = 0;
    double n = strtod(value, &end);
    if (errno || end == value || *end || !isfinite(n) || n < 0.0 || n > 1.0)
        error("%s must be a finite fraction between 0 and 1.", name);
    return n;
}

static int parse_args(int argc, char **argv)
{
    int i;
    const char *tags = NULL;
    const char *thread = NULL;
    const char *chunk = NULL;
    const char *file_thread = NULL;
    const char *map_qual = NULL;
    const char *vague_edge = NULL;
    const char *gene_min_bp = NULL, *gene_minor_bp = NULL;
    const char *gene_min_frac = NULL, *gene_minor_frac = NULL;
    
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        // common options
        if (strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-@") == 0) var = &file_thread;
        else if (strcmp(a, "-chunk") == 0) var = &chunk;
        else if (strcmp(a, "-q") == 0) var = &map_qual;
        else if (strcmp(a, "-gene-only") == 0 || strcmp(a, "--gene-only") == 0) {
            args.gene_only = 1;
            continue;
        }
        else if (strcmp(a, "-gene-min-bp") == 0) var = &gene_min_bp;
        else if (strcmp(a, "-gene-minor-bp") == 0) var = &gene_minor_bp;
        else if (strcmp(a, "-gene-min-frac") == 0) var = &gene_min_frac;
        else if (strcmp(a, "-gene-minor-frac") == 0) var = &gene_minor_frac;
        else if (strcmp(a, "-sam") == 0) {
            args.input_sam = 1;
            continue;
        }
        else if (strcmp(a, "-debug") == 0) {
            args.debug_mode = 1;
            continue;
        }
        else if (strcmp(a, "-ignore-strand") == 0 || strcmp(a, "-is") == 0) {
            args.ignore_strand = 1;
            continue;
        }
        else if (strcmp(a, "-splice-consider") == 0 || strcmp(a, "-splice") == 0) {
            args.exonintron_flag = 1;
            continue;
        }
        else if (strcmp(a, "-intron") == 0 || strcmp(a, "-velo") == 0) {
            args.intron_flag = 1;
            continue;
        }
        else if (strcmp(a, "-as") == 0) {
            args.antisense = 1;
            continue;
        }
        else if (strcmp(a, "-exon") == 0) {
            args.exon_level = 1;
            continue;
        }
        else if (strcmp(a, "-flatten") == 0) {
            args.flatten_flag = 1;
            continue;
        }
        else if (strcmp(a, "-psi") == 0) {
            args.psi = 1;
            continue;
        }
        else if (strcmp(a, "-rev") == 0) {
            args.reverse_trans = 1;
            continue;
        }
        // group options
        else if (strcmp(a, "-group") == 0) var = &args.group_tag;
        
        // chrom binding
        else if (strcmp(a, "-chr-species") == 0) var = &args.chr_spec_fname;
        else if (strcmp(a, "-btag") == 0) var = &args.btag;
        
        // bed options
        else if (strcmp(a, "-bed") == 0) var = &args.bed_fname;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;

        // gtf options
        else if (strcmp(a, "-gtf") == 0) var = &args.gtf_fname;
        else if (strcmp(a, "-tags") == 0) var = &tags;

        else if (strcmp(a, "-vcf") == 0) var = &args.vcf_fname;
        else if (strcmp(a, "-vtag") == 0) var = &args.vtag;
        else if (strcmp(a, "-ref-alt") == 0) {
            args.ref_alt = 1;
            continue;
        }
        else if (strcmp(a, "-vcf-ss") == 0) {
            // args.vcf_ss = 1;
            warnings("Strand sensitive mode is enabled by default. -vcf-ss is removed since v1.3. Use -is to disable.");
            continue;
        }
        else if (strcmp(a, "-phased") == 0) {
            args.phased = 1;
            continue;
        }
        else if (strcmp(a, "-ctag") == 0) var = &args.ctag;
        
        else if (strcmp(a, "-anno-only") == 0) {
            args.anno_only = 1;
            continue;
        }

        else if (strcmp(a, "-tss") == 0) {
            args.tss_mode = 1;
            continue;
        }

        else if (strcmp(a, "-vague-edge") == 0) var = &vague_edge;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (a[0] == '-' && a[1] != '\0') error("Unknown argument, %s",a);
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }        
        error("Unknown argument: %s", a);
    }
    
    if (args.bed_fname == NULL && args.gtf_fname == NULL && args.chr_spec_fname == NULL && args.vcf_fname == NULL) 
        error("-bed or -gtf or -chr-species or -vcf must be set.");

    if (args.tss_mode == 1 && args.ctag == NULL) error("-ctag must be set if -tss enable.");
    
    if (args.output_fname == NULL) error("-o must be set.");

    if (gene_min_bp) args.gene_min_bp = gene_only_int_arg(gene_min_bp, "-gene-min-bp", 1);
    if (gene_minor_bp) args.gene_minor_bp = gene_only_int_arg(gene_minor_bp, "-gene-minor-bp", 0);
    if (gene_min_frac) args.gene_min_frac = gene_only_fraction_arg(gene_min_frac, "-gene-min-frac");
    if (gene_minor_frac) args.gene_minor_frac = gene_only_fraction_arg(gene_minor_frac, "-gene-minor-frac");
    if (!args.gene_only && (gene_min_bp || gene_minor_bp || gene_min_frac || gene_minor_frac))
        error("Gene coverage thresholds require -gene-only.");
    if (args.gene_only) {
        if (!args.gtf_fname) error("-gene-only requires -gtf.");
        if (args.anno_only) error("-gene-only preserves all records and cannot be combined with -anno-only.");
        if (tags || args.exon_level || args.flatten_flag || args.psi || args.tss_mode || vague_edge || args.antisense)
            error("-gene-only cannot be combined with -tags/-exon/-flatten/-psi/-tss/-vague-edge/-as.");
        if (args.bed_fname || args.vcf_fname || args.chr_spec_fname)
            error("-gene-only is a dedicated GTF mode; BED, VCF and species annotation must be run separately.");
    }

    if (args.exon_level == 0 && args.flatten_flag == 1) error("-flatten only used with -exon.");
    if (args.exon_level == 0 && args.psi == 1) error("-psi only used with -exon.");
    
    if (args.input_fname == NULL && !isatty(fileno(stdin))) args.input_fname = "-";
    if (args.input_fname == NULL) error("No input bam.");
    
    if (thread) args.n_thread = str2int((char*)thread);
    if (chunk) args.chunk_size = str2int((char*)chunk);
    if (map_qual) args.map_qual = str2int((char*)map_qual);
    if (args.map_qual < 0) args.map_qual = 0;

    if (vague_edge) args.vague_edge = str2int((char*)vague_edge);
    if (args.ignore_strand) {
        args.vcf_ss = 0;
    }

sam_file:
    if (args.input_sam == 1) {
        args.fp_sam = strcmp(args.input_fname, "-") ? gzopen(args.input_fname,"r") : gzdopen(fileno(stdin), "r");
        if (args.fp_sam == NULL) error("%s : %s.", args.input_fname, strerror(errno));
        args.ks = ks_init(args.fp_sam);
        
        kstring_t str = {0,0,0}; // cache first record
        args.hdr = sam_parse_header(args.ks, &str);
        if (args.hdr == NULL) error("Failed to parse header. %s", args.input_fname);
        if (str.s[0] != '@')
            args.preload_record = strndup(str.s, str.l);    
        free(str.s);
    } else {
        args.fp  = hts_open(args.input_fname, "r");
        CHECK_EMPTY(args.fp, "%s : %s.", args.input_fname, strerror(errno));
        htsFormat type = *hts_get_format(args.fp);
        if (type.format == sam) {
            warnings("Input is SAM file; enable -sam now..");
            sam_close(args.fp);
            args.fp = NULL;
            args.input_sam = 1;
            goto sam_file;
        }
        if (type.format != bam && type.format != cram)
            error("Unsupported input format, only support BAM/SAM/CRAM format.");
        args.hdr = sam_hdr_read(args.fp);
        CHECK_EMPTY(args.hdr, "Failed to open header.");
        hts_set_threads(args.fp, args.n_thread);
    }
    
    if (args.bed_fname) {
        if (args.tag) {
            memcpy(PK_tag, args.tag, 2*sizeof(char));
        }
        args.B = bed_read(args.bed_fname);
        if (args.B == 0 || args.B->n == 0) error("Bed is empty.");
    }

    if (args.vcf_fname) {
        if (args.vtag) {
            memcpy(VR_tag, args.vtag, 2*sizeof(char));
        }
        args.V = bed_read_vcf(args.vcf_fname);
        if (args.V == NULL || args.V->n == 0) error("VCF is empty.");
    }
    
    if (args.gtf_fname) {
        args.G = gtf_read_lite(args.gtf_fname);
        if (args.G == NULL) error("GTF is empty.");
        if (tags) {
            kstring_t str = {0,0,0};
            kputs(tags, &str);
            if (str.l != 17 && str.l != 11) error("-tags required 4 or 6 tag names.");
            int n;
            int *s = ksplit(&str, ',', &n);
            if (n != 4 && n != 6) error("-tags required 4 or 6 tag names.");
            
            memcpy(TX_tag, str.s+s[0], 2*sizeof(char));
            // memcpy(AT_tag, str.s+s[1], 2*sizeof(char));
            memcpy(GN_tag, str.s+s[1], 2*sizeof(char));
            memcpy(GX_tag, str.s+s[2], 2*sizeof(char));
            memcpy(RE_tag, str.s+s[3], 2*sizeof(char));
            if (args.exon_level && n == 6) {
                memcpy(EX_tag, str.s + s[4], 2*sizeof(char));
                memcpy(JC_tag, str.s + s[5], 2*sizeof(char));
            }
            
            free(str.s);
            free(s);
        }

        if (args.flatten_flag) {
            struct bed_spec *bed = gtf2bed(args.G, NULL, 3, 1, 0);
            args.flatten = bed_spec_flatten(bed, 1);
            bed_spec_destroy(bed);
            bed_build_index(args.flatten);
        }
    }

    if (args.chr_spec_fname) {
        if (args.btag) {
            memcpy(SP_tag, args.btag, 2*sizeof(char));
        }
        args.chr_binding = chr_binding(args.chr_spec_fname, args.hdr);
    }

    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        CHECK_EMPTY(args.fp_report, "%s : %s", args.report_fname, strerror(errno));
    }
    else args.fp_report =stderr;

    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");

    hts_set_threads(args.out,args.n_thread);
    // set compress level from 6 to 2, save ~1x runtime, but will also increase ~0.5x file size
    if (args.out->is_bgzf)
        args.out->fp.bgzf->compress_level = 2;
    
    args.group_stat = dict_init();
    int idx;
    idx = dict_push(args.group_stat, "__ALL__");

    dict_set_value(args.group_stat);
    
    struct read_stat *stat = malloc(sizeof(*stat));
    memset(stat, 0, sizeof(*stat));
    dict_assign_value(args.group_stat, idx, stat);
    
    return 0;
}

struct pair {
    int start; // 1 based
    int end; // 1 based
};
struct isoform {
    int n, m;
    struct pair *p;
};
struct isoform *bend_sam_isoform(bam1_t *b)
{
    struct isoform *S = malloc(sizeof(*S));
    memset(S, 0, sizeof(*S));
    int i;
    int start = b->core.pos;
    int l = 0;
    for (i = 0; i < b->core.n_cigar; ++i) {
        int cig = bam_cigar_op(bam_get_cigar(b)[i]);
        int ncig = bam_cigar_oplen(bam_get_cigar(b)[i]);
        if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
            l += ncig;
        }
        else if (cig == BAM_CDEL) {
            l += ncig;
        }
        else if (cig == BAM_CREF_SKIP) {
            if (S->n == S->m) {
                S->m = S->m == 0 ? 2 : S->m*2;
                S->p = realloc(S->p, (S->m)*sizeof(struct pair));
            }

            S->p[S->n].start = start +1; // 0 based to 1 based
            S->p[S->n].end = start + l;
            // reset block
            start = start + l + ncig;
            l = 0;
            S->n++;
        }
    }

    if (S->n == S->m) {
        S->m = S->m == 0 ? 2 : S->m*2;
        S->p = realloc(S->p, (S->m)*sizeof(struct pair));
    }
    
    S->p[S->n].start = start +1; // 0 based to 1 based
    S->p[S->n].end = start + l;
    S->n++;
    return S;    
}

struct ret_dat {
    struct bam_pool *p;
    struct dict *group_stat;
    uint64_t reads_input;
    uint64_t reads_pass_qc;
};

void gtf_anno_destroy(struct gtf_anno_type *ann)
{
    int i, j;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *a = &ann->a[i];
        for (j = 0; j < a->n; ++j) {
            struct trans_type *t = &a->a[j];
            if (t->m_exon) free(t->exon);            
            if (t->m_exclude) free(t->exl);
        }
        if (a->m) free(a->a);
        if (a->m_flatten) free(a->flatten);
    }
    free(ann->a);
    free(ann);
}
static void gtf_anno_print(struct gtf_anno_type *ann, struct gtf_spec const *G)
{
    fprintf(stderr, "Type : %s\n", exon_type_name(ann->type));
    int i;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g = &ann->a[i];
        fprintf(stderr, "Gene : %s, %s\n", dict_name(G->gene_name, g->gene_name),  exon_type_name(g->type)); 
        int j;
        for (j = 0; j < g->n; ++j)
            fprintf(stderr, "  Trans : %s, %s\n", dict_name(G->transcript_id, g->a[j].trans_id), exon_type_name(g->a[j].type));       
    }
}

static int cmpfunc(const void *_a, const void *_b)
{
    const struct trans_type *a = (const struct trans_type*) _a;
    const struct trans_type *b = (const struct trans_type*) _b;

    return (a->type > b->type) - (a->type < b->type);
}
static int cmpfunc1(const void *_a, const void *_b)
{
    const struct gene_type *a = (const struct gene_type*) _a;
    const struct gene_type *b = (const struct gene_type*) _b;

    return (a->type > b->type) - (a->type < b->type);
}

// type_splice > type_exon > type_exon_intron > type_intron_retain > type_intron > type_exclude > type_anitisense > type_antisense_intron > type_ambiguous > type_intergenic
// 
static void gene_most_likely_type(struct gene_type *g)
{
    int i;
    qsort(g->a, g->n, sizeof(struct trans_type), cmpfunc);    
    for (i = 0; i < g->n; ++i) {
        struct trans_type *t = &g->a[i];
        if (t->type == type_unknown) continue;
        else if (t->type == type_splice) {
            g->type = type_splice; // same with exon
            break;
        }
        else if (t->type == type_exon) {
            g->type = t->type;
            break; // exon is already the best hit
        }
        else if (t->type == type_exon_intron) {
            g->type = t->type;
            break;
        }
        else if (t->type == type_intron_retain) {
            g->type = t->type;
            break;
        }
        else if (t->type == type_intron) {
            g->type = t->type;
            break;
        }
        else if (t->type == type_exclude) {
            // only label transcript, but not gene
        }
        else if (t->type == type_antisense) {
            if (g->type == type_unknown) g->type = t->type;
            break;
        }
        else if (t->type == type_antisense_intron) {
            if (g->type == type_unknown) g->type = t->type;
            break;
        }
        else if (t->type == type_ambiguous) {
            if (g->type == type_unknown) g->type = t->type;
            break;
        }
        else error("Unknown gene annotation type. %d", t->type);
    }
}

static void gtf_anno_most_likely_type(struct gtf_anno_type *ann)
{
    int i;
    for (i = 0; i < ann->n; ++i) { // refresh type of each gene first
        struct gene_type *g = &ann->a[i];
        gene_most_likely_type(g);
    }

    qsort(ann->a, ann->n, sizeof(struct gene_type), cmpfunc1);
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g = &ann->a[i];
        if (g->type == type_unknown) continue;
        else if (g->type == type_exon) {
            ann->type = g->type;
            break; // exon is already the best hit
        }
        else if (g->type == type_splice) {
            ann->type = g->type; // same with exon
            break;
        }
        else if (g->type == type_exon_intron) {
            ann->type = g->type;
            break;
        }
        else if (g->type == type_intron_retain) {
            ann->type = g->type;
            break;
        }
        else if (g->type == type_intron) {
            //if (ann->type != type_exon && ann->type != type_splice && ann->type != type_exon_intron)
            ann->type = g->type;
            break;
        }
        else if (g->type == type_antisense) {
            if (ann->type == type_unknown) ann->type = g->type;
            else if (ann->type == type_antisense_intron) ann->type = g->type;
        }
        else if (g->type == type_antisense_intron) {
            if (ann->type == type_unknown) ann->type = g->type;
        }
        else if (g->type == type_ambiguous) {
            if (ann->type == type_unknown) ann->type = g->type; // better than unknown
        }
        else error("Unknown gene annotation type. %d", g->type);
    }
}

static struct gtf *query_exon_id(struct gtf const *g, int id)
{
    int i;
    int j = 0;
    for (i = 0; i < g->n_gtf; ++i) {
        struct gtf *g0 = g->gtf[i];
        if (g0->type != feature_exon) continue;
        j++;
        if (id == j) return g0;
    }
    return NULL;
}
//static enum exon_type
static struct gtf *query_exon(int start, int end, struct gtf const *G, int *exon, enum exon_type *exon_type, int vague_edge)
{
    assert(G->type == feature_transcript);
    int j = 0;
    int i;
    int e1 = 0;
    int e2 = -1;
    for (i = 0; i < G->n_gtf; ++i) {
        // from v0.4, transcript and exon in GTF_Spec struct will be sorted by coordinate
        //struct gtf_lite *g0 = G->strand == 0 ? &G->son[i] : &G->son[G->n_son-i-1];
        struct gtf *g0 = G->gtf[i];
        if (g0->type != feature_exon) continue;
        j++;
        if (start >= g0->end) continue; // check next exon
        if (start >= g0->start - vague_edge && start <= g0->end + vague_edge) {
            e1 = j;
        }
        if (end >= g0->start - vague_edge && end <= g0->end + vague_edge) {
            e2 = j;
            if (e1 == e2) {
                *exon = e1<<2| (start==g0->start)<<1 | (end == g0->end);
                *exon_type = type_exon;
                return g0;
            }
            if (e1 < e2) {
                *exon = e1<<2;//| (start==g0->start)<<1 | (end == g0->end);
                *exon_type = type_intron_retain;
                return NULL;
            }
            assert(1); // should not come here
        }
        
        if (end <= g0->start) {
            if (i==0) {
                *exon_type = type_unknown; // out of range
            }
            if (e1 == 0) {
                *exon_type = type_intron;

            } else {
                *exon = e1<<2 | (start==g0->start)<<1 | (end == g0->end);
                *exon_type = type_exon_intron;
            }
            return g0;
        }

        if (start < g0->start && end > g0->start) {
            *exon = j<<2 | (start==g0->start)<<1 | (end == g0->end);
            *exon_type = type_exon_intron;
            return g0;
        }
    }

    *exon_type = type_unknown; // out of range
    return NULL;
}
static void query_exon_inside(int start, int end, struct gtf const *G, int *ex1, int *ex2)
{
    assert(G->type == feature_transcript);
    *ex1 = -1;
    *ex2 = -1;
    int i;
    int j = 0;
    int k = 0;
    for (i = 0; i < G->n_gtf; ++i) {
        struct gtf *g0 = G->gtf[i];
        if (g0->type != feature_exon) continue;
        j++;
        if (start >= g0->end) continue; // check next exon
        if (end <= g0->start) break;
        if (start < g0->start && end > g0->end) {
            if (*ex1 == -1) *ex1 = j;
            k++;
        }
    }
    *ex2 = *ex1 + k;
}
// for each transcript, return a type of alignment record
static struct trans_type *gtf_anno_core(struct isoform *S, struct gtf const *g, int antisense, int vague)
{
    struct trans_type *tp = malloc(sizeof(*tp));
    memset(tp, 0, sizeof(*tp));
    tp->trans_id = g->transcript_id;
    // tp->type = type_unknown;

    if (args.psi && S->n >1 && antisense == 0) {
        int i;
        for (i = 1; i < S->n; ++i) {
            struct pair *p1 = &S->p[i-1];
            struct pair *p2 = &S->p[i];
            int ex1, ex2;
            query_exon_inside(p1->end, p2->start, g, &ex1, &ex2);
            if (ex1 > 0) {
                tp->type = type_exclude;
                int j;
                for (j = ex1; j < ex2; ++j) {
                    if (tp->n_exclude == tp->m_exclude) {
                        tp->m_exclude = tp->m_exclude == 0 ? 1 : tp->m_exclude *2;
                        tp->exl = realloc(tp->exl, sizeof(void*)*tp->m_exclude);
                    }
                    tp->exl[tp->n_exclude++] = query_exon_id(g, j);
                }                
            }
        }
        if (tp->type == type_exclude) return tp;
    }
    
    int exon;
    int last_exon = -1;

    // linear search
    int i;
    for (i = 0; i < S->n; ++i) {
        struct pair *p = &S->p[i];
        enum exon_type t0;
        struct gtf *e = query_exon(p->start, p->end, g, &exon, &t0, vague);

        if (t0 == type_unknown) {
            if (tp->type != type_unknown) tp->type = type_ambiguous;
            // at least some part of read cover this transcript
            if (i > 0) tp->type = type_ambiguous;
            break; // not covered this exon
        }
        // allow vague alignments when exon edges are not exactly matched
        else if (t0 == type_exon || t0 == type_exon_intron) {
            if (tp->type == type_unknown) {
                tp->type = t0;
                last_exon = exon;
                if (tp->n_exon == tp->m_exon) {
                    tp->m_exon = tp->m_exon == 0 ? 1 : tp->m_exon*2;
                    tp->exon = realloc(tp->exon, sizeof(void*)*(tp->m_exon));
                }
                tp->exon[tp->n_exon++] = e;
                continue;
            }
            else if (tp->type == type_exon || tp->type == type_splice || tp->type == type_exon_intron) {
                // || tp->type == type_exclude) {
                assert(last_exon != -1);
                int e2 = exon>>2;
                int e1 = last_exon>>2;
                
                if (e2-e1>1) {  // check the exon number
                    tp->type = type_ambiguous;
                    break;
                }

                if ((last_exon & 0x1) && (exon & 0x2)) { // check the edge
                    if (tp->type == type_exon) tp->type = type_splice;
                    last_exon = exon;
                    if (tp->n_exon == tp->m_exon) {
                        tp->m_exon = tp->m_exon == 0 ? 2 : tp->m_exon*2;
                        tp->exon = realloc(tp->exon, sizeof(void*)*(tp->m_exon));
                    }
                    tp->exon[tp->n_exon++] = e;
                    continue;
                }

                tp->type = type_ambiguous;
                break;
            }
        }
        else if (t0 == type_intron) {
            if (tp->type == type_unknown) { // intron
                tp->type = t0;
                break;
            }
            else if (tp->type == type_exon || tp->type == type_splice || tp->type == type_exon_intron) {
                // junction reads that partial mapped to exons but others mapped to intron 
                tp->type = type_ambiguous;
                break;
            }
            else {
                // should not come here
                assert(0);
            }
        }
        else if (t0 == type_intron_retain) {
            tp->type = type_intron_retain;
            break;
        }
    }
    if (antisense == 1) {
        if (tp->type == type_exon)  tp->type = type_antisense;
        else if (tp->type == type_splice)  tp->type = type_antisense;
        else if (tp->type == type_intron_retain)  tp->type = type_antisense;
        else if (tp->type == type_intron) tp->type = type_antisense_intron;
        else if (tp->type == type_exon_intron) tp->type = type_antisense_intron;
        else if (tp->type == type_ambiguous) tp->type = type_antisense; // ?
        else if (tp->type == type_exclude) tp->type = type_antisense;
        //else if (tp->type == type_unknown)
    }

    if (tp->type != type_exon && tp->type != type_splice && tp->type != type_exon_intron) {
        if (tp->n_exon > 0) {
            free(tp->exon);
            tp->n_exon = 0;
            tp->m_exon = 0;
        }
    }

    /* if (tp->type != type_exclude && tp->n_exclude > 0) { */
    /*     assert("should not come here"); */
    /*     free(tp->exl); */
    /*     tp->n_exclude = 0; */
    /*     tp->m_exclude = 0; */
    /* } */
    return tp;  
}

// add new trans node to the tree
void gtf_anno_push(struct trans_type *a, struct gtf_anno_type *ann, int gene_id, int gene_name)
{
    if (ann->n == ann->m) {
        ann->m = ann->m == 0 ? 2 : ann->m*2;
        ann->a = realloc(ann->a, ann->m*sizeof(struct gene_type));
    }
    int i;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g0 = &ann->a[i];
        if (gene_id == g0->gene_id) {
            if (g0->m == g0->n) {
                g0->m += 5;
                g0->a = realloc(g0->a, sizeof(struct trans_type)*g0->m);
            }
            memcpy(&g0->a[g0->n], a, sizeof(struct trans_type));
            g0->n++;
            return;
        }
    }

    struct gene_type *g = &ann->a[ann->n];
    ann->n++;
    memset(g, 0, sizeof(*g));
    g->gene_id = gene_id;
    g->gene_name = gene_name;
    g->type = type_unknown;
    g->m = 2;
    g->a = malloc(sizeof(struct trans_type)*g->m);
    memcpy(&g->a[g->n], a, sizeof(struct trans_type));
    g->n++;
}
// return 1 if annotate more than one gene, otherwise return 0
int gtf_anno_string(bam1_t *b, struct gtf_anno_type *ann, struct gtf_spec const *G)
{
    int ret = 0;
    // 
    if (ann->type == type_unknown) return ret;
    else if (ann->type == type_intron && args.intron_flag == 0) return ret;
    else if (ann->type == type_exon_intron && args.exonintron_flag == 0 && args.intron_flag == 0) return ret; // 
    else if (ann->type == type_ambiguous) return ret;
    else if (ann->type == type_antisense && args.antisense == 0) return ret;
    else if (ann->type == type_antisense_intron && args.antisense == 0) return ret;
    else if (ann->type == type_exclude) return ret;
    // only exon or splice come here
    kstring_t gene_name = {0,0,0};
    kstring_t gene_id   = {0,0,0};
    kstring_t trans_id  = {0,0,0};
    kstring_t tmp = {0,0,0};
    struct dict *exons = NULL;
    struct dict *juncs = NULL;
    struct dict *exl = NULL;
    struct dict *flatten = NULL;
    
    if (args.exon_level) {
        exons = dict_init();
        juncs = dict_init();
    }

    if (args.psi) {
        exl = dict_init();
    }

    if (args.flatten_flag) {
        flatten = dict_init();
    }
    
    int i;
    for (i = 0; i < ann->n; ++i) {
        struct gene_type *g = &ann->a[i];

        if (g->type == ann->type) { // from v0.12, antisense also be annotated
            char *gene = NULL;
            char *id = NULL;
            if (g->gene_name != -1) gene = dict_name(G->gene_name, g->gene_name);
            if (g->gene_id != -1) id = dict_name(G->gene_id, g->gene_id);
            
            if (gene_name.l) {
                kputc(';', &gene_name);
                kputc(';', &gene_id);
                kputc(';', &trans_id);
                ret = 1;
            }

            if (gene == NULL && id == NULL) error("No gene name or gene id in gtf? %s", (char*)b->data);
            if (gene == NULL) gene = id;
            if (id == NULL)  id = gene;
            // debug_print("gene, %s, id, %s", gene, id);
            kputs(gene, &gene_name);
            kputs(id, &gene_id);
            int j;
            int n_trans = 0;
            for (j = 0; j < g->n; ++j) {
                struct trans_type *t = &g->a[j];
                if (t->type == g->type) {                    
                    if (n_trans) kputc(',', &trans_id);
                    char *trans = dict_name(G->transcript_id, t->trans_id);
                    kputs(trans, &trans_id);
                    n_trans++;

                    if (args.exon_level) {
                        int k;
                        for (k = 0; k < t->n_exon; ++k) {
                            tmp.l = 0;
                            struct gtf *e = t->exon[k];
                            ksprintf(&tmp,"%s:%d-%d/", dict_name(args.G->name, e->seqname), e->start, e->end);
                            if (args.ignore_strand) {
                                kputs(gene, &tmp);
                            } else {
                                ksprintf(&tmp, "%c/%s", "+-"[e->strand], gene);
                            }
                            kputs("", &tmp);
                            dict_push(exons,tmp.s);
                        }          
                        
                        // junction
                        if (t->type == type_splice && t->n_exon > 1) {
                            for (k = 0; k < t->n_exon -1; ++k) {
                                tmp.l = 0;
                                struct gtf *e1 = t->exon[k];
                                struct gtf *e2 = t->exon[k+1];
                                ksprintf(&tmp,"%s:%d-%d/", dict_name(args.G->name, e1->seqname), e1->end, e2->start);
                                if (args.ignore_strand) {
                                    kputs(gene, &tmp);
                                } else {
                                    ksprintf(&tmp, "%c/%s", "+-"[e1->strand], gene);
                                }
                                kputs("", &tmp);
                                dict_push(juncs,tmp.s);
                            }
                        }
                        if (t->m_exon) free(t->exon);
                        t->n_exon = 0;
                        t->m_exon = 0;
                    }
                }

                if (args.psi && t->type == type_exclude) {
                    int k;
                    for (k = 0; k < t->n_exclude; ++k) {
                        tmp.l = 0;
                        struct gtf *e = t->exl[k];
                        ksprintf(&tmp,"%s:%d-%d/", dict_name(args.G->name, e->seqname), e->start, e->end);
                        if (args.ignore_strand) {
                            kputs(gene, &tmp);
                        } else {
                            ksprintf(&tmp, "%c/%s", "+-"[e->strand], gene);
                        }
                        kputs("", &tmp);
                        dict_push(exl,tmp.s);
                    }
                    if (t->m_exclude) free(t->exl);
                    t->m_exclude = 0;
                    t->n_exclude = 0;
                }
            }

            if (args.flatten_flag) {
                int k;
                for (k = 0; k < g->n_flatten; ++k) {
                    tmp.l = 0;
                    struct bed *bed = g->flatten[k];
                    ksprintf(&tmp,"%s:%d-%d/", bed_seqname(args.flatten, bed->seqname), bed->start+1, bed->end);
                    if (args.ignore_strand) {
                        kputs(gene, &tmp);
                    } else {
                        ksprintf(&tmp, "%c/%s", "+-"[bed->strand], gene);
                    }
                    kputs("", &tmp);
                    dict_push(flatten,tmp.s);
                }          
            }

        }
    }
    
    if (gene_name.l) {
        if (ann->type == type_antisense) {
            bam_aux_append(b, AT_tag, 'Z', gene_name.l+1, (uint8_t*)gene_name.s);
        } else {
            bam_aux_append(b, GX_tag, 'Z', gene_id.l+1, (uint8_t*)gene_id.s);
            bam_aux_append(b, GN_tag, 'Z', gene_name.l+1, (uint8_t*)gene_name.s);
            bam_aux_append(b, TX_tag, 'Z', trans_id.l+1, (uint8_t*)trans_id.s);
            if (args.exon_level && dict_size(exons)>0) {
                tmp.l = 0;
                int k;
                for (k = 0; k < dict_size(exons); ++k) {
                    if (tmp.l) kputc(',', &tmp);
                    kputs(dict_name(exons, k), &tmp);  
                }            
                bam_aux_append(b, EX_tag, 'Z', tmp.l+1, (uint8_t*)tmp.s);
                
                tmp.l = 0;
                for (k = 0; k < dict_size(juncs); ++k) {
                    if (tmp.l) kputc(',', &tmp);
                    kputs(dict_name(juncs, k), &tmp);  
                }
                if (tmp.l > 0) {
                    bam_aux_append(b, JC_tag, 'Z', tmp.l+1, (uint8_t*)tmp.s);
                }
            }
            
            if (args.psi && dict_size(exl)>0) {
                tmp.l = 0;
                int k;
                for (k = 0; k < dict_size(exl); ++k) {
                    if (tmp.l) kputc(',', &tmp);
                    kputs(dict_name(exl, k), &tmp);  
                }
                bam_aux_append(b, ER_tag, 'Z', tmp.l+1, (uint8_t*)tmp.s);
            }
            
            if (args.flatten_flag && dict_size(flatten) > 0) {
                tmp.l = 0;
                int k;
                for (k = 0; k < dict_size(flatten); ++k) {
                    if (tmp.l) kputc(',', &tmp);
                    kputs(dict_name(flatten, k), &tmp);  
                }
                //debug_print("flatten : %s", tmp.s);
                bam_aux_append(b, FL_tag, 'Z', tmp.l+1, (uint8_t*)tmp.s);
            }
        }
        free(gene_id.s);
        free(gene_name.s);
        free(trans_id.s);
    }

    if (tmp.m) free(tmp.s);
    if (args.exon_level) {
        dict_destroy(exons);
        dict_destroy(juncs);
    }

    if (args.psi) dict_destroy(exl);
    if (args.flatten_flag) dict_destroy(flatten);
    return ret;
}

/* Gene-only annotations are alignment evidence. Reads with several records
 * must be resolved by read_id before UMI grouping; these functions do not
 * classify fusions or count independent molecules. All coordinates below are
 * zero-based, half-open, irrespective of the GTF's one-based representation. */
struct gene_only_hit {
    const char *id;
    const char *name;
    const struct gtf *gene;
    int64_t bp;
};

static int gene_only_hit_compare(const void *left, const void *right)
{
    const struct gene_only_hit *a = left, *b = right;
    return strcmp(a->id, b->id);
}

static void gene_only_string(bam1_t *b, const char *tag, const char *value)
{
    if (bam_aux_append(b, tag, 'Z', strlen(value) + 1, (const uint8_t *)value) < 0)
        error("Failed to write gene-only tag %.2s for %s.", tag, bam_get_qname(b));
}

static void gene_only_clear_tags(bam1_t *b)
{
    const char *tags[] = {"GX", "GN", "TX", "RE", "AT", "EX", "JC", "FL", "ER",
                          "ZA", "ZG", "ZN", "ZB", "ZQ", "ZF", "ZC", "ZR"};
    size_t i;
    for (i = 0; i < sizeof(tags) / sizeof(tags[0]); ++i) {
        uint8_t *data;
        while ((data = bam_aux_get(b, tags[i])) != NULL)
            if (bam_aux_del(b, data) < 0) error("Failed to remove stale annotation tag.");
    }
}

static int64_t gene_only_matched_bp(const bam1_t *b)
{
    int64_t total = 0;
    uint32_t i;
    const uint32_t *cigar = bam_get_cigar(b);
    for (i = 0; i < b->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
            total += bam_cigar_oplen(cigar[i]);
    }
    return total;
}

static int64_t gene_only_overlap_region(const bam1_t *b, int64_t start, int64_t end)
{
    if (end <= start) return 0;
    int64_t ref = b->core.pos, total = 0;
    const uint32_t *cigar = bam_get_cigar(b);
    uint32_t i;
    for (i = 0; i < b->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int64_t len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int64_t lo = ref > start ? ref : start;
            int64_t hi = ref + len < end ? ref + len : end;
            if (hi > lo) total += hi - lo;
            ref += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            ref += len;
        }
    }
    return total;
}

static int64_t gene_only_overlap(const bam1_t *b, const struct gtf *gene)
{
    return gene_only_overlap_region(b, (int64_t)gene->start - 1, gene->end);
}

struct gene_only_interval {
    int64_t start, end;
};

static int gene_only_interval_compare(const void *left, const void *right)
{
    const struct gene_only_interval *a = left, *b = right;
    if (a->start != b->start) return (a->start > b->start) - (a->start < b->start);
    return (a->end > b->end) - (a->end < b->end);
}

/* Support outside the UNION of all competing gene bodies. Pairwise exclusion
 * would incorrectly accept a top gene whose support is collectively explained
 * by several short overlapping genes. Count only M/=/X within merged intervals. */
static int64_t gene_only_independent_bp(const bam1_t *b,
                                        const struct gene_only_hit *hits, int n, int best)
{
    const int64_t start = (int64_t)hits[best].gene->start - 1;
    const int64_t end = hits[best].gene->end;
    struct gene_only_interval *shared = malloc(n * sizeof(*shared));
    if (!shared) error("Out of memory calculating shared gene support.");
    int i, count = 0;
    for (i = 0; i < n; ++i) {
        if (i == best) continue;
        int64_t lo = (int64_t)hits[i].gene->start - 1, hi = hits[i].gene->end;
        if (lo < start) lo = start;
        if (hi > end) hi = end;
        if (hi > lo) shared[count++] = (struct gene_only_interval){lo, hi};
    }
    if (count > 1) qsort(shared, count, sizeof(*shared), gene_only_interval_compare);
    int merged = 0;
    for (i = 0; i < count; ++i) {
        if (merged && shared[i].start <= shared[merged - 1].end) {
            if (shared[i].end > shared[merged - 1].end) shared[merged - 1].end = shared[i].end;
        } else shared[merged++] = shared[i];
    }
    int64_t shared_bp = 0;
    for (i = 0; i < merged; ++i)
        shared_bp += gene_only_overlap_region(b, shared[i].start, shared[i].end);
    free(shared);
    return hits[best].bp - shared_bp;
}

static int gene_only_strand_matches(const bam1_t *b, const struct gtf *gene)
{
    if (args.ignore_strand) return 1;
    int reverse = !!(b->core.flag & BAM_FREVERSE);
    if (args.reverse_trans) reverse = !reverse;
    if (b->core.flag & BAM_FREAD2) reverse = !reverse;
    return gene->strand == GTF_STRAND_UNK ||
        gene->strand == (reverse ? GTF_STRAND_REV : GTF_STRAND_FWD);
}

static void bam_gene_only_anno(bam1_t *b, const struct gtf_spec *G,
                               const bam_hdr_t *h, struct read_stat *stat)
{
    gene_only_clear_tags(b);
    const int64_t total = gene_only_matched_bp(b);
    if (total > INT32_MAX) error("Aligned length exceeds ZQ int32 capacity for %s.", bam_get_qname(b));
    if (bam_aux_update_int(b, "ZQ", (int32_t)total) < 0)
        error("Failed to write ZQ for %s.", bam_get_qname(b));

    enum gene_only_status status = GO_UNASSIGNED;
    if ((b->core.flag & BAM_FUNMAP) || b->core.tid < 0 || b->core.tid >= h->n_targets) {
        status = GO_UNMAPPED;
        goto done;
    }

    /* The existing GTF index bins one-based starts. Expand the query end by
     * one base so a gene starting exactly at a bin boundary is not missed;
     * the exact half-open M overlap below removes all extra candidates. */
    int query_end = bam_endpos(b);
    if (query_end < INT_MAX) query_end++;
    struct region_itr *itr = gtf_query(G, h->target_name[b->core.tid], b->core.pos, query_end);
    struct gene_only_hit *hits = NULL;
    int n = 0, i;
    if (itr && itr->n > 0) {
        hits = malloc(itr->n * sizeof(*hits));
        if (!hits) error("Out of memory allocating gene candidates.");
        for (i = 0; i < itr->n; ++i) {
            const struct gtf *g = itr->rets[i];
            /* A display name is not a substitute for an absent gene_id. */
            if (g->gene_id < 0 || !gene_only_strand_matches(b, g)) continue;
            int64_t bp = gene_only_overlap(b, g);
            if (bp <= 0) continue;
            hits[n].id = dict_name(G->gene_id, g->gene_id);
            hits[n].name = g->gene_name < 0 ? hits[n].id : dict_name(G->gene_name, g->gene_name);
            hits[n].gene = g;
            hits[n++].bp = bp;
        }
    }
    if (itr) region_itr_destroy(itr);

    if (n > 1) qsort(hits, n, sizeof(*hits), gene_only_hit_compare);
    /* gtf_read_lite consolidates each gene. Defensive duplicate suppression
     * prevents duplicated GTF records from inflating overlap evidence. */
    int kept = 0;
    for (i = 0; i < n; ++i) {
        if (kept && strcmp(hits[kept - 1].id, hits[i].id) == 0) {
            if (hits[i].bp > hits[kept - 1].bp) hits[kept - 1] = hits[i];
        } else hits[kept++] = hits[i];
    }
    n = kept;
    int best = -1, supported = 0;
    kstring_t ids = {0,0,0}, names = {0,0,0}, counts = {0,0,0};
    for (i = 0; i < n; ++i) {
        if (i) { kputc(';', &ids); kputc(';', &names); kputc(';', &counts); }
        kputs(hits[i].id, &ids);
        kputs(hits[i].name, &names);
        ksprintf(&counts, "%lld", (long long)hits[i].bp);
        if (best < 0 || hits[i].bp > hits[best].bp) best = i;
        if (hits[i].bp >= args.gene_min_bp) supported++;
    }
    if (n) {
        gene_only_string(b, "ZG", ids.s);
        gene_only_string(b, "ZN", names.s);
        gene_only_string(b, "ZB", counts.s);
    }
    free(ids.s); free(names.s); free(counts.s);

    if (b->core.flag & BAM_FQCFAIL) status = GO_QC_FAIL;
    else if (b->core.qual < args.map_qual) status = GO_LOW_MAPQ;
    else if (n) {
        int unique = total > 0 && hits[best].bp >= args.gene_min_bp &&
                     (double)hits[best].bp / total >= args.gene_min_frac;
        int64_t independent_bp = -1;
        for (i = 0; i < n && unique; ++i) {
            if (i == best) continue;
            if (hits[i].bp > args.gene_minor_bp) {
                unique = 0;
            } else if ((double)hits[i].bp / total > args.gene_minor_frac) {
                int64_t lo = (int64_t)hits[i].gene->start - 1;
                int64_t hi = hits[i].gene->end;
                if (lo < (int64_t)hits[best].gene->start - 1) lo = (int64_t)hits[best].gene->start - 1;
                if (hi > hits[best].gene->end) hi = hits[best].gene->end;
                if (gene_only_overlap_region(b, lo, hi) != hits[i].bp) {
                    unique = 0;
                } else {
                    if (independent_bp < 0)
                        independent_bp = gene_only_independent_bp(b, hits, n, best);
                    if (independent_bp < args.gene_min_bp) unique = 0;
                }
            }
        }
        if (unique) {
            status = GO_UNIQUE;
            gene_only_string(b, "GX", hits[best].id);
            gene_only_string(b, "GN", hits[best].name);
        } else if (supported >= 2) {
            status = GO_AMBIGUOUS;
            kstring_t key = {0,0,0}, display = {0,0,0};
            kputs("AMBIG:", &key);
            int j = 0;
            for (i = 0; i < n; ++i) {
                if (hits[i].bp < args.gene_min_bp) continue;
                if (j++) { kputc('|', &key); kputc('|', &display); }
                kputs(hits[i].id, &key);
                kputs(hits[i].name, &display);
            }
            gene_only_string(b, "GX", key.s);
            gene_only_string(b, "GN", display.s);
            free(key.s); free(display.s);
        } else status = GO_INSUFFICIENT;
    }
    free(hits);

done:
    gene_only_string(b, "ZA", gene_only_status_names[status]);
    stat->gene_only_counts[status]++;
}

struct gtf_anno_type *bam_gtf_anno_core(bam1_t *b, struct gtf_spec const *G, bam_hdr_t *h, int vague_edge)
{
    bam1_core_t *c;
    c = &b->core;
    
    char *name = h->target_name[c->tid];
    int endpos = bam_endpos(b);

    if (c->tid <= -1 || c->tid > h->n_targets || (c->flag & BAM_FUNMAP)) return NULL;
    
    struct gtf_anno_type *ann = malloc(sizeof(*ann));
    memset(ann, 0, sizeof(*ann));
    ann->type = type_unknown;

    struct region_itr *itr = gtf_query(G, name, c->pos, endpos);

    // non-overlap, intergenic
    if (itr == NULL || itr->n == 0) {
        ann->type = type_intergenic;
        return ann; // no hit
    }

    // exon == splice > intron > antisense
    // see online manual for details

    struct isoform *S = bend_sam_isoform(b);
    int i;
    for (i = 0; i < itr->n; ++i) {
        int antisense = 0; // DO NOT CHANGE HERE
        
        struct gtf const *g0 = (struct gtf*)itr->rets[i];
        // check if fully enclosed in the gene region
        if (g0->start - vague_edge > c->pos+1 || endpos  > g0->end + vague_edge) continue; 

        if (args.ignore_strand == 0) {
            int bam_strand = b->core.flag & BAM_FREVERSE;
            if (args.reverse_trans) bam_strand = bam_strand ? 0 : 1;
            
            if (bam_strand) {
                if (g0->strand == GTF_STRAND_FWD) {
                    antisense = 1; 
                }
            }
            else {
                if (g0->strand == GTF_STRAND_REV) {
                    antisense = 1;
                }
            }

            if (c->flag & BAM_FREAD2) {
                antisense = antisense ? 0 : 1;
            }
        }

        // if (antisense == 1 && args.antisense == 0) continue;
        // from v0.12, also annotate detailed coverage information for antisense
        int j;
        for (j = 0; j < g0->n_gtf; ++j) {
            struct gtf const *g1 = g0->gtf[j];
            if (g1->type != feature_transcript) continue;
            struct trans_type *a = gtf_anno_core(S, g1, antisense, vague_edge);
            gtf_anno_push(a, ann, g1->gene_id, g1->gene_name);
            free(a);
        }
    }

    region_itr_destroy(itr);
    
    // stat type
    gtf_anno_most_likely_type(ann);
    
    if (args.flatten_flag &&
        (ann->type == type_exon || ann->type == type_splice || (args.exonintron_flag && ann->type == type_exon_intron))) {
        int i, j;
        for (i = 0; i < ann->n; ++i) {
            struct gene_type *a = &ann->a[i];
            if (a->type != ann->type) continue;
            for (j = 0; j < S->n; ++j) {
                struct pair *p = &S->p[j];
                struct region_itr *itr = bed_query(args.flatten, name, p->start-1, p->end, BED_STRAND_IGN);
                assert(itr);
                int k;
                for (k = 0; k < itr->n; ++k) {
                    struct bed *bed = (struct bed*)itr->rets[k];
                    // bed->start is 0 based
                    if (bed->start >= p->end || bed->end < p->start) continue; // not covered
                    char *gene = dict_name(args.flatten->name, bed->name);
                    if (a->gene_name != dict_query(args.G->gene_name, gene)) continue;

                    if (a->n_flatten == a->m_flatten) {
                        a->m_flatten = a->m_flatten == 0 ? 2 : a->m_flatten*2;
                        a->flatten = realloc(a->flatten, a->m_flatten*sizeof(void*));
                    }
                    a->flatten[a->n_flatten++] = bed;
                }
                region_itr_destroy(itr);    
            }
        }
    }

    if (ann->type == type_unknown) {
        ann->type = type_intergenic; // not fully convered
    }
    
    free(S->p); free(S);

    if (args.debug_mode) {
        fprintf(stderr, "%s   ", b->data);
        gtf_anno_print(ann, G);
    }
    
    return ann;
}
int bam_gtf_anno(bam1_t *b, struct gtf_spec const *G, struct read_stat *stat)
{
    // cleanup all exist tags
    uint8_t *data;
    if ((data = bam_aux_get(b, TX_tag)) != NULL) bam_aux_del(b, data);
    // if ((data = bam_aux_get(b, AT_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, GN_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, GX_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, RE_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, EX_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, JC_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, FL_tag)) != NULL) bam_aux_del(b, data);
    if ((data = bam_aux_get(b, ER_tag)) != NULL) bam_aux_del(b, data);

    struct gtf_anno_type *ann = bam_gtf_anno_core(b, G, args.hdr, args.vague_edge);

    bam_aux_append(b, RE_tag, 'A', 1, (uint8_t*)RE_tag_name(ann->type));

    // in default, not annotate gene name for Antisense
    int overlap = gtf_anno_string(b, ann, G);

    if (overlap) stat->reads_anno_genes++;

    if (ann->type == type_exon) stat->reads_in_exon++;
    else if (ann->type == type_splice) stat->reads_in_exon++; // reads cover two exomes
    else if (ann->type == type_intron) stat->reads_in_intron++;
    else if (ann->type == type_exon_intron) stat->reads_in_exonintron++;
    else if (ann->type == type_intron_retain) stat->reads_intron_retain++;
    else if (ann->type == type_exclude) stat->reads_exclude++; // new transcript
    else if (ann->type == type_ambiguous) stat->reads_ambiguous++; // new transcript
    else if (ann->type == type_intergenic) stat->reads_in_intergenic++;
    else if (ann->type == type_antisense) stat->reads_antisense++;
    else if (ann->type == type_antisense_intron) stat->reads_antisenseintron++;
    else error("Unknown type? %s", exon_type_name(ann->type));
    
    if (args.tss_mode == 1) {
        if ((data = bam_aux_get(b, args.ctag)) != NULL) bam_aux_del(b, data);
        if (ann->type == type_exon || ann->type == type_splice) {
            kstring_t str = {0,0,0};
            int i;
            for (i = 0; i < ann->n; ++i) {
                struct gene_type *g = &ann->a[i];
                int j;
                for (j = 0; j < g->n; ++j) {
                    struct trans_type *tx = &g->a[j];
                    if (tx->type != type_exon && tx->type != type_splice) continue;
                    struct gtf *tx_gtf = dict_query_value(G->transcript_id, tx->trans_id);
                    if (tx_gtf->strand == GTF_STRAND_FWD) {
                        if (b->core.pos+1 == tx_gtf->start) {
                            char *gene_name =  dict_name(G->gene_name, tx_gtf->gene_name);
                            if (str.l >0) kputc(';', &str);
                            ksprintf(&str, "%s:%d/+", gene_name, tx_gtf->start);
                            break;
                        }                            
                    }
                    else {
                        int endpos = bam_endpos(b);
                        if (endpos == tx_gtf->end) {
                            char *gene_name =  dict_name(G->gene_name, tx_gtf->gene_name);
                            if (str.l >0) kputc(';', &str);
                            ksprintf(&str, "%s:%d/-", gene_name, tx_gtf->end);
                            break;
                        }      
                    }
                }
            }
            if (str.l) {
                bam_aux_append(b, args.ctag, 'Z', str.l+1, (uint8_t*)str.s);
                free(str.s);
                stat->reads_tss++;
            }
        }
    }

    int ret = ann->type == type_intergenic ? 0 : 1;
    gtf_anno_destroy(ann);
    
    return ret;
}

int bam_bed_anno(bam1_t *b, struct bed_spec const *B, struct read_stat *stat)
{
    bam_hdr_t *h = args.hdr;
    
    bam1_core_t *c;
    c = &b->core;
    
    // cleanup exist tag
    uint8_t *data;
    if ((data = bam_aux_get(b, PK_tag)) != NULL) bam_aux_del(b, data);
    
    char *name = h->target_name[c->tid];
    // int endpos = bam_endpos(b);

    int i, j;
    kstring_t temp = {0,0,0};
    int read_in_peak = 0;
    int read_diff_strand = 0;
    struct dict *val = dict_init();
    struct isoform *isf = bend_sam_isoform(b);
    
    for (j = 0; j < isf->n; ++j) {
        struct pair *s = &isf->p[j];
        struct region_itr *itr = bed_query(args.B, name, s->start, s->end, BED_STRAND_IGN);
        if (itr == NULL) continue; // query failed
        if (itr->n == 0) continue; // no hit
        for (i = 0; i < itr->n; ++i) {
            struct bed *bed = (struct bed*)itr->rets[i];
            // bed->start is 0 based
            if (bed->start >= s->end || bed->end < s->start) continue; // not covered
            
            temp.l = 0;
            
            if (bed->name == -1) {          
                ksprintf(&temp, "%s:%d-%d", dict_name(B->seqname, bed->seqname), bed->start, bed->end);
                if (bed->strand == 0) kputs("/+", &temp);
                else if (bed->strand == 1) kputs("/-", &temp);
            } else {
                kputs(dict_name(B->name, bed->name), &temp);
            }
            if (bed->strand == BED_STRAND_UNK) read_in_peak = 1;
            // stat->reads_in_region++;
            else {
                int bam_strand = c->flag & BAM_FREVERSE;
                if (args.reverse_trans) bam_strand = bam_strand ? 0 : 1;
                
                //if (c->flag & BAM_FREVERSE) {
                if (bam_strand) {
                    if (bed->strand == BED_STRAND_REV) read_in_peak = 1;
                    // stat->reads_in_region++;               
                    else {
                        read_diff_strand = 1;
                        // stat->reads_in_region_diff_strand++;
                        temp.l = 0; // reset
                    }
                }
                else {
                    if (bed->strand == BED_STRAND_FWD) read_in_peak = 1;
                    // stat->reads_in_region++;               
                    else {
                        read_diff_strand = 1;
                        // stat->reads_in_region_diff_strand++;
                        temp.l = 0;
                    }
                }
            }
            
            if (temp.l) dict_push(val, temp.s);
        }
        region_itr_destroy(itr);
    }

    if (temp.m) free(temp.s);

    if (read_in_peak) stat->reads_in_region++;
    else if (read_diff_strand) stat->reads_in_region_diff_strand++;
    free(isf->p);
    free(isf);
    
    if (dict_size(val)) {
        kstring_t str = {0,0,0};
        int i;
        for (i = 0; i < dict_size(val); ++i) {
            if (i) kputc(';', &str);
            kputs(dict_name(val, i), &str);
        }

        bam_aux_append(b, PK_tag, 'Z', str.l+1, (uint8_t*)str.s);
        free(str.s);
        dict_destroy(val);
        return 1;
    }
    dict_destroy(val);
    return 0;
}

extern int bam_vcf_anno(bam1_t *b, bam_hdr_t *h, struct bed_spec const *B, const char *vtag, int ref_alt, int vcf_ss, int phased);
extern int sam_safe_check(kstring_t *str);
extern int parse_name_str(kstring_t *s);
void *run_it(void *_d)
{
    bam_hdr_t *h = args.hdr;
    struct ret_dat *dat = malloc(sizeof(struct ret_dat));
    memset(dat, 0, sizeof(*dat));

    if (args.input_sam) {
        struct kstring_pool *d = (struct kstring_pool *)_d;
        struct bam_pool *p = bam_pool_init(d->n);
        int i;
        for (i = 0; i < d->n; ++i) {
            if (!args.gene_only) parse_name_str(&d->str[i]);
            if (sam_safe_check(&d->str[i])) {
                warnings("Failed to parse %s", d->str[i].s);
                continue;
            }
            if (sam_parse1(&d->str[i], h, &p->bam[p->n])) {
                warnings ("Failed to parse SAM., %s", bam_get_qname(&p->bam[i]));
                continue;
            }
            p->n++;
        }

        kstring_pool_destroy(d);

        dat->p = p;
        
    } else {
        dat->p = (struct bam_pool*)_d;
    }
    
    dat->group_stat = dict_init();

    int idx;
    idx = dict_push(dat->group_stat, "__ALL__");

    dict_set_value(dat->group_stat);
    
    struct read_stat *stat = malloc(sizeof(*stat));
    memset(stat, 0, sizeof(*stat));
    dict_assign_value(dat->group_stat, idx, stat);
    
    int i;
    
    for (i = 0; i < dat->p->n; ++i) {
        int ann = 0;
        bam1_t *b = &dat->p->bam[i];

        if (args.gene_only) stat = dict_query_value(dat->group_stat, 0);
    
        if (args.group_tag) {
            uint8_t *data;
            if ((data = bam_aux_get(b, args.group_tag)) != NULL) {
                char *tag_val = (char*)(data+1);
                int idx = dict_query(dat->group_stat, tag_val);
                if (idx == -1) {
                    idx = dict_push(dat->group_stat, tag_val);
                    struct read_stat *s = malloc(sizeof(*s));
                    memset(s, 0, sizeof(*s));
                    dict_assign_value(dat->group_stat, idx, s);
                }
                stat = dict_query_value(dat->group_stat, idx);
            }
        }

        bam1_core_t *c;
        c = &b->core;
        if (args.gene_only) {
            dat->reads_input++;
            if (!(c->flag & BAM_FUNMAP) && c->tid >= 0 && c->tid < h->n_targets &&
                c->qual >= args.map_qual && !(c->flag & BAM_FQCFAIL))
                dat->reads_pass_qc++;
            bam_gene_only_anno(b, args.G, h, stat);
            continue;
        }
        // move the mapping quality check above of unmap and secondary alignment check. because both of these two jump to check -anno-only
        if (c->qual < args.map_qual) { // if -q set, filter low quailty directly
            b->core.flag |= BAM_FQCFAIL;
            continue;
        } //goto check_continue;

        // secondary alignment
        if (c->flag & BAM_FSECONDARY) goto check_continue;
        
        dat->reads_input++;
        
        // QC
        if (c->tid <= -1 || c->tid > h->n_targets || (c->flag & BAM_FUNMAP)) goto check_continue;

        dat->reads_pass_qc++;

        if (args.G) 
            if (bam_gtf_anno(b, args.G, stat)) ann = 1;

        if (args.B)
            if (bam_bed_anno(b, args.B, stat)) ann = 1;

        if (args.V)
            if (bam_vcf_anno(b, args.hdr, args.V, VR_tag, args.ref_alt, args.vcf_ss, args.phased)) ann = 1;

        if (args.chr_binding) {
            char *v = args.chr_binding[b->core.tid];
            if (v != NULL) {
                bam_aux_append(b, SP_tag, 'Z', strlen(v)+1, (uint8_t*)v);
                ann=1;
            }
        }

      check_continue:
        if (args.anno_only && ann == 0) { // if only export annotated reads, intergenic reads will be filter
            b->core.flag |= BAM_FQCFAIL;
        } 
    }
    return dat;
}

static void write_out(void *_d)
{
    struct ret_dat *dat = (struct ret_dat *)_d;
    int i;
   
    for (i = 0; i < dat->p->n; ++i) {
        if (!args.gene_only && (dat->p->bam[i].core.flag & BAM_FQCFAIL)) continue; // gene-only preserves records
        if (sam_write1(args.out, args.hdr, &dat->p->bam[i]) == -1)
            error("Failed to write SAM.");
    }
   
    args.reads_input   += dat->reads_input;
    args.reads_pass_qc += dat->reads_pass_qc;

    for (i = 0; i < dict_size(dat->group_stat); ++i) {
        int idx = dict_query(args.group_stat, dict_name(dat->group_stat, i));
        if (idx == -1) {
            idx = dict_push(args.group_stat, dict_name(dat->group_stat, i));
            struct read_stat *s = malloc(sizeof(*s));
            memset(s, 0, sizeof(*s));
            dict_assign_value(args.group_stat, idx, s);
        }

        struct read_stat *s0 = dict_query_value(args.group_stat, idx);
        struct read_stat *s1 = dict_query_value(dat->group_stat, i);
        s0->reads_in_region += s1->reads_in_region;
        s0->reads_in_region_diff_strand += s1->reads_in_region_diff_strand;
        s0->reads_in_intergenic += s1->reads_in_intergenic;
        s0->reads_in_exon += s1->reads_in_exon;
        s0->reads_in_intron += s1->reads_in_intron;
        s0->reads_intron_retain += s1->reads_intron_retain;
        s0->reads_antisense += s1->reads_antisense;
        s0->reads_antisenseintron += s1->reads_antisenseintron;
        s0->reads_ambiguous += s1->reads_ambiguous;
        s0->reads_in_exonintron += s1->reads_in_exonintron;
        s0->reads_tss += s1->reads_tss;
        s0->reads_anno_genes += s1->reads_anno_genes;
        int j;
        for (j = 0; j < GO_STATUS_COUNT; ++j)
            s0->gene_only_counts[j] += s1->gene_only_counts[j];
    }
    bam_pool_destory(dat->p);

    // free assign memory manually
    for (i = 0; i < dict_size(dat->group_stat); ++i) {
        void *v = dict_query_value(dat->group_stat, i);
        if (v) free(v);
    }

    dict_destroy(dat->group_stat);
    free(dat);
}
void write_report()
{
    if (args.gene_only) {
        uint64_t counts[GO_STATUS_COUNT] = {0};
        int i, j;
        for (i = 0; i < dict_size(args.group_stat); ++i) {
            const struct read_stat *s = dict_query_value(args.group_stat, i);
            for (j = 0; j < GO_STATUS_COUNT; ++j) counts[j] += s->gene_only_counts[j];
        }
        fprintf(args.fp_report, "Metric,Value\n");
        fprintf(args.fp_report, "Mode,gene_only_alignment_evidence\n");
        fprintf(args.fp_report, "Counting_unit,alignment_record_not_read_or_UMI\n");
        fprintf(args.fp_report, "Alignment_records,%llu\n", (unsigned long long)args.reads_input);
        fprintf(args.fp_report, "Mapped_alignment_records_passing_MAPQ_and_QC,%llu\n",
                (unsigned long long)args.reads_pass_qc);
        for (j = 0; j < GO_STATUS_COUNT; ++j)
            fprintf(args.fp_report, "%s,%llu\n", gene_only_status_names[j], (unsigned long long)counts[j]);
        return;
    }
    if (dict_size(args.group_stat) == 1) {
        struct read_stat *s0 = (struct read_stat*)dict_query_value(args.group_stat, 0);
        fprintf(args.fp_report, "Reads Mapped to Genome (Map Quality >= %d),%.1f%%\n", args.map_qual, (float)args.reads_pass_qc/args.reads_input*100);
        
        if (args.B) {
            fprintf(args.fp_report, "Reads Mapped to BED regions / Peaks,%.1f%%\n", (float)s0->reads_in_region/args.reads_pass_qc*100);
            if (s0->reads_in_region_diff_strand)
                fprintf(args.fp_report, "Reads Mapped on different strand of BED regions,%.1f%%\n", (float)s0->reads_in_region_diff_strand/args.reads_pass_qc*100);
        }
        if (args.G) {
            fprintf(args.fp_report, "Reads Mapped to Exonic Regions,%.1f%%\n", (float)s0->reads_in_exon/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to Intronic Regions,%.1f%%\n", (float)s0->reads_in_intron/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to both Exonic and Intronic Regions,%.1f%%\n", (float)s0->reads_in_exonintron/args.reads_pass_qc*100);
            if (s0->reads_antisenseintron > 0) {
                fprintf(args.fp_report, "Reads Mapped Antisense to Exon,%.1f%%\n", (float)s0->reads_antisense/args.reads_pass_qc*100);
                fprintf(args.fp_report, "Reads Mapped Antisense to Intron,%.1f%%\n", (float)s0->reads_antisenseintron/args.reads_pass_qc*100);
            } else {
                fprintf(args.fp_report, "Reads Mapped Antisense to Gene,%.1f%%\n", (float)s0->reads_antisense/args.reads_pass_qc*100);
            }
            
            fprintf(args.fp_report, "Reads Mapped to Intergenic Regions,%.1f%%\n", (float)s0->reads_in_intergenic/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to Gene but Failed to Interpret Type,%.1f%%\n", (float)s0->reads_ambiguous/args.reads_pass_qc*100);
            fprintf(args.fp_report, "Reads Mapped to Overlapped genes,%.1f%%\n", (float)s0->reads_anno_genes/args.reads_pass_qc*100);
            if (s0->reads_tss>0) {
                fprintf(args.fp_report, "Reads map start from TSS,%.1f%%\n", (float)s0->reads_tss/args.reads_pass_qc*100);
            }
        }
    }
    else {
    }
}
static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    if (args.fp) sam_close(args.fp);
    if (args.fp_sam) {
        gzclose(args.fp_sam);
        ks_destroy(args.ks);
    }
    sam_close(args.out);
    int i;
    for (i = 0; i < dict_size(args.group_stat); ++i) {
        void *v = dict_query_value(args.group_stat, i);
        if (v) free(v);
    }
    dict_destroy(args.group_stat);
    if (args.B) bed_spec_destroy(args.B);
    if (args.G) gtf_destroy(args.G);
    if (args.V) bed_spec_var_destroy(args.V);
    if (args.flatten_flag) bed_spec_destroy(args.flatten);
    if (args.fp_report != stderr) fclose(args.fp_report);
}

extern int anno_usage();
static void *read_chunk()
{
    if (args.input_sam) {
        struct kstring_pool *b =  kstring_pool_read(args.ks, args.chunk_size);
        if (b == NULL) return NULL;
        if (b->n == 0) {
            free(b->str);
            free(b);
            return NULL;
        }
        return b;
    }
    
    struct bam_pool *b = bam_pool_create();
    bam_read_pool((struct bam_pool*)b, args.fp, args.hdr, args.chunk_size);
    if (b == NULL) return NULL;
    if (b->n == 0) {
        free(b->bam);
        free(b);
        return NULL;
    }
    return b;
}
void process_bam()
{
    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        void *b0 = read_chunk();
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b0, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
                write_out(d);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }
    
    hts_tpool_process_flush(q);
    
    while ((r = hts_tpool_next_result(q))) {
        struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
}

void process_bam1()
{
#ifdef _OPENMP
#pragma omp parallel num_threads(args.n_thread)
#endif
    for (;;) {
        void *b = NULL;
#ifdef _OPENMP
#pragma omp critical (read)
#endif
        b = read_chunk();
        if (b == NULL) break;

        b = run_it(b);

#ifdef _OPENMP
#pragma omp critical (write)
#endif
        write_out(b);
    }
}
int bam_anno_attr(int argc, char *argv[])
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return anno_usage();
    
    process_bam1();

    write_report();
    memory_release();    

    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Speed : %d records/sec; Peak RSS: %.3f GB.",
              realtime() - t_real, cputime(), (int)(args.reads_pass_qc/cputime()), peakrss() / 1024.0 / 1024.0 / 1024.0);
    
    return 0;    
}
