# PISA 长读基因注释与 read 级汇总

这个流程为长读 BAM 生成基因／候选组合标签，并准备每条 read 只有一条记录的 UMI-tools 输入。**它目前只做注释，不运行 UMI 纠错、不删除重复 reads、不生成已经去重的 BAM，也不是经过生物学验证的融合检测器。**

当前版本为 `v1.8-gene-only.1`，使用 gene-only v1 规则。参数是可调整的工程起始值，尚未证明对所有长读数据均最优。

实现分成两层：

1. 修改后的 `PISA anno -gene-only` 对每条 alignment 计算基因区域的实际比对碱基支持。它包含 primary、secondary、supplementary 的证据。
2. `scripts/pisa_resolve_gene_reads.py` 汇总同一文库、同一 QNAME 的完整比对集合，给出 read 级标签，并输出完整 BAM 与 primary-only 工作 BAM。

Python 脚本需要 Python 3.10+、`pysam` 和 `samtools`。它用普通 `samtools collate` 将同名记录聚到一起，逐 read 处理；没有使用会过滤记录的 `collate -f`，也不会在内存中保存全 BAM。磁盘上需要容纳中间 BAM、最终 BAM 及分类子集。`--threads` 用于 collate、输出排序和索引；Python 逐 read 的归属判定仍串行执行。

## 运行

先构建修改后的 PISA，使用与原分析一致的完整 GTF。注释输入 BAM 不要求排序或索引；下面将注释结果直接交给 resolver，避免提前排序：

```bash
make -j 8

./PISA anno -gene-only -ignore-strand \
  -gtf genes.gtf \
  -t 4 \
  -o alignment_gene_only.bam \
  -report alignment_gene_only.report.csv \
  input.bam

python3 scripts/pisa_resolve_gene_reads.py \
  --input alignment_gene_only.bam \
  --gtf genes.gtf \
  --output-dir read_gene_results \
  --threads 4
```

**当前 `PISA anno` 的 `-t`（默认 4）同时设置注释工作线程和 HTS 输入、输出线程，不是整个进程的线程总数上限。`-@` 被解析但未用于设置线程，所以当前无效；`samtools -@` 属于 samtools 自己的参数，仍然有效。** `-intron` 在 gene-only 模式没有必要：这个模式使用整个 gene body，包括内含子。`-ignore-strand` 是示例选择；应根据文库方向信息决定，不等于文库一定没有链信息。

多线程 PISA 注释输出不保证坐标有序，即使继承的 header 仍标为有序。若直接建立索引或用 IGV 查看，先执行 `samtools sort`；若交给 resolver，无需预先排序或索引，它会 collate 输入并将所有输出 BAM 坐标排序、建立索引。

输出目录必须不存在或为空，脚本拒绝覆盖已有结果。可用 `--temp-dir /path/to/scratch` 指定临时文件位置。

## Alignment 层的证据

只累计 CIGAR 的 `M/=/X`，`D/N` 只移动参考坐标，`I/S/H` 不提供基因覆盖。使用 gene feature 的区间，不要求整个 alignment 被 gene 完整包含，不检查剪接位点是否属于已知 transcript。

| 标签 | 含义 |
|---|---|
| `ZG:Z` | 有正重叠的 gene_id，以分号分隔 |
| `ZN:Z` | 与 ZG 对齐的 gene_name |
| `ZB:Z` | 与 ZG 对齐的实际支持 bp，以分号分隔 |
| `ZQ:i` | 该 alignment 的 `M/=/X` 总长度 |
| `ZA:Z` | PISA 的初步状态，汇总脚本会替换为 read 级状态 |

PISA 的 gene-only 模式清理旧的基因、转录本、区域注释以及本流程旧标签，避免过时注释残留。它保留输入记录、比对字段、序列、质量、CB、UR 和其他非注释标签，包括带 QC-fail、secondary、supplementary 的记录。传统注释模式则仍沿用原逻辑，另应用C类型分支修复。

当前基因支持的默认判定参数为：

| PISA 参数 | 汇总脚本参数 | 默认值 |
|---|---|---:|
| `-gene-min-bp` | `--min-gene-bp` | 100 bp |
| `-gene-min-frac` | `--min-gene-fraction` | 0.8 |
| `-gene-minor-bp` | `--minor-bp` | 100 bp |
| `-gene-minor-frac` | `--minor-fraction` | 0.2 |

这些是可调整的工程起始值，不是已验证的最优生物学阈值。两层参数应保持一致。PISA 总是输出所有正重叠基因的证据，因此调整汇总阈值时不必改变比对或截断候选基因列表。

在 C 注释层，对通过当前质量检查的 alignment，主基因达到最小支持和比例门槛、且其他候选满足次要支持约束时，写入单基因 GX。否则，若至少两个基因各达到 `-gene-min-bp`，写入排序后的 `AMBIG` GX；该分支不额外要求所有候选的支持并集达到80%。其他支持不足的情况不写 GX。C 层不汇总同名记录，也不生成融合判定。

下面的 Python 汇总另加完整比对集合、query 区间并集、竞争 secondary 和可靠 anchor 检查。因此即使阈值相同，alignment 与 read 级标签也不保证一致；不能把 C 层的 AMBIG 数量当作最终 read 数量。

## Read 层的归属规则

同名记录按 header 中 RG 对应的 LB 区分文库。没有 LB 时按 RG 区分；没有 RG 时视为同一输入文库。不同文库的数据不应在没有可区分的 RG/LB 标识时合并。

每组必须有唯一 primary。主要比对链包含 primary 和不带 secondary 标志的 supplementary。先统一原始 read 的 query 坐标：包括软／硬剪切所占长度，并对反向比对做坐标转换。重叠 query 区间使用并集，不能重复累加。

归属按以下顺序判定：

- 主比对链的实际覆盖对每个基因形成支持区间。共享区支持不能作为两个基因各自的独立支持；每个基因的独立支持还要扣除其他候选在同一 query 区间的支持。
- 至少两个基因各有 `--min-anchor-bp`（默认 100 bp）的可靠独立支持，基因支持区间并集覆盖至少 `--min-gene-fraction` 的全部实际比对区间，且没有显著的未解决第三候选、竞争 secondary 时，才标为 `fusion_candidate`。可靠 anchor 要求 MAPQ ≥ `--min-mapq`（默认 20），且不是 QC-fail；MAPQ 255 表示不可用，不当作高质量。

  这里的独立支持是多个 `M/=/X` 区间的累计长度，不要求连续 100 bp，也没有验证融合断点。80% 的分母是实际比对区间并集，不是包括软剪切的整条 read 长度；整条 read 的比对比例另见 TSV 的 `matched_query_fraction`。
- 唯一主导基因需支持 ≥100 bp、占全部实际比对区间 ≥80%，每个其他候选通常同时满足 ≤100 bp 和 ≤20%。对共享边界设一个例外：其他候选仍须 ≤100 bp；如果其全部支持区间包含在主基因支持区间内，并且主基因扣除**所有其他候选支持的并集**后仍有 ≥100 bp 独立支持，则可超过 20%。例如基因A支持230 bp、基因B的66 bp全部来自共享区，A还剩164 bp独立支持，可以归属A。这个例外使用实际比对支持，不要求融合 anchor 的高 MAPQ 条件；primary MAPQ 的整体准入规则仍然生效。
- 有多个显著候选但无法区分来源时，保存 `AMBIG` 候选集合。基因支持并集覆盖不足时保存 `insufficient`。
- secondary 用作替代解释，不作为另一条 read，也不加入主要比对链的支持总量。与某个主要比对片段的 query 覆盖相互达到默认 80%，且 AS（两边都有时）达到对方的默认 95%，才视为有竞争力。缺少 AS 时采取保守策略。现版本逐条比较替代片段，没有对完整替代 chimeric chain 做联合打分。
- 一个基因内的新剪接、截短或末端越界不会单独导致拒绝。多个候选集合之间不做传递合并，也不使用 UMI 或同组其他 reads 自动改写基因归属。

缺 primary、多 primary、CB/UR 冲突、主比对链总 query 长度矛盾、SA 指向缺失片段、低／不可用 primary MAPQ 会旁路 UMI 输入。SA 完整性比较忽略软／硬剪切的表示差异，但仍要求规范化后的CIGAR相同。**已知限制：某些比对器在SA中写简化CIGAR，即使对应记录存在，也可能被误标为 `incomplete_alignment_set`；本版本尚未修正，应复核这些旁路结果。** 没有SA也不能证明某个区域切片已经包含完整比对集合。

## Read resolver 的最终标签

| 状态 | GX | 其他标签 |
|---|---|---|
| `unique` | 单个 gene_id | `ZA=unique` |
| `ambiguous` | `AMBIG:A\|B` | `ZA=ambiguous` |
| `fusion_candidate` | `FUSION:A\|B` | `ZA=fusion_candidate`，`ZF=fusion_candidate` |
| 不足、未注释或不完整等旁路状态 | 不写 GX | `ZA`、`ZR` 解释原因 |

`GN` 保存对应的显示名称。`ZC` 保存最终 read 级候选 gene_id 列表，`ZR` 保存处理原因。同一个文库、同一个 read 的所有记录获得一致的 read 级标签。原 alignment 的 `ZG/ZN/ZB/ZQ` 保持不变，支持逐条复核。

`AMBIG`、`FUSION` 的 ID 使用排序后的候选基因集合。**FUSION 目前是无序的基因组合 ID，没有编码方向或断点，不应当作精确融合事件 ID。** 可能的 read-through、文库嵌合、错误比对都需要后续分析。

## 输出文件

resolver 输出的所有 BAM 均坐标排序并带 `.bai`：

| 文件 | 内容 |
|---|---|
| `annotated.sorted.bam` | 所有原始 alignment records，新增最终 read 级注释 |
| `umi_input.primary.sorted.bam` | 每个合格 read 一条 primary；没有执行去重 |
| `fusion_candidates.sorted.bam` | 融合候选 read 在输入中存在的全部记录 |
| `ambiguous.sorted.bam` | 歧义 read 在输入中存在的全部记录 |
| `bypass.sorted.bam` | 未进入 UMI 输入的 read 的全部记录 |
| `read_assignments.tsv.gz` | 每个文库、read 一行的候选、支持 bp、状态、CB/UR 和比对记录数 |
| `summary.json` | read 与 alignment 分开的统计、参数、文库统计和运行耗时 |
| `libraries.tsv` | 各文库的 read 组数、合格 UMI 输入 read 数 |
| `collate.log` | 外部 collate 的日志 |

缺少 primary CB/UR 或 primary 带 QC-fail 时，基因注释可以保留，但该 read 不进入 UMI 输入，因此分类 BAM 与 bypass BAM 可能有交集。需要保留记录时以完整 annotated BAM 和 read_id 为准，不能把所有分类 BAM 不加判断地拼接。

若 `summary.json` 的 `umi_input_requires_library_split=true`，**必须按文库拆分 primary 工作 BAM 后再运行 UMI-tools**。聚合工作 BAM 并不表示这些文库可以一起去重，UMI-tools 的 `--per-cell` 不会自动包含 LB。

## 与后续去重的边界

普通单文库的 primary 工作 BAM 可以作为以下分组参数的输入：

```text
--per-gene --gene-tag=GX --per-cell --cell-tag=CB
--extract-umi-method=tag --umi-tag=UR
```

但本程序没有运行 UMI-tools。后续仍需进行 UMI 分组、结构兼容性检查和代表 read 选择，并从完整 BAM 按所选 read_id 取回其全部比对记录。

- UMI-tools 按标签完整字符串分组，不会自动把 `A` 与 `AMBIG:A|B` 合并或消歧。
- 普通 BAM `--per-gene` 模式在染色体切换时结束分组。同一融合／歧义组合的 primary 若落在不同染色体，需要额外处理跨染色体分组；仅换 tag 不会解决。
- 同一基因组合不等于同一断点或剪接结构。直接每组删除到一条可能损失后续 IsoQuant 的模型证据。
- 区域切片中缺失的 primary 或 SA partner 无法凭空恢复，应对完整 BAM 运行。脚本只标记缺失并保留记录，不能将切片结果解释为全量样本的融合或 UMI 统计。

## 回归检查

```bash
PISA_TEST_BINARY="$PWD/PISA" \
  python3 -m unittest discover -s scripts/tests -p 'test_pisa_gene_only*.py' -v
```

合成用例覆盖 gene-only v1 和 read resolver：基因边界、新剪接、CIGAR的M/=/X与D/N区别、共享与嵌套区域、替代比对、正反链及硬剪切分段、输入记录保留、缺失比对、标签冲突，以及primary-only输出。测试是规则实现验证，不代表生物学准确率验证。
