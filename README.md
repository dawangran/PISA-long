# PISA-long

面向长读 RNA 比对的基因注释工具，基于 [PISA](https://github.com/shiquan/PISA) 扩展。根据实际比对片段给 BAM 写入基因 ID、候选集合和支持证据，供后续 read 级汇总及 UMI 分组使用。

当前采用 **gene-only v1**，二进制版本为 `v1.8-gene-only.1`。规则阈值可调，尚未经过全面生物学准确率验证。本工具不执行 UMI 纠错或去重，也不把跨基因标签当作已经确认的融合。

## 安装

需要 C 编译器和 GNU Make。仓库包含构建所需的 HTSlib、zlib 源码；运行 Python 汇总脚本及测试还需要 Python 3.10+、`pysam` 和 PATH 中的 `samtools`。

```bash
git clone https://github.com/dawangran/PISA-long.git
cd PISA-long
make -j 8
python3 -m pip install -r requirements-test.txt
```

构建后可执行文件为仓库根目录的 `./PISA`，项目名称不改变原 PISA 命令名。`samtools` 需通过自己的系统或 Conda 环境安装。上游命令仍保留；原始使用文档见 [PISA User's Guide](https://shiquan.github.io/PISA.html)。

## 工作流程

`./PISA anno -gene-only` 对每条 alignment 计算实际比对碱基与 gene body 的重叠，使用可调的 100 bp / 80% 等规则写入初步 GX。可选的 `scripts/pisa_resolve_gene_reads.py` 再汇总同名 primary、secondary、supplementary，重新判定 read 级标签并准备 primary-only 工作 BAM。

不加 `-gene-only` 时仍运行传统 PISA 注释。C 注释和 Python read 汇总处理不同层级，后者会重算 GX；两者的状态及计数不能混用。

## 快速运行

以下命令从仓库根目录执行，替换 `genes.gtf` 和 `input.bam` 为自己的完整输入。GTF 与 BAM 应使用相同参考基因组和染色体命名。注释输入 BAM 不要求预先排序，也不需要索引。

```bash
./PISA anno -gene-only -ignore-strand \
  -gtf genes.gtf -t 4 \
  -o gene_only.bam -report gene_only.csv input.bam

samtools sort -@ 4 -o gene_only.sorted.bam gene_only.bam
samtools index gene_only.sorted.bam
```

**当前 `PISA anno` 的 `-t`（默认 4）同时设置注释工作线程和 HTS 输入、输出线程，不是整个进程的线程总数上限。`-@` 虽然被解析，但当前实现没有使用其值，不能靠它增加读写线程。** 上面的 `samtools -@` 则是 samtools 自己的有效参数。gene-only 已包括内含子范围，不需要另加 `-intron`。示例的 `-ignore-strand` 表示忽略链方向；应根据文库方向信息决定是否使用，不能只根据 read 的 FLAG 推断原始 RNA 来源。

只统计 CIGAR 的 `M/=/X` 实际比对碱基；`N/D` 跨越的参考区域不算覆盖。新剪接位点本身不会导致拒绝注释。输出保留输入比对记录、CB/UR 及非注释标签，但会替换旧基因和结构注释标签。

多线程注释输出不保证坐标有序，即使继承的 BAM header 仍标为有序，也要先 `samtools sort` 再建立索引或用 IGV 查看。

如需在完整 BAM 上汇总同一 read 的全部比对，可直接使用未排序的 `gene_only.bam`，跳过上面的排序和索引，避免重复处理：

```bash
python3 scripts/pisa_resolve_gene_reads.py \
  --input gene_only.bam --gtf genes.gtf \
  --output-dir read_gene_results --threads 4
```

resolver 内部先按名称 collate，再将所有输出 BAM 坐标排序并建立索引。`--threads` 用于 collate、排序和索引，Python 的逐 read 判定仍串行执行。输出包括完整 `annotated.sorted.bam`、每个合格 read 一条 primary 的 `umi_input.primary.sorted.bam`、分类 BAM、逐 read TSV 和 JSON 汇总；均未执行去重。输出目录必须不存在或为空。详细规则与输出见 [使用说明](scripts/README.pisa_gene_only.md)。

## GX 与后续 UMI 分组

- `GX` 保存 gene_id 或排序后的组合键，例如 `AMBIG:geneA|geneB`；`GN` 保存可读的基因名。
- `ZG/ZN/ZB/ZQ` 保留 alignment 层的原始候选和实际支持量；PISA 写入初步 `ZA`，read resolver 再写入统一的 `ZA/ZR` 状态与原因。
- `AMBIG` 表示仍有多个候选。read resolver 的 `FUSION` 也只表示候选组合，不是已验证融合或断点。
- 后续工具按 GX 的完整字符串分组，不会自动把 `geneA` 与 `AMBIG:geneA|geneB` 合并，也不会自动根据自定义 ZA 排除低置信结果。
- 不同文库应独立处理；同基因、同 CB、同 UR 也不保证剪接结构相同。基因注释不等于长读结构层面的去重。

## 已知限制

1. **固定长度门槛**可能保留较短的次要命中，也可能忽略真实但较短的证据；默认阈值是工程起始值，需要根据数据评估。
2. **resolver 的 SA 完整性检查仍偏严格**：它规范化软/硬剪切后比较 CIGAR；某些比对器在 SA 中使用简化 CIGAR，即使对应记录存在也可能被标为 `incomplete_alignment_set`。应复核这类结果，不能据此认定比对丢失。
3. C 注释是 alignment 级，不判断融合，也不会整合其他比对记录；区域 BAM 切片可能缺少 primary 或 supplementary，应在完整 BAM 上汇总 read 级证据。
4. 标签覆盖率及合成测试通过不等于基因归属准确率。候选支持不足、质量不满足等情况下，v1 可以不写 GX；自定义 ZA 也不会自动保护后续去重。

仓库提供可复现的合成测试，不包含研究样本 BAM、参考 GTF 或私有分析结果。

## 测试

```bash
make test
```

合成回归测试覆盖 gene-only v1 和 read resolver，包括记录保留、CIGAR、共享/嵌套基因、替代及补充比对。测试通过表示实现符合当前规则，不表示这些规则已通过生物学准确率验证。

## 来源与许可

保留 PISA 作者的 MIT 许可及第三方源码中的许可声明，详见 [LICENSE](LICENSE) 和 [上游来源](docs/UPSTREAM.md)。本项目是在 PISA 基础上的独立扩展，新增长读 gene-only 模式、read 汇总与相关测试。

PISA 原论文：[Bioinformatics, btac562](https://doi.org/10.1093/bioinformatics/btac562)。
