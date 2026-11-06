# Docker 使用说明

PISA-long 容器提供 gene-only v1（`v1.8-gene-only.1`），以及 Python 3.12、pysam 和 samtools。`PISA`、`samtools`、`python3`、`pisa_resolve_gene_reads.py` 均在 PATH 中。没有固定的 ENTRYPOINT，镜像名后可以直接指定这些命令。

镜像名称为 [`dawang02/pisa-long`](https://hub.docker.com/r/dawang02/pisa-long)，平台为 **linux/amd64**。优先用 `v1.8-gene-only.1` 版本标签固定分析环境；`latest` 是随发布移动的标签。需要精确复现某次镜像内容时，另外记录所用镜像的 digest。容器不包含参考 GTF、样本 BAM 或 UMI-tools，也不执行 UMI 去重。

## 拉取与准备文件

```bash
docker pull --platform linux/amd64 dawang02/pisa-long:v1.8-gene-only.1
PISA_IMAGE=dawang02/pisa-long:v1.8-gene-only.1
```

切换到自己的数据目录，其中包含 `input.bam` 和 `genes.gtf`，或在下列命令中替换为实际文件名。`-v "$PWD":/data` 将该宿主机目录映射到容器 `/data`；容器中的文件参数必须使用挂载后的路径。BAM 与 GTF 需使用同一参考基因组和染色体命名，注释输入不需要预先排序或建立索引。

`--user "$(id -u):$(id -g)"` 让输出使用当前宿主机用户的 UID/GID。当前目录需要对该用户可写。示例的 `--platform linux/amd64` 指定镜像架构；非 amd64 主机需要其容器运行环境支持该架构。

## 注释 BAM

```bash
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" -v "$PWD":/data -w /data \
  "$PISA_IMAGE" \
  PISA anno -gene-only -ignore-strand -t 4 \
  -gtf /data/genes.gtf -o /data/gene_only.bam \
  -report /data/gene_only.csv /data/input.bam
```

`gene_only.bam` 和 `gene_only.csv` 会出现在宿主机当前目录。`-ignore-strand` 是示例选择，应根据自己的文库方向信息决定是否使用。gene-only 已按整个基因区域统计支持，不需要额外添加 `-intron`。

当前 PISA 的 `-t` 同时设置注释工作线程和 HTS 输入、输出线程，默认值为 4，并不是整个进程的线程总数上限。`PISA anno -@` 目前只解析参数而未使用其值，因此示例不使用它。下面的 `samtools -@` 是 samtools 自己的有效参数。

## 排序和索引

多线程注释输出不保证坐标有序，即使继承的 header 仍声明有序。若要建立索引或在 IGV 中查看，先排序：

```bash
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" -v "$PWD":/data -w /data \
  "$PISA_IMAGE" \
  samtools sort -@ 4 -T /data/pisa_sort_tmp \
  -o /data/gene_only.sorted.bam /data/gene_only.bam

docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" -v "$PWD":/data -w /data \
  "$PISA_IMAGE" samtools index /data/gene_only.sorted.bam
```

`-T` 指定排序临时文件前缀，请为临时文件和结果保留足够的磁盘空间。若后续直接使用下面的 resolver，可以跳过此处排序与索引，避免重复处理。

## 可选：汇总同一 read 的全部比对

```bash
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" -v "$PWD":/data -w /data \
  "$PISA_IMAGE" pisa_resolve_gene_reads.py \
  --input /data/gene_only.bam --gtf /data/genes.gtf \
  --output-dir /data/read_gene_results --threads 4
```

输出目录必须不存在或为空。resolver 内部执行 collate，再将全部输出 BAM 坐标排序并建立索引。`--threads` 仅用于 collate、排序和索引，Python 的逐 read 判定仍是串行的。完整输出为 `read_gene_results/annotated.sorted.bam`；`umi_input.primary.sorted.bam` 每个合格 read 只保留一条 primary，**尚未执行去重**。

这里会根据完整比对集合重新判定 read 级 GX。应使用完整输入 BAM；SA 简化 CIGAR 仍可能使已有片段被误判缺失，详见 [规则与已知限制](../scripts/README.pisa_gene_only.md)。

## 从源码构建

在 PISA-long 仓库根目录执行：

```bash
docker build --platform linux/amd64 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  -t dawang02/pisa-long:v1.8-gene-only.1 \
  -t dawang02/pisa-long:latest .
```

`VCS_REF` 将当前 Git 提交写入镜像来源标签；它表示源码提交，工作区另有未提交修改时还应单独记录这些修改。构建命令只生成本地镜像，不会上传 Docker Hub。

可以读取镜像记录的源码提交：

```bash
docker image inspect dawang02/pisa-long:v1.8-gene-only.1 \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}'
```
