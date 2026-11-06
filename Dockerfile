ARG PYTHON_IMAGE=python:3.12-slim-bookworm@sha256:935c3c5b12357a7e10caa9c1e6ca5692344d97cf229832ef584be2cec9fc03df

# Keep the same Python and system libraries in the build and runtime stages.
FROM ${PYTHON_IMAGE} AS runtime-base
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates libbz2-1.0 libcurl4 libgomp1 liblzma5 samtools \
    && rm -rf /var/lib/apt/lists/*
COPY requirements-test.txt /tmp/pisa-requirements.txt
RUN python -m pip install --no-cache-dir --disable-pip-version-check \
        --only-binary=:all: -r /tmp/pisa-requirements.txt \
    && rm /tmp/pisa-requirements.txt

FROM runtime-base AS builder
ARG BUILD_JOBS=4
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential libbz2-dev libcurl4-openssl-dev liblzma-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /build/pisa-long
COPY Makefile pisa_version.h LICENSE README.md ./
COPY src/ ./src/
COPY third_party/ ./third_party/
COPY scripts/pisa_resolve_gene_reads.py scripts/compare_pisa_gene_annotations.py ./scripts/
COPY scripts/tests/test_pisa_gene_only.py ./scripts/tests/
COPY scripts/README.pisa_gene_only.md ./scripts/
COPY docs/ ./docs/
COPY config/ ./config/

# The runtime image depends on this tested stage: a failed test stops the build.
RUN make -j "${BUILD_JOBS}" \
    && make test PYTHON=python \
    && install -D -m 755 PISA /install/bin/PISA \
    && install -m 755 scripts/pisa_resolve_gene_reads.py \
        scripts/compare_pisa_gene_annotations.py /install/bin/ \
    && mkdir -p /install/share/pisa-long/licenses /install/share/pisa-long/scripts \
    && cp LICENSE /install/share/pisa-long/licenses/PISA-LICENSE \
    && cp LICENSE README.md /install/share/pisa-long/ \
    && cp scripts/README.pisa_gene_only.md /install/share/pisa-long/scripts/ \
    && cp -r docs config /install/share/pisa-long/

# Preserve notices from the bundled static libraries without shipping their
# implementation or the synthetic fixtures in the runtime image.
RUN python -c 'from pathlib import Path; import re; notices = {}; \
    [(notices.setdefault(block, []).append(str(path))) \
     for path in sorted(Path("third_party").rglob("*")) \
     if path.suffix in (".c", ".h") \
     for block in re.findall(r"/\*.*?\*/", path.read_text(errors="replace"), flags=re.S) \
     if re.search(r"copyright|permission|redistribution|license", block, flags=re.I)]; \
    Path("/install/share/pisa-long/licenses/THIRD_PARTY_NOTICES.txt").write_text( \
        "\n\n".join("Sources: " + ", ".join(paths) + "\n" + block \
                    for block, paths in notices.items()) + "\n")'

FROM runtime-base AS runtime
ARG VCS_REF=unknown
ARG BUILD_DATE
LABEL org.opencontainers.image.title="PISA-long" \
      org.opencontainers.image.description="Long-read gene-body annotation and read-level gene reconciliation; gene-only v1." \
      org.opencontainers.image.source="https://github.com/dawangran/PISA-long" \
      org.opencontainers.image.url="https://github.com/dawangran/PISA-long" \
      org.opencontainers.image.version="v1.8-gene-only.1" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}" \
      org.opencontainers.image.licenses="MIT"
COPY --from=builder /install/ /usr/local/
WORKDIR /data

# No fixed entrypoint: PISA, samtools, Python and both helper scripts are on PATH.
# PISA's no-argument help intentionally retains its upstream exit status of 1.
CMD ["PISA"]
