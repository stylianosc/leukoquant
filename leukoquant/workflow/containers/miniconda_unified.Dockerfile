FROM continuumio/miniconda3:25.3.1-1

# Build context must be the workflow/ directory:
#   docker buildx build --platform linux/amd64 \
#     -f containers/miniconda_unified.Dockerfile -t ... .

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends bc time libgomp1 gcc && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY envs/miniconda_unified_env.yaml /tmp/miniconda_unified_env.yaml
RUN conda env update --name base --file /tmp/miniconda_unified_env.yaml && \
    conda clean --all --force-pkgs-dirs -y && \
    pip cache purge && \
    rm -f /tmp/miniconda_unified_env.yaml && \
    find /opt/conda -type f -name '*.pyc' -delete && \
    find /opt/conda -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true && \
    rm -rf /tmp/* /var/tmp/* /root/.cache 2>/dev/null || true

CMD [ "/bin/bash" ]
