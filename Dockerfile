# Use the miniconda container
FROM continuumio/miniconda3:main

# Version argument (set during build)
ARG CONTEXTSV_VERSION

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates && rm -rf /var/lib/apt/lists/*
RUN conda update -y conda

# Install ContextSV and plotting dependencies.
RUN conda config --add channels wglab \
	&& conda config --add channels conda-forge \
	&& conda config --add channels bioconda \
	&& conda create -y -n contextsv python=3.10 \
	&& conda install -y -n contextsv -c wglab -c conda-forge -c bioconda \
	   contextsv=${CONTEXTSV_VERSION} plotly python-kaleido \
	&& conda clean -afy

# Smoke test both commands at build time.
RUN conda run -n contextsv contextsv --help \
	&& conda run -n contextsv contextsv-cnv-plot --help

SHELL ["/bin/bash", "--login", "-c"]

# Default command remains contextsv, but this allows overriding with contextsv-cnv-plot.
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "contextsv"]
CMD ["contextsv"]
