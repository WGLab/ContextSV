# Use the miniconda container
FROM continuumio/miniconda3:main

# Version argument (set during build)
ARG CONTEXTSV_VERSION

WORKDIR /app

RUN apt-get update
RUN conda update conda

# Install ContextSV
RUN conda config --add channels wglab
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda create -n contextsv python=3.9
RUN echo "conda activate contextsv" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
RUN conda install -n contextsv -c wglab -c conda-forge -c jannessp -c bioconda contextsv=${CONTEXTSV_VERSION} && conda clean -afy

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "contextsv", "contextsv"]
