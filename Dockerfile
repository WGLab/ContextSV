# Build and run:
# docker build -t clion/ubuntu/cpp-env:1.0 -f Dockerfile .
# docker run -it clion/ubuntu/cpp-env:1.0 /bin/bash

FROM ubuntu:20.04

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get -y install tzdata

RUN apt-get update \
  && apt-get install -y build-essential \
      gcc \
      g++ \
      gdb \
      clang \
      make \
      ninja-build \
      cmake \
      autoconf \
      automake \
      libtool \
      valgrind \
      locales-all \
      dos2unix \
      rsync \
      tar \
      python \
      python-dev \
      wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev \
  && apt-get clean

# HTSLib
WORKDIR /deps
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2

WORKDIR /deps/htslib-1.3.2
RUN ./configure && \
    make && \
    make install

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/deps/htslib-1.3.2/

WORKDIR /
