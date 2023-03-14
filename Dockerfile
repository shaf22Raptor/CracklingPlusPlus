FROM ubuntu:jammy
WORKDIR /crackling
RUN apt-get update && apt-get install -y \
    cmake \
    build-essential \
    libboost-all-dev \
    vienna-rna \ 
    bowtie2 \
    && rm -rf /var/lib/apt/lists/*

