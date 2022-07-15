FROM ubuntu:22.04

LABEL maintainer="rdemko2332@gmail.com"

WORKDIR /usr/bin/

RUN apt-get -qq update --fix-missing && \
    apt-get install -y \
    wget \
    perl \
    bowtie2=2.4.4-1 \
    bowtie=1.3.1-1 \
    samtools=1.13-4 \
    bedtools=2.30.0+dfsg-2
  
RUN chmod +x *

WORKDIR /work
