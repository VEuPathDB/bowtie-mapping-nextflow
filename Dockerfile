FROM ubuntu:22.04

LABEL maintainer="rdemko2332@gmail.com"

WORKDIR /usr/bin/

RUN apt-get update && \
    apt-get install -y \
    wget \
    perl \
    make zip unzip gcc  build-essential \
    bowtie2=2.4.4-1 \
    samtools=1.13-4 \
    bedtools=2.30.0+dfsg-2 \
  && rm -rf /var/lib/apt/lists/*

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz \
    && tar -xvzf sratoolkit.3.0.0-ubuntu64.tar.gz



RUN mkdir -p /root/.ncbi
RUN printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg
RUN printf '/libs/cloud/report_instance_identity = "true"\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/libs/cloud/accept_aws_charges = "false"\n/libs/cloud/accept_gcp_charges = "false"\n' >> /root/.ncbi/user-settings.mkfg

RUN chmod +x *

RUN cp /root/.ncbi/user-settings.mkfg /usr/bin/user-settings.mkfg

WORKDIR /homer
RUN wget http://homer.ucsd.edu/homer/configureHomer.pl
RUN perl configureHomer.pl -install

ENV PATH=$PATH:/usr/bin/sratoolkit.3.0.0-ubuntu64/bin:/homer/bin

WORKDIR /work
