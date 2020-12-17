FROM python:3.8

# install fortran (required for dropkick python package)
RUN apt-get update; apt-get install -y gfortran=4:8.3.0-1

# install conda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -O ~/anaconda.sh
RUN bash ~/anaconda.sh -b -p /usr/local/anaconda
RUN /bin/bash -c "source /usr/local/anaconda/bin/activate"
ENV PATH="/usr/local/anaconda/bin:${PATH}"
RUN conda init bash

# install python dependencies via conda
COPY qc_pipe_env.yml environment.yaml
RUN conda env create -f environment.yaml

# install OS packages required for singularity
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config

# install go
RUN wget https://dl.google.com/go/go1.15.6.linux-amd64.tar.gz
RUN tar -C /usr/local -xzf go1.15.6.linux-amd64.tar.gz
ENV PATH="/usr/local/go/bin:${PATH}"
RUN go version

# install Singularity
ARG SINGULARITY_VERSION
RUN mkdir -p $(go env GOPATH)/src/github.com/sylabs && \
    cd $(go env GOPATH)/src/github.com/sylabs && \
    wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-${SINGULARITY_VERSION}.tar.gz && \
    tar -xzf singularity-${SINGULARITY_VERSION}.tar.gz && \
    cd ./singularity && \
    ./mconfig
RUN cd $(go env GOPATH)/src/github.com/sylabs/singularity/builddir && \
    make && \
    make install
