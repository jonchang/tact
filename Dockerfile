ARG version=20.04
FROM ubuntu:$version
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    gcc \
    gfortran \
    liblapack3 \
    liblapack-dev \
    libopenblas0 \
    libopenblas-dev \
    locales \
    pypy3 \
    pypy3-dev \
  && rm -rf /var/lib/apt/lists/* \
  && localedef -i en_US -f UTF-8 en_US.UTF-8 \
  && curl -L https://bootstrap.pypa.io/get-pip.py | pypy3 -

COPY . /tact

RUN pypy3 -mpip install ./tact --verbose

RUN apt-get update \
  && apt-get remove -y \
    pypy3-dev \
    liblapack-dev \
    libopenblas-dev \
    gcc \
    gfortran \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*
