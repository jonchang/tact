FROM ubuntu:21.04
ARG DEBIAN_FRONTEND=noninteractive

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

COPY . /tact

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    g++ \
    gcc \
    gfortran \
    liblapack-dev \
    liblapack3 \
    libopenblas-dev \
    libopenblas0 \
    locales \
    pypy3 \
    pypy3-dev \
  && localedef -i en_US -f UTF-8 en_US.UTF-8 \
  && curl -L https://bootstrap.pypa.io/get-pip.py | pypy3 - \
  && pypy3 -mpip install ./tact --verbose \
  && rm -rf tact \
  && apt-get remove -y \
    g++ \
    gcc \
    gfortran \
    liblapack-dev \
    libopenblas-dev \
    pypy3-dev \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /root/.cache \
  && rm -rf /var/cache/* \
  && rm -rf /var/log/*
