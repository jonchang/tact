FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

COPY . /tact

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    g++ \
    gcc \
    gfortran \
    liblapack-dev \
    liblapack3 \
    libopenblas-pthread-dev \
    libopenblas0-pthread \
    libopenblas64-0-pthread \
    libopenblas64-pthread-dev \
    locales \
    locales-all \
    pypy3 \
    pypy3-dev \
    wget \
  && wget -nv https://bootstrap.pypa.io/get-pip.py \
  && pypy3.8 get-pip.py \
  && pypy3.8 -mpip install ./tact \
  && rm -rf tact get-pip.py \
  && apt-get remove -y \
    g++ \
    gcc \
    gfortran \
    liblapack-dev \
    libopenblas-pthread-dev \
    libopenblas64-pthread-dev \
    pypy3-dev \
    wget \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /root/.cache \
  && rm -rf /var/cache/* \
  && rm -rf /var/log/*
