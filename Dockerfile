FROM ubuntu:23.04
ARG DEBIAN_FRONTEND=noninteractive

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH=/root/.local/bin:$PATH

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
    pkg-config \
    python3 \
    python3-venv \
    pypy3 \
    pypy3-dev \
    wget \
  && wget -nv https://bootstrap.pypa.io/get-pip.py \
  && wget -nv -O get-poetry.py https://install.python-poetry.org \
  && python3 get-pip.py \
  && python3 get-poetry.py \
  && cd tact \
  && poetry env use $(which pypy3.9) \
  && poetry install --only main \
  && ln -s $(poetry env info --path)/bin/tact_* /root/.local/bin \
  && cd .. \
  && rm -rf get-pip.py get-poetry.py \
  && apt-get remove -y \
    g++ \
    gcc \
    gfortran \
    liblapack-dev \
    libopenblas-pthread-dev \
    libopenblas64-pthread-dev \
    pkg-config \
    python3 \
    python3-venv \
    pypy3-dev \
    wget \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /var/cache/* \
  && rm -rf /var/log/*
