FROM debian:trixie-slim
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
    libopenblas-dev \
    libopenblas0 \
    locales \
    meson \
    pkg-config \
    pypy3 \
    pypy3-dev \
    pypy3-venv \
    wget \
  && wget -nv https://bootstrap.pypa.io/get-pip.py \
  && wget -nv -O get-poetry.py https://install.python-poetry.org \
  && rm /usr/lib/pypy3.10/EXTERNALLY-MANAGED \
  && pypy3.10 get-pip.py \
  && pypy3.10 get-poetry.py \
  && cd tact \
  && poetry self add poetry-plugin-export \
  && poetry export -f requirements.txt -o requirements.txt --without-hashes --only main \
  && pypy3.10 -mpip install -r requirements.txt --compile . \
  && cd .. \
  && pypy3.10 get-poetry.py --uninstall \
  && rm -rf get-pip.py get-poetry.py \
  && apt-get remove -y \
    g++ \
    gcc \
    gfortran \
    liblapack-dev \
    libopenblas-dev \
    meson \
    pkg-config \
    pypy3-dev \
    pypy3-venv \
    wget \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /var/cache/* \
  && rm -rf /var/log/* \
  && rm -rf ~/.cache \
  && rm -rf ~/.local
