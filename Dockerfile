ARG PYTHON_VERSION=3.15

# Stage 1: build wheel and pin dependencies to poetry.lock versions
FROM python:${PYTHON_VERSION}-slim AS build

ENV PIP_NO_CACHE_DIR=1 \
  PIP_DISABLE_PIP_VERSION_CHECK=1 \
  PIP_ROOT_USER_ACTION=ignore

RUN pip install poetry poetry-plugin-export

COPY . /src
WORKDIR /src
RUN poetry build -f wheel -o /dist \
  && poetry export -f requirements.txt -o /dist/requirements.txt --only main

# Stage 2: install from wheels only to avoid compiling numpy/scipy
FROM python:${PYTHON_VERSION}-slim

ENV PIP_NO_CACHE_DIR=1 \
  PIP_DISABLE_PIP_VERSION_CHECK=1 \
  PIP_ROOT_USER_ACTION=ignore \
  # avoid thread pool penalty from small linalg operations
  OPENBLAS_NUM_THREADS=1

RUN --mount=type=bind,from=build,source=/dist,target=/dist \
  pip install --only-binary=:all: --require-hashes -r /dist/requirements.txt \
  && pip install --no-deps /dist/*.whl

CMD ["bash"]
