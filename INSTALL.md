# Installing TACT

TACT is a Python-based tool that relies on several key scientific packages, including [click](https://click.palletsprojects.com), [DendroPy](https://jeetsukumaran.github.io/DendroPy/), and [NumPy](https://numpy.org), which are automatically installed when you install TACT.

For faster performance, especially with large datasets, using the [PyPy 3 implementation of Python](https://www.pypy.org/) is recommended (when available), as it can significantly speed up analyses.

!!! tip "Which installation method should I use?"
    - **Docker** (recommended): Easiest to install and fastest for large datasets. Works on Windows, Mac, and Linux.
    - **pipx**: Good alternative if you can't use Docker. Also supports PyPy for better performance.
    - **Other methods**: Not recommended or supported.

## Docker (recommended)

If you can use Docker, this is the **recommended installation method**. It's both convenient to install and [fast for large datasets](troubleshooting.md#why-is-tact-so-slow) thanks to PyPy.

Install [Docker Desktop](https://www.docker.com/products/docker-desktop) and run the following to download the TACT image:

```sh
docker pull jonchang/tact:latest
```

Then, run TACT from the container image, giving it access to your current working directory:

```bash
mkdir -p examples
cd examples
curl -LO https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.csv
curl -LO https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.tre
docker run -it -v "$(pwd)":/workdir -w /workdir jonchang/tact tact_build_taxonomic_tree Carangaria.csv --output Carangaria.taxonomy.tre
docker run -it -v "$(pwd)":/workdir -w /workdir jonchang/tact tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre --output Carangaria.tacted
```

Here's a [screencast](https://asciinema.org/a/347571) of how to use the Docker commands:

<script id="asciicast-347571" src="https://asciinema.org/a/347571.js" async></script>

The above Docker image defaults to the latest tagged release. In the rare case you need to use a different version, a full list of tags is available on [Docker Hub](https://hub.docker.com/r/jonchang/tact/tags).

## pipx

[pipx](https://pipx.pypa.io) is a tool that installs Python applications in isolated environments. This is a good alternative if you can't use Docker.

### Basic Installation

First, [install `pipx`](https://pipx.pypa.io/stable/installation/), then run:

```sh
pipx install tact
```

This installs TACT with the standard Python interpreter (CPython). This works fine for small to medium datasets.

### Faster Installation with PyPy (Optional)

If you have PyPy3 installed, you can install a faster version of TACT:

```sh
pipx install --python pypy3 tact
```

!!! warning "PyPy installation issues"
    - This will take much longer to install than the standard version
    - Installation may fail if the proper dependencies (mainly openblas) aren't set up
    - On macOS, you'll need to run: `brew install openblas gcc pypy3 pipx`
    - You may also need to set environment variables:
        - `PKG_CONFIG_PATH`: should point to where `openblas.pc` lives (e.g., `/opt/homebrew/opt/openblas/lib/pkgconfig`)
        - `MACOSX_DEPLOYMENT_TARGET`: should be your macOS version (e.g., `11.0`)

## Verifying Your Installation

After installing, verify that TACT is working correctly:

```sh
tact_add_taxa --version
```

You should see the version number. If you get a "command not found" error, make sure TACT is in your PATH.

!!! tip "Next steps"
    Once TACT is installed, check out the **[Tutorial](tutorial.md)** to learn how to use it, or read the **[Background](background.md)** to understand what TACT does.
