# Installing TACT

TACT is a Python-based tool that relies on several key scientific packages, including [click](https://click.palletsprojects.com), [DendroPy](https://jeetsukumaran.github.io/DendroPy/), and [NumPy](https://numpy.org), which are automatically installed when you install TACT.

!!! tip "Which installation method should I use?"
    - **pipx** (recommended): Fastest, and works on Windows, Mac, and Linux with Python 3.12 or newer.
    - **Docker**: Use this if you can't install Python, or want a fixed environment. It is [a little slower](troubleshooting.md#how-long-will-tact-take) than pipx.
    - **Other methods**: Not recommended or supported.

## pipx (recommended)

[pipx](https://pipx.pypa.io) is a tool that installs Python applications in isolated environments.

First, [install `pipx`](https://pipx.pypa.io/stable/installation/).

Then run:

```sh
pipx install tact
```

## Docker

If you can't install Python, you can run TACT from its Docker image instead.

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

## Verifying your installation

After installing, verify that TACT is working correctly:

```sh
tact_add_taxa --version
```

You should see the version number. If you get a "command not found" error, make sure TACT is in your PATH.

!!! tip "Next steps"
    Once TACT is installed, check out the **[Tutorial](tutorial.md)** to learn how to use it, or read the **[Background](background.md)** to understand what TACT does.
