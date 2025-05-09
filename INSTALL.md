# Installing

TACT requires Python 3. When possible, I recommend using the [PyPy 3 implementation](https://www.pypy.org/) as it can significantly speed up TACT analyses, particularly on large datasets. In addition, TACT depends on the [click](https://click.palletsprojects.com), [DendroPy](https://dendropy.org), [NumPy](https://numpy.org), and [SciPy](https://www.scipy.org/) packages.

## Docker

If you can use Docker, this is the recommended method as it is both convenient to install and [fast for large datasets](https://tact.jonathanchang.org/troubleshooting/#why-is-tact-so-slow) thanks to PyPy. Note that Docker images are only provided for x86_64 processors (`amd64`). Docker images are not automatically built for ARM processors and must either be built manually, or consider installing TACT with one of the other options below.

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

[Install `pipx`](https://pipxproject.github.io/pipx/installation/), then run:

```sh
pipx install tact
```

If you have PyPy3 installed, you can try to install a faster version using:

```sh
pipx install --python pypy3 tact
```

Note that this will take much longer to install and could fail if the proper dependencies (mainly openblas) aren't set up.

On macOS, you'll need to run `brew install openblas gcc pypy3 pipx`. You may also need to set certain environment variables:

* `PKG_CONFIG_PATH`: this should be where `openblas.pc` lives (e.g., `/opt/homebrew/opt/openblas/lib/pkgconfig`)
* `MACOSX_DEPLOYMENT_TARGET`: this should be your macOS version (e.g., `11.0`)

## Other

Other ways of installing TACT, including unpacking the tarball somewhere or directly using `pip`, are neither supported nor recommended.
# 
