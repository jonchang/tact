# Contributing

Development on TACT uses [`poetry`](https://poetry.eustace.io/). Simply clone the repository and install:

```console
$ git clone https://github.com/jonchang/tact.git
$ cd tact
$ poetry install
```

## Releasing

When releasing a new version of TACT, run its tests and bump its revision like so:

```console
$ poetry run pytest  # optionally with --script-launch-mode=subprocess
$ poetry version patch  # or minor, etc.
$ git commit -p
$ git tag VERSION  # (v0.4.0)
$ git push --atomic origin BRANCH_NAME TAG  # (master, v0.4.0)
```

A GitHub Actions workflow will build and publish the new version on PyPI, as well as releasing container images to Docker Hub and GitHub Packages.
