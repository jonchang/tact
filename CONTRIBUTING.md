# Contributing

## Bug reports and feature requests

Open a bug report or a feature request via [GitHub Issues](https://github.com/jonchang/tact/issues/new/choose).

## Developing

Development uses [`poetry`](https://poetry.eustace.io/). Simply clone the repository and install:

```console
$ git clone https://github.com/jonchang/tact.git
$ cd tact
$ poetry install
```

## Pull requests

If you spot a bug, opening [a pull request](https://guides.github.com/activities/hello-world/) to fix it would be appreciated!

On the web version of this documentation, click the edit button (looks like a pencil) to open the text editor to propose changes.

## Testing

Tests use [pytest](https://docs.pytest.org/). Mostly these are integration tests but there are some unit tests as well.

```console
$ poetry run pytest  # optionally with --script-launch-mode=subprocess
```

## Releasing

To release a new version, increment the version using Poetry, then tag the new commit and push to GitHub:

```console
$ poetry version patch  # or minor, etc.
$ git commit -p
$ git tag VERSION  # (v0.4.0)
$ git push --atomic origin BRANCH_NAME TAG  # (master, v0.4.0)
```

A GitHub Actions workflow will build and publish the new version on PyPI, and build container images for Docker Hub and GitHub Packages.
