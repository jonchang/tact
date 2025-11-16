# Contributing

## Bug reports and feature requests

Open a bug report or a feature request via [GitHub Issues](https://github.com/jonchang/tact/issues/new/choose).

## Developing

Development uses [`poetry`](https://python-poetry.org/). Simply clone the repository and install:

```console
$ git clone https://github.com/jonchang/tact.git
$ cd tact
$ poetry install
```

**Important:** Always use `poetry run` when executing Python code or commands in this repository to ensure the correct virtual environment and dependencies are used:

```console
$ poetry run python script.py
$ poetry run pytest
$ poetry run tact_add_taxa --version
```

**Python version:** TACT requires Python 3.10 or higher (for PyPy compatibility).

## Pull requests

If you spot a bug, opening [a pull request](https://guides.github.com/activities/hello-world/) to fix it would be appreciated!

On the web version of this documentation, click the edit button (looks like a pencil) to open the text editor to propose changes.

## Testing

Tests use [pytest](https://docs.pytest.org/). Mostly these are integration tests but there are some unit tests as well.

```console
$ poetry run pytest  # optionally with --script-launch-mode=subprocess
```

Tests are automatically run on all supported Python versions via GitHub Actions CI.

## Releasing

### Pre-Release Tasks

1. **Verify CI passes**: [GitHub Actions](https://github.com/jonchang/tact/actions)

2. **Update outdated packages** if needed. Always check for the presence of relevant compiled wheels, especially for NumPy with PyPy.
   ```bash
   poetry show --outdated
   poetry update [package]
   ```

3. **Update version in `pyproject.toml`**. For details on versioning, see [Semantic Versioning](https://semver.org/).
   ```bash
   poetry version patch # (or minor, major)
   ```

4. **Review all commits since the last version** to ensure nothing is missing.
   ```bash
   git log --oneline $(git describe --tags --abbrev=0)..HEAD
   ```

5. **Update NEWS.md**:
    - Document all user-facing changes, relevant dependency updates, and breaking changes.
    - Note: `docs/news.md` is a symlink to `NEWS.md`, so updating `NEWS.md` automatically updates the docs.

6. **Ensure working directory is clean**:
   ```bash
   git status
   ```

7. **Commit version and NEWS.md updates**:
   ```bash
   git add pyproject.toml NEWS.md
   git commit -m "tact $(poetry version -s)"
   ```

### Releasing

1. **Create and push tag**:
   ```bash
   VERSION=$(poetry version -s)
   git tag v$VERSION
   git push --atomic origin master v$VERSION
   ```

**Note:** The tag format must be `v*.*.*` (e.g., `v0.7.0`) for GitHub Actions to automatically trigger the release workflows.

### Post-Release Verification

After pushing the tag:

1. **GitHub Actions CI/CD**: [Check workflows succeed](https://github.com/jonchang/tact/actions)

2. **PyPI release**: [Check package appears](https://pypi.org/project/tact/) and version numbers are correct

3. **Docker Hub release**: [Check image appears](https://hub.docker.com/r/jonchang/tact) and tags are correct (e.g., `0.7.0`, `0.7`, `0`, `latest`)

4. **GitHub Packages release**: [Check image appears](https://github.com/jonchang/tact/pkgs/container/tact) and tags are correct (e.g., `0.7.0`, `0.7`, `0`, `latest`)

5. **Update GitHub release notes**: [Create a release](https://github.com/jonchang/tact/releases/new) and include highlights from NEWS.md

6. **Announce new version**: Announce on Twitter or other relevant channels

**Important notes:**

- GitHub Actions workflows automatically build and publish the new version to PyPI, and build container images for Docker Hub and GitHub Packages when a tag matching `v*.*.*` is pushed.
- PyPI publishing uses trusted publishing aka OIDC (no manual upload needed).
- Docker images are built using PyPy 3.10, which also serves as the PyPy compatibility test.
- All tests run automatically in CI before the release is published.
