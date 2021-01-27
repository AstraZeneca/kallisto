<div align="center">
<img src="./assets/logo.svg" alt="Kallisto" width="300">
</div>

##

[![Tests](https://github.com/AstraZeneca/kallisto/workflows/Tests/badge.svg)](https://github.com/AstraZeneca/kallisto/actions?workflow=Tests)
[![Codecov](https://codecov.io/gh/AstraZeneca/kallisto/branch/master/graph/badge.svg)](https://codecov.io/gh/AstraZeneca/kallisto)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)

Table of Contents
-----------------

- Full Author List
- Introduction
- Installation
- Reference

Full Author List
----------------

Eike Caldeweyher and Philipp Pracht

Introduction
------------

We developed the `kallisto` program for the efficient and robust calculation of atomic features using molecular geometries either in a ``xmol`` or a ``Turbomole`` format.
Furthermore, several modelling tools are implemented, e.g., to calculate root-mean squared deviations via quaternions (including rotation matrices), sorting of molecular geometries and many more. All features of ``kallisto`` are described in detail within our [documentation](https://app.gitbook.com/@ehjc/s/kallisto/).

Installation
------------

`kallisto` runs on `python3`

Python development setup. Install the `pyenv` python version manager:
```bash
curl https://pyenv.run | bash
```
and add this to the `~/.bashrc` and source it:
```bash
export PATH="~/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
```
Install the latest python versions:
```bash
pyenv install 3.8.2
pyenv install 3.7.7
pyenv local 3.8.2 3.7.7
```

Now we are ready to set up `kallisto`.
Clone the repository:
```bash
git clone git@github.com:f3rmion/kallisto.git
```

Install a python dependency manager. We choose to go with `poetry`:
```bash
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
source ~/.poetry/env
```
or alternatively via `pip`:
```bash
pip install --user poetry
```

Now, if you haven't already done so, change into the cloned `kallisto` directory and
download the dependencies via `poetry`:
```bash
cd kallisto
poetry install
```

Finally install the test automation environment `nox` via `pip``:
```bash
pip install --user --upgrade nox
```

Run `nox` to test the setup.

Reference
---------

tba
