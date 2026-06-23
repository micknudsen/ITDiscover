# ITDiscover

Tool for discovering FLT3 ITDs from amplicon sequencing of AML samples

## Installation

```bash
python -m pip install -e .[dev]
```

For conda-based development:

```bash
conda env create -f environment.yml
conda activate itdiscover-dev
```

## Usage

```bash
itdiscover --help
itdiscover --version
```

Example:

```bash
itdiscover \
  --reference reference.fasta \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --unique-support-report-out unique-support.html
```

## Tests

```bash
pytest
```

## Conda Packaging

Build the conda package locally with:

```bash
conda build conda-recipe/
```

The repository includes a tag-triggered GitHub workflow in `.github/workflows/publish.yml`.
Pushing a tag like `v0.1.0` will:

1. build the conda package,
2. smoke test it,
3. upload it to the `micknudsen` Anaconda channel.

The workflow expects a repository secret named `ANACONDA_API_TOKEN`.

To bump the project version consistently across package metadata:

```bash
python tools/set_version.py 0.1.1
```
