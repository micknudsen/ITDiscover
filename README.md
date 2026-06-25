# ITDiscover

[![Anaconda-Server Badge](https://anaconda.org/MOMA-AUH/itdiscover/badges/version.svg)](https://anaconda.org/MOMA-AUH/itdiscover) [![Anaconda-Server Badge](https://anaconda.org/MOMA-AUH/itdiscover/badges/downloads.svg)](https://anaconda.org/MOMA-AUH/itdiscover)

Tool for discovering FLT3 ITDs from amplicon sequencing of AML samples

## Installation

The recommended way to install **ITDiscover** is via [conda](https://docs.conda.io/), using the `MOMA-AUH` channel:

```bash
conda install MOMA-AUH::itdiscover
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
  --output report.html
```
