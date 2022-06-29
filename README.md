[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/koesterlab/alignoth/CI)](https://github.com/koesterlab/alignoth/actions)
[![codecov](https://codecov.io/gh/koesterlab/alignoth/branch/main/graph/badge.svg?token=G751JNS6PU)](https://codecov.io/gh/koesterlab/alignoth)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/alignoth.svg?style=flat)](http://bioconda.github.io/recipes/alignoth/README.html)
[![Bioconda version](https://img.shields.io/conda/vn/bioconda/alignoth.svg?style=flat)](http://bioconda.github.io/recipes/alignoth/README.html)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/alignoth/README.html)
[![Licence](https://img.shields.io/conda/l/bioconda/alignoth.svg?style=flat)](http://bioconda.github.io/recipes/alignoth/README.html)

# alignoth

A tool for creating alignment plots from bam files. The generated [vega-lite](https://vega.github.io/vega-lite/) plots are written to stdout per default.

## Usage

```alignoth -b path/to/my.bam -r path/to/my/reference.fa -g chr1:200-300 > plot.vl.json```

### Arguments

The following options are available when using bamboo:

| argument              | short | explanation                                                                                                 | default |
|-----------------------|-------|-------------------------------------------------------------------------------------------------------------|---------|
| bam-path              | -b    | The bam file to be visualized.                                                                              |         |
| reference             | -r    | The path to the reference fasta file                                                                        |         |
| region                | -g    | Chromosome and region for the visualization. Example: 2:132424-132924                                       |         |
| highlight             | -h    | Interval that will be highlighted in the visualization. Example: 132400-132500                              |         |
| max-read-depth        | -d    | Set the maximum rows of reads that will be shown in the alignment plots                                     | 500     |
| max-width             | -w    | Set the maximum width of the resulting alignment plot                                                       | 1024    |
| output                | -o    | If present, data and vega-lite specs of the generated plot will be split and written to the given directory |         |
| data-format           | -f    | Sets the output format for the read, reference and highlight data                                           | json    |
| spec-output           |       | If present vega-lite specs will be written to the given file path                                           |         |
| read-data-output      |       | If present read data will be written to the given file path                                                 |         |
| ref-data-output       |       | If present reference data will be written to the given file path                                            |         |
| highlight-data-output |       | If present highlight data will be written to the given file path                                            |         |

## Installation

There a multiple ways to install alignoth:

#### Bioconda

Rust-Bio-Tools is available via [Bioconda](https://bioconda.github.io).
With Bioconda set up, installation is as easy as

    conda install alignoth

#### Cargo

If the [Rust](https://www.rust-lang.org/tools/install) compiler and associated [Cargo](https://github.com/rust-lang/cargo/) are installed, alignoth may be installed via

    cargo install alignoth

#### Source

Download the source code and within the root directory of source run

    cargo install

## Authors

* [Felix Wiegand](https://github.com/fxwiegand)
* [Johannes Köster](https://github.com/johanneskoester) (https://koesterlab.github.io)
