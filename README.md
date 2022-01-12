# MousiPLIER

MousiPLIER (a portmanteau of Mouse and [MultiPLIER](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30119-X)) consists of two parts. The first part is a pipeline to learn MultiPLIER latent variables from 190,111 recount3 mouse RNAseq samples. The second part uses the latent variables to analyze TODO.

## Quickstart
### Python + dependencies
If you don't already have a working conda installation, follow [these instructions](https://docs.conda.io/en/latest/miniconda.html) to set it up.
Once Conda is installed, you can install the dependencies and activate the environment with
```
conda env create -f env.yml
conda activate mousiplier
```

### R + dependencies
If you don't have an R installation, go ahead and install R by following [these instructions](https://cran.r-project.org/bin/linux/ubuntu/).
We used R 4.1.0 in development, and doing the same will reduce the chance of getting different results due to library differences.

Renv, our R package manager of choice, should bootstrap itself when you run R from the `mousiplier` main directory.
If you are run into R package issues after renv sets up, try running `renv::restore()` and see if that fixes it.

### Running the pipeline
Once you have the dependencies installed and your conda environment activated, run the command `snakemake -j 8` from the `mousiplier` directory and Snakemake will take care of the rest.
The whole pipeline takes a week or two to run, so we don't recommend sitting at the computer waiting on it to finish

## Development environment
The pipeline was developed on an ubuntu 18.04 LTS system with 64GB RAM.
The on-disk PLIER portion of the pipeline uses a few hundred GB of disk space for temp files; our dev computer had around 500GB open.

## Repo Structure

| file/directory | Description |
| -------------- | ----------- |
| `data/`        | This directory stores the raw and processed data that get fed into PLIER |
| `src/`         | This directory contains the R and Python code that makes up the PLIER pipeline |
| `env.yml`      | The file tracking Python (conda) dependencies |
| `renv.lock`    | The file tracking R dependencies |
| `Snakefile`    | The [Snakemake](https://snakemake.readthedocs.io/en/stable/) file that allows the whole pipeline to be run with a single command |

### Data
Our analysis uses mouse expression data from the [Recount3](https://rna.recount.bio/) database.
Our pipeline downloads all the mouse expression data, removes the single-cell RNAseq, and TPM normalizes the remaining samples.

For the matrix of prior pathways we use cell type marker genes from the [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/help.jsp) database, and pathways from [Reactome](https://reactome.org/).

### Src
The scripts that make up data processing/PLIER running pipeline live in `src/`.
More information about the specific files can be found in the README within `src/`.
