# Transcription Start Site HMM for GROcap data

[doi:10.1038/ng.3142](http://dx.doi.org/10.1038/ng.3142)

## Dependencies

This code requires R and is only know to build on Linux and OS X.

The R scripts in the 'pipeline' folder require the following R packages from CRAN:
- parallel
- grid
- Roxygen or Roxygen2

Additionally, the following R packages, which can be found in the "rpkg" folder, must be built from source
- bigWig
- QHMM
- grocaptss

## Setup

Start by clonning this repository with it's submodules:
```
git clone --recursive https://github.com/andrelmartins/grocap.tsshmm
```
After cloning this repository, and assuming you have installed the pre-requisite CRAN packages, go into the "rpkg" folder and run ```make```. This will build the required source R packages and install them in your system.

## Datasets

This pipeline requires the following datasets (in bigWig format, one bigWig file per strand):
- GRO-cap TAP+
- GRO-cap TAP-
- Nuclear CAGE PolyA Plus (single base pair resolution pileups, not the smoothed versions found in the UCSC Genome Browser)
- ChromHMM Predictions (for post-processing statistics)

## Running

- define your dataset (see common.R and rpkg/grocaptss/pkg/R/dataset.R)
- compute normalization totals for both the main dataset (TAP+) (example: data/gm12878.totals.norm.Rdata) and background (TAP-) (example: data/gm12878.totals.norm.back.Rdata) (see pipeline/hmm.parse.gm12878.R)
- create your copy of the pipeline scripts pointing to your dataset
- execute the pipeline steps in turn:
    - hmm.parse (example: Rscript hmm.parser.gm12878.R)
    - postproc
    - postproc.stats [optional]
    - classify
