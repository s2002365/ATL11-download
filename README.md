# ATL11 Subset Downloader

Spatially subset and download ICESat-2 ATL11 (land-ice height) granules directly from NASA's Earthdata Cloud over HTTPS -- only the requested variables and rows are transferred, not full granules.

Originally written by Ben Smith (UW); adapted by Bryony Freer (SIO).

## Setup

1. Create the conda environment:

   `conda env create -f environment.yml`

   `conda activate atl11-subset`

2. Set up NASA Earthdata Login credentials, either: 
    - a `~/.netrc` entry: 

        `machine urs.earthdata.nasa.gov`

        `login <your username>`

        `password <your password>`

    - or the `$EARTHDATA_TOKEN` environment variable. 


## Usage

See `ATL11-download-subset.ipynb` for a full walkthrough.

## Files
- `is2_subset.py` -- core module: spatial search, subsetting, and output.
- `ATL11-download-subset.ipynb` -- worked example notebook.
- `environment.yml` -- conda environment specification.
