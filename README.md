# ATL11 Subset Downloader

Spatially subset and download ICESat-2 ATL11 (land-ice height) granules directly from NASA's Earthdata Cloud over HTTPS -- only the requested variables and rows are transferred, not full granules.

Originally written by Ben Smith (UW); adapted by Bryony Freer (SIO).

## Setup

1. Create the conda environment:
    ```bash
   conda env create -f environment.yml
   conda activate atl11-subset
   ```

2. Set up NASA Earthdata Login credentials, either: 
    - a `~/.netrc` entry: 
    ```bash 
    machine urs.earthdata.nasa.gov
    login <your username>
    password <your password>
    ```
    - or the `$EARTHDATA_TOKEN` environment variable. 


## Usage

### Notebook

See `ATL11-download-subset.ipynb` for a full walkthrough.

### Command line

`ATL11-download-subset.py` downloads the same variables for the polygon(s) in any vector file (shapefile, GeoPackage, GeoJSON, ...), in any projection:
```bash
python ATL11-download-subset.py --roi-path roi.shp --out-dir /path/to/output
```

**Defaults:** `--out-dir data/` (relative to your current working directory), `--max-workers` 128

Download only specific reference ground tracks (RGTs), keeping all pair tracks: 
```bash
python ATL11-download-subset.py --roi-path roi.shp --rgt 1269
python ATL11-download-subset.py --roi-path roi.shp --rgt 1025 80 422
```

An `ERROR` line is printed for any requested RGT with no data in the polygon, and the rest of the run continues.

A file with one polygon is written to `<out-dir>/`. A file with several is written one polygon at a time to `<out-dir>/feature_NN/`, numbered by row. To download just one of them: 
```bash
python ATL11-download-subset.py --roi-path rois.gpkg --feature 3
```

Raise the worker ceiling on a machine with spare capacity and a high `ulimit -n`: 
```bash
python ATL11-download-subset.py --roi-path roi.shp --max-workers 256
```

See `python ATL11-download-subset.py --help` for all options.

### Variables

The script downloads the variables in the `VARIABLES` list at the top of `ATL11-download-subset.py`:

- `ref_pt`, `latitude`, `longitude`, `delta_time`
- `h_corr`, `h_corr_sigma`, `h_corr_sigma_systematic`, `quality_summary`
- `cycle_number`
- `ref_surf/geoid_h`, `ref_surf/dem_h`, `ref_surf/x_atc`
- `cycle_stats/tide_ocean`, `cycle_stats/dac`
- `ancillary_data/start_rgt`

To change them, edit that list. The notebook (`variables` cell) and `ATL11-AIS-GZ/ATL11_download_AIS_GZ.py` (`VARIABLES`) each keep their own copy, so update those too if you want them to match.

- Paths are relative to each pair track (`pt1`/`pt2`/`pt3`), e.g. `h_corr` or `ref_surf/dem_h`. Granule-wide values use the full path from the file root, e.g. `ancillary_data/start_rgt`.
- A variable that isn't found is skipped without an error, so check a new path appears in your output files.
- Variables are assumed to have one row per point. Anything else needs adding to `non_row` under `ATL11` in `PRODUCTS` in `is2_subset.py`, which reads it whole, as for `cycle_number`.

Full list of ATL11 variables: https://nsidc.org/sites/default/files/icesat2_atl11_data_dict_v002_0.pdf

## Files
- `is2_subset.py` -- core module: spatial search, subsetting, and output.
- `ATL11-download-subset.ipynb` -- worked example notebook.
- `ATL11-download-subset.py` -- command-line download script for any ROI file.
- `environment.yml` -- conda environment specification.
