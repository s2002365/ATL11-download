# ATL11 Antarctic Grounding-Zone Basin Downloader

Downloads spatially-subsetted ATL11 (ICESat-2 land-ice height) granules for each of the 27 Antarctic drainage-basin grounding-zone regions (50 km buffer around each basin's RINGS 2022 grounding line), one basin at a time. 

Script can be re-run periodically as the ATL11 product is updated with new cycles.

Depends on `is2_subset.py` in the parent directory (`../is2_subset.py`) -- see the main repo README for environment setup and Earthdata credentials.

## Usage

**Run this from inside the `ATL11-AIS-GZ/` directory**. The default `--roi-path`/`--out-dir` are relative to your current working directory, not the script's location:

```bash
cd ATL11-AIS-GZ
```

Run every basin:
```bash
python ATL11_download_AIS_GZ.py
```

**Defaults:** `--roi-path RINGS_GL_GSFC_basins_50km_buffer.gpkg` (the file in this folder) and `--out-dir data/` (created inside this folder).


Run just one basin, or point at a different geopackage/output directory: 
```bash 
python ATL11_download_AIS_GZ.py --basin-id 11
python ATL11_download_AIS_GZ.py --roi-path other_basins.gpkg --out-dir /other/directory/atl11_output
```

See `python ATL11_download_AIS_GZ.py --help` for all options.

## Files 

- `ATL11_download_AIS_GZ.py` -- batch download script, one basin at a time.
- `RINGS_GL_GSFC_basins_50km_buffer.gpkg` -- the 27 basin regions (50 km grounding-zone buffer).

