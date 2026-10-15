#!/usr/bin/env python
"""
Download and save ATL11 subsets for every Antarctic GZ basin region.

Usage:
    python ATL11_download_AIS_GZ.py
    python ATL11_download_AIS_GZ.py --roi-path other_basins.gpkg --out-dir /other/atl11_output
    python ATL11_download_AIS_GZ.py --basin-id 11

--roi-path, --out-dir, and --basin-id are all optional; see
`python ATL11_download_AIS_GZ.py --help` for their defaults.
"""

import argparse
from pathlib import Path
import sys
import geopandas as gpd

sys.path.insert(0, str(Path(__file__).parent.parent))  # add parent directory to path
from is2_subset import read_subset, write_h5_by_granule

VARIABLES = [
    "ref_pt",
    "latitude",
    "longitude",
    "delta_time",
    "h_corr",
    "h_corr_sigma",
    "h_corr_sigma_systematic",
    "quality_summary",
    "cycle_number",
    "ref_surf/geoid_h",
    "ref_surf/dem_h",
    "cycle_stats/tide_ocean",
    "cycle_stats/dac",
    "ancillary_data/start_rgt",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--roi-path', type=Path,
                   default=Path('RINGS_GL_GSFC_basins_50km_buffer.gpkg'),
                   help='Path to the basin regions geopackage (default: %(default)s)')
    p.add_argument('--out-dir', type=Path, default=Path('data'),
                   help='Output directory for subsetted ATL11 files (default: %(default)s)')
    p.add_argument('--basin-id', type=int, default=None,
                   help='Process only this Basin_ID, instead of every basin in the geopackage')
    return p.parse_args()


def main():
    args = parse_args()
    gdf = gpd.read_file(args.roi_path)

    basin_ids = [args.basin_id] if args.basin_id is not None else sorted(gdf["Basin_ID"].unique())
    for basin_id in basin_ids:
        matches = gdf[gdf["Basin_ID"] == basin_id]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one row for Basin_ID={basin_id}, found {len(matches)}")
        polygon = matches.geometry.iloc[0]

        print(f"\n=== Basin {basin_id} ===")
        try:
            results = read_subset(
                polygon=polygon,
                polygon_crs=gdf.crs,
                product="ATL11",
                version="007",
                variables=VARIABLES,
                workers="auto",
                max_workers=128,
                trim=True,
                verbose=True,
            )

            write_h5_by_granule(
                results,
                outdir=args.out_dir / f"basin_{basin_id:02d}",
                flat=True,
            )
        except Exception as e:
            print(f"FAILED basin {basin_id}: {e!r}")
            continue


if __name__ == "__main__":
    main()
