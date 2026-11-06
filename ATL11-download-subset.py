#!/usr/bin/env python
"""
Download and save ATL11 subsets for the polygon(s) in any vector file
(shapefile, GeoPackage, GeoJSON, ...), in any projection.

Usage:
    python ATL11-download-subset.py --roi-path roi.shp
    python ATL11-download-subset.py --roi-path roi.shp --out-dir /other/atl11_output
    python ATL11-download-subset.py --roi-path rois.gpkg --feature 3
    python ATL11-download-subset.py --roi-path roi.shp --rgt 1269
    python ATL11-download-subset.py --roi-path roi.shp --rgt 1025 80 422
    python ATL11-download-subset.py --roi-path roi.shp --max-workers 256

A file with one polygon is written to <out-dir>/; a file with several is
written one polygon at a time to <out-dir>/feature_NN/, numbered by row.

--out-dir, --feature, --rgt, and --max-workers are all optional; see
`python ATL11-download-subset.py --help` for their defaults.
"""

import argparse
from pathlib import Path
import geopandas as gpd

from is2_subset import _granule_rgt, read_subset, write_h5_by_granule

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
    "ref_surf/x_atc",
    "cycle_stats/tide_ocean",
    "cycle_stats/dac",
    "ancillary_data/start_rgt",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--roi-path', type=Path, required=True,
                   help='Path to the ROI vector file (shapefile, GeoPackage, GeoJSON, ...)')
    p.add_argument('--out-dir', type=Path, default=Path('data'),
                   help='Output directory for subsetted ATL11 files (default: %(default)s)')
    p.add_argument('--feature', type=int, default=None,
                   help='Process only this row of the file (counting from 0), instead of every polygon')
    p.add_argument('--rgt', type=int, nargs='+', default=None,
                   help='Download only these reference ground tracks, keeping all pair tracks '
                        '(default: every RGT crossing the ROI)')
    p.add_argument('--max-workers', type=int, default=128,
                   help='Ceiling for workers="auto" -- one worker per granule, up to this many '
                        '(default: %(default)s). Raise this on a machine with a high ulimit -n '
                        'and spare RAM/network capacity.')
    return p.parse_args()


def main():
    args = parse_args()
    gdf = gpd.read_file(args.roi_path)
    if gdf.crs is None:
        raise ValueError(f"{args.roi_path} has no CRS defined")

    rgt_pt = [(rgt, None) for rgt in args.rgt] if args.rgt else None
    features = [args.feature] if args.feature is not None else range(len(gdf))
    for feature in features:
        polygon = gdf.geometry.iloc[feature]
        out_dir = args.out_dir if len(gdf) == 1 else args.out_dir / f"feature_{feature:02d}"

        print(f"\n=== Feature {feature} of {args.roi_path.name} ===")
        try:
            results = read_subset(
                polygon=polygon,
                polygon_crs=gdf.crs,
                product="ATL11",
                version="007",
                variables=VARIABLES,
                rgt_pt=rgt_pt,
                workers="auto",
                max_workers=args.max_workers,
                trim=True,
                verbose=True,
            )

            if args.rgt:
                found = {_granule_rgt(r["granule"], "ATL11") for r in results}
                missing = sorted(set(args.rgt) - found)
                if missing:
                    print(f"ERROR feature {feature}: no data for RGT(s) {missing} -- check the "
                          "RGT number, that it crosses this polygon, and any failed reads above")

            write_h5_by_granule(
                results,
                outdir=out_dir,
                flat=True,
            )
        except Exception as e:
            print(f"FAILED feature {feature}: {e!r}")
            continue


if __name__ == "__main__":
    main()
