#!/usr/bin/env python
"""
Spatially subset ICESat-2 ATL11 (or ATL06) granules straight from NASA
Earthdata Cloud over HTTPS, transferring only the HDF5 chunks that are needed.

For each granule / pair track (ATL11) or beam (ATL06):
  1. read the full latitude array (small),
  2. keep rows whose latitude falls in any region's latitude range,
  3. read longitude for just those rows and test them against each region
     in the region's own CRS,
  4. read the requested variables for the contiguous row range spanning the
     points inside each region, then (optionally) drop rows outside it.

Regions are given as bounding boxes (in lat/lon or any EPSG) and/or polygons
(shapely geometries, coordinate lists, or GeoJSON-like dicts) in any EPSG.

Authentication: an Earthdata Login bearer token is taken from $EARTHDATA_TOKEN,
or obtained with the urs.earthdata.nasa.gov credentials in ~/.netrc.

Example:
    from is2_subset import read_subset, concat
    res = read_subset(bbox=[(-1610e3, -260e3, -1590e3, -240e3)], bbox_crs=3031,
                      variables=['ref_pt', 'h_corr', 'delta_time'])
    pig = concat(res, region=0)
"""

import argparse
import netrc
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path

import fsspec
import h5py
import numpy as np
import pyproj
import requests
import shapely
import shapely.ops
from shapely.geometry import Polygon, box, shape
from shapely.geometry.polygon import orient

CMR_URL = 'https://cmr.earthdata.nasa.gov/search/granules.umm_json'
URS_HOST = 'urs.earthdata.nasa.gov'

PRODUCTS = {
    'ATL11': {
        'version': '007',
        'groups': ['pt1', 'pt2', 'pt3'],
        'lat': 'latitude',
        'lon': 'longitude',
        # variables that are not indexed by ref_pt; read whole
        'non_row': ['cycle_number'],
        'default_vars': ['ref_pt', 'latitude', 'longitude', 'h_corr',
                         'h_corr_sigma', 'delta_time', 'quality_summary',
                         'cycle_number'],
        # ATL11_{rgt}{region}_{cycle range}_{version}_{revision}.h5
        'rgt_re': r'^ATL11_(?P<rgt>\d{4})\d{2}_\d{4}_\d{3}_\d{2}\.h5$',
    },
    'ATL06': {
        'version': '007',
        'groups': [f'gt{n}{s}/land_ice_segments' for n in (1, 2, 3) for s in ('l', 'r')],
        'lat': 'latitude',
        'lon': 'longitude',
        'non_row': [],
        'default_vars': ['segment_id', 'latitude', 'longitude', 'h_li',
                         'h_li_sigma', 'atl06_quality_summary', 'delta_time'],
        # ATL06_{yyyymmddhhmmss}_{rgt}{cycle}{region}_{version}_{revision}.h5
        'rgt_re': r'^ATL06_\d{14}_(?P<rgt>\d{4})\d{4}_\d{3}_\d{2}\.h5$',
    },
}


def _granule_rgt(name, product):
    """Parse the reference ground track number out of a granule filename."""
    m = re.match(PRODUCTS[product]['rgt_re'], name)
    if not m:
        raise ValueError(f'unrecognized {product} granule name: {name}')
    return int(m.group('rgt'))


def _groups_for_pt(cfg, pt):
    """Resolve a short pair-track/beam name (e.g. 'pt1', 'gt1l') to group path(s)."""
    matches = [g for g in cfg['groups'] if g == pt or g.startswith(pt + '/')]
    if not matches:
        raise ValueError(f"no group matches pt {pt!r}; available: {cfg['groups']}")
    return matches


# ---------------------------------------------------------------------------
# regions
# ---------------------------------------------------------------------------

class Region:
    """A polygon in a given CRS, with its lon/lat footprint for CMR searches."""

    def __init__(self, geom, crs, densify=1000.):
        self.geom = geom
        self.crs = pyproj.CRS.from_user_input(crs)
        shapely.prepare(self.geom)
        if self.crs.equals(pyproj.CRS.from_epsg(4326)):
            self.to_xy = None
            ll = geom
        else:
            self.to_xy = pyproj.Transformer.from_crs(4326, self.crs, always_xy=True)
            to_ll = pyproj.Transformer.from_crs(self.crs, 4326, always_xy=True)
            # densify so edges stay accurate after reprojection
            step = densify if self.crs.is_projected else 0.01
            ll = shapely.ops.transform(to_ll.transform, shapely.segmentize(geom, step))
        self.lonlat = ll
        lon0, self.lat_min, lon1, self.lat_max = ll.bounds
        # a footprint spanning >180 deg of longitude wraps the pole or the
        # antimeridian; fall back to a latitude band for the CMR search
        self.wraps = (lon1 - lon0) > 180

    def contains(self, lon, lat):
        if self.to_xy is None:
            x, y = lon, lat
        else:
            x, y = self.to_xy.transform(lon, lat)
        return shapely.contains_xy(self.geom, x, y)

    def cmr_spatial(self):
        if self.wraps:
            return {'bounding_box': f'-180,{self.lat_min},180,{self.lat_max}'}
        # convex hull keeps the CMR query short and is always a superset of
        # the true shape (never excludes a genuinely intersecting granule);
        # the exact polygon test in contains() still does the real trimming
        hull = orient(self.lonlat.convex_hull, sign=1.0)
        return {'polygon': ','.join(f'{x:.6f},{y:.6f}' for x, y in hull.exterior.coords)}


def _as_geometry(poly):
    if isinstance(poly, shapely.Geometry):
        return poly
    if isinstance(poly, dict):
        return shape(poly.get('geometry', poly))
    return Polygon(poly)


def make_regions(bbox=None, bbox_crs=4326, polygon=None, polygon_crs=4326):
    """
    bbox: (xmin, ymin, xmax, ymax) or a list of them, in bbox_crs
          (for EPSG:4326, x is longitude and y is latitude).
    polygon: shapely geometry, coordinate list, GeoJSON geometry/feature,
             or a list of any of these, in polygon_crs.
    Regions are numbered bboxes first, then polygons.
    """
    regions = []
    if bbox is not None:
        boxes = [bbox] if np.isscalar(bbox[0]) else bbox
        regions += [Region(box(*b), bbox_crs) for b in boxes]
    if polygon is not None:
        if isinstance(polygon, dict) and polygon.get('type') == 'FeatureCollection':
            polys = polygon['features']
        elif isinstance(polygon, (shapely.Geometry, dict)):
            polys = [polygon]
        elif not isinstance(polygon[0], (shapely.Geometry, dict)) and np.isscalar(polygon[0][0]):
            polys = [polygon]       # a single coordinate list [(x, y), ...]
        else:
            polys = polygon
        regions += [Region(_as_geometry(p), polygon_crs) for p in polys]
    if not regions:
        raise ValueError('specify at least one bbox or polygon')
    return regions


# ---------------------------------------------------------------------------
# Earthdata access
# ---------------------------------------------------------------------------

def earthdata_token():
    token = os.environ.get('EARTHDATA_TOKEN')
    if token:
        return token
    auth = netrc.netrc().authenticators(URS_HOST)
    if auth is None:
        raise RuntimeError(f'no {URS_HOST} entry in ~/.netrc and $EARTHDATA_TOKEN unset')
    r = requests.post(f'https://{URS_HOST}/api/users/find_or_create_token',
                      auth=(auth[0], auth[2]), timeout=60)
    r.raise_for_status()
    return r.json()['access_token']


def search_granules(regions, product='ATL11', version=None, temporal=None):
    """Return {granule_name: https_url} for granules intersecting any region."""
    version = version or PRODUCTS[product]['version']
    found = {}
    for reg in regions:
        params = {'short_name': product, 'version': version,
                  'provider': 'NSIDC_CPRD', 'page_size': 2000,
                  **reg.cmr_spatial()}
        if temporal:
            params['temporal'] = f'{temporal[0]},{temporal[1]}'
        headers = {}
        while True:
            r = requests.get(CMR_URL, params=params, headers=headers, timeout=120)
            r.raise_for_status()
            for item in r.json()['items']:
                for u in item['umm'].get('RelatedUrls', []):
                    url = u['URL']
                    if (url.startswith('https://') and url.endswith('.h5')
                            and 'nsidc-cumulus-prod-protected' in url):
                        found[url.rsplit('/', 1)[1]] = url
            if 'CMR-Search-After' not in r.headers:
                break
            headers['CMR-Search-After'] = r.headers['CMR-Search-After']
    return dict(sorted(found.items()))


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------

def _read_group(h5, group, cfg, regions, variables, trim):
    """
    Subset one pair-track/beam group to the points inside each region.
    cfg: the product's entry from PRODUCTS (lat/lon field names, non-row
    variables, etc.). Returns a list of per-region dicts with keys region,
    group, rows, n, data, and fill_values.
    """
    if group not in h5 or cfg['lat'] not in h5[group]:
        return []
    g = h5[group]
    lat = g[cfg['lat']][:]
    lat_ok = np.zeros(lat.shape, bool)
    for reg in regions:
        lat_ok |= (lat >= reg.lat_min) & (lat <= reg.lat_max)
    cand = np.flatnonzero(lat_ok)
    if cand.size == 0:
        return []
    # longitude for just the candidate span
    c0, c1 = cand[0], cand[-1] + 1
    lon = g[cfg['lon']][c0:c1]
    lat_c = lat[c0:c1]

    out = []
    for ir, reg in enumerate(regions):
        sub = (lat_c >= reg.lat_min) & (lat_c <= reg.lat_max)
        inside = np.zeros(sub.shape, bool)
        inside[sub] = reg.contains(lon[sub], lat_c[sub])
        idx = np.flatnonzero(inside)
        if idx.size == 0:
            continue
        r0, r1 = c0 + idx[0], c0 + idx[-1] + 1
        keep = inside[idx[0]:idx[-1] + 1]
        data, fill_values = {}, {}
        for v in variables:
            if v not in g:
                continue
            fv = g[v].attrs.get('_FillValue')
            if fv is not None:
                fill_values[v] = fv
            if v in cfg['non_row']:
                data[v] = g[v][:]
                continue
            arr = g[v][r0:r1]
            data[v] = arr[keep] if trim else arr
        if not trim:
            data['inside'] = keep
        out.append({'region': ir, 'group': group, 'rows': (int(r0), int(r1)),
                    'n': int(keep.sum()), 'data': data, 'fill_values': fill_values})
    return out


def read_granule(url, regions, product='ATL11', variables=None, token=None,
                 trim=True, block_size=2**18, groups=None):
    """
    Subset one granule; returns (list of per-(region, group) results, bytes
    transferred).  block_size sets the HTTP range-request size: smaller blocks
    transfer fewer bytes but need more round trips (ATL11 v007 chunks are
    ~900 rows, and 256 kB moved ~1/6 of a granule in tests).
    groups: pair-track/beam groups to read; defaults to all of them.
    variables: paths relative to each pair/beam group (e.g. 'h_corr',
               'ref_surf/dem_h', 'cycle_stats/tide_ocean'), or a path rooted
               outside the pair/beam groups (e.g. 'ancillary_data/start_rgt')
               to read once per granule instead of once per group.
    """
    cfg = PRODUCTS[product]
    groups = cfg['groups'] if groups is None else groups
    # always return coordinates so results can be located
    variables = list(dict.fromkeys([cfg['lat'], cfg['lon']] + list(variables or cfg['default_vars'])))
    fs = fsspec.filesystem('https', client_kwargs={
        'headers': {'Authorization': f'Bearer {token}'}})
    name = url.rsplit('/', 1)[1]
    with fs.open(url, 'rb', block_size=block_size, cache_type='blockcache') as f, \
            h5py.File(f, 'r') as h5:
        # variables rooted at the file top level (outside the pair/beam
        # groups), e.g. 'ancillary_data/start_rgt', are granule-wide: read
        # once here rather than per group.
        root_vars = [v for v in variables
                    if '/' in v and v.split('/', 1)[0] in h5
                    and v.split('/', 1)[0] not in cfg['groups']]
        row_vars = [v for v in variables if v not in root_vars]
        ancillary = ({v.rsplit('/', 1)[-1]: np.asarray(h5[v][()]).squeeze()
                     for v in root_vars} if root_vars else None)
        results = []
        for group in groups:
            for r in _read_group(h5, group, cfg, regions, row_vars, trim):
                r['granule'] = name
                if ancillary:
                    r['ancillary'] = ancillary
                results.append(r)
        nbytes = getattr(f.cache, 'total_requested_bytes', None)
    return results, nbytes


def read_subset(bbox=None, bbox_crs=4326, polygon=None, polygon_crs=4326,
                product='ATL11', version=None, variables=None, temporal=None,
                granules=None, rgt_pt=None, trim=True, workers=4, max_workers=128,
                block_size=2**18, verbose=False):
    """
    Search CMR and subset every granule that intersects the regions.

    bbox, bbox_crs, polygon, polygon_crs: see make_regions().
    product: 'ATL11' or 'ATL06'.
    variables: paths relative to each pair/beam group, e.g. 'h_corr',
               'cycle_stats/h_rms_mean' (ATL11) or 'fit_statistics/dh_fit_dx' (ATL06).
    temporal: optional (start, end) ISO strings.
    granules: optional {name: url} to skip the CMR search.
    rgt_pt: optional list of (rgt, pt) pairs restricting which reference
            ground tracks and pair tracks/beams are read, e.g.
            [(205, 'pt1'), (205, 'pt3'), (811, None)] for ATL11, or
            [(205, 'gt1l')] for ATL06. pt=None keeps all pair tracks/beams
            for that rgt. Granules for RGTs not listed are skipped entirely
            (not even downloaded); spatial trimming to the region(s) still
            applies to whatever rgt/pt combos are kept.
    trim: if False, return the full contiguous row range plus an 'inside' mask.
    workers: number of granules read in parallel (separate processes), or
             'auto' to use one worker per granule (min(len(granules),
             max_workers)) -- the point past which more workers stop
             helping, since reading one granule is itself a serial chain of
             requests.  With workers > 1, calling scripts need an
             `if __name__ == '__main__':` guard (Python's forkserver/spawn
             start methods re-import them); notebooks are fine.  workers=1
             reads serially in-process.
    max_workers: ceiling used when workers='auto', to avoid overwhelming
             file-descriptor limits or triggering server-side throttling.
    block_size: HTTP range-request size in bytes.

    Returns a list of dicts with keys region, granule, group, rows, n, data.
    """
    product = product.upper()
    cfg = PRODUCTS[product]
    regions = make_regions(bbox, bbox_crs, polygon, polygon_crs)
    if granules is None:
        granules = search_granules(regions, product, version, temporal)

    groups_by_name = {}
    if rgt_pt is not None:
        pts_by_rgt = {}
        for rgt, pt in rgt_pt:
            if pt is None:
                pts_by_rgt[rgt] = None  # keep all pair tracks/beams
            else:
                pts_by_rgt.setdefault(rgt, set())
                if pts_by_rgt[rgt] is not None:
                    pts_by_rgt[rgt].add(pt)
        kept = {}
        for name, url in granules.items():
            rgt = _granule_rgt(name, product)
            if rgt not in pts_by_rgt:
                continue
            pts = pts_by_rgt[rgt]
            groups_by_name[name] = (cfg['groups'] if pts is None else
                                    [g for pt in pts for g in _groups_for_pt(cfg, pt)])
            kept[name] = url
        granules = kept

    

    if verbose:
        print(f'{len(granules)} {product} granules intersect {len(regions)} region(s)')

    if workers == 'auto':
        workers = min(len(granules), max_workers)
        print(f'{len(granules)} ATL11 granules; using workers={workers} (auto)')

    token = earthdata_token()
    args = (regions, product, variables, token, trim, block_size)
    results, total = [], 0
    t_start = time.time()
    if workers > 1:
        pool = ProcessPoolExecutor(max_workers=workers)
        futures = {pool.submit(_read_granule_retry, url, *args,
                               groups=groups_by_name.get(name)): name
                   for name, url in granules.items()}
        done = ((futures[f], f.result) for f in as_completed(futures))
    else:
        pool = None
        done = ((name, partial(_read_granule_retry, url, *args,
                               groups=groups_by_name.get(name)))
                for name, url in granules.items())
    try:
        for i, (name, get) in enumerate(done):
            try:
                res, nbytes = get()
            except Exception as e:
                print(f'failed to read {name}: {e!r}', file=sys.stderr)
                continue
            results += res
            total += nbytes or 0
            if verbose:
                print(f'[{i + 1}/{len(granules)}] {name}: '
                      f'{sum(r["n"] for r in res)} pts, {(nbytes or 0) / 1e6:.1f} MB')
    finally:
        if pool is not None:
            pool.shutdown(cancel_futures=True)
    if verbose:
        print(f'{sum(r["n"] for r in results)} points, {total / 1e6:.1f} MB '
              f'transferred in {time.time() - t_start:.0f} s')
    results.sort(key=lambda r: (r['region'], r['granule'], r['group']))
    return results


def _read_granule_retry(url, *args, tries=3, **kwargs):
    """
    Wraps every failure in a plain RuntimeError (and drops exception chaining
    with `from None`) so it's always picklable back from a worker process --
    some underlying errors (e.g. aiohttp's, which carry a CIMultiDictProxy of
    response headers) aren't, and would otherwise be replaced by an opaque
    "can't pickle ..." TypeError that hides the real failure.
    """
    for k in range(tries):
        try:
            return read_granule(url, *args, **kwargs)
        except OSError as e:
            if k == tries - 1:
                raise RuntimeError(f'{type(e).__name__}: {e}') from None
            time.sleep(2 ** k)
        except Exception as e:
            raise RuntimeError(f'{type(e).__name__}: {e}') from None


def concat(results, region=None):
    """Concatenate results (optionally for one region) along the row axis."""
    sel = [r for r in results if region is None or r['region'] == region]
    if not sel:
        return {}
    out = {}
    for v in sel[0]['data']:
        if v in ('cycle_number',):
            continue
        out[v] = np.concatenate([r['data'][v] for r in sel if v in r['data']])
    out['granule'] = np.concatenate([[r['granule']] * r['n'] for r in sel])
    out['group'] = np.concatenate([[r['group']] * r['n'] for r in sel])
    out['region'] = np.concatenate([[r['region']] * r['n'] for r in sel])
    return out


def write_h5(results, filename):
    """Write results as /region_N/<granule>/<group>/<variable>, with any
    granule-wide values (see read_granule's `variables`) as attrs on
    /region_N/<granule>."""
    with h5py.File(filename, 'w') as out:
        for r in results:
            g = out.require_group(f"region_{r['region']}/{r['granule']}/{r['group']}")
            g.attrs['rows'] = r['rows']
            # Re-attach each variable's _FillValue attribute
            for v, arr in r['data'].items():
                ds = g.create_dataset(v, data=arr, compression='gzip')
                if v in r.get('fill_values', {}):
                    ds.attrs['_FillValue'] = r['fill_values'][v]
            # Attach granule-wide values (e.g. 'start_rgt') as attributes
            for k, v in r.get('ancillary', {}).items():
                out[f"region_{r['region']}/{r['granule']}"].attrs[k] = v


def _write_h5_flat(results, filename):
    """
    Write results as /<group>/<variable>, matching the raw ATL11 file layout
    (pair-track groups at the file root, ancillary_data/<var> as a dataset)
    with no region/granule nesting -- for e.g. read_atl11()-style readers
    that expect a genuine ATL11 file. Assumes `results` covers a single
    region and a given group ('pt1' etc.) appears at most once -- true as
    long as read_subset() is called with one region (bbox/polygon) at a time.
    """
    with h5py.File(filename, 'w') as out:
        ancillary = {}
        for r in results:
            ancillary.update(r.get('ancillary', {}))
            g = out.require_group(r['group'])
            g.attrs['rows'] = r['rows']
            for v, arr in r['data'].items():
                ds = g.create_dataset(v, data=arr, compression='gzip')
                if v in r.get('fill_values', {}):
                    ds.attrs['_FillValue'] = r['fill_values'][v]
        for k, v in ancillary.items():
            out.require_group('ancillary_data').create_dataset(k, data=np.atleast_1d(v))


def write_h5_by_granule(results, outdir, flat=False):
    """
    Write one <outdir>/<granule>.h5 per granule.

    flat: if True, write each file via _write_h5_flat() (matching the raw
          ATL11 layout, e.g. for read_atl11()) instead of write_h5() (which
          nests under /region_N/<granule>/<group>/).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    writer = _write_h5_flat if flat else write_h5
    granules = dict.fromkeys(r['granule'] for r in results)
    for granule in granules:
        sel = [r for r in results if r['granule'] == granule]
        writer(sel, outdir / granule)


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--product', default='ATL11', choices=sorted(PRODUCTS))
    p.add_argument('--version', default=None)
    p.add_argument('--bbox', nargs=4, type=float, action='append',
                   metavar=('XMIN', 'YMIN', 'XMAX', 'YMAX'),
                   help='bounding box in --bbox-crs; repeat for several boxes')
    p.add_argument('--bbox-crs', default='EPSG:4326',
                   help='CRS of --bbox, e.g. EPSG:3031 (default EPSG:4326, lon/lat)')
    p.add_argument('--polygon-file', action='append',
                   help='GeoJSON file of polygon(s); repeatable')
    p.add_argument('--polygon-crs', default='EPSG:4326',
                   help='CRS of the polygon coordinates (default EPSG:4326)')
    p.add_argument('--vars', nargs='+', default=None,
                   help='variables relative to each pair/beam group')
    p.add_argument('--temporal', nargs=2, metavar=('START', 'END'))
    p.add_argument('--no-trim', action='store_true',
                   help='keep the full row range and write an "inside" mask')
    p.add_argument('--workers', type=int, default=4,
                   help='granules read in parallel (default 4)')
    p.add_argument('--block-size', type=int, default=2**18,
                   help='HTTP range-request size in bytes (default 262144); '
                        'smaller moves fewer bytes but makes more requests')
    p.add_argument('--out', required=True, help='output HDF5 file')
    a = p.parse_args(argv)

    polys = None
    if a.polygon_file:
        import json
        polys = []
        for fn in a.polygon_file:
            with open(fn) as fh:
                gj = json.load(fh)
            polys += gj['features'] if gj.get('type') == 'FeatureCollection' else [gj]
    results = read_subset(bbox=a.bbox, bbox_crs=a.bbox_crs, polygon=polys,
                          polygon_crs=a.polygon_crs, product=a.product,
                          version=a.version, variables=a.vars, temporal=a.temporal,
                          trim=not a.no_trim, workers=a.workers,
                          block_size=a.block_size,
                          verbose=True)
    write_h5(results, a.out)
    print(f'wrote {sum(r["n"] for r in results)} points to {a.out}')


if __name__ == '__main__':
    sys.exit(main())
