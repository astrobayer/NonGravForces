#!/usr/bin/env python3
"""Independent oracle for CERES EBAF-4.2.1 TOA fluxes and the Earth radiation irradiance they imply, pure Python.

Written 2026-09-15 by Claude Opus 5 (Anthropic): written from a specification by an agent that
did not open the MATLAB toolbox (T9). No MATLAB file, no test and no toolbox source was read;
the netCDF reader, the interpolation and the cap quadrature are this file's own code.
Standard library only (struct, math, csv, datetime, argparse, os, re, sys, time) - no numpy,
no netCDF4, no scipy - so it runs on Python 3.12 (Ubuntu) and 3.13 (Windows).

INPUT
  data/ceres/{rsut,rlut,rsdt}_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc
  netCDF-3 classic files (magic "CDF\\x01"; "CDF\\x02" 64-bit offsets is also accepted), read
  by the reader below from the Unidata "NetCDF Classic Format Specification": big-endian header
  numrecs, dim_list, gatt_list, var_list; names and attribute values padded to 4 bytes; each
  variable carries its begin offset; record variables (first dimension of length 0 = the record
  dimension) are read record by record at begin + k * recsize. The dimension order of each
  variable is taken from its dimension ids, never assumed; C order (last dimension fastest).
  Grids are re-indexed by COORDINATE VALUE onto canonical centres
      lat_i = -89.5 + (i - 1), i = 1..180,   lon_j = 0.5 + (j - 1), j = 1..360,
  so a file stored north-to-south or with another dimension order gives the same result.
  A value equal to _FillValue / missing_value, or non-finite, stops the script (the model below
  has no rule for missing data). The global attributes tracking_id and variable_id are printed.

OUTPUT 1: ceres_cells.csv  (proves axis order and orientation)
  columns  variable,month,lat_deg,lon_deg,value,time_days
  For variable in [rsut, rlut, rsdt], month m in {1, 2} (1-based, FILE order) and the cells
  (lat, lon) = (-89.5, 0.5), (-60.5, 120.5), (-0.5, 359.5), (0.5, 0.5), (23.5, 10.5),
  (45.5, 200.5), (70.5, 300.5), (89.5, 179.5): the stored float32 as a Python float (%.17g).
  One more row per variable and month with lat_deg = lon_deg = "global": the area-weighted mean
      sum(F * w) / sum(w),  w = sin(lat + 0.5 deg) - sin(lat - 0.5 deg)
  (the spherical area of a 1-degree band, the same for every longitude), sums by math.fsum.
  time_days is the file's time value of that month (days since 2000-03-01 00:00:00).

OUTPUT 2: ceres_oracle.csv  (Earth radiation irradiance at a satellite, integrated over the cap)
  THE CONTINUOUS MODEL
    a = 6378137.0 m, f = 1/298.257223563, TOA reference height h = 20000.0 m, R = a + h;
    a SPHERE of radius R emits. A point on it has unit direction u = (x, y, z), Earth-fixed:
        phi_d = atan2(z, (1 - f)^2 sqrt(x^2 + y^2))    geodetic latitude [deg]
        lam   = atan2(y, x) mapped to [0, 360)          longitude [deg]
    Time: t = UTC instant as days since 2000-03-01 00:00:00 UTC (proleptic Gregorian, no leap
    seconds). t1 < t2 the two file times, w_t = (t - t1) / (t2 - t1), 0 <= w_t <= 1 required
    (the script stops rather than extrapolate).  F = (1 - w_t) F_month1 + w_t F_month2.
    Space: bilinear between cell centres; phi_d clamped to [-89.5, 89.5] (constant beyond the
    outermost centres); longitude PERIODIC (the cell between 359.5 and 360.5 = 0.5 uses j = 360
    and j = 1).
    Albedo A = F_rsut / F_rsdt where F_rsdt > 0 (both interpolated in space and time), else 0.
    Longwave exitance L = F_rlut (interpolated).
    Satellite at r (Earth-fixed, m, |r| > R), Sun unit vector s (Earth-fixed), irradiance S.
    Surface element dA at p = R u:  rho = r - p, d = |rho|, mu = (u . rho) / d  (seen iff mu > 0),
    mu_s = max(0, u . s):
        dE_alb = A S mu_s mu / (pi d^2) dA,        dE_lw = L mu / (pi d^2) dA
    Direction from the satellite toward the element: e = -rho / d.
        E_alb_vec = integral_{mu > 0} dE_alb e,     E_lw_vec = integral_{mu > 0} dE_lw e
        E_alb_sum = integral_{mu > 0} dE_alb,       E_lw_sum = integral_{mu > 0} dE_lw
  QUADRATURE
    u = cos(alpha) r_hat + sin(alpha) (cos(psi) e1 + sin(psi) e2), dA = R^2 sin(alpha) dalpha dpsi,
    alpha in [0, alpha_max], alpha_max = acos(R / |r|) (the tangent horizon, where mu = 0).
    e1 is the coordinate axis least aligned with r_hat, made orthogonal to r_hat and normalised;
    e2 = r_hat x e1 (safe when r_hat is near the z axis).
    alpha: COMPOSITE 2-POINT GAUSS-LEGENDRE on N_alpha equal panels of [0, alpha_max]. On each psi
      line the panel that holds the terminator (u . s = 0, i.e. b cos(alpha) + c sin(alpha) = 0 with
      b = r_hat . s, c = (cos psi e1 + sin psi e2) . s; at most one root since alpha_max < pi/2) is
      split at the root into two 2-point panels, so the mu_s kink is a panel edge. The bilinear
      kinks on the cell lines are NOT aligned with panels.
    psi: uniform midpoint rule with N_psi nodes psi_j = (j + 1/2) 2 pi / N_psi.
    Why this rule (measured on these cases while writing the script): the composite midpoint rule
      in alpha carried a smooth geometric error of ~3e-5 at N_alpha = 256 (it converged at O(h^2)
      to the uniform-sphere closed forms 2 L (1 - sqrt(1 - (R/|r|)^2)) and L (R/|r|)^2), and the
      thin lit sliver of case 5 did not converge below ~3e-3 without the terminator split.
    d and mu depend on alpha only and are computed without cancellation:
        d^2 = (|r| - R)^2 + 2 |r| R (1 - cos alpha),   u . rho = (|r| - R) - |r| (1 - cos alpha),
        1 - cos alpha = 2 sin^2(alpha / 2),   e = (R u - r) / d.
    Resolution: start at (N_alpha panels, N_psi) = (--n-alpha, --n-psi), then double BOTH until
        rel_err = |value_2N - value_N| / |value_2N| < --tol (default 1e-4)
    for E_alb_sum and E_lw_sum (0/0, i.e. both runs exactly 0, is reported as rel_err = 0).
    The CSV holds the 2N (finer) values; N_alpha (PANELS, 2 nodes each, before any terminator
    split) and N_psi (nodes) are the counts OF THAT FINER run; rel_err compares it with the run at
    half those counts.
  DATA NOTE (printed): in these EBAF files rsut is NOT zero in the polar-night cells where rsdt is
    exactly zero, so A = rsut / rsdt exceeds 1 near the polar-night edge and diverges where the
    interpolated rsdt tends to 0. The script prints, per case, the smallest interpolated rsdt and
    the largest A met at a lit node and the count of lit nodes with rsdt <= 0, so a consumer can
    see whether a case touches that region.
  CASES (lat/lon of r and s are SPHERICAL angles: v = (cos lat cos lon, cos lat sin lon, sin lat))
    1  2020-01-01T00:00:00Z  r (0, 10)      |r| = a + 490 km   s (-23.0, 10.0)   S = 1407.0
    2  2020-01-01T00:00:00Z  r (45, 200)    |r| = a + 490 km   s (-23.0, 250.0)  S = 1407.0
    3  2019-12-20T12:00:00Z  r (-70, 300)   |r| = a + 300 km   s (-23.4, 280.0)  S = 1407.5
    4  2020-01-10T06:00:00Z  r (80, 30)     |r| = a + 490 km   s (-22.0, 120.0)  S = 1407.3
    5  2020-01-10T06:00:00Z  r (88.5, 0)    |r| = a + 490 km   s (-22.0, 0.0)    S = 1407.3
    6  2020-01-01T00:00:00Z  r (10, 359.7)  |r| = a + 490 km   s (-23.0, 0.0)    S = 1407.0
    7  first file time t1    r (25, 15)     |r| = a + 490 km   s (-23.3, 15.0)   S = 1407.9
    8  2020-01-15T00:00:00Z  r (-30, 140)   |r| = a + 490 km   s (60.0, 0.0)     S = 1407.2
  columns  caseId,epoch_utc,days_since_20000301,w_t,rx_m,ry_m,rz_m,sx,sy,sz,S_W_m2,
           E_alb_x,E_alb_y,E_alb_z,E_lw_x,E_lw_y,E_lw_z,E_alb_sum,E_lw_sum,
           rel_err_alb,rel_err_lw,N_alpha,N_psi          (irradiances in W m-2, %.17g)

RUN
  python tools/pyxcheck/ceres_check.py [--data DIR] [--out DIR] [--n-alpha 64] [--n-psi 256]
                                       [--tol 1e-4] [--max-doublings 5]
  Defaults: --data <repo>/data/ceres, --out the folder of this script.
"""
import argparse
import csv
import math
import os
import re
import struct
import sys
import time
from datetime import datetime, timedelta, timezone

# ----------------------------------------------------------------------------------------------
# Constants of the specification
# ----------------------------------------------------------------------------------------------
A_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
H_TOA = 20000.0
R_EMIT = A_WGS84 + H_TOA
VARIABLES = ("rsut", "rlut", "rsdt")
FILE_PATTERN = "%s_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc"
NLAT = 180
NLON = 360
REF_EPOCH = datetime(2000, 3, 1, 0, 0, 0, tzinfo=timezone.utc)
CELLS = ((-89.5, 0.5), (-60.5, 120.5), (-0.5, 359.5), (0.5, 0.5), (23.5, 10.5),
         (45.5, 200.5), (70.5, 300.5), (89.5, 179.5))

# (caseId, epoch or None = first file time, r lat, r lon, altitude m, s lat, s lon, S)
CASES = (
    (1, "2020-01-01T00:00:00Z", 0.0, 10.0, 490e3, -23.0, 10.0, 1407.0),
    (2, "2020-01-01T00:00:00Z", 45.0, 200.0, 490e3, -23.0, 250.0, 1407.0),
    (3, "2019-12-20T12:00:00Z", -70.0, 300.0, 300e3, -23.4, 280.0, 1407.5),
    (4, "2020-01-10T06:00:00Z", 80.0, 30.0, 490e3, -22.0, 120.0, 1407.3),
    (5, "2020-01-10T06:00:00Z", 88.5, 0.0, 490e3, -22.0, 0.0, 1407.3),
    (6, "2020-01-01T00:00:00Z", 10.0, 359.7, 490e3, -23.0, 0.0, 1407.0),
    (7, None, 25.0, 15.0, 490e3, -23.3, 15.0, 1407.9),
    (8, "2020-01-15T00:00:00Z", -30.0, 140.0, 490e3, 60.0, 0.0, 1407.2),
)

# ----------------------------------------------------------------------------------------------
# netCDF-3 classic reader
# ----------------------------------------------------------------------------------------------
NC_DIMENSION = 10
NC_VARIABLE = 11
NC_ATTRIBUTE = 12
NC_CHAR = 2
NC_TYPES = {1: ("b", 1), 2: ("c", 1), 3: ("h", 2), 4: ("i", 4), 5: ("f", 4), 6: ("d", 8)}


def pad4(n):
    return (n + 3) & ~3


class NetCDF3:
    """Minimal netCDF-3 classic (CDF-1) and 64-bit offset (CDF-2) reader."""

    def __init__(self, path):
        with open(path, "rb") as fh:
            self.buf = fh.read()
        self.path = path
        if self.buf[:3] != b"CDF":
            raise ValueError("%s: not a netCDF-3 file (magic %r)" % (path, self.buf[:4]))
        version = self.buf[3]
        if version == 1:
            self.offset_size = 4
        elif version == 2:
            self.offset_size = 8
        else:
            raise ValueError("%s: netCDF-3 version byte %d is not supported" % (path, version))
        self.version = version
        self.pos = 4
        numrecs = self._uint()
        self.streaming = numrecs == 0xFFFFFFFF
        self.numrecs = 0 if self.streaming else numrecs
        self.dims = self._dim_list()            # list of (name, length); length 0 = record dim
        self.gatts = self._att_list()
        self.vars = self._var_list()            # dict name -> record
        self.header_end = self.pos
        self._record_layout()

    # -- primitive readers ------------------------------------------------------------------
    def _uint(self):
        (v,) = struct.unpack_from(">I", self.buf, self.pos)
        self.pos += 4
        return v

    def _offset(self):
        if self.offset_size == 4:
            (v,) = struct.unpack_from(">I", self.buf, self.pos)
        else:
            (v,) = struct.unpack_from(">Q", self.buf, self.pos)
        self.pos += self.offset_size
        return v

    def _name(self):
        n = self._uint()
        s = self.buf[self.pos:self.pos + n].decode("utf-8")
        self.pos += pad4(n)
        return s

    def _values(self, nc_type, nelems):
        if nc_type not in NC_TYPES:
            raise ValueError("%s: unknown nc_type %d" % (self.path, nc_type))
        code, size = NC_TYPES[nc_type]
        nbytes = nelems * size
        raw = self.buf[self.pos:self.pos + nbytes]
        self.pos += pad4(nbytes)
        if nc_type == NC_CHAR:
            return raw.rstrip(b"\x00").decode("utf-8", errors="replace")
        return struct.unpack(">%d%s" % (nelems, code), raw)

    # -- header lists -----------------------------------------------------------------------
    def _dim_list(self):
        tag = self._uint()
        n = self._uint()
        if tag == 0 and n == 0:
            return []
        if tag != NC_DIMENSION:
            raise ValueError("%s: bad dim_list tag %d" % (self.path, tag))
        dims = []
        for _ in range(n):
            name = self._name()
            length = self._uint()
            dims.append((name, length))
        return dims

    def _att_list(self):
        tag = self._uint()
        n = self._uint()
        if tag == 0 and n == 0:
            return {}
        if tag != NC_ATTRIBUTE:
            raise ValueError("%s: bad att_list tag %d" % (self.path, tag))
        atts = {}
        for _ in range(n):
            name = self._name()
            nc_type = self._uint()
            nelems = self._uint()
            atts[name] = self._values(nc_type, nelems)
        return atts

    def _var_list(self):
        tag = self._uint()
        n = self._uint()
        if tag == 0 and n == 0:
            return {}
        if tag != NC_VARIABLE:
            raise ValueError("%s: bad var_list tag %d" % (self.path, tag))
        variables = {}
        order = []
        for _ in range(n):
            name = self._name()
            ndims = self._uint()
            dimids = [self._uint() for _ in range(ndims)]
            atts = self._att_list()
            nc_type = self._uint()
            vsize = self._uint()
            begin = self._offset()
            for d in dimids:
                if d >= len(self.dims):
                    raise ValueError("%s: variable %s has dimension id %d out of range" % (self.path, name, d))
            is_record = ndims > 0 and self.dims[dimids[0]][1] == 0
            variables[name] = {"name": name, "dimids": dimids, "atts": atts, "nc_type": nc_type,
                               "vsize": vsize, "begin": begin, "is_record": is_record}
            order.append(name)
        self.var_order = order
        return variables

    def _record_layout(self):
        rec = [self.vars[k] for k in self.var_order if self.vars[k]["is_record"]]
        if not rec:
            self.recsize = 0
            return
        if len(rec) == 1:
            v = rec[0]
            self.recsize = self._count(v, per_record=True) * NC_TYPES[v["nc_type"]][1]
        else:
            self.recsize = sum(v["vsize"] for v in rec)
        if self.streaming:
            first = min(v["begin"] for v in rec)
            self.numrecs = (len(self.buf) - first) // self.recsize

    # -- public ------------------------------------------------------------------------------
    def dim_names(self, var):
        return [self.dims[d][0] for d in self.vars[var]["dimids"]]

    def shape(self, var):
        return [self.numrecs if self.dims[d][1] == 0 else self.dims[d][1] for d in self.vars[var]["dimids"]]

    def _count(self, v, per_record=False):
        n = 1
        for k, d in enumerate(v["dimids"]):
            if k == 0 and v["is_record"]:
                if per_record:
                    continue
                n *= self.numrecs
            else:
                n *= self.dims[d][1]
        return n

    def read(self, var):
        """Values of a variable as a flat tuple of Python numbers in C order (a str for NC_CHAR)."""
        v = self.vars[var]
        code, size = NC_TYPES[v["nc_type"]]
        if not v["is_record"]:
            n = self._count(v)
            end = v["begin"] + n * size
            if end > len(self.buf):
                raise ValueError("%s: variable %s runs past the end of the file" % (self.path, var))
            if v["nc_type"] == NC_CHAR:
                return self.buf[v["begin"]:end].decode("utf-8", errors="replace")
            return struct.unpack_from(">%d%s" % (n, code), self.buf, v["begin"])
        per = self._count(v, per_record=True)
        out = []
        for k in range(self.numrecs):
            start = v["begin"] + k * self.recsize
            if start + per * size > len(self.buf):
                raise ValueError("%s: record %d of %s runs past the end of the file" % (self.path, k, var))
            out.extend(struct.unpack_from(">%d%s" % (per, code), self.buf, start))
        return tuple(out)


# ----------------------------------------------------------------------------------------------
# CERES fields
# ----------------------------------------------------------------------------------------------
def fail(msg):
    raise SystemExit("ceres_check: " + msg)


def canonical_index(values, first, count, what):
    """File index -> canonical index for 1-degree centres first, first + 1, ..."""
    idx = []
    for v in values:
        k = int(round(v - first))
        if not (0 <= k < count) or abs(v - (first + k)) > 1e-9:
            fail("%s value %r is not a 1-degree cell centre" % (what, v))
        idx.append(k)
    if sorted(idx) != list(range(count)):
        fail("%s axis does not cover every cell centre exactly once" % what)
    return idx


def parse_time_units(units):
    m = re.match(r"^\s*days\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2})(?:[ T](\d{1,2}):(\d{2})(?::(\d{2}(?:\.\d*)?))?)?\s*$", units)
    if not m:
        fail("time units %r are not 'days since YYYY-MM-DD hh:mm:ss'" % units)
    y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
    hh = int(m.group(4) or 0)
    mi = int(m.group(5) or 0)
    ss = float(m.group(6) or 0.0)
    ref = datetime(y, mo, d, hh, mi, tzinfo=timezone.utc) + timedelta(seconds=ss)
    if ref != REF_EPOCH:
        fail("time reference %s differs from the specified 2000-03-01 00:00:00" % ref.isoformat())
    return ref


def load_field(path, varname):
    """Read one CERES file; return dict with time (file order), per-month canonical grids, header ids."""
    nc = NetCDF3(path)
    tracking_id = nc.gatts.get("tracking_id", "<absent>")
    variable_id = nc.gatts.get("variable_id", "<absent>")
    print("%s: netCDF-3 CDF-%d, variable_id = %s, tracking_id = %s"
          % (os.path.basename(path), nc.version, variable_id, tracking_id))
    if variable_id != varname:
        fail("%s: variable_id %r, expected %r" % (path, variable_id, varname))
    for need in ("time", "lat", "lon", varname):
        if need not in nc.vars:
            fail("%s: variable %s missing" % (path, need))

    tatts = nc.vars["time"]["atts"]
    parse_time_units(tatts.get("units", ""))
    calendar = tatts.get("calendar", "standard")
    if calendar not in ("standard", "gregorian", "proleptic_gregorian"):
        fail("%s: calendar %r is not the standard calendar" % (path, calendar))
    times = [float(t) for t in nc.read("time")]
    lats = [float(v) for v in nc.read("lat")]
    lons = [float(v) for v in nc.read("lon")]
    if nc.dim_names("time") != ["time"] or nc.dim_names("lat") != ["lat"] or nc.dim_names("lon") != ["lon"]:
        fail("%s: coordinate variables are not one-dimensional on their own dimension" % path)

    v = nc.vars[varname]
    names = nc.dim_names(varname)
    shape = nc.shape(varname)
    if sorted(names) != ["lat", "lon", "time"]:
        fail("%s: %s has dimensions %r, expected a permutation of (time, lat, lon)" % (path, varname, names))
    if NC_TYPES[v["nc_type"]][0] not in ("f", "d"):
        fail("%s: %s is not a floating-point variable" % (path, varname))
    print("    %s dims (C order) = %s, shape = %s, type = %s, record variable = %s"
          % (varname, names, shape, {5: "float", 6: "double"}[v["nc_type"]], v["is_record"]))
    ntime = shape[names.index("time")]
    if ntime != len(times) or shape[names.index("lat")] != len(lats) or shape[names.index("lon")] != len(lons):
        fail("%s: %s shape does not match its coordinate variables" % (path, varname))

    data = nc.read(varname)
    fills = []
    for key in ("_FillValue", "missing_value"):
        if key in v["atts"]:
            fills.extend(float(x) for x in v["atts"][key])
    bad = sum(1 for x in data if not math.isfinite(x) or x in fills)
    if bad:
        fail("%s: %d values of %s are fill or non-finite; the model has no rule for them" % (path, bad, varname))

    stride = [1, 1, 1]
    for k in (1, 0):
        stride[k] = stride[k + 1] * shape[k + 1]
    st, sy, sx = stride[names.index("time")], stride[names.index("lat")], stride[names.index("lon")]
    ci_lat = canonical_index(lats, -89.5, NLAT, "lat")
    ci_lon = canonical_index(lons, 0.5, NLON, "lon")
    grids = []
    for it in range(ntime):
        g = [0.0] * (NLAT * NLON)
        for iy in range(NLAT):
            row = ci_lat[iy] * NLON
            base = it * st + iy * sy
            for ix in range(NLON):
                g[row + ci_lon[ix]] = float(data[base + ix * sx])
        grids.append(g)
    return {"time": times, "lat": lats, "lon": lons, "grids": grids,
            "tracking_id": tracking_id, "variable_id": variable_id}


def days_since_ref(dt):
    td = dt - REF_EPOCH
    return (td.days * 86400000000 + td.seconds * 1000000 + td.microseconds) / 86400000000


def epoch_string(dt):
    if dt.microsecond:
        return dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def unit_spherical(lat_deg, lon_deg):
    la = math.radians(lat_deg)
    lo = math.radians(lon_deg)
    return (math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la))


# ----------------------------------------------------------------------------------------------
# Cap quadrature
# ----------------------------------------------------------------------------------------------
GL2_X = 1.0 / math.sqrt(3.0)   # 2-point Gauss-Legendre abscissa on [-1, 1]; both weights are 1


def gl2_panel(lo, hi):
    """Nodes and weights of the 2-point Gauss-Legendre rule on [lo, hi]."""
    c = 0.5 * (lo + hi)
    h = 0.5 * (hi - lo)
    return ((c - h * GL2_X, h), (c + h * GL2_X, h))


def terminator_root(b, c, alpha_max):
    """The alpha in (0, alpha_max) where b cos(alpha) + c sin(alpha) = 0, or None.

    b cos + c sin = rho cos(alpha - phi), phi = atan2(c, b); its zeros are pi apart and
    alpha_max < pi / 2, so the open interval holds at most one."""
    if b == 0.0 and c == 0.0:
        return None
    phi = math.atan2(c, b)
    for cand in (phi - 1.5 * math.pi, phi - 0.5 * math.pi, phi + 0.5 * math.pi, phi + 1.5 * math.pi):
        if 0.0 < cand < alpha_max:
            return cand
    return None


def cap_integral(r, s, S, g_rsut, g_rsdt, g_rlut, n_alpha, n_psi):
    """Composite 2-point Gauss-Legendre (alpha, split at the terminator) x uniform midpoint (psi).

    Returns (E_alb_vec, E_lw_vec, E_alb_sum, E_lw_sum, diagnostics)."""
    rn = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    if not rn > R_EMIT:
        fail("satellite radius %.3f m is not above the emitting sphere" % rn)
    rh = (r[0] / rn, r[1] / rn, r[2] / rn)
    k = min(range(3), key=lambda i: abs(rh[i]))
    ax = [0.0, 0.0, 0.0]
    ax[k] = 1.0
    dot = ax[0] * rh[0] + ax[1] * rh[1] + ax[2] * rh[2]
    e1 = [ax[i] - dot * rh[i] for i in range(3)]
    n1 = math.sqrt(e1[0] ** 2 + e1[1] ** 2 + e1[2] ** 2)
    e1 = (e1[0] / n1, e1[1] / n1, e1[2] / n1)
    e2 = (rh[1] * e1[2] - rh[2] * e1[1], rh[2] * e1[0] - rh[0] * e1[2], rh[0] * e1[1] - rh[1] * e1[0])
    rhx, rhy, rhz = rh
    sx, sy, sz = s
    b_sun = rhx * sx + rhy * sy + rhz * sz

    alpha_max = math.acos(R_EMIT / rn)
    ha = alpha_max / n_alpha
    hp = 2.0 * math.pi / n_psi
    h0 = rn - R_EMIT

    def node(a, wa):
        """(cos a, sin a, w, w R/d, w |r|/d) with w = R^2 sin(a) mu / (pi d^2) * wa * hpsi."""
        sh = math.sin(0.5 * a)
        omc = 2.0 * sh * sh                          # 1 - cos(alpha), no cancellation
        d2 = h0 * h0 + 2.0 * rn * R_EMIT * omc
        d = math.sqrt(d2)
        mu = (h0 - rn * omc) / d
        if mu <= 0.0:
            return (math.cos(a), math.sin(a), 0.0, 0.0, 0.0)
        w = R_EMIT * R_EMIT * math.sin(a) * mu / (math.pi * d2) * wa * hp
        return (math.cos(a), math.sin(a), w, w * R_EMIT / d, w * rn / d)

    std = [node(a, wa) for i in range(n_alpha) for (a, wa) in gl2_panel(i * ha, (i + 1) * ha)]

    kk = (1.0 - F_WGS84) ** 2
    deg = 180.0 / math.pi
    atan2 = math.atan2
    sqrt = math.sqrt
    floor = math.floor
    # accumulators: sum of value*w, value*w|r|/d, value*wR/d*(x, y, z)
    l_w = l_p = l_x = l_y = l_z = 0.0
    a_w = a_p = a_x = a_y = a_z = 0.0
    lit_dark = 0            # lit nodes (mu_s > 0) where interpolated rsdt <= 0, so A := 0
    lit_min_rsdt = math.inf
    lit_max_A = 0.0
    n_split = 0
    for j in range(n_psi):
        cp = math.cos((j + 0.5) * hp)
        sp = math.sin((j + 0.5) * hp)
        vx = cp * e1[0] + sp * e2[0]
        vy = cp * e1[1] + sp * e2[1]
        vz = cp * e1[2] + sp * e2[2]
        root = terminator_root(b_sun, vx * sx + vy * sy + vz * sz, alpha_max)
        if root is None:
            nodes = std
        else:
            n_split += 1
            kp = min(int(root / ha), n_alpha - 1)
            sub = [node(a, wa) for (a, wa) in gl2_panel(kp * ha, root) + gl2_panel(root, (kp + 1) * ha)]
            nodes = std[:2 * kp] + sub + std[2 * kp + 2:]
        for ca, sa, w, wq, wp in nodes:
            x = ca * rhx + sa * vx
            y = ca * rhy + sa * vy
            z = ca * rhz + sa * vz
            lat = atan2(z, kk * sqrt(x * x + y * y)) * deg
            lon = atan2(y, x) * deg
            if lon < 0.0:
                lon += 360.0
            if lat < -89.5:
                lat = -89.5
            elif lat > 89.5:
                lat = 89.5
            yi = lat + 89.5
            i0 = int(yi)
            if i0 > NLAT - 2:
                i0 = NLAT - 2
            fy = yi - i0
            xj = lon - 0.5
            j0 = floor(xj)
            fx = xj - j0
            j0 %= NLON
            j1 = (j0 + 1) % NLON
            r0 = i0 * NLON
            r1 = r0 + NLON
            c00 = r0 + j0
            c01 = r0 + j1
            c10 = r1 + j0
            c11 = r1 + j1
            gx = 1.0 - fx
            gy = 1.0 - fy
            val = gy * (gx * g_rlut[c00] + fx * g_rlut[c01]) + fy * (gx * g_rlut[c10] + fx * g_rlut[c11])
            l_w += val * w
            l_p += val * wp
            vq = val * wq
            l_x += vq * x
            l_y += vq * y
            l_z += vq * z
            us = x * sx + y * sy + z * sz
            if us > 0.0:
                dt = gy * (gx * g_rsdt[c00] + fx * g_rsdt[c01]) + fy * (gx * g_rsdt[c10] + fx * g_rsdt[c11])
                if dt > 0.0:
                    ut = gy * (gx * g_rsut[c00] + fx * g_rsut[c01]) + fy * (gx * g_rsut[c10] + fx * g_rsut[c11])
                    alb = ut / dt
                    if dt < lit_min_rsdt:
                        lit_min_rsdt = dt
                    if alb > lit_max_A:
                        lit_max_A = alb
                    val = alb * us
                    a_w += val * w
                    a_p += val * wp
                    vq = val * wq
                    a_x += vq * x
                    a_y += vq * y
                    a_z += vq * z
                else:
                    lit_dark += 1
    # e = (R u - r) / d = (R / d) u - (|r| / d) r_hat
    e_lw = (l_x - l_p * rhx, l_y - l_p * rhy, l_z - l_p * rhz)
    e_alb = (S * (a_x - a_p * rhx), S * (a_y - a_p * rhy), S * (a_z - a_p * rhz))
    diag = {"lit_dark": lit_dark, "lit_min_rsdt": lit_min_rsdt, "lit_max_A": lit_max_A,
            "psi_lines_split": n_split, "alpha_max_deg": math.degrees(alpha_max)}
    return e_alb, e_lw, S * a_w, l_w, diag


def rel_err(fine, coarse):
    if fine == 0.0 and coarse == 0.0:
        return 0.0
    if fine == 0.0:
        return math.inf
    return abs(fine - coarse) / abs(fine)


def fmt(x):
    return "%.17g" % x


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", default=os.path.join(root, "data", "ceres"), help="folder of the three .nc files")
    parser.add_argument("--out", default=here, help="output folder for the two CSV files")
    parser.add_argument("--n-alpha", type=int, default=128, help="alpha panels of the coarse run")
    parser.add_argument("--n-psi", type=int, default=512, help="psi nodes of the coarse run")
    parser.add_argument("--tol", type=float, default=1e-4, help="target rel_err for E_alb_sum and E_lw_sum")
    parser.add_argument("--max-doublings", type=int, default=3, help="refinements allowed beyond the first pair")
    parser.add_argument("--cases", default="", help="comma-separated caseIds to run (default: all)")
    args = parser.parse_args(argv)
    t_start = time.perf_counter()

    fields = {}
    for name in VARIABLES:
        fields[name] = load_field(os.path.join(args.data, FILE_PATTERN % name), name)
    ref = fields[VARIABLES[0]]
    for name in VARIABLES[1:]:
        for key in ("time", "lat", "lon"):
            if fields[name][key] != ref[key]:
                fail("%s coordinate of %s differs from %s" % (key, name, VARIABLES[0]))
    times = ref["time"]
    if len(times) != 2:
        fail("expected 2 time records, found %d" % len(times))
    for m, t in enumerate(times):
        print("month %d (file order): time = %.17g days = %s" % (m + 1, t, epoch_string(REF_EPOCH + timedelta(days=t))))
    order = sorted(range(2), key=lambda k: times[k])
    t1, t2 = times[order[0]], times[order[1]]
    if not t1 < t2:
        fail("the two file times are not distinct")
    for m in range(2):
        dt = fields["rsdt"]["grids"][m]
        ut = fields["rsut"]["grids"][m]
        zero = sum(1 for x in dt if x == 0.0)
        both = sum(1 for x, y in zip(dt, ut) if x == 0.0 and y != 0.0)
        print("month %d: %d cells with rsdt == 0, of which %d have rsut != 0 (A = rsut/rsdt diverges toward them)"
              % (m + 1, zero, both))

    # ---- output 1: cells and global means -----------------------------------------------------
    band_w = [math.sin(math.radians(-89.0 + i)) - math.sin(math.radians(-90.0 + i)) for i in range(NLAT)]
    wsum = math.fsum(band_w[i] for i in range(NLAT) for _ in range(NLON))
    cells_path = os.path.join(args.out, "ceres_cells.csv")
    with open(cells_path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["variable", "month", "lat_deg", "lon_deg", "value", "time_days"])
        for name in VARIABLES:
            for m in range(2):
                g = fields[name]["grids"][m]
                for lat, lon in CELLS:
                    i = int(round(lat + 89.5))
                    j = int(round(lon - 0.5))
                    w.writerow([name, m + 1, fmt(lat), fmt(lon), fmt(g[i * NLON + j]), fmt(times[m])])
                gm = math.fsum(g[i * NLON + j] * band_w[i] for i in range(NLAT) for j in range(NLON)) / wsum
                w.writerow([name, m + 1, "global", "global", fmt(gm), fmt(times[m])])
                print("    %s month %d global area-weighted mean = %.6f W m-2" % (name, m + 1, gm))
    print("wrote %s" % cells_path)

    # ---- output 2: cap integrals ----------------------------------------------------------------
    wanted = set(int(c) for c in args.cases.split(",") if c.strip()) if args.cases else None
    oracle_path = os.path.join(args.out, "ceres_oracle.csv")
    rows = []
    for case_id, epoch, rlat, rlon, alt, slat, slon, S in CASES:
        if wanted is not None and case_id not in wanted:
            continue
        tc = time.perf_counter()
        if epoch is None:
            t = t1
            dt = REF_EPOCH + timedelta(days=t1)
        else:
            dt = datetime.strptime(epoch, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
            t = days_since_ref(dt)
        if not (t1 <= t <= t2):
            fail("case %d: epoch %s (%.17g d) lies outside the file interval [%.17g, %.17g]; no extrapolation"
                 % (case_id, epoch_string(dt), t, t1, t2))
        wt = (t - t1) / (t2 - t1)
        g1 = order[0]
        g2 = order[1]
        grids = {}
        for name in VARIABLES:
            a = fields[name]["grids"][g1]
            b = fields[name]["grids"][g2]
            grids[name] = [(1.0 - wt) * x + wt * y for x, y in zip(a, b)]
        rmag = A_WGS84 + alt
        u = unit_spherical(rlat, rlon)
        r = (rmag * u[0], rmag * u[1], rmag * u[2])
        s = unit_spherical(slat, slon)

        na, npsi = args.n_alpha, args.n_psi
        coarse = cap_integral(r, s, S, grids["rsut"], grids["rsdt"], grids["rlut"], na, npsi)
        level = 0
        while True:
            fine = cap_integral(r, s, S, grids["rsut"], grids["rsdt"], grids["rlut"], 2 * na, 2 * npsi)
            ea = rel_err(fine[2], coarse[2])
            el = rel_err(fine[3], coarse[3])
            print("    case %d  N_alpha %5d N_psi %5d  E_alb_sum %.10g  E_lw_sum %.10g  rel_err alb %.3g lw %.3g"
                  % (case_id, 2 * na, 2 * npsi, fine[2], fine[3], ea, el))
            if (ea < args.tol and el < args.tol) or level >= args.max_doublings:
                break
            coarse = fine
            na *= 2
            npsi *= 2
            level += 1
        if not (ea < args.tol and el < args.tol):
            print("    WARNING case %d: tolerance %.1e not reached after %d doublings" % (case_id, args.tol, level))
        ealb, elw, salb, slw, diag = fine
        print("    case %d  alpha_max %.4f deg, psi lines split at the terminator %d of %d, lit nodes with rsdt <= 0: %d,"
              " lit min rsdt %.6g W m-2, lit max A %.6g"
              % (case_id, diag["alpha_max_deg"], diag["psi_lines_split"], 2 * npsi, diag["lit_dark"],
                 diag["lit_min_rsdt"], diag["lit_max_A"]))
        rows.append([str(case_id), epoch_string(dt), fmt(t), fmt(wt)] + [fmt(x) for x in r] +
                    [fmt(x) for x in s] + [fmt(S)] + [fmt(x) for x in ealb] + [fmt(x) for x in elw] +
                    [fmt(salb), fmt(slw), fmt(ea), fmt(el), str(2 * na), str(2 * npsi)])
        print("case %d  %s  w_t %.6f  E_alb %.6f W m-2  E_lw %.6f W m-2  (%.1f s)"
              % (case_id, epoch_string(dt), wt, salb, slw, time.perf_counter() - tc))

    with open(oracle_path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["caseId", "epoch_utc", "days_since_20000301", "w_t", "rx_m", "ry_m", "rz_m", "sx", "sy", "sz",
                    "S_W_m2", "E_alb_x", "E_alb_y", "E_alb_z", "E_lw_x", "E_lw_y", "E_lw_z", "E_alb_sum",
                    "E_lw_sum", "rel_err_alb", "rel_err_lw", "N_alpha", "N_psi"])
        w.writerows(rows)
    print("wrote %s (%d cases), total %.1f s" % (oracle_path, len(rows), time.perf_counter() - t_start))
    return 0


if __name__ == "__main__":
    sys.exit(main())
