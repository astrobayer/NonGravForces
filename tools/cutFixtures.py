#!/usr/bin/env python3
"""cutFixtures.py -- cut the committed observation excerpts from the full products in DATAPOOL.

    python3 tools/cutFixtures.py --pool <DATAPOOL root> [--out <repository root>]

setupNonGravForces('FetchPool', ...) is what FETCHES (it is the one place of the toolbox that
touches the network, and it owns every URL). This script only READS what setup put into the
pool and writes the excerpts the tests are allowed to see:

  data/spaceweather/dst<yymm>.for.request      Kyoto final Dst, whole months, unchanged
  data/spaceweather/radio_flux_*_observation_2007-09_2008-03_2019-09_2020-03.txt
                                               CLS F30: header unchanged, the daily lines of the
                                               two windows unchanged
  data/samples/itsg/<sat>_<product>_<day>.txt  TU Graz products: the first 600 epochs (10 s), the
                                               GROOPS epoch-count line rewritten, data unchanged
  data/samples/swarm/SW_OPER_DNS<x>POD_2_<first>_<last>.csv
                                               Swarm A/B/C POD densities (VirES HAPI CSV): the JSON
                                               header unchanged, the 30 s lines of the 2014+ days

Every line it keeps is byte-identical to the source (decoded and re-encoded as latin-1). It
prints the size and SHA-256 of every source it read, for the READMEs' digest tables. Standard
library only (CLAUDE.md section 7: no numpy on either tier). Run by
.github/workflows/fetch-data.yml on the self-hosted runner, after setup has filled its pool.

The date filter is the one that cut run 40375328241: year-month-day only, unless the caller
asks for OMNI's year-doy-hour -- "2008    1    1" read as day 1 of 2008 cut all of 2008 once.

Created 2026-10-26 by Claude Opus 5.5 (Anthropic), from tools/fetchData.py of the same PR.
"""
import argparse
import datetime as dt
import gzip
import hashlib
import os
import re
import sys

DAYS = ["2008-01-01", "2020-01-01", "2020-01-02", "2020-01-03"]
DAILY_WINDOWS = [(dt.date(2007, 9, 1), dt.date(2008, 3, 31)),
                 (dt.date(2019, 9, 1), dt.date(2020, 3, 31))]
F30_FILES = ["radio_flux_adjusted_observation.txt", "radio_flux_absolute_observation.txt"]
F30_SUFFIX = "_2007-09_2008-03_2019-09_2020-03.txt"
ITSG_PRODUCTS = ["attitude", "nonConservativeForces", "reducedDynamicOrbit"]
GROOPS_EXCERPT = 600
SWARM_SATS = ["A", "B", "C"]

_DATE_RES = [
    (re.compile(r"^\s*\"?(\d{4})-(\d{2})-(\d{2})"), "ymd"),
    (re.compile(r"^\s*(\d{4})[ /.]+(\d{1,2})[ /.]+(\d{1,2})\b"), "ymd"),   # CLS: "  2008    1    1"
    (re.compile(r"^\s*(\d{4})(\d{2})(\d{2})\b"), "ymd"),
    (re.compile(r"^\s*(\d{4})\s+(\d{1,3})\s+(\d{1,2})\s"), "ydh"),       # OMNI2: year doy hour
]


def line_date(line, kind=None):
    for rx, k in _DATE_RES:
        if k != (kind or "ymd"):
            continue          # year-doy-hour only when asked for: "2008 1 1" is a date otherwise
        m = rx.match(line)
        if not m:
            continue
        a, b, c = (int(x) for x in m.groups())
        try:
            if k == "ymd":
                return dt.date(a, b, c)
            return dt.date(a, 1, 1) + dt.timedelta(days=b - 1)
        except ValueError:
            continue
    return None


def excerpt_by_date(body, windows, kind=None):
    """Header lines (everything before the first dated line) plus the dated lines inside any
    window, each line unchanged. Returns (text, number of data lines)."""
    head, data, started = [], [], False
    for ln in body.decode("latin-1").splitlines(keepends=True):
        d = line_date(ln, kind)
        if d is None:
            if not started:
                head.append(ln)
            continue
        started = True
        if any(lo <= d <= hi for lo, hi in windows):
            data.append(ln)
    return "".join(head + data), len(data)


def cut_groops(body, keep):
    """A GROOPS instrument file: header through the '# ===' rule and, where present, the
    epoch-count line (its count rewritten), then the first `keep` epochs unchanged.
    Returns (text, kept, total) or None when the file is not one."""
    lines = body.decode("latin-1").splitlines(keepends=True)
    rule = next((i for i, ln in enumerate(lines) if ln.startswith("# ====")), None)
    if rule is None or not lines or not lines[0].startswith("groops instrument"):
        return None
    head, data = lines[:rule + 1], lines[rule + 1:]
    if data and re.fullmatch(r"\s*\d+\s*", data[0]):
        count, data = data[0], data[1:]
        kept = data[:keep]
        head.append(re.sub(r"\d+", str(len(kept)), count, count=1))
    else:
        kept = data[:keep]
    return "".join(head + kept), len(kept), len(data)


def satellites(day):
    return ["GRACE-1", "GRACE-2"] if day < "2018" else ["GRACEFO-1", "GRACEFO-2"]


def dst_months(days):
    out = set()
    for day in days:
        y, m = int(day[:4]), int(day[5:7])
        out.add((y, m))
        out.add((y - 1, 12) if m == 1 else (y, m - 1))
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool", required=True, help="DATAPOOL root that setupNonGravForces('FetchPool') filled")
    ap.add_argument("--out", default=".", help="repository root to write the excerpts into")
    ap.add_argument("--days", nargs="*", default=DAYS)
    a = ap.parse_args()
    missing, report = [], []

    def read(path):
        if not os.path.isfile(path):
            missing.append(path)
            return None
        with open(path, "rb") as f:
            b = f.read()
        report.append(f"{len(b):>10}  {hashlib.sha256(b).hexdigest()}  {os.path.relpath(path, a.pool)}")
        return b

    def write(rel, text):
        path = os.path.join(a.out, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as f:
            f.write(text.encode("latin-1") if isinstance(text, str) else text)

    sw = os.path.join(a.pool, "SpaceWeather")
    for y, m in dst_months(a.days):
        name = f"dst{y % 100:02d}{m:02d}.for.request"
        b = read(os.path.join(sw, "Dst", name))
        if b is not None:
            write(os.path.join("data", "spaceweather", name), b)
    for name in F30_FILES:
        b = read(os.path.join(sw, "F30", name))
        if b is not None:
            text, n = excerpt_by_date(b, DAILY_WINDOWS)
            if n != sum((hi - lo).days + 1 for lo, hi in DAILY_WINDOWS):
                missing.append(f"{name}: {n} daily lines in the windows, not one per day")
            write(os.path.join("data", "spaceweather", name[:-4] + F30_SUFFIX), text)
    for day in a.days:
        for sat in satellites(day):
            for product in ITSG_PRODUCTS:
                base = f"{sat}_{product}_{day}.txt"
                b = read(os.path.join(a.pool, "ITSG", "operational", sat, product, day[:4], base + ".gz"))
                if b is None:
                    continue
                cut = cut_groops(gzip.decompress(b), GROOPS_EXCERPT)
                if cut is None:
                    missing.append(f"{base}.gz: not a GROOPS instrument file")
                    continue
                write(os.path.join("data", "samples", "itsg", base), cut[0])

    swarm_days = sorted(d for d in a.days if d >= "2014")
    for sat in SWARM_SATS:
        ds = f"SW_OPER_DNS{sat}POD_2_"
        parts, total = [], 0
        for month in sorted({d[:7] for d in swarm_days}):
            b = read(os.path.join(a.pool, "Swarm", "DNSxPOD", month[:4], f"{ds}{month}.csv"))
            if b is None:
                continue
            want = [dt.date.fromisoformat(d) for d in swarm_days if d.startswith(month)]
            text, n = excerpt_by_date(b, [(d, d) for d in want])
            if n == 0:
                missing.append(f"{ds}{month}.csv: no 30 s line on {', '.join(map(str, want))}")
            parts.append((text, n))
            total += n
        if parts:
            # one file per satellite: the first month's header, every month's dated lines
            head = parts[0][0].splitlines(keepends=True)
            head = [ln for ln in head if line_date(ln) is None]
            body = [ln for text, _ in parts for ln in text.splitlines(keepends=True) if line_date(ln) is not None]
            write(os.path.join("data", "samples", "swarm", f"{ds}{swarm_days[0]}_{swarm_days[-1]}.csv"), "".join(head + body))

    print("sources read (bytes, SHA-256, path below the pool):")
    print("\n".join(report))
    if missing:
        print("NOT CUT:\n  " + "\n  ".join(missing))
        print("-> run setupNonGravForces('FetchPool', \"ifMissing\", 'DataPool', <pool>) first")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
