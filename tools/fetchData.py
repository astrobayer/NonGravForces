#!/usr/bin/env python3
"""fetchData.py -- fetch the public data three open gaps need, and cut excerpts of it.

Run by .github/workflows/fetch-data.yml on a GitHub-hosted runner, because the
development container cannot reach these hosts (2026-10-26). Standard library only,
like every other script in tools/ (CLAUDE.md section 7: no numpy on either tier).

    python3 tools/fetchData.py --out data/incoming --raw <scratch folder>

Targets (HANDOVER Part H and the open issues):
  dst       O-9   hourly Dst: NASA OMNI2 yearly files (which carry Kyoto's Dst) and the
                  Kyoto WDC pages themselves, probed
  f30       O-9   daily 30 cm (1 GHz) solar radio flux: Nobeyama NoRP and CLS, probed
  itsg      N-12  TU Graz operational products per satellite and day: attitude, and for
            O-19  comparison the non-conservative forces and the reduced-dynamic orbit
  classic   N-08  classic GRACE Level-1B for 2008-01-01 (A and B) from GFZ ISDC

This is a STAGING step, not a fixture step: the server formats have not been read by
anyone in this project yet, so every product is committed as (a) the listing that found
it, (b) the SHA-256 and size of the full download, (c) a head sample, and (d) a date
excerpt with every data line unchanged where the format allowed one. A later commit
moves what is wanted into data/samples or data/spaceweather. Nothing is guessed silently:
a URL that fails is recorded with its status, and a reply whose head does not look like
the product is recorded as such and not used ("a status code is not a file", CLAUDE.md).

Full downloads stay in --raw (uploaded as a run artefact); nothing above CAP bytes is
written to --out.

Created 2026-10-26 by Claude Opus 5.5 (Anthropic).
"""
import argparse
import datetime as dt
import hashlib
import html.parser
import io
import json
import os
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile

UA = "NonGravForces-fetchData/1.0 (+https://github.com/astrobayer/NonGravForces)"
CAP = 400_000            # largest file committed to --out, bytes
HEAD_LINES = 60          # lines in a head sample

# Date windows of the excerpts: the GRACE-FO fixture days (2020-01-01..03) and the
# classic fixture day (2008-01-01), with margin for 3-day histories.
HOURLY_WINDOWS = [(dt.date(2007, 12, 25), dt.date(2008, 1, 7)),
                  (dt.date(2019, 12, 25), dt.date(2020, 1, 7))]
# Daily windows match data/spaceweather (2019-09..2020-03) and the same months of 2007/08,
# enough for the 81-day centred means a density model takes.
DAILY_WINDOWS = [(dt.date(2007, 9, 1), dt.date(2008, 3, 31)),
                 (dt.date(2019, 9, 1), dt.date(2020, 3, 31))]

CLASSIC_DAY = "2008-01-01"
CLASSIC_PRODUCTS = ("ACC1B", "SCA1B", "GNV1B", "THR1B", "MAS1B", "AHK1B")
CLASSIC_EXCERPT = 600    # records, as for data/samples (the first 10 min of a 1 Hz product)


# ---------------------------------------------------------------------------- record
class Record:
    """What happened, for PROBE.md (human) and FETCHED.json (machine)."""

    def __init__(self):
        self.lines = []
        self.items = []

    def h(self, text):
        self.lines += ["", "## " + text, ""]

    def say(self, text):
        self.lines.append(text)

    def item(self, **kw):
        self.items.append(kw)
        status = kw.get("status")
        self.say(f"- `{kw.get('url', '')}` -> {status}"
                 + (f", {kw['bytes']} bytes, sha256 `{kw['sha256']}`" if kw.get("sha256") else "")
                 + (f" -- {kw['note']}" if kw.get("note") else ""))


# ---------------------------------------------------------------------------- http
def get(url, max_bytes=2_000_000_000, tries=3, timeout=90):
    """GET url. Returns (status, bytes or None, final url, error text)."""
    err = ""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                chunks, n = [], 0
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    n += len(c)
                    if n > max_bytes:
                        return r.status, None, r.geturl(), f"larger than {max_bytes} bytes, not kept"
                    chunks.append(c)
                return r.status, b"".join(chunks), r.geturl(), ""
        except urllib.error.HTTPError as e:
            if e.code in (404, 403, 401, 410):
                return e.code, None, url, str(e)
            err = str(e)
        except Exception as e:  # network errors: retry
            err = f"{type(e).__name__}: {e}"
        time.sleep(2 * 2 ** k)
    return 0, None, url, err


class _Links(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for k, v in attrs:
                if k == "href" and v:
                    self.hrefs.append(v)


def listing(url, rec, timeout=90):
    """Entries (name, is_dir, absolute url) of an autoindex page, or of an FTP LIST, below url."""
    if not url.endswith("/"):
        url += "/"
    status, body, _, err = get(url, max_bytes=20_000_000, timeout=timeout)
    if body is None:
        rec.item(url=url, status=status or "failed", note=err or "no listing")
        return None
    if url.startswith("ftp://"):
        out = []
        for ln in body.decode("latin-1").splitlines():
            parts = ln.split(None, 8)
            if len(parts) == 9 and parts[8] not in (".", ".."):
                d = parts[0].startswith("d")
                out.append((parts[8] + ("/" if d else ""), d, url + parts[8] + ("/" if d else "")))
        time.sleep(0.2)
        return out
    p = _Links()
    p.feed(body.decode("utf-8", "replace"))
    out, seen = [], set()
    for href in p.hrefs:
        if href.startswith(("?", "#", "mailto:")) or href in ("../", "./", "/"):
            continue
        absu = urllib.parse.urljoin(url, href)
        if not absu.startswith(url) or absu == url:
            continue          # parent, sibling or another host
        name = absu[len(url):]
        if not name or name in seen or "/" in name.rstrip("/"):
            continue
        seen.add(name)
        out.append((urllib.parse.unquote(name), name.endswith("/"), absu))
    time.sleep(0.2)           # be polite to research servers
    return out


# ---------------------------------------------------------------------------- helpers
def sha(b):
    return hashlib.sha256(b).hexdigest()


def is_text(b):
    return b"\x00" not in b[:65536]


def save_raw(raw, name, b):
    path = os.path.join(raw, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(b)
    return path


def save_out(out, name, data, rec):
    """Write a committed file; refuses anything above CAP (it would bloat the history)."""
    b = data.encode("utf-8") if isinstance(data, str) else data
    if len(b) > CAP:
        rec.say(f"  - NOT committed `{name}`: {len(b)} bytes > CAP {CAP}")
        return None
    path = os.path.join(out, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(b)
    rec.say(f"  - committed `{name}` ({len(b)} bytes)")
    return path


def head_sample(b, n=HEAD_LINES):
    if not is_text(b):
        return "(binary) first 256 bytes, hex:\n" + b[:256].hex(" ") + "\n"
    return "".join(b.decode("latin-1").splitlines(keepends=True)[:n])


_DATE_RES = [
    (re.compile(r"^\s*\"?(\d{4})-(\d{2})-(\d{2})"), "ymd"),
    (re.compile(r"^\s*(\d{4})[ /.](\d{1,2})[ /.](\d{1,2})\b"), "ymd"),
    (re.compile(r"^\s*(\d{4})(\d{2})(\d{2})\b"), "ymd"),
    (re.compile(r"^\s*(\d{4})\s+(\d{1,3})\s+(\d{1,2})\s"), "ydh"),   # OMNI2: year doy hour
]


def line_date(line, kind=None):
    for rx, k in _DATE_RES:
        if kind and k != kind:
            continue
        m = rx.match(line)
        if not m:
            continue
        a, b_, c = (int(x) for x in m.groups())
        try:
            if k == "ymd":
                return dt.date(a, b_, c)
            return dt.date(a, 1, 1) + dt.timedelta(days=b_ - 1)
        except ValueError:
            continue
    return None


def excerpt_by_date(b, windows, kind=None):
    """Header lines (everything before the first dated line) plus the dated lines in any
    window, each line unchanged. Returns (text, n data lines, date format seen)."""
    lines = b.decode("latin-1").splitlines(keepends=True)
    head, data, started = [], [], False
    for ln in lines:
        d = line_date(ln, kind)
        if d is None:
            if not started:
                head.append(ln)
            continue
        started = True
        if any(lo <= d <= hi for lo, hi in windows):
            data.append(ln)
    return "".join(head + data), len(data)


# ---------------------------------------------------------------------------- dst (O-9)
def target_dst(out, raw, rec):
    rec.h("Dst (O-9)")
    base = "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/"
    status, b, _, err = get(base + "omni2.text", max_bytes=5_000_000)
    rec.item(url=base + "omni2.text", status=status or err, bytes=len(b or b""),
             sha256=sha(b) if b else "", note="OMNI2 format description")
    if b:
        save_out(out, "omni2/omni2.text", b, rec)
    for year in sorted({w[0].year for w in HOURLY_WINDOWS} | {w[1].year for w in HOURLY_WINDOWS}):
        url = f"{base}omni2_{year}.dat"
        status, b, _, err = get(url, max_bytes=50_000_000)
        ok = bool(b) and is_text(b) and line_date(b.decode("latin-1").splitlines()[0], "ydh") is not None
        rec.item(url=url, status=status or err, bytes=len(b or b""), sha256=sha(b) if b else "",
                 note="" if ok else "head does not look like OMNI2 -> not used")
        if not ok:
            continue
        save_raw(raw, f"omni2/omni2_{year}.dat", b)
        for lo, hi in HOURLY_WINDOWS:
            if lo.year != year and hi.year != year:
                continue
            text, n = excerpt_by_date(b, [(lo, hi)], kind="ydh")
            rec.say(f"  - window {lo}..{hi} from {year}: {n} hourly lines")
            if n:
                save_out(out, f"omni2/omni2_{year}_{lo}_{hi}.dat", text, rec)
    # Kyoto is the primary source of Dst; OMNI copies it. Run 40358650063 found FINAL values in
    # WDC format for both fixture months (provisional/realtime answer 403 for them); the month
    # before each is taken too, for histories across the new year.
    kyoto = "https://wdc.kugi.kyoto-u.ac.jp/dst_final/"
    for ym in ("200712", "200801", "201912", "202001"):
        url = f"{kyoto}{ym}/dst{ym[2:]}.for.request"
        status, b, _, err = get(url, max_bytes=2_000_000, tries=2)
        ok = bool(b) and b.startswith(b"DST" + ym[2:].encode())
        rec.item(url=url, status=status or err, bytes=len(b or b""), sha256=sha(b) if b else "",
                 note="" if ok else "not a WDC Dst record -> not used")
        if ok:
            save_out(out, f"kyoto/dst{ym[2:]}.for.request", b, rec)


# ---------------------------------------------------------------------------- f30 (O-9)
def target_f30(out, raw, rec):
    rec.h("F30, the 30 cm solar radio flux (O-9)")
    norp = "https://solar.nro.nao.ac.jp/norp/data/daily/"
    entries = listing(norp, rec)
    if entries is not None:
        rec.say(f"listing of `{norp}`: " + ", ".join(f"`{e[0]}`" for e in entries[:80]))
    names = [e[0] for e in (entries or []) if not e[1]]
    want = [n for n in names if re.search(r"dailyflux|flux", n, re.I)] or ["TYKW-NoRP_dailyflux.txt"]
    for name in want[:4]:
        url = norp + name
        status, b, _, err = get(url, max_bytes=50_000_000)
        rec.item(url=url, status=status or err, bytes=len(b or b""), sha256=sha(b) if b else "")
        if not b or not is_text(b):
            continue
        save_raw(raw, "norp/" + name, b)
        save_out(out, f"norp/{name}.head.txt", head_sample(b), rec)
        text, n = excerpt_by_date(b, DAILY_WINDOWS)
        rec.say(f"  - daily windows: {n} dated lines")
        if n:
            save_out(out, f"norp/{os.path.splitext(name)[0]}_2007-09_2008-03_2019-09_2020-03.txt", text, rec)
    # CLS publishes F30 as used by DTM2020. Its page links no HTTP data, only an FTP archive
    # (run 40358650063): ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/.
    cls = "ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/"
    entries = listing(cls, rec, timeout=120)
    if entries is not None:
        rec.say(f"listing of `{cls}`: " + ", ".join(f"`{e[0]}`" for e in entries[:80]))
    for name, isdir, url in (entries or [])[:12]:
        if isdir:
            continue
        s2, b2, _, e2 = get(url, max_bytes=100_000_000, tries=2, timeout=120)
        rec.item(url=url, status=s2 or e2, bytes=len(b2 or b""), sha256=sha(b2) if b2 else "")
        if not b2:
            continue
        save_raw(raw, "cls/" + name, b2)
        for m, body in unpack_any(b2, name):
            save_out(out, f"cls/{m}.head.txt", head_sample(body), rec)
            if is_text(body):
                text, n = excerpt_by_date(body, DAILY_WINDOWS)
                rec.say(f"  - `{m}` daily windows: {n} dated lines")
                if n:
                    save_out(out, f"cls/{os.path.splitext(m)[0]}_2007-09_2008-03_2019-09_2020-03.txt", text, rec)


# ---------------------------------------------------------------------------- ITSG (N-12, O-19)
# TU Graz's operational satellite products, found by run 40358650063's crawl (itsg/TREE.txt):
# per satellite and day, attitude, nonConservativeForces and orbits in GROOPS instrument format.
# GRACE-1/-2 are GRACE-A/-B, GRACEFO-1/-2 are GRACE-C/-D.
ITSG = "https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/"
ITSG_DAYS = [("GRACE-1", "2008-01-01"), ("GRACE-2", "2008-01-01"),
             ("GRACEFO-1", "2020-01-01"), ("GRACEFO-1", "2020-01-02"), ("GRACEFO-1", "2020-01-03"),
             ("GRACEFO-2", "2020-01-01")]
ITSG_PRODUCTS = ("attitude", "nonConservativeForces", "reducedDynamicOrbit")
GROOPS_EXCERPT = 600


def cut_groops(body, keep):
    """A GROOPS instrument file: header (through the '# ===' rule and, where present, the
    epoch-count line) unchanged except the count, then the first `keep` epochs unchanged.
    Returns (text, kept, total, how the count was rewritten) or None."""
    if not is_text(body):
        return None
    lines = body.decode("latin-1").splitlines(keepends=True)
    rule = next((i for i, ln in enumerate(lines) if ln.startswith("# ====")), None)
    if rule is None:
        return None
    head, data = lines[:rule + 1], lines[rule + 1:]
    how = "none"
    if data and re.fullmatch(r"\s*\d+\s*", data[0]):
        count, data = data[0], data[1:]
        kept = data[:keep]
        head.append(re.sub(r"\d+", str(len(kept)), count, count=1))
        how = "count line after the rule"
    else:
        kept = data[:keep]
        for i, ln in enumerate(head):
            m = re.fullmatch(r"(\s*-?\d+\s+\d+\s+)(\d+)(\s*)", ln)
            if m and int(m.group(2)) == len(data):
                head[i] = m.group(1) + str(len(kept)).rjust(len(m.group(2))) + m.group(3)
                how = f"third integer of header line {i + 1}"
    return "".join(head + kept), len(kept), len(data), how


def target_itsg(out, raw, rec):
    rec.h("ITSG operational products: attitude, non-conservative forces, orbits (N-12, O-19)")
    for product in ITSG_PRODUCTS:
        for sat, day in ITSG_DAYS:
            folder = f"{ITSG}{sat}/{product}/{day[:4]}/"
            entries = listing(folder, rec)
            if entries is None:
                continue
            hits = [e for e in entries if day in e[0] and not e[1]]
            rec.say(f"`{folder}`: {len(entries)} entries; for {day}: "
                    + (", ".join(f"`{e[0]}`" for e in hits) or "none"))
            for name, _, url in hits[:2]:
                status, b, _, err = get(url, max_bytes=500_000_000)
                rec.item(url=url, status=status or err, bytes=len(b or b""), sha256=sha(b) if b else "")
                if not b:
                    continue
                save_raw(raw, f"itsg/{sat}/{product}/{name}", b)
                for m, body in unpack_any(b, name):
                    rec.say(f"  - `{m}`: {len(body)} bytes, sha256 `{sha(body)}`")
                    save_out(out, f"itsg/{sat}/{product}/{m}.head.txt", head_sample(body), rec)
                    cut = cut_groops(body, GROOPS_EXCERPT)
                    if cut is None:
                        rec.say(f"  - `{m}`: not a GROOPS instrument file -> no excerpt")
                        continue
                    text, k, total, how = cut
                    rec.say(f"  - `{m}`: {k} of {total} epochs kept; count rewritten: {how}")
                    save_out(out, f"itsg/{sat}/{product}/{m}", text, rec)


def unpack_any(b, name):
    """(member name, bytes) for a gz/tar/zip archive, or the file itself."""
    try:
        if zipfile.is_zipfile(io.BytesIO(b)):
            z = zipfile.ZipFile(io.BytesIO(b))
            return [(os.path.basename(i.filename), z.read(i)) for i in z.infolist()[:20] if not i.is_dir()]
        if b[:2] == b"\x1f\x8b" or tarfile.is_tarfile(io.BytesIO(b)):
            try:
                t = tarfile.open(fileobj=io.BytesIO(b), mode="r:*")
                return [(os.path.basename(m.name), t.extractfile(m).read())
                        for m in t.getmembers()[:20] if m.isfile()]
            except tarfile.ReadError:
                import gzip
                return [(re.sub(r"\.gz$", "", name), gzip.decompress(b))]
    except Exception as e:  # keep going: the listing and the digest are still recorded
        return [(name + ".unreadable", repr(e).encode())]
    return [(name, b)]


# ---------------------------------------------------------------------------- classic (N-08)
HEADER_END = re.compile(r"^(# End of YAML header|END OF HEADER)", re.M)
COUNT_RES = [re.compile(r"(num_records:\s*)(\d+)"),
             re.compile(r"(NUMBER OF DATA RECORDS\s*:\s*)(\d+)")]


def cut_l1b(body, keep):
    """Header unchanged except the record count, then the first `keep` data lines unchanged.
    Returns (text, records kept, records in file) or None for a binary product."""
    if not is_text(body):
        return None
    text = body.decode("latin-1")
    m = HEADER_END.search(text)
    if not m:
        return None
    split = text.index("\n", m.end()) + 1
    head, data = text[:split], text[split:].splitlines(keepends=True)
    kept = data[:keep]
    for rx in COUNT_RES:
        head = rx.sub(lambda mm: mm.group(1) + str(len(kept)), head, count=1)
    return head + "".join(kept), len(kept), len(data)


def target_classic(out, raw, rec):
    rec.h(f"Classic GRACE Level-1B, {CLASSIC_DAY} (N-08)")
    # Run 40358650063: the HTTPS listing timed out at 90 s. Wait longer, then try ISDC's FTP.
    tops = []
    for base in ("https://isdc-data.gfz.de/grace/Level-1B/",
                 "ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/"):
        tops = listing(base, rec, timeout=300) or []
        if tops:
            break
    rec.say(f"`{base}`: " + ", ".join(f"`{e[0]}`" for e in tops))
    releases = []
    for centre, cisdir, curl in tops:
        if not cisdir:
            continue
        inst = listing(curl + "INSTRUMENT/", rec, timeout=300)
        if inst is None:
            continue
        rec.say(f"`{curl}INSTRUMENT/`: " + ", ".join(f"`{e[0]}`" for e in inst))
        releases += [(centre.strip("/"), e[0].strip("/"), e[2]) for e in inst if e[1]]
    for centre, rl, rurl in releases:
        year = listing(rurl + CLASSIC_DAY[:4] + "/", rec, timeout=300)
        if year is None:
            continue
        hits = [e for e in year if CLASSIC_DAY in e[0] and not e[1]]
        rec.say(f"`{rurl}{CLASSIC_DAY[:4]}/`: {len(year)} entries; for {CLASSIC_DAY}: "
                + (", ".join(f"`{e[0]}`" for e in hits) or "none"))
        for name, _, url in hits[:2]:
            status, b, _, err = get(url, max_bytes=1_500_000_000, timeout=600)
            rec.item(url=url, status=status or err, bytes=len(b or b""), sha256=sha(b) if b else "")
            if not b:
                continue
            save_raw(raw, f"classic/{centre}_{rl}/{name}", b)
            members = unpack_all(b, name)
            rec.say(f"  - {len(members)} members: "
                    + ", ".join(f"`{m}` ({len(x)} B{'' if is_text(x) else ', binary'})" for m, x in members))
            save_out(out, f"classic/{centre}_{rl}/MEMBERS.txt",
                     "".join(f"{len(x):>12}  {sha(x)}  {'text  ' if is_text(x) else 'binary'}  {m}\n"
                             for m, x in members), rec)
            for m, x in members:
                if not re.match(r"(%s)_%s_[AB]_" % ("|".join(CLASSIC_PRODUCTS), CLASSIC_DAY), m):
                    continue
                small = m.startswith(("THR1B", "MAS1B"))
                cut = cut_l1b(x, 10 ** 9 if small else CLASSIC_EXCERPT)
                if cut is None:
                    hdr = x[:x.find(b"END OF HEADER") + 200] if b"END OF HEADER" in x[:200000] else x[:4096]
                    rec.say(f"  - `{m}`: binary or no header end found -> header only committed")
                    save_out(out, f"classic/{centre}_{rl}/{m}.header.txt", hdr.decode("latin-1"), rec)
                    continue
                text, k, total = cut
                rec.say(f"  - `{m}`: {k} of {total} records kept")
                save_out(out, f"classic/{centre}_{rl}/{m}", text, rec)


def unpack_all(b, name):
    """Every regular member of a (possibly gzipped) tar, recursing one level into nested
    archives (daily L1B tarballs hold per-product gzips on some releases)."""
    first = unpack_any(b, name) if b[:2] == b"\x1f\x8b" or tarfile.is_tarfile(io.BytesIO(b)) \
        or zipfile.is_zipfile(io.BytesIO(b)) else [(name, b)]
    if b[:2] == b"\x1f\x8b" or tarfile.is_tarfile(io.BytesIO(b)):
        try:
            t = tarfile.open(fileobj=io.BytesIO(b), mode="r:*")
            first = [(os.path.basename(m.name), t.extractfile(m).read()) for m in t.getmembers() if m.isfile()]
        except tarfile.ReadError:
            pass
    out = []
    for m, x in first:
        if x[:2] == b"\x1f\x8b" or (len(x) > 512 and tarfile.is_tarfile(io.BytesIO(x))):
            out += unpack_any(x, m)
        else:
            out.append((m, x))
    return out


# ---------------------------------------------------------------------------- main
TARGETS = {"dst": target_dst, "f30": target_f30, "itsg": target_itsg, "classic": target_classic}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="data/incoming")
    ap.add_argument("--raw", default="fetch_raw")
    ap.add_argument("--only", nargs="*", choices=sorted(TARGETS), default=sorted(TARGETS))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.raw, exist_ok=True)
    rec = Record()
    when = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    rec.say(f"# fetchData probe report, {when}")
    rec.say("")
    rec.say("Written by `tools/fetchData.py` on a GitHub-hosted runner. Staging only: see the "
            "script's docstring. Each line is a request and what came back.")
    for name in a.only:
        try:
            TARGETS[name](a.out, a.raw, rec)
        except Exception as e:  # one target's failure must not hide the others' results
            rec.say(f"**target `{name}` raised {type(e).__name__}: {e}**")
    with open(os.path.join(a.out, "PROBE.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(rec.lines) + "\n")
    with open(os.path.join(a.out, "FETCHED.json"), "w", encoding="utf-8") as f:
        json.dump({"generated": when, "items": rec.items}, f, indent=1)
    print("\n".join(rec.lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
