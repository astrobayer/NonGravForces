#!/usr/bin/env python3
"""Static audit of NonGravForces (stdlib only). Exit 1 on any finding.

Mechanisms, not reminders (HANDOVER P-5). Each check names the rule it enforces.
  1. doc/functionlist.csv <-> files in +NonGravForces agree both ways (R1, R3).
  2. Every public unit carries the help template: H1, INPUT/OUTPUT (functions),
     ERRORS, KEYWORDS, Created line with the authoring model (E.3).
  3. Every NonGravForces.fail call uses a taxonomy mnemonic (E.2).
  4. Package code contains no inputParser, addpath, codegen/coder, and no
     printing except in units on the allow-list (E.3, C1).
  5. Every copy of the version string (Contents.m "Version" lines, README.md
     first line, doc/html footers) equals the one composed from the numeric
     fields of toolboxVersion.m (R1, decision 2026-09-06).
  6. EVERY page of doc/html equals what tools/buildDocs.py renders from the
     sources that own the facts: staleness by content, not by stamp (O-11).
"""
import hashlib
import csv, glob, os, re, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(ROOT, "+NonGravForces")
TAXONOMY = ("input", "dep", "data", "numeric", "unsupported", "notImplemented")
PRINT_ALLOWED = {"dependencyReport.m"}
findings = []

def note(msg):
    findings.append(msg)

# 1. catalogue <-> files
rows = list(csv.DictReader(open(os.path.join(ROOT, "doc", "functionlist.csv"), encoding="utf-8")))
catalogued = set()
for r in rows:
    sub = "" if r["package"] == "NonGravForces" else "+" + r["package"]
    path = os.path.join(PKG, sub, r["name"] + ".m")
    catalogued.add(os.path.normpath(path))
    if not os.path.isfile(path):
        note(f"catalogue row without file: {r['package']}.{r['name']}")
on_disk = set()
private = set()
for dp, dn, fn in os.walk(PKG):
    for f in fn:
        if not f.endswith(".m") or f == "Contents.m":
            continue
        p = os.path.normpath(os.path.join(dp, f))
        if os.path.basename(dp) == "private":
            private.add(p)
        else:
            on_disk.add(p)
for p in sorted(on_disk - catalogued):
    note(f"file without catalogue row: {os.path.relpath(p, ROOT)}")

# 2. help template, 3. fail taxonomy, 4. banned constructs (private helpers: 3 and 4 only)
for p in sorted(on_disk | private):
    rel = os.path.relpath(p, ROOT)
    src = open(p, encoding="utf-8").read()
    lines = src.splitlines()
    is_class = lines[0].lstrip().startswith("classdef")
    head = "\n".join(lines[:80])
    if p in on_disk:
        if not re.search(r"^%[A-Z][A-Z0-9_]+ ", head, re.M):
            note(f"{rel}: no H1 line (%NAME description)")
        if not re.search(r"^% KEYWORDS ", head, re.M) and "KEYWORDS" not in src:
            note(f"{rel}: no KEYWORDS line")
        if not re.search(r"Created \d{4}-\d{2}-\d{2} by ", src):
            note(f"{rel}: no 'Created <date> by <model>' line")
        if not is_class:
            if "INPUT ARGUMENTS" not in src and "NAME-VALUE" not in src and "OUTPUT ARGUMENTS" not in src:
                note(f"{rel}: function help lacks INPUT/OUTPUT ARGUMENTS blocks")
        if "ERRORS" not in src and "notImplemented" not in src:
            note(f"{rel}: help lacks an ERRORS block")
    for m in re.finditer(r'NonGravForces\.fail\(\s*"[^"]+"\s*,\s*"([^"]+)"', src):
        if not m.group(1).startswith(TAXONOMY):
            note(f"{rel}: fail mnemonic '{m.group(1)}' outside the taxonomy")
    for bad in ("inputParser", "addpath(", "codegen(", "coder."):
        if bad in src:
            note(f"{rel}: uses banned construct '{bad}'")
    if os.path.basename(p) not in PRINT_ALLOWED:
        code = "\n".join(l for l in lines if not l.lstrip().startswith("%"))
        if re.search(r"\b(disp|fprintf)\s*\(", code):
            note(f"{rel}: library code prints (disp/fprintf) - not on the allow-list")

# 5. one version, composed from the numeric fields, every copy compared
vsrc = open(os.path.join(PKG, "toolboxVersion.m"), encoding="utf-8").read()
m = re.search(r"major\s*=\s*uint16\((\d+)\);\s*minor\s*=\s*uint16\((\d+)\);\s*patch\s*=\s*uint16\((\d+)\);", vsrc)
if not m:
    note("toolboxVersion.m: numeric fields major/minor/patch not found in the expected form")
    version = None
else:
    version = "%s.%s.%s" % m.groups()
    for lit in re.findall(r"'semantic',\s*\"(\d+\.\d+\.\d+)\"", vsrc):
        note(f"toolboxVersion.m: semantic is a literal '{lit}' - it must be composed from major/minor/patch")
    for c in sorted(glob.glob(os.path.join(PKG, "**", "Contents.m"), recursive=True)):
        txt = open(c, encoding="utf-8").read()
        vm = re.search(r"^% Version (\S+)", txt, re.M)
        if not vm:
            note(f"{os.path.relpath(c, ROOT)}: no '% Version x.y.z' line")
        elif vm.group(1) != version:
            note(f"{os.path.relpath(c, ROOT)}: Version {vm.group(1)} != toolboxVersion {version}")
    readme = open(os.path.join(ROOT, "README.md"), encoding="utf-8").readline().strip()
    if readme != f"# NonGravForces {version}":
        note(f"README.md: first line '{readme}' does not name version {version}")
    for h in sorted(glob.glob(os.path.join(ROOT, "doc", "html", "*.html"))):
        txt = open(h, encoding="utf-8").read()
        for v in set(re.findall(r"NonGravForces (\d+\.\d+\.\d+)", txt)):
            if v != version:
                note(f"{os.path.relpath(h, ROOT)}: mentions version {v}, toolbox is {version}")
    # CLAUDE.md is an orientation file, not an authority, and it carries two numbers that are
    # authored elsewhere: the version and the unit count. It claimed since 2026-09-06 that "CI
    # checks ... this sentence's number against it" and CI did NOT - the sentence said 53 while
    # the catalogue held 55 (LP-C). Claiming a check is not having one (R1).
    claude = open(os.path.join(ROOT, "CLAUDE.md"), encoding="utf-8").read()
    cv = re.search(r"Toolbox version (\d+\.\d+\.\d+)", claude)
    if not cv:
        note("CLAUDE.md: no 'Toolbox version x.y.z' in the header")
    elif cv.group(1) != version:
        note(f"CLAUDE.md: Toolbox version {cv.group(1)} != toolboxVersion {version}")
    cn = re.search(r"^(\d+) public units, catalogued in `doc/functionlist.csv`", claude, re.M)
    if not cn:
        note("CLAUDE.md: no '<n> public units, catalogued in `doc/functionlist.csv`' sentence")
    elif int(cn.group(1)) != len(rows):
        note(f"CLAUDE.md: claims {cn.group(1)} public units, the catalogue holds {len(rows)}")

# 6. Every doc/html page is what the generator renders. Until 2026-09-09 only FunctionList.html
#    was covered and the rest carried a hand-edited footer checked for its version string alone:
#    by then Sensors.html credited a class with a filter it deliberately does not build. A page
#    nobody generates is a page that lies eventually (O-11, WP-21).
sys.path.insert(0, os.path.join(ROOT, "tools"))
try:
    import buildDocs
    for name, text in buildDocs.render_all().items():
        path = os.path.join(ROOT, "doc", "html", name)
        current = open(path, encoding="utf-8").read() if os.path.isfile(path) else ""
        if current != text:
            note(f"doc/html/{name} is stale -> python tools/buildDocs.py")
except Exception as exc:                      # a broken generator must be a finding, not a crash
    note(f"tools/buildDocs.py could not render the documentation: {exc}")

# 7. the catalogue must survive a naive CSV read: an unquoted comma inside an h1 silently
#    truncated env.eclipseFactor's purpose into Contents.m and FunctionList.html until 2026-09-07.
for i, raw in enumerate(open(os.path.join(ROOT, "doc", "functionlist.csv"), encoding="utf-8").read().splitlines(), 1):
    if raw and raw.count(",") != 6:
        note(f"doc/functionlist.csv line {i} has {raw.count(',') + 1} comma-separated fields, not 7 "
             f"-> remove the comma from the h1 text (the catalogue is read unquoted)")

# 8. NRLMSISE-00 is no longer vendored (O-34), so there is nothing here to digest. Its integrity
#    is covered by two other mechanisms: setupNonGravForces pins the upstream commit, and
#    buildMex digests whatever source it actually compiled, so a different release marks the
#    binary stale. This check refuses only a PARTIAL tree, which is the one state that would
#    compile into something nobody chose.
mexdir = os.path.join(ROOT, "mex", "nrlmsise00")
if os.path.isdir(mexdir):
    need = ["nrlmsise-00.c", "nrlmsise-00.h", "nrlmsise-00_data.c"]
    absent = [f for f in need if not os.path.isfile(os.path.join(mexdir, f))]
    if absent and len(absent) != len(need):
        note(f"mex/nrlmsise00 exists but is incomplete (missing {', '.join(absent)}) "
             f"-> delete the folder and let setupNonGravForces clone it again")

# 9. NO TOOLBOX CALLS. Constraint C2 allows base MATLAB and the Parallel Computing Toolbox and
#    nothing else, and that constraint has been broken twice by hand in one day: rms() reached
#    both a docstring EXAMPLE and a test before anyone noticed. A reader who copies an example
#    that needs a product the toolbox refuses to need has been misled by the toolbox itself, so
#    the scan covers comments and help text as well as code.
#
#    The list is deliberately CONSERVATIVE: only names that are unambiguous calls into a
#    product. fit, smooth, filter and median are left out because they collide with ordinary
#    variable names and a false positive here costs more than a missed exotic call -- a check
#    people learn to work around is worse than no check.
#
#    A NAME THAT THE FILE ASSIGNS TO IS A VARIABLE, NOT A CALL. Model.m and forceBudget.m both
#    keep a matrix called rms and index it; flagging those would have been the scan teaching
#    people to ignore it on its first run. So each file's assignment targets are collected first
#    and excluded within that file, which is sound: MATLAB cannot mean both at once in one scope.
#
#    A line that legitimately NAMES one of these (this comment, for instance) carries the marker
#    toolboxScan:allow, so every exception is visible in the diff that introduced it.
TOOLBOX_FUNCTIONS = {
    "rms": "Signal Processing", "filtfilt": "Signal Processing", "butter": "Signal Processing",
    "freqz": "Signal Processing", "periodogram": "Signal Processing", "pwelch": "Signal Processing",
    "findpeaks": "Signal Processing", "resample": "Signal Processing", "medfilt1": "Signal Processing",
    "sgolayfilt": "Signal Processing", "spectrogram": "Signal Processing", "hann": "Signal Processing",
    "hamming": "Signal Processing", "blackman": "Signal Processing", "xcorr": "Signal Processing",
    "normpdf": "Statistics", "normcdf": "Statistics", "prctile": "Statistics", "quantile": "Statistics",
    "nanmean": "Statistics", "nanstd": "Statistics", "regress": "Statistics", "fitlm": "Statistics",
    "robustfit": "Statistics", "ttest": "Statistics", "mvnrnd": "Statistics", "boxplot": "Statistics",
    "fmincon": "Optimization", "fminunc": "Optimization", "lsqnonlin": "Optimization",
    "lsqcurvefit": "Optimization", "fsolve": "Optimization", "quadprog": "Optimization",
    "csaps": "Curve Fitting", "cfit": "Curve Fitting",
    "sym": "Symbolic Math", "syms": "Symbolic Math", "vpa": "Symbolic Math",
    "dcm2quat": "Aerospace", "quatmultiply": "Aerospace", "quatrotate": "Aerospace",
    "geoshow": "Mapping", "wgs84Ellipsoid": "Mapping", "geodetic2ecef": "Mapping",
}
scan_roots = [os.path.join(ROOT, "+NonGravForces"), os.path.join(ROOT, "tests"),
              os.path.join(ROOT, "tools"), os.path.join(ROOT, "demos")]
scan_files = [os.path.join(ROOT, f) for f in ("setupNonGravForces.m", "ciGate.m", "runNgfTests.m")
              if os.path.isfile(os.path.join(ROOT, f))]
for r in scan_roots:
    for dirpath, _, filenames in os.walk(r):
        scan_files += [os.path.join(dirpath, f) for f in filenames if f.endswith(".m")]
pattern = re.compile(r"(?<![\w.])(" + "|".join(sorted(TOOLBOX_FUNCTIONS)) + r")\s*\(")
for f in sorted(scan_files):
    rel = os.path.relpath(f, ROOT).replace(os.sep, "/")
    text = open(f, encoding="utf-8").read()
    assigned = set(re.findall(r"(?<![\w.])(" + "|".join(sorted(TOOLBOX_FUNCTIONS)) +
                              r")\s*(?:\([^)\n]*\))?\s*=(?!=)", text))
    for n, line in enumerate(text.splitlines(), 1):
        if "toolboxScan:allow" in line:
            continue
        m = pattern.search(line)
        if m and m.group(1) not in assigned:
            note(f"{rel}:{n} calls {m.group(1)}(), which needs the {TOOLBOX_FUNCTIONS[m.group(1)]} "
                 f"Toolbox -> constraint C2 allows base MATLAB and Parallel Computing only; "
                 f"write it out, or mark the line toolboxScan:allow if it only names the function")

if findings:
    print(f"auditPackage: {len(findings)} finding(s)")
    for f in findings:
        print("  -", f)
    sys.exit(1)
print(f"auditPackage: clean ({len(on_disk)} units, {len(private)} private helpers, {len(rows)} catalogue rows, version {version})")
