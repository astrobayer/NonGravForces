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
  6. doc/html/FunctionList.html equals what tools/buildFunctionList.py renders
     from the catalogue: staleness by content, not by stamp (O-11).
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

# 6. FunctionList.html is what the generator renders
sys.path.insert(0, os.path.join(ROOT, "tools"))
try:
    import buildFunctionList
    rendered = buildFunctionList.render()
    fl = os.path.join(ROOT, "doc", "html", "FunctionList.html")
    current = open(fl, encoding="utf-8").read() if os.path.isfile(fl) else ""
    if current != rendered:
        note("doc/html/FunctionList.html is stale -> python tools/buildFunctionList.py")
except SystemExit as e:
    note(f"tools/buildFunctionList.py: {e}")

# 6b. Only base MATLAB and the Parallel Computing Toolbox may appear on any code path
#     (HANDOVER D.7 C2). The owner's host has Signal Processing and Statistics installed and the
#     CI runner does not, so a call to one of them passes locally and fails on CI: that is exactly
#     what `corr` in a test did on 2026-09-08. Names below are functions with no base equivalent.
TOOLBOX_ONLY = {
    "corr": "Statistics", "partialcorr": "Statistics", "nanmean": "Statistics",
    "nanstd": "Statistics", "nanmedian": "Statistics", "normpdf": "Statistics",
    "normcdf": "Statistics", "norminv": "Statistics", "ttest": "Statistics",
    "ttest2": "Statistics", "regress": "Statistics", "fitlm": "Statistics",
    "pca": "Statistics", "kmeans": "Statistics", "prctile": "Statistics",
    "zscore": "Statistics", "mvnrnd": "Statistics", "boxplot": "Statistics",
    "filtfilt": "Signal", "butter": "Signal", "cheby1": "Signal", "freqz": "Signal",
    "periodogram": "Signal", "pwelch": "Signal", "findpeaks": "Signal",
    "spectrogram": "Signal", "resample": "Signal", "designfilt": "Signal",
    "xcorr": "Signal", "hann": "Signal", "hamming": "Signal",
}
for rel in sorted(glob.glob(os.path.join(ROOT, "+NonGravForces", "**", "*.m"), recursive=True)
                  + glob.glob(os.path.join(ROOT, "tests", "*.m"))
                  + [os.path.join(ROOT, "ciGate.m")]):
    src = open(rel, encoding="utf-8").read()
    body = "\n".join(ln for ln in src.splitlines() if not ln.lstrip().startswith("%"))
    for fn, box in TOOLBOX_ONLY.items():
        if re.search(r"(?<![\w.])" + fn + r"\s*\(", body):
            note(f"{os.path.relpath(rel, ROOT)} calls {fn}(), which is the {box} toolbox -> "
                 f"HANDOVER D.7 C2 allows base MATLAB and PCT only; it is installed on the owner's "
                 f"host and absent from the CI runner")

# 7. the catalogue must survive a naive CSV read: an unquoted comma inside an h1 silently
#    truncated env.eclipseFactor's purpose into Contents.m and FunctionList.html until 2026-09-07.
for i, raw in enumerate(open(os.path.join(ROOT, "doc", "functionlist.csv"), encoding="utf-8").read().splitlines(), 1):
    if raw and raw.count(",") != 6:
        note(f"doc/functionlist.csv line {i} has {raw.count(',') + 1} comma-separated fields, not 7 "
             f"-> remove the comma from the h1 text (the catalogue is read unquoted)")

# 8. the vendored NRLMSISE-00 sources must match the digests their README records, or a silent
#    edit to third-party C would ship unnoticed. The digests were recorded when it was vendored.
mexdir = os.path.join(ROOT, "mex", "nrlmsise00")
readme = os.path.join(mexdir, "README.md")
if os.path.isfile(readme):
    recorded = dict(re.findall(r"\|\s*`([\w.\-]+)`\s*\|\s*([0-9a-f]{64})\s*\|", open(readme, encoding="utf-8").read()))
    if not recorded:
        note("mex/nrlmsise00/README.md records no SHA-256 digests -> the vendored sources are unverifiable")
    for fname, want in recorded.items():
        fpath = os.path.join(mexdir, fname)
        if not os.path.isfile(fpath):
            note(f"mex/nrlmsise00/{fname} is recorded with a digest but missing")
            continue
        got = hashlib.sha256(open(fpath, "rb").read()).hexdigest()
        if got != want:
            note(f"mex/nrlmsise00/{fname} digest {got[:16]} does not match the recorded {want[:16]} "
                 f"-> the vendored third-party source was edited")

if findings:
    print(f"auditPackage: {len(findings)} finding(s)")
    for f in findings:
        print("  -", f)
    sys.exit(1)
print(f"auditPackage: clean ({len(on_disk)} units, {len(private)} private helpers, {len(rows)} catalogue rows, version {version})")
