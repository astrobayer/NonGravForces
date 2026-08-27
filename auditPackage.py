#!/usr/bin/env python3
"""Static audit of NonGravForces (stdlib only). Exit 1 on any finding.

Mechanisms, not reminders (HANDOVER P-5). Each check names the rule it enforces.
  1. doc/functionlist.csv <-> files in +NonGravForces agree both ways (R1, R3).
  2. Every public unit carries the help template: H1, INPUT/OUTPUT (functions),
     ERRORS, KEYWORDS, Created line with the authoring model (E.3).
  3. Every NonGravForces.fail call uses a taxonomy mnemonic (E.2).
  4. Package code contains no inputParser, addpath, codegen/coder, and no
     printing except in units on the allow-list (E.3, C1).
"""
import csv, os, re, sys

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
for dp, dn, fn in os.walk(PKG):
    if os.path.basename(dp) == "private":
        continue
    for f in fn:
        if f.endswith(".m") and f != "Contents.m":
            on_disk.add(os.path.normpath(os.path.join(dp, f)))
for p in sorted(on_disk - catalogued):
    note(f"file without catalogue row: {os.path.relpath(p, ROOT)}")

# 2. help template, 3. fail taxonomy, 4. banned constructs
for p in sorted(on_disk):
    rel = os.path.relpath(p, ROOT)
    src = open(p, encoding="utf-8").read()
    lines = src.splitlines()
    is_class = lines[0].lstrip().startswith("classdef")
    head = "\n".join(lines[:80])
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

if findings:
    print(f"auditPackage: {len(findings)} finding(s)")
    for f in findings:
        print("  -", f)
    sys.exit(1)
print(f"auditPackage: clean ({len(on_disk)} units, {len(rows)} catalogue rows)")
