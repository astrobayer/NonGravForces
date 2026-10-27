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
  6. EVERY page of doc/html, doc/src/generated.tex and doc/src/api_generated.tex
     equal what tools/buildDocs.py renders from the sources that own the facts:
     staleness by content, not by stamp (O-11).
"""
import hashlib
import csv, glob, os, re, sys
# the argument/option auditors of check 12 live in tools/, so that folder must be importable
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools'))

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
    # HANDOVER.md's title names the version too, and nothing compared it: it said 0.11.0 through
    # the whole of 0.12.0 until release 0.13.0 found it (2026-10-15).
    hand = open(os.path.join(ROOT, "HANDOVER.md"), encoding="utf-8").readline().strip()
    hm = re.match(r"# HANDOVER — NonGravForces (\S+) ", hand)
    if not hm or hm.group(1) != version:
        note(f"HANDOVER.md: title '{hand}' does not name version {version}")
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
    # The README carried "51 of 53 catalogued units ... 2 documented stubs" beside "the four
    # remaining stubs" -- two numbers that disagreed with each other AND with the catalogue, in the
    # first thing a reader sees. Same rule, same instrument.
    readme_all = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    rn = re.search(r"^Status: (\d+) catalogued units", readme_all, re.M)
    if not rn:
        note("README.md: no 'Status: <n> catalogued units' sentence")
    elif int(rn.group(1)) != len(rows):
        note(f"README.md: claims {rn.group(1)} catalogued units, the catalogue holds {len(rows)}")
    stubs = [r for r in rows if r["status"] == "stub"]
    if stubs and "no stubs left" in readme_all:
        note(f"README.md says 'no stubs left', the catalogue holds {len(stubs)}: "
             + ", ".join(f"{r['package']}.{r['name']}" for r in stubs))

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
    # The theory document's derivable facts (O-38). Its title page said "version 0.1.0" for three
    # weeks and it offered three options the toolbox refuses, because no instrument compared it
    # with anything. Now its version and every enumerated option list come from generated.tex.
    gen = buildDocs.render_generated_tex()
    cur = open(buildDocs.TEX_GENERATED, encoding="utf-8").read() if os.path.isfile(buildDocs.TEX_GENERATED) else ""
    if cur != gen:
        note("doc/src/generated.tex is stale -> python tools/buildDocs.py")
    # The API chapter is the same mechanism pointed at a bigger fact: every signature, argument
    # table and error identifier in the user guide comes from the catalogue and the help blocks,
    # so a renamed argument or a new option makes THIS file stale and reddens the audit. An API
    # reference maintained by hand is the documentation defect O-11 was, at 56 units instead of 6.
    api = buildDocs.render_api_tex()
    cur = open(buildDocs.TEX_API, encoding="utf-8").read() if os.path.isfile(buildDocs.TEX_API) else ""
    if cur != api:
        note("doc/src/api_generated.tex is stale -> python tools/buildDocs.py")
except Exception as exc:                      # a broken generator must be a finding, not a crash
    note(f"tools/buildDocs.py could not render the documentation: {exc}")

# 6a. HANDOVER Part H has ONE status field, and nothing else carries a status. Until 2026-09-17
#     it had none: status lived in the body's prose, in the Owner column or in the Re-entry column,
#     and eleven of forty rows disagreed with themselves -- nine announcing CLOSED while Re-entry
#     still named a live trigger. A summary of "what is open" was wrong because there was nothing
#     to read. The Status column is now the authority and this refuses the drift back.
handover = open(os.path.join(ROOT, "HANDOVER.md"), encoding="utf-8").read()
part_h = handover[handover.index("## Part H"):handover.index("## Part I")]
if "| # | Item | Status | Owner | Re-entry |" not in part_h:
    note("HANDOVER Part H: the obligations table has no Status column -> status would scatter again (O-26)")
obligations = re.findall(r"^\| (O-\d+) \| (.*?) \| (.*?) \| (.*?) \| (.*?) \|\s*$", part_h, re.M)
if not obligations:
    note("HANDOVER Part H: no obligation rows parsed -> the table shape changed and nothing can read it")
statuses = {}
for num, body, status, owner, reentry in obligations:
    st = status.strip().lower()
    statuses[num] = st
    if st not in ("open", "closed", "void"):
        note(f"HANDOVER {num}: status {status.strip()!r} is not open, closed or void")
    # A row that announces its own closure in CAPS must not be marked open, and vice versa.
    announced = re.match(r"\*\*(CLOSED|VOID)\b", body.strip())
    if announced and st == "open":
        note(f"HANDOVER {num}: the item announces {announced.group(1)} but its status is open")
    if st in ("closed", "void") and reentry.strip() not in ("\u2014", "-", ""):
        note(f"HANDOVER {num}: status {st} but Re-entry still names {reentry.strip()[:40]!r}")
    for col, name in ((owner, "Owner"), (reentry, "Re-entry")):
        if re.fullmatch(r"(?i)\s*(closed|void|half closed)\s*", col):
            note(f"HANDOVER {num}: the {name} column carries a status ({col.strip()!r}) -> the Status column is the authority")
# CLAUDE.md quotes the range; it is derivable, so it is compared.
rng = re.search(r"currently O-1 to O-(\d+) with an owner", claude)
if not rng:
    note("CLAUDE.md: no 'currently O-1 to O-N' sentence naming the obligation range")
elif obligations and int(rng.group(1)) != max(int(n[2:]) for n in statuses):
    note(f"CLAUDE.md: says O-1 to O-{rng.group(1)}, Part H holds up to O-{max(int(n[2:]) for n in statuses)}")

# 6b. NO BUILT PDF IN THE REPOSITORY (O-38). A PDF committed beside its own .tex has nothing
#     keeping the two in step, and this one drifted nine minor versions before anyone looked.
#     doc/src is the deliverable; the CI job "theory" builds the PDF and uploads it as an artefact.
for pdf in sorted(glob.glob(os.path.join(ROOT, "doc", "**", "*.pdf"), recursive=True)):
    rel = os.path.relpath(pdf, ROOT).replace(os.sep, "/")
    if rel.startswith("doc/references/"):
        continue                              # third-party papers, not built from anything here
    note(f"{rel} is a built PDF in the repository -> doc/src is the deliverable and CI builds the "
         f"PDF as an artefact; delete it and keep it out (O-38)")
# The .tex must exist and must pull the generated facts in, or the macros above prove nothing.
tex = os.path.join(ROOT, "doc", "src", "NonGravForces_Theory_and_Implementation.tex")
if not os.path.isfile(tex):
    note("doc/src/NonGravForces_Theory_and_Implementation.tex is missing -> it is the deliverable")
else:
    body = open(tex, encoding="utf-8").read()
    # ANCHORED AND UNCOMMENTED. A plain substring test passed a "% \input{generated}" when this
    # check was fault-injected on 2026-09-16 -- a check that its own injection does not fire is
    # not a check, and finding that out is the only reason to inject one.
    if not re.search(r"^[ \t]*\\input\{generated\}", body, re.M):
        note("the theory document does not \\input{generated} on a live line -> its version and "
             "option lists would be literals again (O-38)")
    for literal in re.findall(r"version (\d+\.\d+\.\d+)", body):
        note(f"the theory document writes version {literal} as a literal -> use \\ngfversion")

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

# 10. THE NATIVE SETUP IS ONE FACT IN THREE PLACES: the C struct `ngfRhsSetup` declares a field,
#     app.DynamicsSource packs it into a MATLAB struct, and mex/ngfRhsGateway.c copies it across BY
#     NAME. Nothing compared them until 2026-09-17, when `eV` and `torr` were added to the first two
#     and not the third: they reached C as uninitialised memory, SESAM's Langmuir constant became
#     garbage, and the native path evaluated a different model from the MATLAB one. The parity test
#     caught it at 6.5e-02 against a 1e-13 bound, three minutes into a gate; this catches it in a
#     second, which is the difference between a check and a post-mortem.
#     Only the scalar constants are compared: they are the block that is declared as one comma list
#     and is therefore the easy one to extend without noticing the gateway.
try:
    hdr = open(os.path.join(ROOT, "mex", "ngfRhs.h"), encoding="utf-8").read()
    gw = open(os.path.join(ROOT, "mex", "ngfRhsGateway.c"), encoding="utf-8").read()
    block = re.search(r"^\s*double cLight,([^;]*);", hdr, re.M)
    if block is None:
        note("mex/ngfRhs.h no longer declares the scalar constant block starting at cLight -> "
             "check 10 cannot compare it with the gateway; re-point it")
    else:
        declared = ["cLight"] + [f.strip() for f in block.group(1).split(",") if f.strip()]
        copied = set(re.findall(r'scalar\(st,\s*"(\w+)"\)', gw))
        missing = [f for f in declared if f not in copied]
        if missing:
            note("mex/ngfRhsGateway.c never copies %s, which mex/ngfRhs.h declares -> the field "
                 "reaches C as uninitialised memory and the native path silently evaluates "
                 "something else; add scalar(st, \"<name>\") beside the others"
                 % ", ".join(missing))
except OSError as exc:
    note(f"check 10 could not read the native setup sources: {exc}")

# 11. EVERY VERSION STRING AND TEST COUNT IN THE README, not only the ones check 5 looks at.
#     Check 5 compares the README's FIRST LINE with toolboxVersion and the unit count with the
#     catalogue, and stopped there. The status paragraph carried a SECOND version string and a test
#     count -- "LP-C is complete: 365/365 green ... version 0.10.0" -- which survived LP-D, LP-E and
#     two version bumps, because the instrument that read the first line had never been pointed at
#     the rest of the file. Found 2026-09-17 while cutting the release it describes.
#     The rule that produced it is this level's own: A FACT WRITTEN DOWN TWICE AND COMPARED NEVER.
#     So every MAJOR.MINOR.PATCH in the README must be the current version, and every "N/M green"
#     must match what tests/last_run.txt actually recorded.
readme = os.path.join(ROOT, "README.md")
if os.path.isfile(readme):
    body = open(readme, encoding="utf-8").read()
    #     A line that deliberately QUOTES a historical version -- the PDF that said 0.1.0 for
    #     three weeks, say -- carries versionScan:allow, so every exception is visible in the
    #     diff that introduced it. Same convention as check 9's toolboxScan:allow, and for the
    #     same reason: a check people learn to work around is worse than no check.
    live = "\n".join(ln for ln in body.splitlines() if "versionScan:allow" not in ln)
    stale = sorted({v for v in re.findall(r"\b\d+\.\d+\.\d+\b", live) if v != version})
    if stale:
        note("README.md names version(s) %s but the toolbox is %s -> one authority per fact; "
             "toolboxVersion.m is it" % (", ".join(stale), version))
    run = os.path.join(ROOT, "tests", "last_run.txt")
    if os.path.isfile(run):
        head = open(run, encoding="utf-8").read(400)
        m = re.search(r"tests\s+(\d+)\s+passed\s+(\d+)", head)
        for claimed in re.findall(r"(\d+)\s*/\s*(\d+)\s+green", body, re.I):
            if m and (claimed[0] != m.group(2) or claimed[1] != m.group(1)):
                note("README.md claims %s/%s green but tests/last_run.txt recorded %s passed of %s "
                     "-> quote the run, or re-run before quoting"
                     % (claimed[0], claimed[1], m.group(2), m.group(1)))

# 12. EVERY DECLARED ARGUMENT IS A DOCUMENTED ARGUMENT, and every option value is explained.
#     Check 2 verifies the INPUT/NAME-VALUE/OUTPUT sections EXIST; it never compared them with what
#     the functions actually take. Run 2026-09-18, that comparison found: four PUBLIC METHODS named
#     nowhere in their class help (Options.toStruct, DynamicsSource.diagnostics and epochAt,
#     MacroModel.silhouetteDefect); three MacroModel name-value arguments whose PROPERTIES were
#     documented while the way to SET them was not; seven Model method name-values a caller can pass
#     and the single point of access described nowhere (torque/evaluate B_body and Dipole_Am2,
#     evaluate T0_K, series Quantities, partials Method); and four qc.screen options carrying no
#     size, type or units. The same shape as everything else this level found.
#     The three auditors live in tools/ and run standalone; this runs them so CI reads them.
for _mod in ("audit_args", "audit_class", "audit_options"):
    try:
        _m = __import__(_mod)
        import io as _io
        import contextlib as _ctx
        _buf = _io.StringIO()
        with _ctx.redirect_stdout(_buf):
            _rc = _m.main()
        if _rc:
            for _line in _buf.getvalue().splitlines():
                if _line.strip() and not _line.startswith(("0 finding", "\n")) and "finding(s) over" not in _line:
                    note(f"{_mod}: {_line.strip()}")
    except Exception as _exc:                      # a broken auditor is a finding, not a crash
        note(f"tools/{_mod}.py could not run: {_exc}")

# 13. EVERY PER-EPOCH COST QUOTED ANYWHERE IS A COST SOMETHING MEASURED.
#     O-29 recorded 88 ms per epoch on 2026-09-08. WP-C3 fixed it to 1.21 ms on 2026-09-16 and
#     recorded that too, 772 lines away in the same file, and nobody closed the obligation row. The
#     row therefore went on instructing readers to repeat 88 ms, and four more files plus one
#     shipped PDF chapter did -- while tests/tNgfModel.m had held the correct value, in an executed
#     test, for three days. Checking that a fact appears ONCE would not have caught this: the second
#     copy was an obligation to change the first.
#     So this compares the numbers themselves. Every "<x> ms per epoch" or "<x> ms each" claim in the
#     repository must be within 25 % of some perEpoch_ms in tools/benchmark_reference.csv, which
#     tools/benchmark.m writes and which is the single authority for this fact.
#     THE EXEMPTION IS DELIBERATELY NOT "the line carries a date". The original O-29 row read
#     "88 ms per epoch in the standard preset (measured 2026-09-08)" -- it WAS dated, and it
#     propagated into six files anyway, because a date says when a number was taken and not whether
#     it is still true. A check that exempted dated lines would not have caught the defect it exists
#     for. The exemption is instead an explicit statement, on the same line, that the figure is no
#     longer current: this project keeps wrong numbers on purpose, as warnings, and that is the only
#     reason to carry one.
#     WHAT IT DOES NOT CHECK: whether benchmark_reference.csv is itself current, nor that it was
#     measured well. Its header carries a date and the code moves; and on 2026-09-18 the first run
#     of the repaired instrument was still wrong, because one warm-up call is not warm (D.9).
#     Re-run tools/benchmark.m when a number here is meant to carry weight.
#     THE AUTHORITY IS tools/benchmark_reference.csv, WHICH IS COMMITTED -- not benchmark_last.txt,
#     which .gitignore excludes. That distinction is the whole check: the ignored log exists only on
#     the machine that wrote it, so on the CI audit job and in every fresh clone it is absent, and a
#     check whose authority is missing everywhere but one laptop protects nothing. (Caught before
#     this check was ever pushed, by asking what the runner does not have that this host does.)
BENCH = os.path.join(ROOT, "tools", "benchmark_reference.csv")
if not os.path.isfile(BENCH):
    note("tools/benchmark_reference.csv is missing -> run tools/benchmark.m and commit it; it is "
         "the authority for the per-epoch costs this repository quotes")
else:
    bench_text = open(BENCH, encoding="utf-8").read()
    measured = []
    for line in bench_text.splitlines():
        if line.startswith("#") or line.startswith("path,"):
            continue
        # rsplit from the RIGHT: the path column legitimately contains commas
        # ("Model.accel, fast preset"), and only the three numeric columns are fixed.
        cells = line.rsplit(",", 3)
        if len(cells) == 4 and cells[3].strip() not in ("", "NaN"):
            try:
                measured.append(float(cells[3]))
            except ValueError:
                pass
    if not measured:
        note("tools/benchmark_reference.csv holds no perEpoch_ms value -> it measured nothing")
    # AND THE AUTHORITY ITSELF MUST DESCRIBE THIS CODE. Its header records the toolbox version it
    # was measured against; comparing that with toolboxVersion.m is the difference between "these
    # numbers agree with a table" and "these numbers agree with a table that describes the current
    # toolbox". Without it an arbitrarily old measurement licenses any prose within 25 %, which is
    # O-29's own defect moved one level up: metadata present, read by nobody.
    stamp = re.search(r"^#\s*generated\s+(\S+\s+\S+)\s+toolbox\s+(\S+)", bench_text, re.M)
    if not stamp:
        note("tools/benchmark_reference.csv has no '# generated ... toolbox <version>' header -> "
             "re-run tools/benchmark.m; the header is what dates the numbers")
    elif version and stamp.group(2) != version:
        note("tools/benchmark_reference.csv was measured against toolbox %s but this tree is %s "
             "-> re-run tools/benchmark.m and commit it, or the per-epoch figures it licenses "
             "describe code that has moved" % (stamp.group(2), version))
    # A CLOSED set, every member of which ASSERTS that the number is not current. The marker must
    # sit on the line itself, so that deleting the number and leaving the excuse behind is not
    # possible. "measured on <date>" is deliberately not a member: see the note above.
    HISTORY = re.compile(r"superseded|stale|no longer|not current|wrong by|was wrong", re.I)
    CLAIM = re.compile(r"([\d.]+)\s*(?:ms|\\milli\\second\})\s*(?:per epoch|each|/epoch)", re.I)
    for dirpath, dirnames, filenames in os.walk(ROOT):
        # Vendored third-party trees are EXCLUDED. setupNonGravForces clones NRLMSISE-00 and
        # HWM14 into mex/, and a "3 ms per epoch" in someone else's comment would be a finding
        # with no honest fix available: editing upstream source, or marking it superseded.
        # This check governs what this repository asserts, not what it fetched.
        dirnames[:] = [d for d in dirnames
                       if d not in (".git", "html", "helpsearch-v4_en",
                                    "nrlmsise00", "hwm14", "msis21", "nrlmsis21", "dtm2020", "node_modules")]
        for fn in filenames:
            if not fn.endswith((".m", ".md", ".tex", ".c", ".h")):
                continue
            fp = os.path.join(dirpath, fn)
            rel = os.path.relpath(fp, ROOT).replace("\\", "/")
            try:
                lines = open(fp, encoding="utf-8").read().splitlines()
            except (UnicodeDecodeError, OSError):
                continue
            for i, line in enumerate(lines, 1):
                if HISTORY.search(line):
                    continue
                for hit in CLAIM.finditer(line):
                    try:
                        val = float(hit.group(1))
                    except ValueError:
                        continue          # "1.2.3 ms per epoch", "... ms per epoch": not a claim
                    if not any(abs(val - mv) <= 0.25 * max(mv, 1e-9) for mv in measured):
                        note("%s:%d quotes %g ms per epoch, which no row of "
                             "tools/benchmark_reference.csv measured (%s) -> re-run tools/benchmark.m "
                             "and quote it, or mark the line as superseded history"
                             % (rel, i, val, ", ".join("%g" % mv for mv in sorted(set(measured)))))

# 14. THE REGRESSION BASELINE WAS VERIFIED AT THIS VERSION (N-10a, ORC-10).
#     tools/regressionBaseline.m calibrates the reference month 2020-01 (GRACE-C and -D) and either
#     freezes it into tools/baseline_reference.json or compares a new run with it, writing the
#     toolbox version into meta.toolbox only when every number stayed within k sigma. The month
#     takes ~1.5 h, so it runs at a release on Tier A and not in CI -- and an instrument nobody is
#     made to run is the O-29 disease. This makes it unavoidable: a tree whose version is not the
#     one the baseline was last verified at is refused, exactly as a stale benchmark is (check 13).
#     WHAT IT DOES NOT CHECK: the numbers themselves (that is the tool's comparison), nor that the
#     file came from a real run rather than an edit.
import json
BASE = os.path.join(ROOT, "tools", "baseline_reference.json")
if not os.path.isfile(BASE):
    note("tools/baseline_reference.json is missing -> run tools/regressionBaseline('Mode', \"freeze\") "
         "on the DATAPOOL and commit it (N-10a)")
else:
    try:
        bmeta = json.load(open(BASE, encoding="utf-8")).get("meta", {})
    except (ValueError, OSError) as e:
        bmeta = {}
        note("tools/baseline_reference.json is not readable JSON (%s) -> regenerate it" % e)
    if version and bmeta.get("toolbox") != version:
        note("tools/baseline_reference.json was last verified at toolbox %s but this tree is %s "
             "-> run tools/regressionBaseline on the DATAPOOL (Mode \"compare\") before the release"
             % (bmeta.get("toolbox"), version))

if findings:
    print(f"auditPackage: {len(findings)} finding(s)")
    for f in findings:
        print("  -", f)
    sys.exit(1)
print(f"auditPackage: clean ({len(on_disk)} units, {len(private)} private helpers, {len(rows)} catalogue rows, version {version})")
