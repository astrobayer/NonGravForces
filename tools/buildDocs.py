#!/usr/bin/env python3
"""Render every page of doc/html from the sources that own the facts (stdlib only).

WP-21. Until 2026-09-09 only FunctionList.html was generated and the other five
pages carried a hand-edited footer that was checked for its version string and
nothing else. That is O-11, and it had already failed: Sensors.html described
Electrostatic as having a "CRN filter" -- which the delivered class explicitly
does NOT construct -- and called all three sensors future work months after two
of them shipped. A page nobody generates is a page that lies eventually.

So every derivable fact is now derived, and each from its single authority
(rule R1):

    the catalogue          doc/functionlist.csv    -> the function reference,
                                                      the sensor table, the
                                                      delivered/stub counts,
                                                      which units the API
                                                      chapter documents
    the help blocks        +NonGravForces/**/*.m   -> every signature, argument
                                                      table and error
                                                      identifier in the API
                                                      chapter
    the error taxonomy     +NonGravForces/fail.m   -> the mnemonic table
    the version and date   toolboxVersion.m        -> every footer
    the test classes       tests/tNgf*.m           -> the test catalogue: every
                                                      test method, its category,
                                                      the unit it is credited
                                                      to, its comment, and each
                                                      recorded bound with its
                                                      stated reason
    the last gate run      tests/last_run.txt      -> the measured value beside
                                                      each bound, with the run
                                                      it came from
    the pages themselves   this file               -> helptoc.xml

The prose that is NOT derivable lives here, in this generator, because a
generator that reads prose from somewhere else is two files to keep in step
instead of one.

WHY PYTHON AND NOT tools/buildDocs.m, which is what the level prompt named:
this has to run where MATLAB does not. The audit job is the fast gate, it has
no MATLAB licence, and a whole session was spent on Tier B this week with the
bridge down. Help blocks are plain text; parsing them needs no MATLAB. The
deviation is recorded in HANDOVER Part I.

Usage:  python tools/buildDocs.py            # writes every page
        python tools/buildDocs.py --check    # exit 1 if any page is stale
        python tools/buildDocs.py --dump-api # print the parsed API, for eyeballing
"""
import csv, html, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "doc", "functionlist.csv")
HTML = os.path.join(ROOT, "doc", "html")
PKG = os.path.join(ROOT, "+NonGravForces")
FAIL_FILE = os.path.join(ROOT, "+NonGravForces", "fail.m")
VERSION_FILE = os.path.join(ROOT, "+NonGravForces", "toolboxVersion.m")
OPTIONS_FILE = os.path.join(ROOT, "+NonGravForces", "Options.m")
TEX_GENERATED = os.path.join(ROOT, "doc", "src", "generated.tex")
TEX_API = os.path.join(ROOT, "doc", "src", "api_generated.tex")
TEX_TESTS = os.path.join(ROOT, "doc", "src", "tests_generated.tex")
TESTS = os.path.join(ROOT, "tests")
LAST_RUN = os.path.join(ROOT, "tests", "last_run.txt")

STYLE = ("body{font-family:Helvetica,Arial,sans-serif;max-width:900px;margin:2em auto;line-height:1.45}"
         "code,pre{background:#f4f4f4;padding:2px 4px}pre{padding:8px;overflow-x:auto}"
         "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 8px;text-align:left}"
         "h1{border-bottom:2px solid #555}")


def read_version():
    """Version, date and model as toolboxVersion.m states them (numeric fields are the authority)."""
    src = open(VERSION_FILE, encoding="utf-8").read()
    m = re.search(r"major\s*=\s*uint16\((\d+)\);\s*minor\s*=\s*uint16\((\d+)\);\s*patch\s*=\s*uint16\((\d+)\);", src)
    d = re.search(r"'date',\s*\"(\d{4}-\d{2}-\d{2})\"", src)
    g = re.search(r"'model',\s*\"([^\"]+)\"", src)
    if not (m and d and g):
        sys.exit("buildDocs: cannot read major/minor/patch, date and model from toolboxVersion.m")
    return "%s.%s.%s" % m.groups(), d.group(1), g.group(1)


def read_taxonomy():
    """The mnemonic classes and their meanings, from the help block of fail.m -- its authority."""
    src = open(FAIL_FILE, encoding="utf-8").read()
    rows = re.findall(r"^%\s{5}(\w+\*?)\s{2,}(.+?)\s*$", src, re.M)
    if len(rows) < 5:
        sys.exit("buildDocs: cannot read the mnemonic taxonomy from fail.m")
    return rows


def read_option_sets():
    """The admissible value of every enumerated option, from Options.validateValue's `sets` struct.

    That struct is the single authority the class itself validates against, so a document built
    from it cannot claim an option the toolbox would refuse. The theory document claimed three
    that it would (`sesam`, `langmuir`, `erp.anisotropy`) and two backends that had been removed,
    for three weeks, because nothing compared the two (LP-C, O-38).
    """
    src = open(OPTIONS_FILE, encoding="utf-8").read()
    block = re.search(r"sets\s*=\s*struct\((.*?)\);", src, re.S)
    if not block:
        sys.exit("buildDocs: cannot find Options.validateValue's `sets` struct")
    rows = re.findall(r"'(\w+)',\s*\[([^\]]*)\]", block.group(1))
    if len(rows) < 8:
        sys.exit("buildDocs: read only %d option sets from Options.m, expected the full table" % len(rows))
    out = []
    for key, values in rows:
        out.append((key, re.findall(r'"([^"]+)"', values)))
    return out


def tex_macro_name(key):
    """A LaTeX control sequence is letters only: drag_density -> \\ngfoptdragdensity."""
    return "ngfopt" + re.sub(r"[^A-Za-z]", "", key)


def render_generated_tex():
    """doc/src/generated.tex: every fact the theory document must not state twice.

    The .tex \\input{generated} and uses these macros, so the version on its title page and every
    option list it prints come from the same authorities the code does. Written as a separate file
    rather than into the .tex so that the prose stays hand-written and reviewable.
    """
    version, date, model = read_version()
    day, month, year = date[8:10], int(date[5:7]), date[0:4]
    months = ["January", "February", "March", "April", "May", "June",
              "July", "August", "September", "October", "November", "December"]
    lines = ["% GENERATED by tools/buildDocs.py from the authorities named beside each macro.",
             "% DO NOT EDIT: auditPackage.py compares this file by content and refuses a stale one.",
             "% Regenerate with: python tools/buildDocs.py",
             "\\newcommand{\\ngfversion}{%s}                 %% +NonGravForces/toolboxVersion.m" % version,
             "\\newcommand{\\ngfdate}{%d %s %s}%s%% toolboxVersion.m .date"
             % (int(day), months[month - 1], year, " " * 10),
             "\\newcommand{\\ngfmodel}{%s}%s%% toolboxVersion.m .model" % (model, " " * 6),
             # The .tex said "51 units" from LP-5 until 2026-09-17, four levels after it stopped being
             # true, because nothing compared it with anything. Same defect as CLAUDE.md's count (F8).
             "\\newcommand{\\ngfunitcount}{%d}%s%% rows of doc/functionlist.csv"
             % (len(list(csv.DictReader(open(CSV, encoding="utf-8")))), " " * 6)]
    lines.append("%% The admissible option values, from Options.validateValue's `sets` struct:")
    for key, values in read_option_sets():
        pretty = ", ".join("\\texttt{%s}" % v for v in values)
        lines.append("\\newcommand{\\%s}{%s}" % (tex_macro_name(key), pretty))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------------------------
# The API reference (doc/src/api_generated.tex)
#
# Same rule as everything above: one authority per fact. The catalogue says WHICH units exist and
# gives each one's H1; each unit's own help block says what it takes, returns and refuses. Nothing
# about an argument is written here -- this file only re-renders what `help NonGravForces.x.y`
# already prints, so a signature cannot drift from the code the way the theory document's option
# lists did (O-38) or Sensors.html's filter did (O-11).
#
# PARSER NOTES. The help tables are laid out for a human reader in a fixed-width font, not for a
# machine, and every rule below was earned by a false positive on a real file (the same list that
# tools/audit_args.py carries, plus what splitting the COLUMNS -- which that auditor never had to
# do -- costs on top):
#   - the Name column is not a fixed width and many rows use a SINGLE space after the name
#     ("HeatCapacity_J_m2K 1x2 double J/(m^2 K)"), so a cell split on 2-space runs alone loses the
#     size, type and units of such a row into its name;
#   - "Name" is both a column header AND a real name-value option of geom.boxWing, so the header
#     row is found by the whole line carrying Name AND Size AND Type, never by the token alone;
#   - some rows group names ("Bias, Scale, Whitening, Window"), which is terser than the template
#     but IS documentation, so the group is kept as written;
#   - a section heading is ALONE on its line and at indent one, so a row named F, T or GM does not
#     truncate the table, and an indented all-caps paragraph heading (accommodation's MODELS) is
#     not mistaken for a section;
#   - prose lines occur inside tables ("All arrays must be compatible for implicit expansion"): a
#     real row's second column is a size or a type, so anything else is kept as a note under the
#     table instead of becoming a row with a sentence for a name;
#   - names may be dotted (const's c.S0) and a continuation row may carry only a further key
#     (Model's METHOD NAME-VALUE table continues 'torque' with 'evaluate'), which is merged into
#     the row above rather than emitted as a row with four empty cells.
# Sections that are not Name/Size/Type tables at all -- METHODS, PROPERTIES, Sensor's ABSTRACT
# METHODS -- are column-aligned text, so they are reproduced verbatim rather than forced into a
# grid they were never written as. ERRORS is a list of identifiers and is rendered as one.

LAYERS = [("L0", "core", "constants, options, provenance, uniform errors and the MEX build"),
          ("L1", "io", "Level-1B products, the mass model, space weather and external grids"),
          ("L2", "env", "the environment: density, wind, shadow, Earth radiance, magnetic field"),
          ("L3", "geom", "panel geometry, attitude and occlusion"),
          ("L4", "kernel", "per-panel forces, their derivatives and their reductions"),
          ("L5", "Model", "the single point of access"),
          ("L6", "sensor", "instrument models: true dynamics mapped to what one instrument reads"),
          ("L8", "app", "applications built from the layers below")]

# The sections of a help block that describe the API. ACCURACY, REFERENCES and EXAMPLE are
# documentation too, but they are prose about the unit, not its interface, and `help` is where a
# reader meets them.
API_SECTIONS = ("INPUT ARGUMENTS", "NAME-VALUE ARGUMENTS", "NAME-VALUE",
                "CONSTRUCTOR ARGUMENTS", "CONSTRUCTOR INPUTS", "CONSTRUCTOR NAME-VALUE ARGUMENTS",
                "METHOD NAME-VALUE ARGUMENTS", "PROPERTIES", "METHODS", "ABSTRACT METHODS",
                "OUTPUT ARGUMENTS", "OUTPUT", "ERRORS", "WARNINGS")
VERBATIM_SECTIONS = ("METHODS", "ABSTRACT METHODS", "PROPERTIES")

# A heading is alone on its line, at indent one, all capitals, with an optional parenthetical
# ("PROPERTIES (read-only after construction; change through set())"). Three characters minimum:
# a prose line opening "O(1) bilinear read ..." or "T(n+1, m+1, :) with ..." is not a section.
HEADING = re.compile(r"^ ([A-Z][A-Z0-9 \-]{2,}?)\s*(\(.*\))?\s*$")
# A call, for telling the syntax lines from the first paragraph of the description. THE
# PARENTHESES ARE OPTIONAL: `c = NonGravForces.const` is how five units are invoked, and demanding
# a "(" silently dropped the whole syntax block of const, toolboxVersion, dependencyReport,
# buildMex and Options -- a block is only taken when EVERY line of it is a call, so one
# argument-less line rejected the four beside it.
CALL = re.compile(r"^\s*(?:\[[^\]]*\]\s*=\s*|[A-Za-z]\w*(?:\([^)]*\))?\s*=\s*)?"
                  r"[A-Za-z]\w*(?:\.\w+)*\s*(\(|;|%|$)")
# A real row's second cell, as tools/audit_args.py defines it: a size or a type, never a word.
SIZE_OR_TYPE = re.compile(r"^(\d+|[A-Za-z])x[\dA-Za-z]|^any$|^1x1\|1xK$|^--|^\(|^struct$|^double$"
                          r"|^string$|^logical$|^table$|^cell$|^datetime$|^char$|^function_handle$"
                          r"|^duration$|^timetable$|^NonGravForces\.|^-$")
IDENTIFIER = re.compile(r"^[A-Za-z][\w<>.]*:[\w<>*.]+")


def unit_path(row):
    """Where the catalogue says a unit lives. Check 1 of auditPackage.py proves the file is there."""
    sub = "" if row["package"] == "NonGravForces" else "+" + row["package"]
    return os.path.join(PKG, sub, row["name"] + ".m")


def qualified(row):
    """NonGravForces.kernel.accommodation, or NonGravForces.Model for the root package."""
    if row["package"] == "NonGravForces":
        return "NonGravForces." + row["name"]
    return "NonGravForces.%s.%s" % (row["package"], row["name"])


def help_block(src):
    """MATLAB's own help: the comment block after the function or classdef line, '%' removed.

    Indentation is preserved, because the column layout of the argument tables is the only thing
    that says which cell a piece of text belongs to.
    """
    lines = src.splitlines()
    i = 0
    while i < len(lines) and not re.match(r"^\s*(function|classdef)\b", lines[i]):
        i += 1
    if i == len(lines):
        return []
    while i < len(lines) and lines[i].rstrip().endswith("..."):
        i += 1
    out = []
    for line in lines[i + 1:]:
        s = line.lstrip()
        if s.startswith("%"):
            out.append(s[1:].rstrip())
        elif s:
            break
        else:
            break
    return out


def help_sections(help_lines):
    """(lead lines, [(title, parenthetical, body lines)]) in the order the help states them."""
    lead, out, cur = [], [], None
    for line in help_lines:
        m = HEADING.match(line)
        if m:
            cur = (m.group(1).strip(), (m.group(2) or "").strip(), [])
            out.append(cur)
        elif cur is None:
            lead.append(line)
        else:
            cur[2].append(line)
    return lead, out


def syntax_lines(lead):
    """The call lines under the H1: every paragraph of them, up to the first paragraph of prose."""
    blocks, cur = [], []
    for line in lead[1:]:
        if line.strip():
            cur.append(line)
        elif cur:
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)
    out = []
    for block in blocks:
        if not all_calls(block):
            break
        out.extend(block)
    return dedent(out)


def all_calls(block):
    """Is every line of this paragraph a call? A MATLAB CONTINUATION LINE IS PART OF THE ONE ABOVE.

    kernel.partials wraps its first syntax line with `...`, and judging the wrapped half on its own
    ("'T_K', T, 'Vrel_body', v, ...") rejected the paragraph and left that unit with no syntax at
    all. The continuation is still printed as the help wrote it; it is only not judged twice.
    """
    continued = False
    for line in block:
        if not continued and not CALL.match(line):
            return False
        continued = line.rstrip().endswith("...")
    return True


def dedent(lines):
    """Drop the common indentation; the table columns keep their relative positions."""
    body = [l for l in lines if l.strip()]
    if not body:
        return []
    cut = min(len(l) - len(l.lstrip()) for l in body)
    return [l[cut:] if l.strip() else "" for l in lines]


def split_cells(line):
    """[(column, text)] where the columns of a help table are separated by two or more spaces."""
    out, pos = [], 0
    for m in re.finditer(r"\s{2,}", line):
        chunk = line[pos:m.start()]
        if chunk.strip():
            out.append((pos + len(chunk) - len(chunk.lstrip()), chunk.strip()))
        pos = m.end()
    chunk = line[pos:]
    if chunk.strip():
        out.append((pos + len(chunk) - len(chunk.lstrip()), chunk.strip()))
    return out


COLUMN_WORD = re.compile(r"\b(Method|Name|Size|Type|Units|Description)\b")


def is_header(line):
    """The header row, by the WHOLE line and never by one token.

    "Name" is also a real name-value option of geom.boxWing, so a header carries Name AND
    Description AND at least one of Size, Type and Units. qc.eventTable's table has no Size
    column at all, so demanding all three would have left it to be guessed at.
    """
    words = set(COLUMN_WORD.findall(line))
    return "Name" in words and "Description" in words and bool(words & {"Size", "Type", "Units"})


def header_columns(line):
    """The header's column names and the character column each one starts in.

    NOT split on two-space runs: geom.graceFO writes its header with single spaces
    ("Name       Size Type   Units Description") and a two-space split reads three columns
    called Name, "Size Type" and "Units Description", which then swallow every row.
    """
    found = [(m.start(), m.group()) for m in COLUMN_WORD.finditer(line)]
    return [c for _, c in found], [p for p, _ in found]


SIZE_NAME = re.compile(r"^(\d+|[A-Za-z])x[\dA-Za-z]|^any$|^expanded$|^\(|^1x1\|1xK$")
TYPE_NAME = re.compile(r"^(double|string|struct|logical|table|cell|datetime|char|duration|timetable"
                       r"|object|function_handle|MException|NonGravForces\.[\w.]+|[A-Z]\w*)"
                       r"(\|[\w.]+)*$")


def synth_columns(cells):
    """Name the columns of a table whose help block omits the header row, BY WHAT THEY HOLD.

    qc.eventTable's second output column is a type ("table") and geom.graceFO's is a size ("1x1"),
    so numbering them Name, Size, Type, ... in order mislabels one of the two. The first column is
    the name and the last the description; each column between is read.
    """
    labels = ["Name"]
    for c in cells[1:-1]:
        if SIZE_NAME.match(c):
            labels.append("Size")
        elif TYPE_NAME.match(c):
            labels.append("Type")
        else:
            labels.append("Units")
    labels.append("Description")
    order = {"Name": 0, "Size": 1, "Type": 2, "Units": 3, "Description": 4}
    if [order[l] for l in labels] != sorted(order[l] for l in labels) or len(set(labels)) != len(labels):
        labels = ["Name", "Size", "Type", "Units"][:len(cells) - 1] + ["Description"]
    return labels


def table_header(lines):
    """(columns, column starts, index of the header line) for a help table, or None.

    The header is optional in practice: core.provenance and core.dependencyReport write their
    NAME-VALUE rows straight under the heading. A table without one is still a table, so the
    columns are synthesised from its widest row -- which is also where its column positions come
    from, and those are what every other row is read against.
    """
    for i, line in enumerate(lines):
        if is_header(line):
            cols, starts = header_columns(line)
            # THE ROWS, NOT THE HEADER, SAY WHERE THE COLUMNS ARE. env.wind's header is written to
            # one alignment and its rows to another, four characters wider; read against the
            # header, every Units cell of that table fell into the description. The header still
            # owns the NAMES. The dashes row is skipped: it follows the header, not the rows.
            for row in lines[i + 1:]:
                if not row.strip():
                    break
                if set(row.strip()) <= set("- "):
                    continue
                cells = split_cells(row)
                if len(cells) == len(cols):
                    starts = [p for p, _ in cells]
                break
            return cols, starts, i
    best = None
    for line in lines:
        cells = split_cells(line)
        if len(cells) >= 3 and SIZE_OR_TYPE.match(cells[1][1]) and (best is None or len(cells) > len(best)):
            best = cells
    if best is not None:
        n = min(len(best), 5)
        return synth_columns([t for _, t in best[:n]]), [p for p, _ in best[:n]], -1
    # geom.graceFO writes its options with single spaces throughout, so there are only two cells:
    # the name, size and type in one, the description in the other.
    for line in lines:
        cells = split_cells(line)
        if len(cells) != 2:
            continue
        toks = [(m.start(), m.group()) for m in re.finditer(r"\S+", cells[0][1])]
        if len(toks) < 2 or not SIZE_OR_TYPE.match(toks[1][1]):
            continue
        n = min(len(toks), 4)
        starts = [cells[0][0] + p for p, _ in toks[:n]] + [cells[1][0]]
        return synth_columns([t for _, t in toks[:n]] + [cells[1][1]]), starts, -1
    return None


def parse_table(lines):
    """{columns, rows, lead, notes} for a help table, or None when the section is not one.

    Each TOKEN is assigned to the column it starts in, not each two-space-separated cell: a units
    cell as wide as its column leaves a single space before the description (const's
    "W/m^2/K^4 Stefan-Boltzmann, ...") and cell-wise reading buries the description in the units.
    The one layout column positions cannot read is a row written with single spaces throughout
    (boxWing's "HeatCapacity_J_m2K 1x2 double J/(m^2 K)"), which lands a column to the right of
    where it belongs; that row is rebuilt IN ORDER, and only when the token behind the name is a
    size, so a sentence is never taken apart this way.
    """
    hdr = table_header(lines)
    if hdr is None:
        return None
    columns, starts, hi = hdr
    lower = [c.lower() for c in columns]
    name_col = lower.index("name") if "name" in lower else 0
    size_col = lower.index("size") if "size" in lower else None
    type_col = lower.index("type") if "type" in lower else None

    def place(pos):
        col = 0
        for j, s in enumerate(starts):
            if s <= pos + 1:
                col = j
        return col

    rows, notes = [], []
    for line in lines[hi + 1:]:
        if not line.strip():
            break                                           # a blank line ends the table
        if set(line.strip()) <= set("- "):
            continue                                        # the dashes under the header
        toks = [(m.start(), m.group()) for m in re.finditer(r"\S+", line)]
        cells = split_cells(line)
        row = [""] * len(columns)
        if len(starts) > 1 and toks[0][0] >= starts[1] - 1:
            # Nothing in the first column: a wrapped cell. It wraps INTO ITS OWN COLUMN when it
            # stays in one (kernel.partials continues "N per unit of" with "the parameter"), and
            # into the description when it crosses a boundary, which is how a description that is
            # a space or two out of alignment cannot land in the units.
            for pos, text in toks:
                row[place(pos)] = (row[place(pos)] + " " + text).strip()
            if len([c for c in row if c]) > 1:
                row = [""] * len(columns)
                row[-1] = " ".join(t for _, t in toks)
            if rows:
                for i, c in enumerate(row):
                    if c:
                        rows[-1][i] = (rows[-1][i] + " " + c).strip()
            else:
                notes.append(" ".join(t for _, t in toks))
            continue
        for pos, text in toks:
            row[place(pos)] = (row[place(pos)] + " " + text).strip()
        first = cells[0][1].split()
        if (size_col is not None and not SIZE_OR_TYPE.match(row[size_col] or "?")
                and abs(cells[0][0] - starts[name_col]) <= 1
                and len(first) >= 3 and SIZE_OR_TYPE.match(first[1])):
            row = [""] * len(columns)                       # rebuild the single-spaced row in order
            row[name_col] = first[0]
            rest = first[1:]
            for col in range(name_col + 1, len(columns) - 1):
                if rest:
                    row[col] = rest.pop(0)
            if rest:
                row[len(columns) - 2] = (row[len(columns) - 2] + " " + " ".join(rest)).strip()
            row[-1] = " ".join(t for _, t in cells[1:])
        # A REAL ROW'S SECOND COLUMN IS A SIZE OR A TYPE. Prose inside a table ("All arrays must be
        # compatible for implicit expansion") otherwise becomes a row whose name is a sentence.
        real = row[name_col] and ((size_col is not None and SIZE_OR_TYPE.match(row[size_col] or "?"))
                                  or (type_col is not None and SIZE_OR_TYPE.match(row[type_col] or "?")))
        head = [c for c in row[:name_col + 1] if c]
        middle = [c for c in row[name_col + 1:-1] if c]
        if real:
            rows.append(row)
        elif rows and not middle and all(" " not in c for c in head):
            for i, c in enumerate(row[:-1]):                # 'evaluate' under 'torque': same row
                if c:
                    rows[-1][i] = (rows[-1][i] + ", " + c).strip(", ")
            rows[-1][-1] = (rows[-1][-1] + " " + row[-1]).strip()
        else:
            notes.append(" ".join(c for c in row if c))     # prose inside the table
    return {"columns": columns, "rows": rows, "lead": lines[:max(hi, 0)], "notes": notes}


def parse_errors(lines):
    """([(identifier, meaning)], [prose lines]) from an ERRORS block."""
    items, notes = [], []
    for line in lines:
        if not line.strip():
            continue
        parts = split_cells(line)
        first = parts[0][1]
        if IDENTIFIER.match(first):
            rest = [p for _, p in parts[1:]]
            # an identifier as wide as its column leaves ONE space before its meaning
            # (kernel.accommodation's inputSurfaceMass), which a two-space split cannot see. NOT
            # when the cell is a LIST of identifiers ("referenceSystem:*, gracefoReading:*"),
            # which the same cut would have broken in half.
            if " " in first:
                head, _, tail = first.partition(" ")
                if not head.endswith(",") and not IDENTIFIER.match(tail.strip()):
                    first = head
                    rest.insert(0, tail.strip())
            items.append([first, " ".join(rest)])
        elif items and parts[0][0] > 4:
            items[-1][1] = (items[-1][1] + " " + " ".join(p for _, p in parts)).strip()
        else:
            notes.append(line.strip())
    return items, notes


def chunks(lines):
    """The paragraphs of a section body, blank-line separated."""
    out, cur = [], []
    for line in lines:
        if line.strip():
            cur.append(line)
        elif cur:
            out.append(cur)
            cur = []
    if cur:
        out.append(cur)
    return out


def is_structured(chunk):
    """Column-aligned text (a method signature list, a property list), not a prose paragraph.

    Three or more consecutive spaces inside a line is the mark: prose does not carry them, and a
    line written to line up under a heading always does. ONE such line is enough, because the
    continuation lines of an aligned block carry no interior run at all -- counting a majority
    read Options' METHODS block, four fifths of it continuation, as a paragraph.
    """
    return any(re.search(r"\S {3,}\S", l) for l in chunk)


# --- LaTeX -----------------------------------------------------------------------------------

TEX_SPECIAL = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&",
               "#": r"\#", "_": r"\_", "%": r"\%", "^": r"\textasciicircum{}",
               "~": r"\textasciitilde{}", "|": r"\textbar{}", "<": r"\textless{}",
               ">": r"\textgreater{}"}
# The five non-ASCII characters the help blocks actually contain (W\"oske, \S, em dash), as macros
# that need no input encoding. Counted from the sources, not guessed.
TEX_UNICODE = {"\u00f6": r'\"{o}', "\u00e4": r'\"{a}', "\u00fc": r'\"{u}',
               "\u00a7": r"\S{}", "\u2014": "---"}
VERB_UNICODE = {"\u00f6": "oe", "\u00e4": "ae", "\u00fc": "ue", "\u00a7": "sect.", "\u2014": "--"}


def tex(s):
    """Escape for LaTeX text: the specials first, so the macros the next step inserts survive."""
    out = "".join(TEX_SPECIAL.get(c, c) for c in s)
    for ch, macro in TEX_UNICODE.items():
        out = out.replace(ch, macro)
    return out


def tex_tight(s):
    """Escape, and allow a line break after a dot, colon, slash, bar or underscore.

    A p{} column can only break at a space, and the narrow columns hold names that have none:
    "NonGravForces.geom.MacroModel" as a type is 29 characters in a 2.2 cm column. \\allowbreak
    gives it somewhere to break instead of running into the column beside it.
    """
    out = tex(s)
    for mark in (".", ":", "/", "\\_", "\\textbar{}"):
        out = out.replace(mark, mark + "\\allowbreak{}")
    return out


def verb(s):
    """For a listings environment, which is verbatim: transliterate instead of escaping."""
    for ch, plain in VERB_UNICODE.items():
        s = s.replace(ch, plain)
    return "".join(c if 32 <= ord(c) < 127 else "?" for c in s)


def tex_listing(lines):
    return ["\\begin{lstlisting}"] + [verb(l) for l in lines] + ["\\end{lstlisting}"]


def column_widths(columns):
    """Fixed widths in cm: the roles that are always short get a fixed share, Description the rest.

    A longtable of l columns would be set by its widest cell and overflow the page on
    NonGravForces.geom.MacroModel as a type; p{} columns cannot.
    """
    fixed = {"name": 3.0, "method": 1.8, "size": 1.5, "type": 2.2, "units": 1.5}
    widths = [fixed.get(c.lower(), 2.0) for c in columns[:-1]]
    rest = 14.6 - sum(widths)
    if rest < 4.0:                                          # squeeze the fixed columns, not the text
        scale = (14.6 - 4.0) / sum(widths)
        widths = [w * scale for w in widths]
        rest = 4.0
    return widths + [rest]


def tex_table(table):
    """A longtable that breaks across pages, with the header repeated and Description wrapped."""
    cols, rows = table["columns"], table["rows"]
    if not rows:
        return []
    widths = column_widths(cols)
    # \footnotesize in the COLUMN SPEC, not as a group around the environment: a longtable is
    # typeset across pages and wrapping one in a group is the documented way to confuse it.
    spec = "".join(">{\\footnotesize\\raggedright\\arraybackslash%s}p{%.2fcm}"
                   % ("\\ttfamily" if i == 0 else "", w) for i, w in enumerate(widths))
    head = " & ".join("{\\footnotesize\\textbf{%s}}" % tex(c) for c in cols) + " \\\\"
    out = ["\\begin{longtable}{@{}" + spec + "@{}}",
           "\\toprule", head, "\\midrule", "\\endfirsthead",
           "\\toprule", head, "\\midrule", "\\endhead",
           "\\midrule", "\\endfoot", "\\bottomrule", "\\endlastfoot"]
    for row in rows:
        cells = [tex_tight(c) for c in row[:-1]] + [tex(row[-1])]
        out.append(" & ".join(cells) + " \\\\")
    out += ["\\end{longtable}"]
    return out


def tex_section_body(title, lines):
    """One help section as LaTeX: a table where the help drew one, a list of identifiers for
    ERRORS, verbatim for the column-aligned blocks, and prose for the rest."""
    out = []
    if title in ("ERRORS", "WARNINGS"):
        items, notes = parse_errors(lines)
        if items:
            out.append("\\begin{description}[leftmargin=1.2em,style=nextline,itemsep=0pt]")
            for ident, meaning in items:
                out.append("\\item[\\texttt{%s}] %s" % (tex_tight(ident), tex(meaning) or "---"))
            out.append("\\end{description}")
        for n in notes:
            out.append(tex(n))
            out.append("")
        return out
    for i, chunk in enumerate(chunks(lines)):
        verbatim = title in VERBATIM_SECTIONS
        table = None if verbatim else parse_table(chunk)
        # The FIRST paragraph of a METHODS or PROPERTIES block is the list itself, whatever it
        # looks like: geom.MacroModel's is a single signature line, which the alignment test alone
        # would have set as a paragraph of prose.
        if verbatim and i == 0:
            out += tex_listing(dedent(chunk))
        elif table is not None:
            for l in table["lead"]:
                out.append(tex(l.strip()))
            if table["lead"]:
                out.append("")
            out += tex_table(table)
            for n in table["notes"]:
                out += ["", "\\noindent " + tex(n)]
        elif is_structured(chunk):
            out += tex_listing(dedent(chunk))
        else:
            out.append(tex(" ".join(l.strip() for l in chunk)))
            out.append("")
    return out


def render_api_tex():
    """doc/src/api_generated.tex: every catalogued unit, from the catalogue and its own help.

    \\input-able from a report-class document: sections are layers, subsections are units, and
    nothing here opens a document or a chapter.
    """
    version = read_version()[0]
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    known = {lvl for lvl, _, _ in LAYERS}
    unknown = sorted({r["level"] for r in rows} - known)
    if unknown:
        sys.exit("buildDocs: doc/functionlist.csv holds level(s) %s that LAYERS does not name"
                 % ", ".join(unknown))
    out = ["% GENERATED by tools/buildDocs.py from doc/functionlist.csv (which units exist, and",
           "% each one's H1) and the help block of each unit (everything else).",
           "% DO NOT EDIT: auditPackage.py compares this file by content and refuses a stale one.",
           "% Regenerate with: python tools/buildDocs.py",
           "%",
           "% \\input-able from a report-class document; it opens no document and no chapter.",
           "",
           "\\noindent This reference is generated from the catalogue \\texttt{doc/functionlist.csv} "
           "and from the help block of each unit, by \\texttt{tools/buildDocs.py}. It is therefore "
           "what \\texttt{help NonGravForces.\\textless{}unit\\textgreater{}} prints, and cannot "
           "drift from it: \\texttt{auditPackage.py} re-renders this file and refuses a stale one. "
           "%d units, toolbox version %s." % (len(rows), version),
           ""]
    for level, pkg, blurb in LAYERS:
        units = [r for r in rows if r["level"] == level]
        if not units:
            continue
        out += ["\\section{%s \\texttt{%s} --- %s}" % (level, tex(pkg), tex(blurb)), ""]
        for row in units:
            out += render_unit(row)
    return "\n".join(out).rstrip() + "\n"


def render_unit(row):
    """One catalogued unit as a subsection: H1, syntax, its tables, its errors."""
    name = qualified(row)
    src = open(unit_path(row), encoding="utf-8").read()
    help_lines = help_block(src)
    lead, sections = help_sections(help_lines)
    m = re.match(r"^([A-Z][A-Z0-9_]*)\s+(.*)$", lead[0]) if lead else None
    h1 = m.group(2) if m else row["h1"]
    out = ["\\subsection{\\texttt{%s}}" % tex(name),
           "\\label{api:%s.%s}" % (row["package"], row["name"]),
           "\\noindent %s" % tex(h1), ""]
    syn = syntax_lines(lead)
    if syn:
        out += tex_listing(syn)
        out.append("")
    for title, paren, body in sections:
        if title not in API_SECTIONS:
            continue
        rendered = tex_section_body(title, body)
        if not [l for l in rendered if l.strip()]:
            continue
        head = title + (" " + paren if paren else "")
        out += ["\\subsubsection*{%s}" % tex(head)]
        out += rendered
        out.append("")
    return out


# ---------------------------------------------------------------------------------------------
# The test catalogue (doc/src/tests_generated.tex)
#
# Added 2026-10-31 on the owner's request for a PDF that describes every test. Everything is read
# from the test classes themselves (method, TestTags block, the comment under the function line,
# each verifyBoundAndRecord call's label, bound, units and reason) and from tests/last_run.txt (the
# measured value). The unit a method is credited to follows auditTestCategories.m's rule exactly:
# the LONGEST catalogued name the method name begins with, case-insensitively.

def _matlab_args(call):
    """Top-level comma-separated arguments of a call's argument text (quotes and brackets respected)."""
    args, depth, cur, quote = [], 0, [], None
    for i, c in enumerate(call):
        if quote:
            cur.append(c)
            if c == quote:
                quote = None
            continue
        if c in ('"', "'") and (c == '"' or not cur or not (cur[-1].isalnum() or cur[-1] in ")]}.'_")):
            quote = c
            cur.append(c)
            continue
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        if c == "," and depth == 0:
            args.append("".join(cur).strip())
            cur = []
            continue
        cur.append(c)
    if "".join(cur).strip():
        args.append("".join(cur).strip())
    return args


def _literal(arg):
    """The text of a string literal argument, or of a concatenation of literals; None otherwise."""
    parts = re.findall(r'"((?:[^"]|"")*)"|\'((?:[^\']|\'\')*)\'', arg)
    rest = re.sub(r'"(?:[^"]|"")*"|\'(?:[^\']|\'\')*\'', "", arg)
    if not parts or re.sub(r"[\s+\[\],]", "", rest):
        return None
    return "".join((a or b).replace('""', '"').replace("''", "'") for a, b in parts)


def _bound_calls(body):
    """Every verifyBoundAndRecord call in a method body: (label, bound, units, note)."""
    text = re.sub(r"\.\.\.[^\n]*\n\s*", " ", body)
    out = []
    for m in re.finditer(r"verifyBoundAndRecord\(", text):
        i, depth, quote = m.end(), 1, None
        while i < len(text) and depth:
            c = text[i]
            if quote:
                if c == quote:
                    quote = None
            elif c == '"':
                quote = c
            elif c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
            i += 1
        a = _matlab_args(text[m.end():i - 1]) + ["", "", "", "", ""]
        out.append((_literal(a[0]), a[2].strip(), _literal(a[3]) or "", _literal(a[4]) or ""))
    return out


def _last_run():
    """(header line, {label: measured}) from tests/last_run.txt; the first record of a label wins."""
    if not os.path.isfile(LAST_RUN):
        return "", {}
    lines = open(LAST_RUN, encoding="utf-8").read().split("\n")
    header = " ".join(l.strip() for l in lines[:2] if l.strip())
    measured, inside = {}, False
    for l in lines:
        if l.startswith("== "):
            inside = l.startswith("== recorded bounds")
            continue
        m = re.match(r"^(.*?)\s+([-+0-9.eEInfNa]+)\s+<=\s+[-+0-9.eEInfNa]+\s", l) if inside else None
        if m and m.group(1).strip() not in measured:
            measured[m.group(1).strip()] = m.group(2)
    return header, measured


def _test_classes():
    """[(class, H1, [(method, tags, comment, bounds)])] for every tests/tNgf*.m, in file order."""
    classes = []
    for path in sorted(f for f in os.listdir(TESTS) if re.match(r"tNgf\w*\.m$", f)):
        lines = open(os.path.join(TESTS, path), encoding="utf-8").read().split("\n")
        h1 = ""
        for l in lines[1:6]:
            if l.startswith("%"):
                h1 = re.sub(r"^%\s*[A-Za-z0-9_]+\s+", "", l).strip()
                break
        methods, tags, intest, i = [], [], False, 0
        while i < len(lines):
            l = lines[i]
            mb = re.match(r"^    methods\b\s*(\((.*)\))?", l)
            if mb:
                attrs = mb.group(2) or ""
                intest = re.search(r"\bTest\b", attrs) is not None
                tags = re.findall(r"'(\w+)'", attrs.split("TestTags", 1)[1]) if "TestTags" in attrs else []
            mf = re.match(r"^        function\s+(?:\[?[\w, ]*\]?\s*=\s*)?(\w+)\s*\(", l)
            if mf and intest:
                j, comment = i + 1, []
                while j < len(lines) and lines[j].strip().startswith("%"):
                    comment.append(lines[j].strip().lstrip("%").strip())
                    j += 1
                k = j
                while k < len(lines) and lines[k] != "        end":
                    k += 1
                body = "\n".join(lines[i:k])
                text = " ".join(comment)
                text = re.sub(r"^[A-Z][A-Z0-9_]+\s+", "", text)      # a %NAME H1 prefix
                methods.append((mf.group(1), tags, text, _bound_calls(body)))
                i = k
            i += 1
        classes.append((path[:-2], h1, methods))
    return classes


def render_tests_tex():
    """doc/src/tests_generated.tex: every test method, its category, unit, purpose and bounds."""
    version = read_version()[0]
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    units = [r for r in rows if r["status"] not in ("stub", "interface")]
    header, measured = _last_run()
    classes = _test_classes()
    nMethods = sum(len(m) for _, _, m in classes)
    nBounds = sum(len(b) for _, _, ms in classes for _, _, _, b in ms)
    out = ["% GENERATED by tools/buildDocs.py from tests/tNgf*.m (methods, tags, comments, recorded",
           "% bounds) and tests/last_run.txt (measured values).",
           "% DO NOT EDIT: auditPackage.py compares this file by content and refuses a stale one.",
           "% Regenerate with: python tools/buildDocs.py",
           "",
           "\\noindent This catalogue is generated from the test classes by \\texttt{tools/buildDocs.py}: "
           "%d test methods in %d classes, %d recorded bounds, toolbox version %s. For each method it "
           "gives the category its \\texttt{TestTags} block assigns, the unit \\texttt{auditTestCategories} "
           "credits it to (the longest catalogued name the method name begins with), the comment under its "
           "function line, and every bound it records through \\texttt{verifyBoundAndRecord}: label, bound, "
           "units and the reason stated for the bound's size. The measured value beside a bound is the one "
           "\\texttt{tests/last\\_run.txt} holds, from the run: \\emph{%s}; a dash means that run did not "
           "record the label (a newer test, or a label built at run time). A method whose label is computed "
           "at run time is shown as \\emph{(computed label)}. The seven categories and the rules behind them "
           "are those of the section on testing." % (nMethods, len(classes), nBounds, tex(version), tex(header) or "none"),
           ""]
    for cls, h1, methods in classes:
        out += ["\\subsection{\\texttt{%s}}" % tex(cls)]
        if h1:
            out += ["\\noindent %s" % tex(h1), ""]
        if not methods:
            out += ["\\noindent (no test methods)", ""]
            continue
        out.append("\\begin{description}[leftmargin=1em,style=nextline,itemsep=2pt,font=\\normalfont\\ttfamily\\small]")
        for name, tags, comment, bounds in methods:
            owner, best = "", 0
            for r in units:
                if name.lower().startswith(r["name"].lower()) and len(r["name"]) > best:
                    owner, best = r["package"] + "." + r["name"] if r["package"] != "NonGravForces" else r["name"], len(r["name"])
            meta = ", ".join(tags) if tags else "untagged"
            meta += " \\textperiodcentered\\ " + ("\\texttt{%s}" % tex_tight(owner) if owner else "no unit")
            item = "\\item[%s] {\\small\\emph{%s}}" % (tex_tight(name), meta)
            if comment:
                item += " --- " + tex(comment if len(comment) <= 700 else comment[:697].rsplit(" ", 1)[0] + " ...")
            out.append(item)
            if bounds:
                out.append("\\begin{itemize}[leftmargin=1.2em,itemsep=0pt,topsep=1pt]\\small")
                for label, bound, units_, note in bounds:
                    shown = tex(label) if label is not None else "\\emph{(computed label)}"
                    val = measured.get(label, "---") if label is not None else "---"
                    line = "\\item %s: measured %s, bound %s" % (shown, tex(val), tex(bound))
                    if units_ and units_ != "-":
                        line += " %s" % tex(units_)
                    if note:
                        line += " --- \\emph{%s}" % tex(note)
                    out.append(line)
                out.append("\\end{itemize}")
        out += ["\\end{description}", ""]
    return "\n".join(out).rstrip() + "\n"


def page(title, body, version, date, model):
    """The shared chrome. One place for the head, the style and the footer."""
    return ("<!DOCTYPE html><html><head><meta charset='utf-8'><title>" + html.escape(title, quote=False) +
            "</title><style>" + STYLE + "</style></head><body><h1>" + html.escape(title, quote=False) + "</h1>\n" +
            body +
            "<hr><p><small>NonGravForces %s &middot; %s &middot; generated by %s from "
            "<code>tools/buildDocs.py</code></small></p></body></html>\n"
            % (version, date, html.escape(model, quote=False)))


def render_all():
    """Every page, keyed by file name. Facts are derived; only prose is written out here."""
    version, date, model = read_version()
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    done = [r for r in rows if r["status"] in ("done", "interface")]
    stubs = [r for r in rows if r["status"] == "stub"]
    out = {}

    # --- function reference: the catalogue, rendered
    p = ["<p>Generated from <code>doc/functionlist.csv</code>, the single authority for the catalogue, by "
         "<code>tools/buildDocs.py</code>. %d of %d units are delivered; %d are documented stubs that raise "
         "<code>notImplemented</code>.</p><table><tr><th>package</th><th>name</th><th>kind</th>"
         "<th>level</th><th>WP</th><th>status</th><th>purpose</th></tr>" % (len(done), len(rows), len(stubs))]
    for r in rows:
        p.append("<tr><td>%s</td><td><code>%s</code></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % tuple(
            html.escape(r[k], quote=False) for k in ("package", "name", "kind", "level", "wp", "status", "h1")))
    p.append("</table>")
    out["FunctionList.html"] = page("Function Reference", "".join(p), version, date, model)

    # --- sensors: the table comes from the catalogue, so it cannot describe a class that changed
    p = ["<p>The physics core returns the true non-gravitational acceleration, torque, angular rate and angular "
         "acceleration. A <code>NonGravForces.sensor.Sensor</code> subclass maps these to what one instrument "
         "reads and supplies partials with respect to its calibration parameters. Calibration and quality "
         "control are written against the abstract interface, so a new instrument is one new class.</p>"
         "<p>This table is generated from the catalogue. It used to be written by hand, and by 2026-09-09 it "
         "described a filter one of these classes did not apply. Since N-04 (2026-10-14) "
         "<code>Electrostatic</code> applies the Level-1B CRN filter by default, with the kernel taken "
         "from <code>timeSeriesPreparation.crnKernel</code> and never built here.</p>"
         "<table><tr><th>class</th><th>status</th><th>what it models</th></tr>"]
    for r in rows:
        if r["package"] == "sensor":
            p.append("<tr><td><code>%s</code></td><td>%s</td><td>%s</td></tr>" % tuple(
                html.escape(r[k], quote=False) for k in ("name", "status", "h1")))
    p.append("</table><p>An optomechanical or MEMS instrument is a further subclass with its own noise model; "
             "the core does not change for it.</p>")
    out["Sensors.html"] = page("Sensor Layer", "".join(p), version, date, model)

    # --- error handling: the taxonomy comes from fail.m
    p = ["<p>Every error passes through <code>NonGravForces.fail(unit, mnemonic, message, ...)</code>. The "
         "identifier is always <code>NonGravForces:&lt;unit&gt;:&lt;mnemonic&gt;</code>; the message always ends "
         "with a correction hint after <code>-&gt;</code>. The classes below are read from the help block of "
         "<code>fail.m</code>, which owns them.</p><table><tr><th>class</th><th>meaning</th></tr>"]
    for name, meaning in read_taxonomy():
        p.append("<tr><td><code>%s</code></td><td>%s</td></tr>" % (html.escape(name, quote=False),
                                                                   html.escape(meaning, quote=False)))
    p.append("</table><p>Principles: refuse, report or proceed &mdash; never degrade silently; no silent fallback "
             "between fidelity tiers; unknown option fields error; dependency errors propagate unwrapped; "
             "warnings are collected into <code>info.provenance.warnings</code>; library code never prints.</p>")
    out["ErrorHandling.html"] = page("Error Handling", "".join(p), version, date, model)

    # --- architecture: the layer list comes from the packages that exist in the catalogue
    pkgs = []
    for r in rows:
        if r["package"] not in pkgs:
            pkgs.append(r["package"])
    p = ["<p>Eight layers, each depending only on those below it (L7, quality control, moved to "
         "<code>timeSeriesPreparation.qc</code> and was deleted in 0.14.0). The kernel returns per-panel forces and never "
         "sums internally; the environment computes no geometry; <code>NonGravForces.Model</code> is the single "
         "point of access and composes the rest.</p><table><tr><th>package</th><th>units</th></tr>"]
    for pk in pkgs:
        names = ", ".join("<code>%s</code>" % html.escape(r["name"], quote=False) for r in rows if r["package"] == pk)
        p.append("<tr><td>%s</td><td>%s</td></tr>" % (html.escape(pk, quote=False), names))
    p.append("</table><p>Two evaluation paths are kept side by side by design: the vectorised MATLAB path "
             "(<code>Model.accel</code>, <code>Model.series</code>) for an arc, and the pointwise compiled path "
             "(<code>kernel.singleEpoch</code> over <code>mex/ngfKernel.c</code>) for an integrator's right-hand "
             "side. Neither is a fidelity tier of the other and neither is justified by speed; they must agree, "
             "and a test measures that they do.</p>"
             "<p>Constraint C2: base MATLAB and the Parallel Computing Toolbox only. "
             "<code>auditPackage.py</code> refuses calls into any other product, in code and in help text.</p>")
    out["Architecture.html"] = page("Architecture", "".join(p), version, date, model)

    # --- getting started (rewritten 2026-10-31: the old page named one setup call of six and "two MEX gateways")
    p = ["<p>This page is the short path from a fresh clone to a first evaluation. The user guide "
         "(<code>doc/src/NonGravForces_UserGuide.tex</code>) is the long one, and the theory document "
         "(<code>doc/src/NonGravForces_Theory_and_Implementation.tex</code>) describes every model and every "
         "test. Neither PDF is committed: the CI job <code>theory</code> builds both on every push and attaches "
         "them to the run as the artefacts <code>userGuide</code> and <code>theoryDocument</code>; locally, "
         "<code>pdflatex</code> in <code>doc/src</code> builds them.</p>"
         "<h2>1. What you need</h2>"
         "<p>MATLAB R2026a, base MATLAB only (the Parallel Computing Toolbox is used when present); the MinGW-w64 "
         "C compiler support package for the MEX targets; gfortran only for the opt-in DTM2020 backend. Five "
         "sibling toolboxes, set up in this order:</p>"
         "<pre>run('...\\referenceSystem\\setupReferenceSystem.m')                 %% time, EOP, frames, Sun/Moon\n"
         "run('...\\gracefoReading\\setupGraceFOReading.m')                   %% the L1B readers\n"
         "run('...\\odeSuite\\gjToolboxSetup.m')                              %% app.DynamicsSource (integrator)\n"
         "run('...\\timeSeriesPreparation\\setupTimeSeriesPreparation.m')     %% CRN filter, QC, Arc\n"
         "run('...\\lsqvce\\lsqvce_setup.m')                                  %% the calibration estimator</pre>"
         "<h2>2. Set up and build</h2>"
         "<pre>cd('...\\NonGravForces'); setupNonGravForces    %% path, dependency report, pinned fetches\n"
         "NonGravForces.buildMex                          %% msis, msis21, hwm14, dtm2020, rhs, kernel</pre>"
         "<p>Setup is the only part of the toolbox that touches the network. Third-party model sources and data "
         "(NRLMSISE-00, NRLMSIS 2.1, HWM14, DTM2020, IGRF, CERES, the DE440 kernel) are not redistributed here: "
         "setup fetches each at a pinned commit or SHA-256 the first time it runs and says so. "
         "<code>'FetchMsis'</code>, <code>'FetchData'</code> and <code>'FetchKernels'</code> take "
         "<code>\"never\"</code> to leave the network alone. With a DATAPOOL (<code>'DataPool'</code> or the "
         "environment variable) setup also fetches the observation data of the reference days. "
         "<code>buildMex</code> refuses rather than fetching: an installer may reach the network, a build step "
         "may not.</p>"
         "<h2>3. Verify</h2>"
         "<pre>runNgfTests                      %% the suite; 'ExcludeSlow', false adds the speed tests\n"
         "matlab -batch ciGate             %% what CI runs: the suite, the category audit, the verdict</pre>"
         "<h2>4. A first evaluation</h2>"
         "<p>On the committed excerpts of GRACE-C, 2020-01-01. The <code>standard</code> preset multiplies "
         "NRLMSISE-00 by an observed daily density scale and refuses to evaluate without one:</p>"
         "<pre>f   = @(varargin) fullfile(pwd, 'data', varargin{:});\n"
         "eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));\n"
         "sw  = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));\n"
         "cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));\n"
         "swarm = string(f('samples', 'swarm', \"SW_OPER_DNS\" + [\"A\" \"B\" \"C\"] + \"POD_2_2020-01-01_2020-01-03.csv\"));\n"
         "S   = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));\n"
         "[ep, r, v] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);\n"
         "[eq, q]    = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));\n"
         "att = NonGravForces.geom.AttitudeProvider(eq, q);\n"
         "M   = NonGravForces.Model(NonGravForces.geom.graceFO(\"C\"), NonGravForces.Options(\"standard\"), ...\n"
         "          'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);\n"
         "[a, info] = M.accel(ep, r, v, att.at(ep));     %% 3xK m/s^2 in the body frame; info.drag, .radiation, .thermal</pre>"
         "<h2>5. The demos</h2>"
         "<p>Each is a working sequence of public calls, and each is EXECUTED BY THE GATE "
         "(<code>tests/tNgfDemos.m</code>), so an API that moves breaks the demo in CI rather than in front of "
         "a reader:</p>"
         "<pre>demoForceBudget         %% per-effect budget, two presets side by side\n"
         "demoDensityInversion    %% recover a planted density scale\n"
         "demoGraceFOCalibration  %% recover a planted bias and scale\n"
         "demoCAI                 %% cold-atom accelerometer, and why it needs the attitude rate</pre>"
         "<p>Then <code>doc NonGravForces</code> for the function reference, and <code>help "
         "NonGravForces.&lt;unit&gt;</code> for any unit: every help block follows one template (syntax, "
         "arguments with units and defaults, outputs, accuracy, errors, references, an example).</p>"
         "<p>The toolbox is MIT licensed. See <code>THIRD-PARTY-NOTICES.md</code> for what that does not "
         "cover.</p>"]
    out["GettingStarted.html"] = page("Getting Started", "".join(p), version, date, model)

    # --- testing
    p = ["<p>Every public unit carries seven test categories: contract, reference, precision, vectorization, "
         "robustness, metamorphic and speed. <code>tests/auditTestCategories.m</code> refuses a unit that is "
         "short of one, and the count is reported by the gate alongside the pass count.</p>"
         "<p>Bounds are measured, never declared. A passing bound leaves its number in "
         "<code>tests/last_run.txt</code> so a later change can be compared against it, and a failing bound is "
         "decomposed rather than widened. Where cancellation governs, the bound is computed from "
         "<code>eps</code> of the quantities involved instead of written as a literal.</p>"
         "<pre>matlab -batch ciGate            %% the full gate\n"
         "bash mex/test/runCChecks.sh     %% the C side, no MATLAB needed\n"
         "python auditPackage.py          %% catalogue, help template, taxonomy, C2 scan</pre>"]
    out["Testing.html"] = page("Testing", "".join(p), version, date, model)

    # --- index
    p = ["<p>A MATLAB toolbox computing non-gravitational forces and torques on GRACE-FO and GRACE satellites: gas-surface "
         "interaction, solar and Earth radiation pressure, thermal recoil, and the sensors that measure them. "
         "%d of %d catalogued units are delivered; %d are documented stubs.</p><ul>"
         % (len(done), len(rows), len(stubs))]
    for f, t in (("GettingStarted.html", "Getting Started"), ("Architecture.html", "Architecture"),
                 ("FunctionList.html", "Function Reference"), ("Sensors.html", "Sensor Layer"),
                 ("ErrorHandling.html", "Error Handling"), ("Testing.html", "Testing")):
        p.append("<li><a href='%s'>%s</a></li>" % (f, t))
    p.append("</ul>")
    out["index.html"] = page("NonGravForces", "".join(p), version, date, model)

    # --- helptoc: generated too, so a page cannot exist without being reachable
    toc = ["<?xml version='1.0' encoding='utf-8'?>\n<toc version='2.0'>\n",
           "<tocitem target='index.html'>NonGravForces\n"]
    for f, t in (("GettingStarted.html", "Getting Started"), ("Architecture.html", "Architecture"),
                 ("FunctionList.html", "Function Reference"), ("Sensors.html", "Sensor Layer"),
                 ("ErrorHandling.html", "Error Handling"), ("Testing.html", "Testing")):
        toc.append("<tocitem target='%s'>%s</tocitem>\n" % (f, t))
    toc.append("</tocitem>\n</toc>\n")
    out["helptoc.xml"] = "".join(toc)
    return out


def dump_api():
    """What the parser made of every help block, for reading before trusting. Prints only."""
    for row in csv.DictReader(open(CSV, encoding="utf-8")):
        src = open(unit_path(row), encoding="utf-8").read()
        lead, sections = help_sections(help_block(src))
        print("=" * 100)
        print(qualified(row), " [%s %s]" % (row["level"], row["kind"]))
        print("  syntax: %d line(s)" % len(syntax_lines(lead)))
        for l in syntax_lines(lead):
            print("     | " + l)
        for title, paren, body in sections:
            if title not in API_SECTIONS:
                continue
            print("  -- %s %s" % (title, paren))
            if title in ("ERRORS", "WARNINGS"):
                items, notes = parse_errors(body)
                for ident, meaning in items:
                    print("     ! %-55s %s" % (ident, meaning[:80]))
                for n in notes:
                    print("     ~ NOTE " + n[:100])
                continue
            for i, chunk in enumerate(chunks(body)):
                verbatim = title in VERBATIM_SECTIONS
                table = None if verbatim else parse_table(chunk)
                if verbatim and i == 0:
                    print("     [verbatim] " + chunk[0].strip()[:90])
                    continue
                if table is None:
                    print("     %s %s" % ("[verbatim]" if is_structured(chunk) else "[prose]",
                                          chunk[0].strip()[:90]))
                    continue
                print("     columns: %s" % table["columns"])
                for r in table["rows"]:
                    print("     row %s" % " | ".join(c[:46] for c in r))
                for n in table["notes"]:
                    print("     ~ NOTE " + n[:100])
    return 0


def main(argv):
    if "--dump-api" in argv:
        return dump_api()
    pages = render_all()
    targets = [(os.path.join(HTML, name), name, text) for name, text in sorted(pages.items())]
    targets.append((TEX_GENERATED, "src/generated.tex", render_generated_tex()))
    targets.append((TEX_API, "src/api_generated.tex", render_api_tex()))
    targets.append((TEX_TESTS, "src/tests_generated.tex", render_tests_tex()))
    stale = []
    for path, name, text in targets:
        current = open(path, encoding="utf-8").read() if os.path.isfile(path) else ""
        if current != text:
            stale.append(name)
    if "--check" in argv:
        if stale:
            print("buildDocs: stale -> %s (run python tools/buildDocs.py)" % ", ".join(stale))
            return 1
        print("buildDocs: all %d generated files current" % len(targets))
        return 0
    for path, _, text in targets:
        open(path, "w", encoding="utf-8", newline="\n").write(text)
    print("buildDocs: wrote %d files (%s)" % (len(targets), ", ".join(sorted(stale)) if stale else "no change"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
