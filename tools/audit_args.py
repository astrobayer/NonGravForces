"""Compare every DECLARED argument of every public unit against its DOCUMENTED one.

auditPackage check 2 verifies that INPUT/NAME-VALUE/OUTPUT sections EXIST. It has never verified
that they describe the arguments the function actually takes, which is the same shape of gap this
level kept finding: a fact written down twice and compared never.

PARSER NOTES, each earned by a false positive on the first run:
  - the name column is not a fixed width and many rows use a SINGLE space after the name
    ("MsisSource 1x1 string"), so a >=2-space assumption loses real rows;
  - "Name" is both the header of the first column AND a real option of geom.boxWing, so the header
    is identified by the whole line (Name + Size + Type), never by the token alone;
  - several rows group names ("Bias, Scale, Whitening, Window  as qc.residuals"), which is terser
    than the template but IS documentation, so commas are split rather than flagged.
A checker that reports what is already there teaches people to ignore it.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, '+NonGravForces')

SECTION = {'in': 'INPUT ARGUMENTS', 'nv': 'NAME-VALUE ARGUMENTS', 'out': 'OUTPUT ARGUMENTS'}


def help_names(src, title):
    # the heading may carry a parenthetical ('INPUT ARGUMENTS (all 1xN or Nx1, same N)')
    # The TITLE may carry a parenthetical; the section END may not. A heading is alone on its line,
    # and relaxing that made rows named F, T or GM look like new sections and truncated the table.
    m = re.search(r'^%\s*' + title + r'(?:\s*\([^)]*\))?\s*$(.*?)(?=^%\s*[A-Z][A-Z \-]{3,}\s*$|\Z)',
                  src, re.M | re.S)
    if not m:
        return None
    names = []
    for line in m.group(1).splitlines():
        body = line[1:] if line.startswith('%') else line
        if re.search(r'\bName\b', body) and re.search(r'\bSize\b', body) and re.search(r'\bType\b', body):
            continue                                        # the header row, whatever it is called
        if set(body.strip()) <= set('- '):
            continue                                        # the dashes row
        g = re.match(r'^\s{2,5}([A-Za-z][\w.]*(?:\s*,\s*[A-Za-z][\w.]*)*)\s+(\S+)', body)
        if g:
            # a PROSE row inside the table ("All arrays must be compatible ...") is not an argument:
            # a real row's second column is a size or a type, never an ordinary word.
            second = g.group(2)
            if not re.match(r'^(\d+|[A-Za-z])x[\dA-Za-z]|^any$|^1x1\|1xK$|^--|^\(|^struct$|^double$|^string$|^logical$|^table$|^cell$', second):
                continue
            names.extend(n.strip() for n in g.group(1).split(','))
    return names


def declared(src):
    sig = re.search(r'^\s*function\s+(?:(?:\[([^\]]*)\]|(\w+))\s*=\s*)?(\w+)\s*\(([^)]*)\)', src, re.M)   # output clause optional as a whole (see audit_class)
    if not sig:
        return None
    outs = [o.strip() for o in sig.group(1).split(',')] if sig.group(1) else ([sig.group(2)] if sig.group(2) else [])
    outs = [o for o in outs if o]
    ins = [a.strip() for a in sig.group(4).split(',') if a.strip()]
    block = re.search(r'^arguments\b(.*?)^end\s*$', src, re.M | re.S)
    nvs, pos, containers = [], [], set()
    if block:
        for line in block.group(1).splitlines():
            line = line.split('%')[0].strip()
            if not line:
                continue
            g = re.match(r'^(\w+)\.(\w+)', line)
            if g:
                containers.add(g.group(1))
                nvs.append(g.group(2))
                continue
            g = re.match(r'^(\w+)', line)
            if g:
                pos.append(g.group(1))
    if not pos:
        pos = list(ins)
    pos = [p for p in pos if p not in containers]            # only the NAME-VALUE container is dropped
    return outs, pos, nvs


def main():
    findings, files = [], []
    for dirpath, _d, names in os.walk(PKG):
        if os.path.basename(dirpath) == 'private':
            continue
        files += [os.path.join(dirpath, n) for n in sorted(names) if n.endswith('.m') and n != 'Contents.m']

    for path in files:
        rel = os.path.relpath(path, ROOT)
        src = open(path, encoding='utf-8').read()
        if src.lstrip().startswith('classdef'):
            continue
        d = declared(src)
        if d is None:
            continue
        outs, pos, nvs = d
        for kind, decl in (('in', pos), ('nv', nvs), ('out', outs)):
            doc = help_names(src, SECTION[kind])
            if doc is None:
                if decl:
                    findings.append((rel, f'no {SECTION[kind]} block, but declares {", ".join(decl)}'))
                continue
            # varargin is documented by the names it carries, not by the word 'varargin'
            if decl == ['varargin'] and doc:
                continue
            # an output documented by its FIELDS (c.S0, c.kB, ...) is documented
            if kind == 'out' and len(decl) == 1 and any(n == decl[0] for n in doc + [d.split('.')[0] for d in doc]):
                continue
            miss = [x for x in decl if x not in doc]
            extra = [x for x in doc if x not in decl and x not in ('varargin', 'args')]
            if miss:
                findings.append((rel, f'{SECTION[kind]}: declared but NOT documented: {", ".join(miss)}'))
            if extra:
                findings.append((rel, f'{SECTION[kind]}: documented but NOT declared: {", ".join(extra)}'))

    for rel, msg in findings:
        print(f'{rel}: {msg}')
    print(f'\n{len(findings)} finding(s) over {len(files)} file(s)')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
