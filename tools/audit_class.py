"""Audit CLASSES against the convention they actually use.

The first attempt at this applied the FUNCTION convention to class methods and reported 86
findings over 40 methods -- every one of them false. Class methods carry an H1 line only; the
argument tables live in the classdef's own help, as CONSTRUCTOR ARGUMENTS, NAME-VALUE ARGUMENTS
and a METHODS block that gives each method's inputs and outputs. sensor.Electrostatic documents
`dyn`'s fields, `y`'s fields, the (3K)x15 shape of the Jacobian and the frame of every quantity,
which is more than the function template asks for.

So this checks the convention that exists:
  - every PUBLIC method (bar the constructor, accessors and local helpers) is named in METHODS;
  - a constructor taking arguments has a CONSTRUCTOR ARGUMENTS or INPUT ARGUMENTS table;
  - a class declaring name-value arguments anywhere has a NAME-VALUE ARGUMENTS table;
  - every documented method still exists.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, '+NonGravForces')


def class_help(src):
    out = []
    for line in src.splitlines()[1:]:
        if line.lstrip().startswith('%'):
            out.append(line.lstrip()[1:])
        elif line.strip():
            break
    return '\n'.join(out)


def public_methods(src):
    names, lines = [], src.splitlines()
    for i, line in enumerate(lines):
        m = re.match(r'^(\s*)methods\b(.*)$', line)
        if not m:
            continue
        indent, attrs = m.group(1), m.group(2)
        if re.search(r'Access\s*=\s*(private|protected)|Hidden', attrs):
            continue
        for j in range(i + 1, len(lines)):
            if re.match(r'^' + indent + r'end\s*$', lines[j]):
                break
            # The output clause is optional AS A WHOLE: with `(?:...|\w+)?\s*=?` the name of a method
            # without outputs was split, `function save(obj)` reading as a method `e` (2026-10-14).
            g = re.match(r'^\s*function\s+(?:(?:\[[^\]]*\]|\w+)\s*=\s*)?(\w+)\s*\(', lines[j])
            if g and not g.group(1).startswith('local'):
                names.append(g.group(1))
    return names


def main():
    findings, nclass = [], 0
    for dirpath, _d, names in os.walk(PKG):
        if os.path.basename(dirpath) == 'private':
            continue
        for n in sorted(names):
            if not n.endswith('.m'):
                continue
            path = os.path.join(dirpath, n)
            src = open(path, encoding='utf-8').read()
            if not src.lstrip().startswith('classdef'):
                continue
            nclass += 1
            rel, cls = os.path.relpath(path, ROOT), n[:-2]
            doc = class_help(src)
            meths = public_methods(src)
            documented = doc
            for mname in meths:
                if mname == cls:
                    # NOT "is there a table with this heading": the classes spell it CONSTRUCTOR
                    # ARGUMENTS, CONSTRUCTOR INPUTS, or document the arguments in prose and in the
                    # syntax lines (Options, MacroModel). What matters to a reader is whether each
                    # NAME is described somewhere in the class help, so that is what is checked.
                    blk = re.search(r'^\s*arguments\b(.*?)^\s*end\s*$',
                                    src[src.index('function'):], re.M | re.S)
                    if blk:
                        for ln in blk.group(1).splitlines():
                            ln = ln.split('%')[0].strip()
                            g = re.match(r'^(\w+)\.(\w+)', ln)
                            nm = g.group(2) if g else None
                            if nm is None:
                                g2 = re.match(r'^(\w+)', ln)
                                nm = g2.group(1) if g2 else None
                            if nm is None or nm in ('obj', 'varargin', 'options', 'nv', 'opt'):
                                continue
                            if not re.search(r'\b' + re.escape(nm) + r'\b', doc):
                                findings.append((rel, f'constructor argument {nm} is described nowhere in the class help'))
                    continue
                if not re.search(r'\b' + re.escape(mname) + r'\b', documented):
                    findings.append((rel, f'public method {mname}() is named nowhere in the class help'))
                    continue
                # ... and its NAME-VALUE arguments, which a caller can pass and the help may not
                # mention: Model.partials takes WithRespectTo, Method and Frame, and only the first
                # appeared in the class help until 2026-09-18.
                # FROM THE function LINE, not the first textual match: `M.partials(` appears in the
                # class help long before the definition, and starting there read a DIFFERENT
                # method's arguments block and reported everything clean.
                fdef = re.search(r'^\s*function\b[^\n]*\b' + re.escape(mname) + r'\s*\(', src, re.M)
                if not fdef:
                    continue
                seg = src[fdef.end():]
                blk = re.search(r'^\s*arguments\b(.*?)^\s*end\s*$', seg, re.M | re.S)
                if blk:
                    for ln in blk.group(1).splitlines():
                        ln = ln.split('%')[0].strip()
                        g = re.match(r'^(\w+)\.(\w+)', ln)
                        if g and not re.search(r'\b' + re.escape(g.group(2)) + r'\b', documented):
                            findings.append((rel, f'{mname}(): name-value {g.group(2)} is named nowhere in the class help'))
    for rel, msg in findings:
        print(f'{rel}: {msg}')
    print(f'\n{len(findings)} finding(s) over {nclass} class(es)')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
