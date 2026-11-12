"""makeTutorialsTex.py - the tutorials chapter of the user guide PDF (first-time users).

Builds doc/src/tutorials/tutorials.tex from the tutorial scripts tutorials/ngfTutorial*.m and the record
of their last MATLAB run, doc/src/tutorials/<name>/run.json (written by tools/runTutorials.m): per cell the
explanation (the cell's leading comment block), the code, the output it printed and the figures it created.

The PDF must never show output of code that is no longer in the script: a record whose CODE (the cell with
full-line comments and blank lines removed) differs from the script's is refused. Edits to the explanations
need no new run.

  python3 tools/makeTutorialsTex.py          write doc/src/tutorials/tutorials.tex
  python3 tools/makeTutorialsTex.py --check  exit 1 if a record is stale or missing, or if tutorials.tex
                                             is not what this script writes (auditPackage.py check 15)
  python3 tools/makeTutorialsTex.py --allow-missing
                                             write it from the records there are (a draft build while the
                                             records are made on the runner, .github/workflows/tutorials.yml)

Ported 2026-11-12 by Claude (Anthropic) for Matthias Weigelt from lgdToolbox's docsrc/make_tutorials_tex.py
(its R-055); standard library only, like every tool the audit job runs.
"""
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TUT = os.path.join(ROOT, 'doc', 'src', 'tutorials')
OUT = os.path.join(TUT, 'tutorials.tex')

ASCII = {'\u00b5': 'u', '\u03bc': 'u', '\u00d7': 'x', '\u2212': '-', '\u00b0': ' deg', '\u2013': '-',
         '\u2014': '-', '\u00b1': '+-', '\u2264': '<=', '\u2265': '>=', '\u00a0': ' '}


def tutorials():
    return sorted(os.path.splitext(os.path.basename(p))[0]
                  for p in glob.glob(os.path.join(ROOT, 'tutorials', 'ngfTutorial*.m')))


def split_cells(path):
    """The rule of runTutorials.m: a cell starts at a line beginning with %%."""
    lines = open(path, encoding='utf-8').read().replace('\r', '').split('\n')
    starts = [i for i, ln in enumerate(lines) if ln.startswith('%%')]
    cells = []
    for k, i in enumerate(starts):
        j = starts[k + 1] if k + 1 < len(starts) else len(lines)
        block = lines[i:j]
        while block and block[-1].strip() == '':
            block.pop()
        cells.append({'title': block[0][2:].strip(), 'code': '\n'.join(block)})
    return cells


def code_only(text):
    return [ln.rstrip() for ln in text.replace('\r', '').split('\n')
            if ln.strip() and not ln.strip().startswith('%')]


def parts(cell):
    """(prose lines, code lines): the leading comment block after the title, the rest."""
    lines = cell['code'].split('\n')[1:]
    k = 0
    while k < len(lines) and lines[k].startswith('%'):
        k += 1
    prose = [re.sub(r'^% ?', '', ln) for ln in lines[:k]]
    code = lines[k:]
    while code and code[0].strip() == '':
        code.pop(0)
    while code and code[-1].strip() == '':
        code.pop()
    return prose, code


def esc(s):
    s = ''.join(ASCII.get(ch, ch) for ch in s)
    rep = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '#': r'\#', '_': r'\_', '{': r'\{',
           '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}', '<': r'\textless{}',
           '>': r'\textgreater{}'}
    return ''.join(rep.get(ch, ch) for ch in s)


def inline(s):
    """Publish markup to LaTeX: |code|, *bold*, $math$ (kept), the rest escaped."""
    out = []
    for tok in re.split(r'(\$[^$]+\$|\|[^|]+\||\*[^*\s][^*]*\*)', s):
        if not tok:
            continue
        if tok.startswith('$') and tok.endswith('$'):
            out.append(tok)
        elif tok.startswith('|') and tok.endswith('|') and len(tok) > 2:
            out.append(r'\texttt{' + esc(tok[1:-1]) + '}')
        elif tok.startswith('*') and tok.endswith('*') and len(tok) > 2:
            out.append(r'\textbf{' + esc(tok[1:-1]) + '}')
        else:
            # m/s^2, 1e-9, (R/r)^(n+1): a power written with ^ becomes a superscript
            t = esc(tok)
            t = re.sub(r'"([^"]*)"', r"``\1''", t)    # "prior" -> ``prior''
            t = re.sub(r'\\textasciicircum\{\}\(([^()]*)\)', r'\\textsuperscript{\1}', t)
            t = re.sub(r'\\textasciicircum\{\}(-?[0-9A-Za-z.]+)', r'\\textsuperscript{\1}', t)
            out.append(t)
    return ''.join(out)


def prose_tex(lines):
    """Paragraphs, bullet lists (* item, continued by indented lines) and $$ display math $$."""
    lines = [ln for ln in lines if not re.match(r'^(Requirements:|NonGravForces v\d)', ln.strip())]
    out, para, items, math = [], [], [], None

    def flush():
        nonlocal para, items
        if para:
            out.append(inline(' '.join(p.strip() for p in para)) + '\n')
            para = []
        if items:
            out.append('\\begin{itemize}\\setlength\\itemsep{0pt}\n' +
                       ''.join('\\item ' + inline(it) + '\n' for it in items) + '\\end{itemize}\n')
            items = []
    for ln in lines:
        st = ln.strip()
        if math is not None:
            math.append(st)
            if st.endswith('$$'):
                out.append('\\[' + ' '.join(math)[2:-2] + '\\]\n')
                math = None
            continue
        if st.startswith('$$'):
            flush()
            if st.endswith('$$') and len(st) > 4:
                out.append('\\[' + st[2:-2] + '\\]\n')
            else:
                math = [st]
            continue
        if st == '':
            flush()
            continue
        if st.startswith('* '):
            if para:
                flush()
            items.append(st[2:])
            continue
        if items and ln.startswith('  '):
            items[-1] += ' ' + st
            continue
        if items:
            flush()
        para.append(st)
    flush()
    return '\n'.join(out)


def clean_output(txt):
    txt = re.sub(r'</?strong>', '', txt)
    txt = re.sub(r'<a href="[^"]*">([^<]*)</a>', r'\1', txt)
    txt = ''.join(ASCII.get(ch, ch) if ord(ch) > 127 else ch for ch in txt)
    txt = ''.join(ch if ord(ch) < 128 else '?' for ch in txt)
    lines = txt.replace('\r', '').split('\n')
    while lines and lines[0].strip() == '':
        lines.pop(0)
    while lines and lines[-1].strip() == '':
        lines.pop()
    # MATLAB's display puts blank lines around every table and struct: keep at most one
    out = []
    for ln in lines:
        if ln.strip() == '' and out and out[-1].strip() == '':
            continue
        out.append(ln.rstrip())
    return '\n'.join(out)


def as_list(x):
    if x is None:
        return []
    return [x] if isinstance(x, str) else list(x)


def build():
    problems, tex = [], []
    names = tutorials()
    tex.append('% generated by tools/makeTutorialsTex.py from tutorials/ngfTutorial*.m and\n'
               '% doc/src/tutorials/*/run.json (tools/runTutorials.m) - do not edit\n')
    for name in names:
        cells = split_cells(os.path.join(ROOT, 'tutorials', name + '.m'))
        rp = os.path.join(TUT, name, 'run.json')
        if not os.path.isfile(rp):
            problems.append(f'{name}: no run record ({os.path.relpath(rp, ROOT)}); run runTutorials("{name}") or dispatch release-checks.yml with what tutorials')
            continue
        run = json.load(open(rp, encoding='utf-8'))
        rcells = run['cells'] if isinstance(run['cells'], list) else [run['cells']]
        if len(rcells) != len(cells):
            problems.append(f'{name}: {len(cells)} cells in the script, {len(rcells)} in the record')
            continue
        for k, (c, r) in enumerate(zip(cells, rcells), 1):
            if code_only(c['code']) != code_only(r['code']):
                problems.append(f'{name}: cell {k} ("{c["title"]}") changed since the run')
        m = run['meta']
        num = int(re.search(r'ngfTutorial(\d+)', name).group(1))
        title = re.sub(r'^Tutorial \d+ - ', '', cells[0]['title'])
        tex.append('\\clearpage\\phantomsection\n'
                   f'\\addcontentsline{{toc}}{{section}}{{Tutorial {num}: {inline(title)}}}\n'
                   f'\\section*{{Tutorial {num}: {inline(title)}}}\\label{{tut:{num:02d}}}\n')
        prose, _ = parts(cells[0])
        tex.append(prose_tex(prose))
        tex.append(f'\\noindent{{\\small\\itshape Script \\texttt{{tutorials/{esc(name)}.m}}; '
                   f'the output below was recorded with NonGravForces {esc(m["toolbox"])}, MATLAB '
                   f'{esc(m["matlab"])}, {esc(m["date"])} (\\texttt{{tools/runTutorials.m}}).}}\n')
        for k, (c, r) in enumerate(zip(cells, rcells)):
            if k == 0:
                continue
            prose, code = parts(c)
            tex.append(f'\\subsection*{{{inline(c["title"])}}}\n')
            if prose:
                tex.append(prose_tex(prose))
            if code:
                tex.append('\\begin{lstlisting}[style=tutcode]\n' + '\n'.join(code) + '\n\\end{lstlisting}\n')
            out = clean_output(r.get('output') or '')
            if out:
                tex.append('\\begin{lstlisting}[style=tutout]\n' + out + '\n\\end{lstlisting}\n')
            for f in as_list(r.get('figures')):
                fp = os.path.join(TUT, name, f)
                if not os.path.isfile(fp):
                    problems.append(f'{name}: figure {f} missing')
                tex.append('\\begin{center}\\includegraphics[width=0.92\\linewidth]'
                           f'{{tutorials/{name}/{f}}}\\end{{center}}\n')
    return '\n'.join(tex), problems


def main():
    check = '--check' in sys.argv
    tex, problems = build()
    for p in problems:
        print('STALE:', p)
    if check:
        cur = open(OUT, encoding='utf-8').read() if os.path.isfile(OUT) else ''
        if cur != tex:
            print('STALE: doc/src/tutorials/tutorials.tex differs from what tools/makeTutorialsTex.py writes')
            problems.append('tex')
        if problems:
            sys.exit(1)
        print(f'tutorials: {len(tutorials())} scripts, records and tutorials.tex in sync')
        return
    if problems and '--allow-missing' not in sys.argv:
        sys.exit(1)
    os.makedirs(TUT, exist_ok=True)
    open(OUT, 'w', encoding='utf-8').write(tex)
    print('wrote', os.path.relpath(OUT, ROOT))


if __name__ == '__main__':
    main()
