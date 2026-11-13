function report = runTutorials(names, options)
%RUNTUTORIALS Run the first-time tutorials cell by cell and record their code, output and figures.
%
%   report = runTutorials
%   report = runTutorials(names)
%   report = runTutorials(names, 'OutDir', folder, 'FigureSize', [w h])
%
%   The tutorials tutorials/ngfTutorial*.m are documentation in three places: the scripts
%   themselves, their published pages in the Help browser (tools/publishTutorials.m) and the
%   tutorials chapter of the user guide PDF. A PDF cannot run MATLAB, so this tool runs each
%   tutorial ONCE, cell by cell (a cell starts at a line beginning with %%), and records per cell
%   the code that ran, what it printed (evalc) and every figure it created (PNG at 200 dpi: a
%   built PDF is never committed, auditPackage.py check 6b, so the figures are raster).
%   tools/makeTutorialsTex.py builds the chapter from that record and REFUSES a record whose code
%   differs from the script (full-line comments and blank lines excepted), so the PDF never shows
%   output of code that is no longer there.
%
%   Each tutorial runs in a workspace of its own (a local function's; the caller's workspace is
%   neither seen nor changed), with every figure closed first and figures invisible. A cell that
%   errors stops the run; the record of that tutorial is then not written.
%
%   Lives in tools/: a documentation instrument, not part of the scientific surface, not gated
%   (the gate RUNS the tutorials through tests/tNgfTutorials.m; this tool records them).
%   .github/workflows/tutorials.yml runs it on the CI runner and commits the record.
%
% INPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   names  1xK   string  -      tutorial names without .m; default: every tutorials/ngfTutorial*.m
%
% NAME-VALUE ARGUMENTS
%   OutDir      1x1  string  -    default doc/src/tutorials
%   FigureSize  1x2  double  px   figure width and height; default [760 420]
%
% OUTPUT ARGUMENTS
%   report  1xK struct  .name .cells .seconds .folder; files per tutorial in <OutDir>/<name>/:
%                       run.json (meta: toolbox, MATLAB, date, computer; cells: title, code,
%                       output, figures, seconds) and cellNN_figM.png
%
% ERRORS
%   Errors of the tutorials propagate.
%
% Created 2026-11-12 by Claude (Anthropic) for Matthias Weigelt, after lgdToolbox's
% docsrc/runTutorials.m (R-055).
arguments
    names (1,:) string = string.empty(1, 0)
    options.OutDir (1,1) string = ""
    options.FigureSize (1,2) double {mustBePositive} = [760 420]
end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'tutorials'));
if isempty(names)
    d = dir(fullfile(root, 'tutorials', 'ngfTutorial*.m'));
    names = sort(string(erase({d.name}, ".m")));
end
if options.OutDir == ""
    options.OutDir = fullfile(root, 'doc', 'src', 'tutorials');
end
try
    feature('hotlinks', false);
catch
end
set(groot, 'DefaultFigureVisible', 'off', 'DefaultFigurePosition', [100 100 options.FigureSize]);
restore = onCleanup(@() set(groot, 'DefaultFigureVisible', 'remove', 'DefaultFigurePosition', 'remove'));
v = NonGravForces.toolboxVersion();
report = struct('name', {}, 'cells', {}, 'seconds', {}, 'folder', {});
for name = names
    out = fullfile(options.OutDir, name);
    cells = localSplitCells(fullfile(root, 'tutorials', name + ".m"));
    close all force
    total = tic;
    [rec, figs] = localRunCells(cells);
    % written only after the whole tutorial ran
    if isfolder(out), rmdir(out, 's'); end
    mkdir(out);
    for j = 1:numel(figs)
        exportgraphics(figs{j}{1}, fullfile(out, figs{j}{2}), 'Resolution', 200);
    end
    meta = struct('toolbox', string(v.semantic), 'matlab', string(version), 'date', ...
        string(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm ''UTC''')), ...
        'computer', string(computer));
    fid = fopen(fullfile(out, 'run.json'), 'w', 'n', 'UTF-8');
    fprintf(fid, '%s', jsonencode(struct('name', name, 'meta', meta, 'cells', rec), 'PrettyPrint', true));
    fclose(fid);
    report(end + 1) = struct('name', name, 'cells', numel(cells), 'seconds', toc(total), 'folder', out); %#ok<AGROW>
    close all force
end
clear restore
end

function [ngfTut__rec, ngfTut__figs] = localRunCells(ngfTut__cells)
%LOCALRUNCELLS The cells one after the other in THIS function's workspace: the tutorial's variables
%persist from cell to cell and nothing of the caller's workspace is seen or changed. The function's
%own variables carry a prefix no tutorial uses.
ngfTut__rec = struct('title', {}, 'code', {}, 'output', {}, 'figures', {}, 'seconds', {});
ngfTut__figs = {};
for ngfTut__k = 1:numel(ngfTut__cells)
    ngfTut__before = findall(groot, 'Type', 'figure');
    ngfTut__t0 = tic;
    ngfTut__txt = evalc(ngfTut__cells(ngfTut__k).code);
    ngfTut__sec = toc(ngfTut__t0);
    ngfTut__after = findall(groot, 'Type', 'figure');
    ngfTut__new = ngfTut__after(~ismember(ngfTut__after, ngfTut__before));
    if ~isempty(ngfTut__new)
        [~, ngfTut__order] = sort([ngfTut__new.Number]);
        ngfTut__new = ngfTut__new(ngfTut__order);
    end
    ngfTut__files = strings(1, numel(ngfTut__new));
    for ngfTut__j = 1:numel(ngfTut__new)
        ngfTut__files(ngfTut__j) = sprintf('cell%02d_fig%d.png', ngfTut__k, ngfTut__j);
        ngfTut__figs{end + 1} = {ngfTut__new(ngfTut__j), ngfTut__files(ngfTut__j)}; %#ok<AGROW>
    end
    ngfTut__rec(ngfTut__k) = struct('title', ngfTut__cells(ngfTut__k).title, ...
        'code', ngfTut__cells(ngfTut__k).code, 'output', string(ngfTut__txt), ...
        'figures', ngfTut__files, 'seconds', ngfTut__sec);
end
end

function cells = localSplitCells(file)
%LOCALSPLITCELLS Cells of a script: a cell starts at a line beginning with %% (its title follows on that
%line). Line endings normalised to LF, trailing blank lines dropped. tools/makeTutorialsTex.py splits
%with the same rule.
lines = splitlines(string(fileread(file)));
lines = erase(lines, char(13));
start = find(startsWith(lines, "%%"));
cells = struct('title', {}, 'code', {});
for k = 1:numel(start)
    last = numel(lines);
    if k < numel(start), last = start(k + 1) - 1; end
    block = lines(start(k):last);
    while ~isempty(block) && strtrim(block(end)) == ""
        block(end) = [];
    end
    cells(k).title = strtrim(extractAfter(block(1), 2));
    cells(k).code = strjoin(block, newline);
end
end
