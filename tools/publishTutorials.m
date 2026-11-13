function files = publishTutorials(names)
%PUBLISHTUTORIALS Run and publish the first-time tutorials into doc/html/tutorials for the Help browser.
%
%   files = publishTutorials
%   files = publishTutorials(["ngfTutorial01FirstSteps" "ngfTutorial06ForceModel"])
%
%   Publishing RUNS each tutorial, so a published page is proof that it executed on this MATLAB
%   release; a broken tutorial fails here instead of in front of a user. The pages go to
%   doc/html/tutorials/<name>.html with their figures beside them; doc/html/Tutorials.html (written
%   by tools/buildDocs.py) links them, and doc/html/helptoc.xml lists them under "Tutorials".
%   .github/workflows/tutorials.yml runs it on the CI runner together with tools/runTutorials.m and
%   commits both outputs.
%
%   Lives in tools/: a documentation instrument, not part of the scientific surface, not gated.
%
% INPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   names  1xK   string  -      tutorial names without .m; default: every tutorials/ngfTutorial*.m
%
% OUTPUT ARGUMENTS
%   files  1xK string  paths of the published pages
%
% ERRORS
%   Errors of the tutorials propagate (publish stops at the first one).
%
% Created 2026-11-12 by Claude (Anthropic) for Matthias Weigelt, after lgdToolbox's
% docsrc/publishExamples.m (R-055).
arguments
    names (1,:) string = string.empty(1, 0)
end
root = fileparts(fileparts(mfilename('fullpath')));
src = fullfile(root, 'tutorials');
addpath(src);
if isempty(names)
    d = dir(fullfile(src, 'ngfTutorial*.m'));
    names = sort(string(erase({d.name}, ".m")));
end
out = fullfile(root, 'doc', 'html', 'tutorials');
if ~isfolder(out), mkdir(out); end
set(groot, 'DefaultFigureVisible', 'off');
restore = onCleanup(@() set(groot, 'DefaultFigureVisible', 'remove'));
files = strings(1, numel(names));
for k = 1:numel(names)
    close all force
    files(k) = string(publish(fullfile(src, names(k) + ".m"), 'format', 'html', 'outputDir', out, ...
        'catchError', false, 'imageFormat', 'png'));
    close all force
end
clear restore
end
