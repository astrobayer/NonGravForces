function info = setupNonGravForces(options)
%SETUPNONGRAVFORCES Put the NonGravForces toolbox on the path and check dependencies.
%
%   setupNonGravForces
%   info = setupNonGravForces('Verbose', false)
%
%   Adds the toolbox root, tests/ and tools/ to the MATLAB path (no savepath),
%   registers the documentation for the Help browser, and reports which of the
%   four ecosystem dependencies are reachable. It never adds a dependency to
%   the path itself: each dependency ships its own setup function and owns its
%   own path.
%
% INPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   (none)
%
% NAME-VALUE ARGUMENTS
%   Verbose   1x1   logical  -      Print the dependency table (default true).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   info      1x1   struct   -      .root (char), .version (char),
%                                   .dependencies (struct array with fields
%                                   name, required, present, version, setup),
%                                   .mex (1xT NonGravForces.buildMex('Check') report,
%                                   one element per MEX target: msis, kernel).
%
% ERRORS
%   None raised for absent dependencies: absence is REPORTED (rule R4/R5).
%
% EXAMPLE
%   info = setupNonGravForces;
%   info.dependencies(1)
%
% See also NonGravForces.dependencyReport, runNgfTests.
%
% KEYWORDS setup path dependencies
% Created 2026-08-27 by Claude Fable 5 (Anthropic) for Matthias Weigelt.
arguments
    options.Verbose (1,1) logical = true
end
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root,'tests'), fullfile(root,'tools'));
% Help-browser registration: info.xml at the root is picked up by doc.
try
    builddocsearchdb(fullfile(root,'doc','html'));
catch
    % Search database is optional; the pages are readable without it.
end
info.root    = root;
info.version = NonGravForces.toolboxVersion('string');
info.dependencies = NonGravForces.dependencyReport('Verbose', options.Verbose);
info.mex = NonGravForces.buildMex('Check', true);      % 1xT since 0.6.1: msis, kernel
if options.Verbose
    for m = 1:numel(info.mex)
        if info.mex(m).upToDate, state = 'built, current';
        elseif info.mex(m).present, state = 'built, STALE (sources changed): run NonGravForces.buildMex';
        else, state = 'NOT BUILT: run NonGravForces.buildMex (the caller of this gateway refuses until then)';
        end
        [~, gw] = fileparts(info.mex(m).file);
        fprintf('MEX  %-14s %s\n', gw, state);
    end
end
end
