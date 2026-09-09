function info = setupNonGravForces(options)
%SETUPNONGRAVFORCES Put the NonGravForces toolbox on the path and check dependencies.
%
%   setupNonGravForces
%   info = setupNonGravForces('Verbose', false)
%
%   Adds the toolbox root, tests/, tools/ and demos/ to the MATLAB path (no savepath),
%   registers the documentation for the Help browser, and reports which of the
%   ecosystem dependencies are reachable. It never adds a dependency to the path
%   itself: each dependency ships its own setup function and owns its own path.
%
%   IT ALSO FETCHES THE NRLMSISE-00 C SOURCE WHEN THAT IS MISSING, and this is
%   the one place in the toolbox that touches the network. The source is not
%   redistributed here -- its package carries no such grant
%   (THIRD-PARTY-NOTICES.md, obligation O-34) -- so it is cloned from its
%   maintainer into mex/nrlmsise00 the first time this runs. An installer is the
%   right place for that; NonGravForces.buildMex is not, because a BUILD step
%   that reaches the network without being asked is a build step nobody can
%   audit. The clone prints one line saying what it did and from where, whatever
%   'Verbose' says, because a network action has to be visible.
%
%   Nothing is fetched when the folder is already there, when NGF_NRLMSISE00
%   points at a copy, or with 'FetchMsis', "never". A failed clone is REPORTED,
%   not raised: setup reports and buildMex refuses, which keeps one error for
%   one cause.
%
% INPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   (none)
%
% NAME-VALUE ARGUMENTS
%   Verbose    1x1  logical  -      Print the dependency table (default true).
%   FetchMsis  1x1  string   -      "ifMissing" (default) clones the NRLMSISE-00 source only
%                                   when it is absent; "never" leaves the network alone;
%                                   "force" re-clones and moves to the pinned commit even if a
%                                   copy is there.
%
% OUTPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   info      1x1   struct   -      .root (char), .version (char),
%                                   .dependencies (struct array with fields
%                                   name, required, present, version, setup),
%                                   .mex (1xT NonGravForces.buildMex('Check') report,
%                                   one element per MEX target: msis, kernel),
%                                   .msis (.dir, .action "present"|"fetched"|"skipped"|"failed",
%                                   .url, .commit, .message).
%
% ERRORS
%   None raised for absent dependencies: absence is REPORTED (rule R4/R5).
%
% EXAMPLE
%   info = setupNonGravForces;
%   info.dependencies(1)
%
% See also NonGravForces.dependencyReport, NonGravForces.buildMex, runNgfTests.
%
% KEYWORDS setup path dependencies
% Created 2026-08-27 by Claude Fable 5 (Anthropic) for Matthias Weigelt.
arguments
    options.Verbose (1,1) logical = true
    options.FetchMsis (1,1) string {mustBeMember(options.FetchMsis, ["ifMissing","never","force"])} = "ifMissing"
end
% THE PIN LIVES HERE AND NOWHERE ELSE (rule R1). The owner's ruling of 2026-09-09 is to track
% upstream and take its bug fixes; pinning a COMMIT rather than a branch is what makes that
% deliberate instead of at the mercy of someone else's push. Moving it is an act with a run
% behind it. bc9a2fe is "gtd7: fix output at 32.5 km" and carries, over the 2007 release this
% toolbox used until 2026-09-09, fixes for Argon densities and a memory corruption.
MSIS_URL = "https://git.linta.de/~brodo/nrlmsise-00.git";
MSIS_COMMIT = "bc9a2feba4344e74201281e563332688a4d09cc3";
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root,'tests'), fullfile(root,'tools'), fullfile(root,'demos'));
% Help-browser registration: info.xml at the root is picked up by doc.
try
    builddocsearchdb(fullfile(root,'doc','html'));
catch
    % Search database is optional; the pages are readable without it.
end
info.root    = root;
info.version = NonGravForces.toolboxVersion('string');
info.dependencies = NonGravForces.dependencyReport('Verbose', options.Verbose);
info.msis = localFetchMsis(root, MSIS_URL, MSIS_COMMIT, options.FetchMsis);
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
    fprintf('MSIS %-14s %s (%s)\n', 'nrlmsise00', info.msis.action, info.msis.dir);
end
end

function m = localFetchMsis(root, url, commit, mode)
%LOCALFETCHMSIS Obtain the NRLMSISE-00 C source if it is not already here.
% Absence is reported rather than raised (R4/R5): buildMex is the unit that refuses, and one
% cause should produce one error and not two.
m = struct('dir', "", 'action', "skipped", 'url', url, 'commit', commit, 'message', "");
env = string(getenv('NGF_NRLMSISE00'));
if strlength(env) > 0
    m.dir = env;
    m.action = "skipped";
    m.message = "NGF_NRLMSISE00 points at a copy, so nothing is fetched";
    return
end
dst = string(fullfile(root, 'mex', 'nrlmsise00'));
m.dir = dst;
here = isfile(fullfile(dst, 'nrlmsise-00.c'));
if mode == "never" || (here && mode ~= "force")
    if here, m.action = "present"; else, m.action = "skipped"; end
    return
end
if here && mode == "force"
    rmdir(dst, 's');
end
% One line, whatever Verbose says: a network action has to be visible.
fprintf('setupNonGravForces: fetching NRLMSISE-00 from %s (commit %s)\n', url, extractBefore(commit, 8));
% No --depth: this server speaks dumb HTTP, which does not support shallow clones (measured
% 2026-09-09). The repository is a few hundred kilobytes, so a full clone costs nothing.
[st, out] = system(sprintf('git clone --quiet "%s" "%s"', url, dst));
if st == 0
    [st, out] = system(sprintf('git -C "%s" checkout --quiet %s', dst, commit));
end
if st ~= 0
    if isfolder(dst), rmdir(dst, 's'); end        % never leave a half-clone for buildMex to digest
    m.action = "failed";
    m.message = "git failed: " + string(strtrim(out)) + " -> clone it by hand: git clone " + url + " mex/nrlmsise00";
    fprintf('setupNonGravForces: %s\n', m.message);
    return
end
m.action = "fetched";
end
