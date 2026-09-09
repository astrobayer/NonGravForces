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
%   IT ALSO FETCHES WHAT THE TOOLBOX NEEDS AND DOES NOT SHIP, and it is the one
%   place in the toolbox that touches the network: the NRLMSISE-00 C source, the
%   IERS Earth-orientation file and the CelesTrak space-weather file. Everything
%   it fetches is optional and only when missing; what it cannot fetch, it names.
%
%   GRACE-FO LEVEL-1B IS NOT FETCHABLE HERE. PO.DAAC requires an Earthdata login
%   and no anonymous mirror was established (probed 2026-09-09), so setup reports
%   where the files come from and what to call them, and never pretends. Beyond
%   that, supplying data is the user's responsibility -- the owner's ruling of
%   2026-09-09. The source is not
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
%   FetchData  1x1  string   -      Same three settings for the Earth-orientation and space-weather
%                                   files. GRACE-FO Level-1B is NOT downloadable without an
%                                   Earthdata account, so it is reported and never fetched.
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
    options.FetchData (1,1) string {mustBeMember(options.FetchData, ["ifMissing","never","force"])} = "ifMissing"
end
% THE PIN LIVES HERE AND NOWHERE ELSE (rule R1). The owner's ruling of 2026-09-09 is to track
% upstream and take its bug fixes; pinning a COMMIT rather than a branch is what makes that
% deliberate instead of at the mercy of someone else's push. Moving it is an act with a run
% behind it. bc9a2fe is "gtd7: fix output at 32.5 km" and carries, over the 2007 release this
% toolbox used until 2026-09-09, fixes for Argon densities and a memory corruption.
MSIS_URL = "https://git.linta.de/~brodo/nrlmsise-00.git";
MSIS_COMMIT = "bc9a2feba4344e74201281e563332688a4d09cc3";
% The two data sources that can be fetched WITHOUT an account, probed 2026-09-09. Both were
% downloaded end to end before these lines were written; CelesTrak answered 503 once and 200 on
% the retry, which is why localFetchData retries. GRACE-FO Level-1B is NOT here: it needs an
% Earthdata login at PO.DAAC and no anonymous mirror was established, so setup reports where to
% get it instead of pretending it can.
EOP_URL = "https://datacenter.iers.org/data/9/finals2000A.all";
SW_URL = "https://celestrak.org/SpaceData/SW-All.csv";
% GRACE-FO Level-1B IS anonymously available after all, at GFZ's ISDC -- the owner supplied the
% path on 2026-09-09 after my probe went to the wrong host and concluded it needed an Earthdata
% login. Verified the same day: the directory lists one archive per day and the 2020-01-01 one
% holds every product under exactly the filenames this toolbox's fixtures use. It is 147 MB, so
% it is fetched only when the sample files are missing and only the needed members are kept.
L1B_BASE = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04";
L1B_DAY = "2020-01-01";
L1B_MEMBERS = ["GNV1B_2020-01-01_C_04.txt", "SCA1B_2020-01-01_C_04.txt"];
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
info.data = localFetchData(root, EOP_URL, SW_URL, L1B_BASE, L1B_DAY, L1B_MEMBERS, options.FetchData);
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
    for d = 1:numel(info.data.items)
        fprintf('DATA %-14s %s\n', info.data.items(d).name, info.data.items(d).action);
    end
end
end

function d = localFetchData(root, eopUrl, swUrl, l1bBase, l1bDay, l1bMembers, mode)
%LOCALFETCHDATA The two files that can be had without an account, and an honest word about the third.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
targets = {"eop", fullfile(root, 'data', 'eop'), "finals2000A.all", eopUrl, "*.all"; ...
           "spaceweather", fullfile(root, 'data', 'spaceweather'), "SW-All.csv", swUrl, "*.csv"};
for t = 1:size(targets, 1)
    folder = targets{t, 2};
    present = ~isempty(dir(fullfile(folder, targets{t, 5})));
    rec = struct('name', targets{t, 1}, 'action', "present", 'file', "", 'source', targets{t, 4});
    if mode == "never" || (present && mode ~= "force")
        if ~present, rec.action = "missing (FetchData never)"; end
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if ~isfolder(folder), mkdir(folder); end
    dst = fullfile(folder, char(targets{t, 3}));
    fprintf('setupNonGravForces: fetching %s from %s\n', targets{t, 1}, targets{t, 4});
    ok = false;
    for attempt = 1:2                     % CelesTrak answered 503 once and 200 on the retry
        try
            websave(dst, char(targets{t, 4}));
            ok = true;
            break
        catch err
            rec.action = "failed: " + string(err.message);
        end
    end
    if ok
        rec.action = "fetched";
        rec.file = string(dst);
    else
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, targets{t, 4});
    end
    items(end+1) = rec; %#ok<AGROW>
end
% GRACE-FO Level-1B: one archive per day, 147 MB, so only fetched when the members are absent
% and only the members are kept. The archive holds far more than this toolbox reads --
% ACT1B, MAG1B, MAS1B, THR1B among them -- which is worth knowing before writing a reader.
samples = fullfile(root, 'data', 'samples');
url = l1bBase + "/" + extractBefore(l1bDay, 5) + "/gracefo_1B_" + l1bDay + "_RL04.ascii.noLRI.tgz";
missing = l1bMembers(~isfile(fullfile(samples, l1bMembers)));
rec = struct('name', "samples", 'action', "present", 'file', "", 'source', url);
if mode == "never" || (isempty(missing) && mode ~= "force")
    if ~isempty(missing), rec.action = "missing (FetchData never)"; end
else
    if ~isfolder(samples), mkdir(samples); end
    tmp = [tempname '.tgz'];
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching GRACE-FO Level-1B %s (about 147 MB) from %s\n', l1bDay, url);
    try
        websave(tmp, char(url));
        untar(tmp, samples, cellstr(l1bMembers));
        got = l1bMembers(isfile(fullfile(samples, l1bMembers)));
        if numel(got) == numel(l1bMembers)
            rec.action = "fetched";
        else
            rec.action = "incomplete: the archive did not hold " + strjoin(setdiff(l1bMembers, got), ", ");
        end
    catch err
        rec.action = "failed: " + string(err.message);
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, url);
    end
end
items(end+1) = rec;
d = struct('items', items, 'mode', mode);
end

function localQuietDelete(f)
%LOCALQUIETDELETE Remove a temporary file without turning a cleanup into an error.
if isfile(f)
    try, delete(f); catch, end
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
