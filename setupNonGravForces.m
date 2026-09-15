function info = setupNonGravForces(options)
%SETUPNONGRAVFORCES Put the NonGravForces toolbox on the path, fetch what it does not ship, report the dependencies.
%
%   setupNonGravForces
%   info = setupNonGravForces('Verbose', false, 'FetchData', "never")
%
%   Adds the toolbox root, tests/, tools/ and demos/ to the MATLAB path (no savepath),
%   registers the documentation for the Help browser, and reports which of the
%   ecosystem dependencies are reachable. It never adds a dependency to the path
%   itself: each dependency ships its own setup function and owns its own path.
%
%   IT IS THE ONE PLACE IN THE TOOLBOX THAT TOUCHES THE NETWORK. What it fetches,
%   and why it is fetched rather than shipped:
%
%     model sources ('FetchMsis')
%       NRLMSISE-00 C       git clone at a pinned commit into mex/nrlmsise00 (no
%                           redistribution grant, O-34)
%       NRLMSIS 2.1 Fortran NRL's tarball, pinned by SHA-256, unpacked into
%                           mex/nrlmsis21 (Open Source Academic Research License:
%                           research and non-profit use only; the statement is printed
%                           when it is fetched and recorded in THIRD-PARTY-NOTICES.md)
%     data ('FetchData')
%       finals2000A.all     IERS Earth orientation, into data/eop
%       SW-All.csv          CelesTrak space weather, into data/spaceweather
%       igrf14coeffs.txt    IAGA V-MOD IGRF-14 coefficients, into data/igrf
%       GRACE-FO Level-1B   GNV1B, SCA1B and MAG1B of 2020-01-01 (GRACE-C), full day, from the
%                           archive of GFZ's ISDC, into data/l1b -- NOT data/samples, whose
%                           committed ten-minute excerpts carry the same file names and are
%                           the tests' oracles
%
%   PRESENCE IS DECIDED BY THE EXACT FILE, AND A DOWNLOAD BY ITS CONTENT. Until
%   2026-09-15 a glob decided presence, so the committed seven-month excerpts
%   (finals2000A_2019-09_2020-03.all, ...) counted as the full products and the
%   full products were never fetched; and whatever a server returned was kept, an
%   HTML error page included (HANDOVER LP-C, finding F15; CLAUDE.md: "a status code
%   is not a file"). Every fetched file is now checked by its header or digest
%   before it replaces anything, and a failed check is REPORTED, never raised.
%   The Level-1B members are read out of the archive by this function's own tar
%   reader: MATLAB's untar takes no member list (its three-argument call, used
%   until 2026-09-15, raised MATLAB:narginchk:tooManyInputs on every attempt, and
%   the failure was reported as a download failure).
%
%   A build step never fetches: NonGravForces.buildMex refuses when a source is
%   missing and names this function, because a build that reaches the network
%   without being asked is a build nobody can audit. Every network action prints
%   one line whatever 'Verbose' says.
%
% INPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   (none)
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type     Units  Description
%   ----       ----  ----     -----  -----------
%   Verbose    1x1   logical  -      Print the dependency and fetch table (default true).
%   FetchMsis  1x1   string   -      "ifMissing" (default) fetches the model sources only when absent;
%                                    "never" leaves the network alone; "force" fetches again.
%   FetchData  1x1   string   -      The same three settings for the data files.
%
% OUTPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   info      1x1   struct   -      .root (char), .version (char),
%                                   .dependencies (struct array: name, required, present, version, setup),
%                                   .mex (1xT NonGravForces.buildMex('Check') report, one element per target),
%                                   .msis (.dir, .action "present"|"fetched"|"skipped"|"failed", .url, .commit, .message),
%                                   .sources (struct array: name, action, dir, source, message),
%                                   .data (.items struct array: name, action, file, source; .mode;
%                                   the Level-1B item's file lists the members, one per line)
%
% ERRORS
%   None raised for absent dependencies or failed fetches: absence is REPORTED (rule R4/R5).
%
% EXAMPLE
%   info = setupNonGravForces('FetchData', "never");
%   [info.sources.name; info.sources.action]
%
% See also NonGravForces.dependencyReport, NonGravForces.buildMex, runNgfTests.
%
% KEYWORDS setup path dependencies fetch
% Created 2026-08-27 by Claude Fable 5 (Anthropic) for Matthias Weigelt; content-checked fetches and NRLMSIS 2.1 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    options.Verbose (1,1) logical = true
    options.FetchMsis (1,1) string {mustBeMember(options.FetchMsis, ["ifMissing","never","force"])} = "ifMissing"
    options.FetchData (1,1) string {mustBeMember(options.FetchData, ["ifMissing","never","force"])} = "ifMissing"
end
% THE PINS AND URLS LIVE HERE AND NOWHERE ELSE (rule R1). Moving a pin is an act with a run behind it.
% NRLMSISE-00: bc9a2fe is "gtd7: fix output at 32.5 km" (Argon densities and a memory corruption fixed
% since the 2007 release this toolbox used until 2026-09-09).
MSIS_URL = "https://git.linta.de/~brodo/nrlmsise-00.git";
MSIS_COMMIT = "bc9a2feba4344e74201281e563332688a4d09cc3";
% NRLMSIS 2.1: NRL's release tarball of 2022-04-04, digest measured 2026-09-15 when it was first fetched.
MSIS21_URL = "https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/nrlmsis2.1.tar.gz";
MSIS21_SHA256 = "41e47b29f795d36a5cc252b2858aa2a384c4a7323ace3d48d3ea2f2b37a1a6a8";
EOP_URL = "https://datacenter.iers.org/data/9/finals2000A.all";
SW_URL = "https://celestrak.org/SpaceData/SW-All.csv";
% IGRF-14, IAGA V-MOD. Verified 2026-09-09: byte-identical to the SHA-256 data/igrf/README.md records.
IGRF_URL = "https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt";
% GRACE-FO Level-1B is anonymously available at GFZ's ISDC (owner, 2026-09-09): one archive per day.
L1B_BASE = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04";
L1B_DAY = "2020-01-01";
L1B_MEMBERS = ["GNV1B_2020-01-01_C_04.txt", "SCA1B_2020-01-01_C_04.txt", "MAG1B_2020-01-01_C_04.txt"];
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
msis21 = localFetchArchive("nrlmsis21", fullfile(root, 'mex', 'nrlmsis21'), "msis_calc.F90", MSIS21_URL, MSIS21_SHA256, ...
    options.FetchMsis, "This software incorporates the MSIS(R) empirical atmospheric model software designed and provided by NRL. " + ...
    "Use is governed by the Open Source Academic Research License Agreement contained in the file nrlmsis2.1_license.txt (research, academic and non-profit use only).");
info.sources = [struct('name', "nrlmsise00", 'action', info.msis.action, 'dir', info.msis.dir, 'source', MSIS_URL, 'message', info.msis.message), msis21];
info.data = localFetchData(root, EOP_URL, SW_URL, IGRF_URL, L1B_BASE, L1B_DAY, L1B_MEMBERS, options.FetchData);
info.mex = NonGravForces.buildMex('Check', true);      % 1xT: one element per target
if options.Verbose
    for m = 1:numel(info.mex)
        if info.mex(m).upToDate, state = 'built, current';
        elseif info.mex(m).present, state = 'built, STALE (sources changed): run NonGravForces.buildMex';
        else, state = 'NOT BUILT: run NonGravForces.buildMex (the caller of this gateway refuses until then)';
        end
        [~, gw] = fileparts(info.mex(m).file);
        fprintf('MEX  %-14s %s\n', gw, state);
    end
    for s = 1:numel(info.sources)
        fprintf('SRC  %-14s %s (%s)\n', info.sources(s).name, info.sources(s).action, info.sources(s).dir);
    end
    for d = 1:numel(info.data.items)
        fprintf('DATA %-14s %s\n', info.data.items(d).name, info.data.items(d).action);
    end
end
end

function d = localFetchData(root, eopUrl, swUrl, igrfUrl, l1bBase, l1bDay, l1bMembers, mode)
%LOCALFETCHDATA The data files that can be had without an account, each checked by its content.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
% name, folder, exact file, url, content check (a pattern the file's head must match)
targets = {"eop", fullfile(root, 'data', 'eop'), "finals2000A.all", eopUrl, '^\s*\d{1,2}\s*\d{1,2}\s*\d{1,2}\s+\d{5}\.\d{2}\s+[IP]'; ...
           "spaceweather", fullfile(root, 'data', 'spaceweather'), "SW-All.csv", swUrl, '^DATE,BSRN,ND,KP1'; ...
           "igrf", fullfile(root, 'data', 'igrf'), "igrf14coeffs.txt", igrfUrl, 'IGRF'};
for t = 1:size(targets, 1)
    folder = targets{t, 2};
    dst = fullfile(folder, char(targets{t, 3}));
    present = isfile(dst);
    rec = struct('name', targets{t, 1}, 'action', "present", 'file', string(dst), 'source', targets{t, 4});
    if mode == "never" || (present && mode ~= "force")
        if ~present, rec.action = "missing (FetchData never)"; rec.file = ""; end
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if ~isfolder(folder), mkdir(folder); end
    tmp = [tempname '.download'];
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching %s from %s\n', targets{t, 1}, targets{t, 4});
    ok = false;
    for attempt = 1:2                     % CelesTrak answered 503 once and 200 on the retry
        try
            websave(tmp, char(targets{t, 4}));
            ok = true;
            break
        catch err
            rec.action = "failed: " + string(err.message);
        end
    end
    if ok && ~localHeadMatches(tmp, targets{t, 5})
        ok = false;
        rec.action = "failed: the server's reply is not the expected product (its first lines do not match) -> it was discarded";
    end
    if ok
        movefile(tmp, dst, 'f');
        rec.action = "fetched";
    else
        rec.file = "";
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, targets{t, 4});
    end
    items(end+1) = rec; %#ok<AGROW>
    clear cleanup
end
% GRACE-FO Level-1B: one archive per day, 147 MB, fetched only when a member is absent, members only.
l1b = fullfile(root, 'data', 'l1b');
url = l1bBase + "/" + extractBefore(l1bDay, 5) + "/gracefo_1B_" + l1bDay + "_RL04.ascii.noLRI.tgz";
missing = l1bMembers(~isfile(fullfile(l1b, l1bMembers)));
rec = struct('name', "l1b", 'action', "present", 'file', strjoin(fullfile(l1b, l1bMembers), newline), 'source', url);
if mode == "never" || (isempty(missing) && mode ~= "force")
    if ~isempty(missing), rec.action = "missing (FetchData never)"; rec.file = ""; end
else
    if ~isfolder(l1b), mkdir(l1b); end
    tmp = [tempname '.tgz']; unpack = tempname;
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching GRACE-FO Level-1B %s (about 147 MB) from %s\n', l1bDay, url);
    try
        websave(tmp, char(url));
        mkdir(unpack);
        tarFile = gunzip(tmp, unpack);
        localTarMembers(tarFile{1}, l1bMembers, unpack);
        delete(tarFile{1});
        got = l1bMembers(isfile(fullfile(unpack, l1bMembers)));
        got = got(arrayfun(@(m) localHeadMatches(fullfile(unpack, m), '^header:'), got));
        for m = got, movefile(fullfile(unpack, m), fullfile(l1b, m), 'f'); end
        if numel(got) == numel(l1bMembers)
            rec.action = "fetched";
        else
            rec.action = "incomplete: the archive did not hold readable " + strjoin(setdiff(l1bMembers, got), ", ");
            rec.file = strjoin(fullfile(l1b, got), newline);
        end
    catch err
        rec.action = "failed: " + string(err.message);
        rec.file = "";
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, url);
    end
    if isfolder(unpack), rmdir(unpack, 's'); end
end
items(end+1) = rec;
d = struct('items', items, 'mode', mode);
end

function localTarMembers(tarFile, names, outDir)
%LOCALTARMEMBERS Write the members of a POSIX (ustar or GNU) tar whose base names are in names.
% A tar is a sequence of 512-byte headers, each followed by its data padded to 512 bytes; an
% all-zero header ends it. GNU long names arrive as a preceding 'L' entry. Only regular files
% ('0' or NUL) are written, under their base name, so a folder inside the archive does not matter.
fid = fopen(tarFile, 'r');
if fid < 0
    error('setupNonGravForces:tar', 'cannot open %s', tarFile);
end
closer = onCleanup(@() fclose(fid));
longName = "";
while true
    h = fread(fid, 512, '*uint8')';
    if numel(h) < 512 || all(h == 0), break; end
    name = string(deblank(char(h(1:100))));
    prefix = string(deblank(char(h(346:500))));
    if strlength(prefix) > 0, name = prefix + "/" + name; end
    digits = char(h(125:136));
    digits = digits(digits >= '0' & digits <= '7');       % octal, NUL- or blank-terminated, sometimes blank-padded
    sz = 0;
    if ~isempty(digits), sz = base2dec(digits, 8); end
    type = char(h(157));
    blocks = ceil(sz / 512) * 512;
    if type == 'L'
        data = fread(fid, blocks, '*uint8')';
        longName = string(deblank(char(data(1:sz))));
        continue
    end
    if strlength(longName) > 0, name = longName; longName = ""; end
    [~, base, ext] = fileparts(name);
    if (type == '0' || type == char(0)) && any(base + ext == names)
        data = fread(fid, sz, '*uint8');
        out = fopen(fullfile(outDir, base + ext), 'w');
        fwrite(out, data, 'uint8');
        fclose(out);
        fseek(fid, blocks - sz, 'cof');
    else
        fseek(fid, blocks, 'cof');
    end
end
end

function tf = localHeadMatches(file, pattern)
%LOCALHEADMATCHES The first ten lines of a file contain the pattern: the content check a status code is not.
tf = false;
fid = fopen(file, 'r');
if fid < 0, return; end
closer = onCleanup(@() fclose(fid));
for k = 1:10
    line = fgetl(fid);
    if ~ischar(line), return; end
    if ~isempty(regexp(line, pattern, 'once'))
        tf = true;
        return
    end
end
end

function localQuietDelete(f)
%LOCALQUIETDELETE Remove a temporary file without turning a cleanup into an error.
if isfile(f)
    try, delete(f); catch, end
end
end

function rec = localFetchArchive(name, dst, marker, url, sha256, mode, notice)
%LOCALFETCHARCHIVE A pinned tarball of model source, unpacked into dst after its digest matched.
rec = struct('name', name, 'action', "present", 'dir', string(dst), 'source', url, 'message', "");
here = isfile(fullfile(dst, marker));
if mode == "never" || (here && mode ~= "force")
    if ~here, rec.action = "skipped"; end
    return
end
fprintf('setupNonGravForces: fetching %s from %s\n', name, url);
tmp = [tempname '.tar.gz'];
cleanup = onCleanup(@() localQuietDelete(tmp));
try
    websave(tmp, char(url));
catch err
    rec.action = "failed"; rec.message = "download failed: " + string(err.message) + " -> fetch " + url + " by hand and unpack it into " + dst;
    fprintf('setupNonGravForces: %s\n', rec.message);
    return
end
got = localSha256(tmp);
if got ~= sha256
    rec.action = "failed";
    rec.message = "the archive's SHA-256 is " + got + ", the pin is " + sha256 + " -> the release changed or the download is not the archive; nothing was unpacked";
    fprintf('setupNonGravForces: %s\n', rec.message);
    return
end
if isfolder(dst), rmdir(dst, 's'); end
untar(tmp, dst);
fprintf('setupNonGravForces: %s\n', notice);
rec.action = "fetched";
rec.message = notice;
end

function h = localSha256(file)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file, 'r');
closer = onCleanup(@() fclose(fid));
md.update(fread(fid, Inf, '*uint8'));
h = string(lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
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
