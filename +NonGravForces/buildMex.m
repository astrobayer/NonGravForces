function report = buildMex(options)
%BUILDMEX Compile the MEX gateways in mex/ into their package private folders, or report their state.
%
%   report = NonGravForces.buildMex                          % build every target
%   report = NonGravForces.buildMex('Check', true)            % report only, nothing compiled
%   report = NonGravForces.buildMex('Target', "kernel")       % one target, scalar report
%
%   The two places C is sanctioned (HANDOVER C3), and no others:
%
%     msis    the NRLMSISE-00 C source plus the marshalling gateway
%             mex/ngfMsisGateway.c, built into +env/private/msisGateway, where
%             NonGravForces.env.density can call it and nothing else can.
%
%             THE NRLMSISE-00 SOURCE IS NOT IN THIS REPOSITORY. Its package
%             carries no redistribution grant (THIRD-PARTY-NOTICES.md), so it is
%             obtained from its maintainer instead of vendored:
%
%                 git clone https://git.linta.de/~brodo/nrlmsise-00.git mex/nrlmsise00
%
%             setupNonGravForces does that clone for you the first time it runs,
%             and it owns the URL and the pinned commit. Three places are
%             searched here, in order: the 'MsisSource' option, the
%             NGF_NRLMSISE00 environment variable, then mex/nrlmsise00. If none
%             holds the source this refuses with depMissing and names setup;
%             it does NOT fetch anything itself, because a build step that
%             reaches the network without being asked is a build step nobody
%             can audit -- an installer may, a build may not.
%     kernel  the single-epoch kernel mex/ngfKernel.c and its gateway
%             mex/ngfKernelGateway.c (LP-6b, WP-17b), built into
%             +kernel/private/kernelGateway, where
%             NonGravForces.kernel.singleEpoch can call it and nothing else can.
%
%   Both binaries are platform-specific and are never committed; run this once
%   per machine and after every change to mex/.
%
%   COMPILER. MATLAB's selected C compiler is used. When none is selected but
%   the MinGW-w64 support package is installed, its location (MW_MINGW64_LOC,
%   %ProgramData%\MATLAB\SupportPackages\<release>\3P.instrset\mingw_w64.instrset,
%   the literal C:\ProgramData path when the variable is absent, or the support
%   package root) is placed in MW_MINGW64_LOC for this session and "mex -setup C"
%   is run, because a fresh MATLAB session forgets that variable (measured
%   2026-09-07, HANDOVER O-14). Without any compiler the function refuses with
%   depCompiler.
%
%   STALENESS IS DECIDED BY CONTENT. A SHA-256 digest of the sources is written
%   beside the binary; 'Check' compares it with the sources on disk, so an
%   edited C file makes the binary report stale even though it still loads.
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type     Units  Description
%   ----     ----  ----     -----  -----------
%   Check    1x1   logical  -      true: report only (default false).
%   Target   1x1   string   -      "all" (default), "msis" or "kernel". "all"
%                                  returns one report per target, in that order.
%   MsisSource 1x1 string   -      Directory holding the NRLMSISE-00 C source. Default: the
%                                  NGF_NRLMSISE00 environment variable if set, else mex/nrlmsise00.
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   report   1xT   struct  -      One element per target (T = 1 for a named
%                                 Target, 2 for "all"), each with
%                                 .target (string), .file (string),
%                                 .present (logical),
%                                 .upToDate (logical: present and digest
%                                 matches the sources), .sourceDigest (string,
%                                 SHA-256 of the sources), .binaryDigest
%                                 (string recorded at build time, "" if none),
%                                 .compiler (string, "" when not built now),
%                                 .built (logical: compiled in this call),
%                                 .sources (1xN string),
%                                 .sourcesPresent (logical: every source file
%                                 exists; false is a REPORT under 'Check' and an
%                                 error otherwise),
%                                 .sourceDir (string, where the source was found,
%                                 "" when it was not).
%                                 The digest covers whatever source was actually
%                                 used, so a different NRLMSISE-00 release marks
%                                 the binary stale and forces a rebuild -- which
%                                 is the mechanism that makes tracking upstream
%                                 safe rather than merely convenient.
%
% ERRORS
%   NonGravForces:buildMex:depCompiler   no C compiler selected and no MinGW-w64
%                                        support package found
%   NonGravForces:buildMex:dataMissing   a source file of mex/ is missing
%   MATLAB:mex:*                         compiler errors propagate unwrapped
%
% EXAMPLE
%   r = NonGravForces.buildMex('Check', true); [r.target; r.upToDate]
%
% See also NonGravForces.env.density, NonGravForces.kernel.singleEpoch, setupNonGravForces.
%
% KEYWORDS mex build compiler nrlmsise kernel
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); second target 2026-09-08 by Claude Opus 5 (Anthropic).
arguments
    options.Check (1,1) logical = false
    options.Target (1,1) string {mustBeMember(options.Target, ["all","msis","kernel"])} = "all"
    options.MsisSource (1,1) string = ""
end
unit = "buildMex";
pkgRoot = fileparts(mfilename('fullpath'));                 % +NonGravForces
root = fileparts(pkgRoot);
mexDir = fullfile(root, 'mex');
msisDir = localMsisDir(mexDir, options.MsisSource);
% One row per target: everything that differs between them lives here and
% nothing below branches on the name (rule R1 applied to the build recipe).
targets = struct( ...
    'name',    {"msis", "kernel"}, ...
    'output',  {'msisGateway', 'kernelGateway'}, ...
    'outDir',  {fullfile(pkgRoot, '+env', 'private'), fullfile(pkgRoot, '+kernel', 'private')}, ...
    'include', {char(msisDir), mexDir}, ...
    'sources', {string({fullfile(mexDir, 'ngfMsisGateway.c'), ...
                        char(fullfile(msisDir, 'nrlmsise-00.c')), ...
                        char(fullfile(msisDir, 'nrlmsise-00_data.c')), ...
                        char(fullfile(msisDir, 'nrlmsise-00.h'))}), ...
                string({fullfile(mexDir, 'ngfKernelGateway.c'), ...
                        fullfile(mexDir, 'ngfKernel.c'), ...
                        fullfile(mexDir, 'ngfKernel.h')})});
if options.Target ~= "all"
    targets = targets([targets.name] == options.Target);
end
cc = [];
report = repmat(struct('target', "", 'file', "", 'sources', {strings(1, 0)}, 'sourceDigest', "", ...
                       'built', false, 'compiler', "", 'sourcesPresent', false, 'sourceDir', "", ...
                       'present', false, 'binaryDigest', "", 'upToDate', false), 1, numel(targets));
for t = 1:numel(targets)
    tg = targets(t);
    missing = tg.sources(~isfile(tg.sources));
    rep.sourcesPresent = isempty(missing);
    rep.sourceDir = "";
    if tg.name == "msis", rep.sourceDir = msisDir; end
    if ~rep.sourcesPresent && ~options.Check
        % Under 'Check' a missing source is REPORTED: "report its state" is what Check is for,
        % and a machine that has not fetched NRLMSISE-00 yet must still be able to ask.
        if tg.name == "msis"
            NonGravForces.fail(unit, "depMissing", ...
                "the NRLMSISE-00 C source is not at %s (missing: %s) -> it is not redistributed with this toolbox (THIRD-PARTY-NOTICES.md). Run setupNonGravForces, which fetches it at the pinned commit; or clone it by hand with\n    git clone https://git.linta.de/~brodo/nrlmsise-00.git mex/nrlmsise00\nor point NGF_NRLMSISE00 or 'MsisSource' at an existing copy", ...
                msisDir, strjoin(missing, ", "));
        end
        NonGravForces.fail(unit, "dataMissing", "source file(s) missing: %s -> restore the mex/ folder from the repository", strjoin(missing, ", "));
    end
    rep.target = tg.name;
    rep.file = string(fullfile(tg.outDir, string(tg.output) + "." + mexext));
    rep.sources = tg.sources;
    if rep.sourcesPresent
        rep.sourceDigest = localDigest(tg.sources);
    else
        rep.sourceDigest = "";           % nothing to digest, and "" can match no recorded digest
    end
    digestFile = fullfile(tg.outDir, [tg.output '.digest.txt']);
    rep.built = false;
    rep.compiler = "";
    if options.Check
        report(t) = localState(rep, digestFile);
        continue
    end
    if isempty(cc)
        cc = localCompiler(unit);                            % selected once for the whole call
    end
    rep.compiler = string(cc(1).Name) + " " + string(cc(1).Version);
    if ~isfolder(tg.outDir), mkdir(tg.outDir); end
    csrc = cellstr(tg.sources(endsWith(tg.sources, ".c")));   % headers are dependencies, not mex arguments
    args = [{'-silent', '-O', '-DNDEBUG', ['-I' char(tg.include)], '-outdir', char(tg.outDir), ...
             '-output', tg.output}, csrc];
    mex(args{:});
    writelines(rep.sourceDigest, digestFile);
    rep.built = true;
    report(t) = localState(rep, digestFile);
end
end

function d = localMsisDir(mexDir, override)
%LOCALMSISDIR Where the NRLMSISE-00 source is expected: option, environment, then the default.
% Resolution order is fixed and reported, so a machine with two copies cannot be surprised by
% which one was compiled: buildMex puts the answer in report.sourceDir.
if strlength(override) > 0
    d = override;
    return
end
env = string(getenv('NGF_NRLMSISE00'));
if strlength(env) > 0
    d = env;
    return
end
d = string(fullfile(mexDir, 'nrlmsise00'));
end

function cc = localCompiler(unit)
%LOCALCOMPILER The selected C compiler, selecting MinGW-w64 for this session if none is (O-14).
cc = mex.getCompilerConfigurations('C', 'Selected');
if ~isempty(cc)
    return
end
% Where the MinGW-w64 support package may live, most specific first. Environment variables can be
% absent in a MATLAB started by a service (measured 2026-09-07: ProgramData empty in the MCP session),
% so the standard Windows location is also tried literally.
rel = "R" + version('-release');
bases = [string(getenv('ProgramData')), string(getenv('ALLUSERSPROFILE')), "C:\ProgramData"];
bases = bases(strlength(bases) > 0);
spr = string(matlabshared.supportpkg.getSupportPackageRoot());
spr = spr(strlength(spr) > 0);
candidates = [string(getenv('MW_MINGW64_LOC')), ...
    fullfile(bases, 'MATLAB', 'SupportPackages', rel, '3P.instrset', 'mingw_w64.instrset'), ...
    fullfile(spr, '3P.instrset', 'mingw_w64.instrset')];
candidates = candidates(strlength(candidates) > 0);     % an unset MW_MINGW64_LOC is ""
candidates = unique(candidates, 'stable');
hit = candidates(arrayfun(@(c) isfile(fullfile(c, 'bin', 'gcc.exe')), candidates));
if isempty(hit)
    NonGravForces.fail(unit, "depCompiler", ...
        "no C compiler is selected and no MinGW-w64 support package was found (looked in %s) -> install MATLAB Support for MinGW-w64 C, or run mex -setup C for another compiler", ...
        strjoin(candidates, "; "));
end
setenv('MW_MINGW64_LOC', char(hit(1)));
evalc('mex -setup C');                                   % selects MinGW for this session; output swallowed (no printing from library code)
cc = mex.getCompilerConfigurations('C', 'Selected');
if isempty(cc)
    NonGravForces.fail(unit, "depCompiler", "MinGW-w64 found at %s but mex -setup C did not select it -> run mex -setup C interactively", hit(1));
end
end

function out = localState(rep, digestFile)
%LOCALSTATE Fill the three state fields, in the field order the report array was preallocated with.
% Assignment into an element of a struct array needs the same fields in the same order, so the
% record is rebuilt here rather than extended; sources is wrapped in a cell or struct() would
% return one element per source file.
present = isfile(rep.file);
if isfile(digestFile)
    binaryDigest = string(strtrim(fileread(digestFile)));
else
    binaryDigest = "";
end
out = struct('target', rep.target, 'file', rep.file, 'sources', {rep.sources}, ...
             'sourceDigest', rep.sourceDigest, 'built', rep.built, 'compiler', rep.compiler, ...
             'sourcesPresent', rep.sourcesPresent, 'sourceDir', rep.sourceDir, ...
             'present', present, 'binaryDigest', binaryDigest, ...
             'upToDate', present && rep.sourcesPresent && binaryDigest == rep.sourceDigest);
end

function h = localDigest(files)
md = java.security.MessageDigest.getInstance('SHA-256');
for k = 1:numel(files)
    md.update(uint8(fileread(files(k))));
end
h = string(lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end
