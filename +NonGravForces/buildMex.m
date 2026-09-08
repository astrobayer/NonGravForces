function report = buildMex(options)
%BUILDMEX Compile the MEX gateways in mex/ into their package private folders, or report their state.
%
%   report = NonGravForces.buildMex                          % build every target
%   report = NonGravForces.buildMex('Check', true)            % report only, nothing compiled
%   report = NonGravForces.buildMex('Target', "kernel")       % one target, scalar report
%
%   The two places C is sanctioned (HANDOVER C3), and no others:
%
%     msis    the vendored, unmodified NRLMSISE-00 source in mex/nrlmsise00
%             plus the marshalling gateway mex/ngfMsisGateway.c, built into
%             +env/private/msisGateway, where NonGravForces.env.density can
%             call it and nothing else can.
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
%                                 .sources (1xN string).
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
end
unit = "buildMex";
pkgRoot = fileparts(mfilename('fullpath'));                 % +NonGravForces
root = fileparts(pkgRoot);
mexDir = fullfile(root, 'mex');
% One row per target: everything that differs between them lives here and
% nothing below branches on the name (rule R1 applied to the build recipe).
targets = struct( ...
    'name',    {"msis", "kernel"}, ...
    'output',  {'msisGateway', 'kernelGateway'}, ...
    'outDir',  {fullfile(pkgRoot, '+env', 'private'), fullfile(pkgRoot, '+kernel', 'private')}, ...
    'include', {fullfile(mexDir, 'nrlmsise00'), mexDir}, ...
    'sources', {string({fullfile(mexDir, 'ngfMsisGateway.c'), ...
                        fullfile(mexDir, 'nrlmsise00', 'nrlmsise-00.c'), ...
                        fullfile(mexDir, 'nrlmsise00', 'nrlmsise-00_data.c'), ...
                        fullfile(mexDir, 'nrlmsise00', 'nrlmsise-00.h')}), ...
                string({fullfile(mexDir, 'ngfKernelGateway.c'), ...
                        fullfile(mexDir, 'ngfKernel.c'), ...
                        fullfile(mexDir, 'ngfKernel.h')})});
if options.Target ~= "all"
    targets = targets([targets.name] == options.Target);
end
cc = [];
report = repmat(struct('target', "", 'file', "", 'sources', {strings(1, 0)}, 'sourceDigest', "", ...
                       'built', false, 'compiler', "", 'present', false, 'binaryDigest', "", ...
                       'upToDate', false), 1, numel(targets));
for t = 1:numel(targets)
    tg = targets(t);
    missing = tg.sources(~isfile(tg.sources));
    if ~isempty(missing)
        NonGravForces.fail(unit, "dataMissing", "source file(s) missing: %s -> restore the mex/ folder from the repository", strjoin(missing, ", "));
    end
    rep.target = tg.name;
    rep.file = string(fullfile(tg.outDir, string(tg.output) + "." + mexext));
    rep.sources = tg.sources;
    rep.sourceDigest = localDigest(tg.sources);
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
             'present', present, 'binaryDigest', binaryDigest, ...
             'upToDate', present && binaryDigest == rep.sourceDigest);
end

function h = localDigest(files)
md = java.security.MessageDigest.getInstance('SHA-256');
for k = 1:numel(files)
    md.update(uint8(fileread(files(k))));
end
h = string(lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end
