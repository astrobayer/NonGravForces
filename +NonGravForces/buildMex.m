function report = buildMex(options)
%BUILDMEX Compile the NRLMSISE-00 MEX gateway (mex/) into +NonGravForces/+env/private, or report its state.
%
%   report = NonGravForces.buildMex                 % build (or rebuild) now
%   report = NonGravForces.buildMex('Check', true)   % report only, nothing compiled
%
%   The only sanctioned C in this toolbox besides LP-6b (HANDOVER C3): the
%   vendored, unmodified NRLMSISE-00 source in mex/nrlmsise00 plus the
%   hand-written marshalling gateway mex/ngfMsisGateway.c. The binary lands
%   in the env package's private folder, where NonGravForces.env.density can
%   call it and nothing else can. It is platform-specific and is never
%   committed; run this once per machine and after every change to mex/.
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
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   report   1x1   struct  -      .file (string), .present (logical),
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
%   r = NonGravForces.buildMex('Check', true); r.upToDate
%
% See also NonGravForces.env.density, setupNonGravForces.
%
% KEYWORDS mex build compiler nrlmsise
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    options.Check (1,1) logical = false
end
unit = "buildMex";
pkgRoot = fileparts(mfilename('fullpath'));                 % +NonGravForces
root = fileparts(pkgRoot);
mexDir = fullfile(root, 'mex');
outDir = fullfile(pkgRoot, '+env', 'private');
sources = string({fullfile(mexDir, 'ngfMsisGateway.c'), ...
                  fullfile(mexDir, 'nrlmsise00', 'nrlmsise-00.c'), ...
                  fullfile(mexDir, 'nrlmsise00', 'nrlmsise-00_data.c'), ...
                  fullfile(mexDir, 'nrlmsise00', 'nrlmsise-00.h')});
missing = sources(~isfile(sources));
if ~isempty(missing)
    NonGravForces.fail(unit, "dataMissing", "source file(s) missing: %s -> restore the mex/ folder from the repository", strjoin(missing, ", "));
end
report.file = string(fullfile(outDir, "msisGateway." + mexext));
report.sources = sources;
report.sourceDigest = localDigest(sources);
digestFile = fullfile(outDir, 'msisGateway.digest.txt');
report.built = false;
report.compiler = "";
if options.Check
    report = localState(report, digestFile);
    return
end
% --- compiler
cc = mex.getCompilerConfigurations('C', 'Selected');
if isempty(cc)
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
report.compiler = string(cc(1).Name) + " " + string(cc(1).Version);
if ~isfolder(outDir), mkdir(outDir); end
mex('-silent', '-O', '-DNDEBUG', ['-I' fullfile(mexDir, 'nrlmsise00')], '-outdir', outDir, '-output', 'msisGateway', ...
    char(sources(1)), char(sources(2)), char(sources(3)));
writelines(report.sourceDigest, digestFile);
report.built = true;
report = localState(report, digestFile);
end

function report = localState(report, digestFile)
report.present = isfile(report.file);
if isfile(digestFile)
    report.binaryDigest = string(strtrim(fileread(digestFile)));
else
    report.binaryDigest = "";
end
report.upToDate = report.present && report.binaryDigest == report.sourceDigest;
end

function h = localDigest(files)
md = java.security.MessageDigest.getInstance('SHA-256');
for k = 1:numel(files)
    md.update(uint8(fileread(files(k))));
end
h = string(lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end
