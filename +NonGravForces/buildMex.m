function report = buildMex(options)
%BUILDMEX Compile the MEX gateways in mex/ into their package private folders, or report their state.
%
%   report = NonGravForces.buildMex                          % build every target
%   report = NonGravForces.buildMex('Check', true)            % report only, nothing compiled
%   report = NonGravForces.buildMex('Target', "kernel")       % one target, scalar report
%
%   The six places compiled code is sanctioned (HANDOVER C3, LP-C WP-C6, WP-C9, WP-C11 and DTM2020, Part I
%   2026-10-27), and no others:
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
%     msis21  the UNMODIFIED NRLMSIS 2.1 Fortran of NRL, a bind(C) shim of this
%             toolbox's own (mex/ngfMsis21Shim.F90) and the C gateway
%             mex/ngfMsis21Gateway.c, built into +env/private/msis21Gateway with the
%             model's parameter file msis21.parm and NRL's licence agreement copied
%             beside it (the agreement travels with every copy, 4c). The NRL package is
%             not redistributed either (Open Source Academic Research License):
%             setupNonGravForces fetches NRL's tarball, checks its SHA-256 and unpacks
%             it into mex/nrlmsis21; 'Msis21Source' and NGF_NRLMSIS21 point elsewhere.
%             The Fortran is compiled with the gfortran of the MinGW-w64 support
%             package (the compiler MATLAB's mex uses for C), in double precision
%             (-cpp -DDBLE, which the package names as the setting that reproduces its
%             reference output exactly), and libgfortran is linked STATICALLY so the
%             binary needs no DLL from the compiler at run time.
%     hwm14   the UNMODIFIED HWM14 Fortran of NRL (Drob et al. 2015), a bind(C) shim
%             (mex/ngfHwm14Shim.F90) and the gateway mex/ngfHwm14Gateway.c, built into
%             +env/private/hwm14Gateway with the model's three data files beside it.
%             setupNonGravForces fetches the files by SHA-256 into mex/hwm14;
%             'Hwm14Source' and NGF_HWM14 point elsewhere. Same compiler and static
%             runtime as msis21; the model itself is single precision.
%     dtm2020 CNES's UNMODIFIED DTM2020 operational subroutines (F10.7, Kp) as SWAMI's MCM
%             repository ships them, a bind(C) shim (mex/ngfDtm2020Shim.F90) and the
%             gateway mex/ngfDtm2020Gateway.c, built into +env/private/dtm2020Gateway with
%             the coefficient file and CNES's licence beside it. Not redistributed (CNES:
%             academic, non-commercial use only, no modification): setupNonGravForces
%             fetches the files by SHA-256 into mex/dtm2020; 'Dtm2020Source' and NGF_DTM2020
%             point elsewhere. Same compiler and static runtime as msis21, compiled in
%             double precision (-fdefault-real-8), as CNES's benchmark program asks.
%     rhs     the WHOLE non-gravitational chain for one state in C (mex/ngfRhs.c),
%             built into +app/private/rhsGateway, where NonGravForces.app.DynamicsSource
%             can reach it and nothing else can. It links the kernel, NRLMSISE-00,
%             NRLMSIS 2.1 and HWM14, and it compiles against gjToolbox's frozen
%             contract csrc/gj_native.h, which is NOT copied here: 'GjToolboxSource',
%             NGF_GJTOOLBOX or the folder of an installed gjToolbox is read (rule R1).
%             Without gjToolbox this target is REPORTED as unbuildable and the other
%             targets still build, because integrating is optional and computing forces
%             is not; asking for it by name refuses with depMissing instead.
%     kernel  the single-epoch kernel mex/ngfKernel.c and its gateway
%             mex/ngfKernelGateway.c (LP-6b, WP-17b), built into
%             +kernel/private/kernelGateway, where
%             NonGravForces.kernel.singleEpoch can call it and nothing else can.
%
%   The binaries are platform-specific and are never committed; run this once
%   per machine and after every change to mex/.
%
%   COMPILER. MATLAB's selected C compiler is used. When none is selected but
%   the MinGW-w64 support package is installed, its location (MW_MINGW64_LOC,
%   %ProgramData%\MATLAB\SupportPackages\<release>\3P.instrset\mingw_w64.instrset,
%   the literal C:\ProgramData path when the variable is absent, or the support
%   package root) is placed in MW_MINGW64_LOC for this session and "mex -setup C"
%   is run, because a fresh MATLAB session forgets that variable (measured
%   2026-09-07, HANDOVER O-14). Without any compiler the function refuses with
%   depCompiler; so do the Fortran targets (msis21, hwm14, dtm2020) when the selected C
%   compiler has no gfortran beside it (a non-MinGW compiler).
%
%   STALENESS IS DECIDED BY CONTENT. A SHA-256 digest of the sources is written
%   beside the binary; 'Check' compares it with the sources on disk, so an
%   edited C file makes the binary report stale even though it still loads.
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type     Units  Description
%   ----     ----  ----     -----  -----------
%   Check    1x1   logical  -      true: report only (default false).
%   Target   1x1   string   -      "all" (default), "msis", "msis21", "hwm14", "dtm2020", "rhs" or "kernel". "all"
%                                  returns one report per target, in that order.
%   MsisSource 1x1 string   -      Directory holding the NRLMSISE-00 C source. Default: the
%                                  NGF_NRLMSISE00 environment variable if set, else mex/nrlmsise00.
%   Msis21Source 1x1 string -      Directory holding the NRLMSIS 2.1 package. Default: the
%                                  NGF_NRLMSIS21 environment variable if set, else mex/nrlmsis21.
%   Hwm14Source 1x1 string  -      Directory holding the HWM14 files. Default: the NGF_HWM14
%                                  environment variable if set, else mex/hwm14.
%   Dtm2020Source 1x1 string -     Directory holding the DTM2020 files. Default: the NGF_DTM2020
%                                  environment variable if set, else mex/dtm2020.
%   GjToolboxSource 1x1 string -   Directory holding gj_native.h (odeSuite's csrc). Default: the
%                                  NGF_GJTOOLBOX environment variable, else the csrc folder of the
%                                  gjToolbox on the path.
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   report   1xT   struct  -      One element per target (T = 1 for a named
%                                 Target, 6 for "all"), each with
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
%                                 .sourceDir (string, where the third-party source
%                                 was looked for; gjToolbox's csrc for rhs; "" for
%                                 the kernel, which has none, and for an unbuildable rhs).
%                                 The digest covers whatever source was actually
%                                 used, so a different NRLMSISE-00 release marks
%                                 the binary stale and forces a rebuild -- which
%                                 is the mechanism that makes tracking upstream
%                                 safe rather than merely convenient.
%
% ERRORS
%   NonGravForces:buildMex:depCompiler   no C compiler selected and no MinGW-w64
%                                        support package found; or no gfortran beside
%                                        the selected compiler (msis21, hwm14, dtm2020), or it
%                                        does not run, or its static runtime is missing, or
%                                        mex -setup C did not select the MinGW-w64 found
%   NonGravForces:buildMex:depMissing    the NRLMSISE-00, NRLMSIS 2.1, HWM14 or DTM2020 source is absent,
%                                        or Target "rhs" without gjToolbox or its gj_native.h
%   NonGravForces:buildMex:dataMissing   a source file of mex/ is missing
%   NonGravForces:buildMex:numericBuild  gfortran refused a file (its output is in the message)
%   MATLAB:mex:*                         compiler errors propagate unwrapped
%
% ACCURACY
%   Not applicable: a build step. Staleness is decided by a SHA-256 over the source
%   bytes (C, Fortran and data files), so any change to them reports the binary stale.
%
% REFERENCES
%   HANDOVER C3, LP-C WP-C6/C9/C11, O-14 (MinGW location), O-34 (NRLMSISE-00 not vendored);
%   THIRD-PARTY-NOTICES.md.
%
% EXAMPLE
%   r = NonGravForces.buildMex('Check', true); [r.target; r.upToDate]
%
% See also NonGravForces.env.density, NonGravForces.kernel.singleEpoch, setupNonGravForces.
%
% KEYWORDS mex build compiler nrlmsise nrlmsis dtm2020 gfortran kernel
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); second target 2026-09-08, msis21 2026-09-15 by Claude Opus 5 (Anthropic); dtm2020 2026-10-27 by Claude (Anthropic).
arguments
    options.Check (1,1) logical = false
    options.Target (1,1) string {mustBeMember(options.Target, ["all","msis","msis21","hwm14","dtm2020","rhs","kernel"])} = "all"
    options.MsisSource (1,1) string = ""
    options.Msis21Source (1,1) string = ""
    options.Hwm14Source (1,1) string = ""
    options.Dtm2020Source (1,1) string = ""
    options.GjToolboxSource (1,1) string = ""
end
unit = "buildMex";
pkgRoot = fileparts(mfilename('fullpath'));                 % +NonGravForces
root = fileparts(pkgRoot);
mexDir = fullfile(root, 'mex');
msisDir = localSourceDir(options.MsisSource, 'NGF_NRLMSISE00', fullfile(mexDir, 'nrlmsise00'));
msis21Dir = localSourceDir(options.Msis21Source, 'NGF_NRLMSIS21', fullfile(mexDir, 'nrlmsis21'));
hwm14Dir = localSourceDir(options.Hwm14Source, 'NGF_HWM14', fullfile(mexDir, 'hwm14'));
dtm2020Dir = localSourceDir(options.Dtm2020Source, 'NGF_DTM2020', fullfile(mexDir, 'dtm2020'));
gjDir = localGjDir(options.GjToolboxSource);
hwm14Fortran = [string(fullfile(hwm14Dir, 'hwm14.f90')), string(fullfile(mexDir, 'ngfHwm14Shim.F90'))];
% NRL's compile order: each file uses the modules of the files before it (readme.txt).
msis21Fortran = fullfile(msis21Dir, ["msis_constants" "msis_utils" "msis_init" "msis_gfn" "msis_tfn" "msis_dfn" "msis_calc"] + ".F90");
% CNES's subroutine file first (the shim calls its external procedures), then the shim.
dtm2020Fortran = [string(fullfile(dtm2020Dir, 'dtm2020_F107_Kp-subr_MCM.f90')), string(fullfile(mexDir, 'ngfDtm2020Shim.F90'))];
% One row per target: everything that differs between them lives here and
% nothing below branches on the name (rule R1 applied to the build recipe).
% fortran: compiled first, in order, with fflags, and linked into the gateway; data: copied beside the binary.
% fflags: -DDBLE selects NRLMSIS 2.1's double-precision kind (rp); DTM2020 has no kind parameter, and
% CNES's benchmark program asks for a double-precision compile, which is -fdefault-real-8.
% Both take part in the digest, the parameter file included: other coefficients are another model.
targets = struct( ...
    'name',    {"msis", "msis21", "hwm14", "dtm2020", "rhs", "kernel"}, ...
    'output',  {'msisGateway', 'msis21Gateway', 'hwm14Gateway', 'dtm2020Gateway', 'rhsGateway', 'kernelGateway'}, ...
    'outDir',  {fullfile(pkgRoot, '+env', 'private'), fullfile(pkgRoot, '+env', 'private'), fullfile(pkgRoot, '+env', 'private'), ...
                fullfile(pkgRoot, '+env', 'private'), fullfile(pkgRoot, '+app', 'private'), fullfile(pkgRoot, '+kernel', 'private')}, ...
    'include', {string(msisDir), string(mexDir), string(mexDir), string(mexDir), [string(mexDir), string(msisDir), gjDir], string(mexDir)}, ...
    'sourceDir', {msisDir, msis21Dir, hwm14Dir, dtm2020Dir, gjDir, ""}, ...
    'fortran', {strings(1, 0), [msis21Fortran, string(fullfile(mexDir, 'ngfMsis21Shim.F90'))], hwm14Fortran, dtm2020Fortran, ...
                [msis21Fortran, string(fullfile(mexDir, 'ngfMsis21Shim.F90')), hwm14Fortran], strings(1, 0)}, ...
    'fflags',  {"", "-DDBLE", "-DDBLE", "-fdefault-real-8", "-DDBLE", ""}, ...
    'data',    {strings(1, 0), fullfile(msis21Dir, ["msis21.parm" "nrlmsis2.1_license..txt"]), ...
                fullfile(hwm14Dir, ["hwm123114.bin" "dwm07b104i.dat" "gd2qd.dat"]), ...
                fullfile(dtm2020Dir, ["DTM_2020_F107_Kp.dat" "LICENSE_MCM_DTM.txt"]), strings(1, 0), strings(1, 0)}, ...
    'sources', {string({fullfile(mexDir, 'ngfMsisGateway.c'), ...
                        char(fullfile(msisDir, 'nrlmsise-00.c')), ...
                        char(fullfile(msisDir, 'nrlmsise-00_data.c')), ...
                        char(fullfile(msisDir, 'nrlmsise-00.h'))}), ...
                string(fullfile(mexDir, 'ngfMsis21Gateway.c')), ...
                string({fullfile(mexDir, 'ngfHwm14Gateway.c'), fullfile(mexDir, 'ngfHwm14Load.c')}), ...
                string(fullfile(mexDir, 'ngfDtm2020Gateway.c')), ...
                string({fullfile(mexDir, 'ngfRhsGateway.c'), fullfile(mexDir, 'ngfRhs.c'), fullfile(mexDir, 'ngfKernel.c'), ...
                        fullfile(mexDir, 'ngfHwm14Load.c'), char(fullfile(msisDir, 'nrlmsise-00.c')), ...
                        char(fullfile(msisDir, 'nrlmsise-00_data.c')), fullfile(mexDir, 'ngfRhs.h'), ...
                        fullfile(mexDir, 'ngfKernel.h'), char(fullfile(gjDir, 'gj_native.h'))}), ...
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
    if tg.name == "rhs" && strlength(gjDir) == 0
        % Integrating is optional and computing forces is not: without gjToolbox's contract header
        % this target is REPORTED, and "all" builds the rest (asking for it by name refuses below).
        if options.Target == "rhs"
            % Say which of the two it is, MEASURED. The message said "gjToolbox is not on the
            % path" unconditionally until 2026-09-16, and on a host where it WAS on the path and
            % only the folder walk was wrong, that sentence sent the reader to the wrong problem.
            if strlength(string(which('gjToolbox.DynamicsSource'))) == 0
                NonGravForces.fail(unit, "depMissing", ...
                    "the rhs target compiles against gjToolbox's csrc/gj_native.h and gjToolbox is not on the path -> run gjToolboxSetup, or point 'GjToolboxSource' or NGF_GJTOOLBOX at odeSuite's csrc folder");
            end
            NonGravForces.fail(unit, "depMissing", ...
                "gjToolbox is on the path (%s) but no gj_native.h was found beside it -> point 'GjToolboxSource' or NGF_GJTOOLBOX at odeSuite's csrc folder", ...
                which('gjToolbox.DynamicsSource'));
        end
        rep.target = tg.name;
        rep.file = string(fullfile(tg.outDir, string(tg.output) + "." + mexext));
        rep.sources = strings(1, 0); rep.sourceDigest = ""; rep.sourcesPresent = false;
        rep.sourceDir = ""; rep.built = false; rep.compiler = "";
        report(t) = localState(rep, fullfile(tg.outDir, [tg.output '.digest.txt']));
        continue
    end
    allFiles = [tg.sources, tg.fortran, tg.data];
    missing = allFiles(~isfile(allFiles));
    rep.sourcesPresent = isempty(missing);
    rep.sourceDir = tg.sourceDir;
    if ~rep.sourcesPresent && ~options.Check
        % Under 'Check' a missing source is REPORTED: "report its state" is what Check is for,
        % and a machine that has not fetched the models yet must still be able to ask.
        if strlength(tg.sourceDir) > 0 && any(startsWith(missing, tg.sourceDir))
            NonGravForces.fail(unit, "depMissing", ...
                "the third-party source of target %s is not at %s (missing: %s) -> it is not redistributed with this toolbox (THIRD-PARTY-NOTICES.md). Run setupNonGravForces, which fetches it at its pin; or point %s at an existing copy", ...
                tg.name, tg.sourceDir, strjoin(missing, ", "), localSourceHint(tg.name));
        end
        NonGravForces.fail(unit, "dataMissing", "source file(s) missing: %s -> restore the mex/ folder from the repository", strjoin(missing, ", "));
    end
    rep.target = tg.name;
    rep.file = string(fullfile(tg.outDir, string(tg.output) + "." + mexext));
    rep.sources = allFiles;
    if rep.sourcesPresent
        rep.sourceDigest = localDigest(allFiles);
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
    extra = {};
    if ~isempty(tg.fortran)
        objDir = tempname;
        mkdir(objDir);
        cleanObj = onCleanup(@() rmdir(objDir, 's'));
        extra = localFortran(unit, cc, tg.name, tg.fortran, tg.fflags, objDir);
    end
    inc = cellfun(@(d) ['-I' char(d)], cellstr(tg.include), 'UniformOutput', false);
    args = [{'-silent', '-O', '-DNDEBUG'}, inc, {'-outdir', char(tg.outDir), '-output', tg.output}, csrc, extra];
    mex(args{:});
    for f = tg.data
        copyfile(f, tg.outDir, 'f');                         % 'f': NRL ships the parameter file read-only
    end
    clear cleanObj
    writelines(rep.sourceDigest, digestFile);
    rep.built = true;
    report(t) = localState(rep, digestFile);
end
end

function d = localSourceDir(override, variable, default)
%LOCALSOURCEDIR Where a third-party source is expected: option, environment, then the default.
% Resolution order is fixed and reported, so a machine with two copies cannot be surprised by
% which one was compiled: buildMex puts the answer in report.sourceDir.
if strlength(override) > 0
    d = override;
    return
end
env = string(getenv(variable));
if strlength(env) > 0
    d = env;
    return
end
d = string(default);
end

function d = localGjDir(override)
%LOCALGJDIR odeSuite's csrc folder, which owns gj_native.h: option, environment, then the
% gjToolbox on the path. "" when there is none.
%   The layout this walks: <odeSuite>/toolbox/+gjToolbox/DynamicsSource.m, so three fileparts
%   reach <odeSuite> and csrc sits directly in it. A stray '..' here sent the search to the
%   PARENT of odeSuite until 2026-09-16, and every run on this host had set NGF_GJTOOLBOX, so
%   the branch CI actually depends on was never once exercised (LP-C). tNgfEnv/buildMexContract
%   now clears the environment variable and resolves through this branch on purpose.
d = "";
if strlength(override) > 0
    d = override;
elseif strlength(string(getenv('NGF_GJTOOLBOX'))) > 0
    d = string(getenv('NGF_GJTOOLBOX'));
else
    w = string(which('gjToolbox.DynamicsSource'));
    if strlength(w) > 0
        d = string(fullfile(fileparts(fileparts(fileparts(w))), 'csrc'));
    end
end
if strlength(d) > 0 && ~isfile(fullfile(d, 'gj_native.h'))
    d = "";
end
end

function h = localSourceHint(name)
%LOCALSOURCEHINT The option and environment variable that relocate a target's third-party source.
if name == "msis21"
    h = "NGF_NRLMSIS21 or 'Msis21Source'";
elseif name == "hwm14"
    h = "NGF_HWM14 or 'Hwm14Source'";
elseif name == "dtm2020"
    h = "NGF_DTM2020 or 'Dtm2020Source'";
else
    h = "NGF_NRLMSISE00 or 'MsisSource'";
end
end

function objs = localFortran(unit, cc, name, files, fflags, objDir)
%LOCALFORTRAN Compile Fortran files in the given order with the gfortran beside the selected C compiler.
% Returns the mex arguments that link them: the objects, then the static runtime archives. The
% compiler's bin folder is put in front of PATH for the call only, because gfortran runs the
% assembler and linker it finds there.
root = string(cc(1).Location);
gf = fullfile(root, 'bin', 'gfortran.exe');
if ~isfile(gf)
    NonGravForces.fail(unit, "depCompiler", "the %s target needs gfortran beside the selected C compiler and there is none at %s -> select the MinGW-w64 C compiler (mex -setup C) or install its support package", name, gf);
end
[st, ver] = system(sprintf('"%s" -dumpversion', gf));
if st ~= 0
    NonGravForces.fail(unit, "depCompiler", "gfortran at %s does not run: %s -> repair the MinGW-w64 support package", gf, strtrim(ver));
end
libDir = fullfile(root, 'lib', 'gcc', 'x86_64-w64-mingw32', strtrim(ver));
archives = fullfile(libDir, ["libgfortran.a" "libquadmath.a"]);
if ~all(isfile(archives))
    NonGravForces.fail(unit, "depCompiler", "the static gfortran runtime is not at %s -> repair the MinGW-w64 support package", libDir);
end
objs = cell(1, numel(files));
for k = 1:numel(files)
    [~, name] = fileparts(files(k));
    o = fullfile(objDir, name + ".o");
    % -O2 and no -ffast-math: the package's reference output is reproducible only without it.
    cmd = sprintf('set "PATH=%s;%%PATH%%" && "%s" -O2 -cpp %s -c "%s" -J "%s" -o "%s" 2>&1', ...
                  fullfile(root, 'bin'), gf, fflags, files(k), objDir, o);
    [st, out] = system(cmd);
    if st ~= 0
        NonGravForces.fail(unit, "numericBuild", "gfortran refused %s:\n%s -> the source is not the pinned release of target %s, or the compiler is broken", files(k), strtrim(out), name);
    end
    objs{k} = char(o);
end
objs = [objs, cellstr(archives)];
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
%LOCALDIGEST SHA-256 over the BYTES of the files, in order (the parameter file is binary).
md = java.security.MessageDigest.getInstance('SHA-256');
for k = 1:numel(files)
    fid = fopen(files(k), 'r');
    bytes = fread(fid, Inf, '*uint8');
    fclose(fid);
    md.update(bytes);
end
h = string(lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end
