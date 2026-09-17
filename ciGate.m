function ciGate(options)
%CIGATE The ONE CI entry point for NonGravForces: install dependencies, run the suite, judge.
%
%   ciGate                       % what .github/workflows/ci.yml runs
%   ciGate('DependencyRoot', p)  % folder holding referenceSystem/, gracefoReading/ and odeSuite/
%
%   Usable locally as a postflight: matlab -batch ciGate. It is never a bare
%   runNgfTests call: runNgfTests RETURNS a result and errors only with
%   'Strict'; the verdict is made here, in one place, and a red run raises so
%   that matlab -batch exits non-zero.
%
% NAME-VALUE ARGUMENTS
%   DependencyRoot  1x1 string  Folder that contains referenceSystem/,
%                               gracefoReading/ and odeSuite/ as siblings.
%                               Default: the parent of this repository, which
%                               is the CI checkout layout and the owner's
%                               local layout.
%   ExcludeSlow     1x1 logical Skip the 'speed' tag (default false in CI:
%                               the speed tests are ratios, not host times).
%
% ERRORS
%   NonGravForces:ciGate:depMissing   a required dependency setup was not found, MICE is absent,
%                                     gjToolbox's stepper did not build, or a MEX target did not
%   NonGravForces:ciGate:dataMissing  a data fetch failed its SHA-256 check
%   NonGravForces:ciGate:red          the gate is not GREEN; message lists why
%   NonGravForces:buildMex:*          no compiler / missing sources (O-14)
%
% KEYWORDS ci gate
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    options.DependencyRoot (1,1) string = string(fileparts(fileparts(mfilename('fullpath'))))
    options.ExcludeSlow (1,1) logical = false
end
here = fileparts(mfilename('fullpath'));
% odeSuite is OPTIONAL for the toolbox (NonGravForces computes forces without an integrator)
% and REQUIRED for the gate: since LP-C WP-C9 the app layer serves gjToolbox's frozen interface
% and the rhs MEX compiles against its csrc/gj_native.h, so a gate without it would leave
% app.DynamicsSource and the whole native evaluator unrun and still report GREEN.
deps = ["referenceSystem/setupReferenceSystem.m", "gracefoReading/setupGraceFOReading.m", ...
        "odeSuite/gjToolboxSetup.m"];
for k = 1:numel(deps)
    f = fullfile(options.DependencyRoot, deps(k));
    if ~isfile(f)
        error("NonGravForces:ciGate:depMissing", "%s not found -> check out the dependency next to this repository", f);
    end
    run(f);
end
cd(here);
% MICE, for the same reason and with the same rule as gjToolbox: optional for the toolbox,
% required by the gate. tNgfEnv and tNgfIntegration evaluate sun.source rs.de440, which is the
% precise preset's OWN DEFAULT, and it cannot be evaluated without NAIF's toolkit. It is reached
% here only because it sits on this machine's MATLAB path, so the gate names it out loud instead
% of depending on it silently -- the same trap that cost two green gates on 2026-09-16.
if exist('cspice_tkvrsn', 'file') == 0
    error("NonGravForces:ciGate:depMissing", ...
        "NAIF MICE is not on the MATLAB path -> the suite evaluates sun.source rs.de440, the precise preset's default, and cannot without it; install MICE or add it to the path");
end
fprintf("ciGate: MICE %s\n", string(cspice_tkvrsn('TOOLKIT')));
% gjToolbox's NATIVE STEPPING KERNEL is an optional build of that toolbox, and this gate asserts
% through it: tNgfIntegration propagates one arc on Binding "native" and requires it to agree with
% the vectorized binding. A fresh CI checkout carries no gj_stepper (the binary is not committed,
% rightly), so the gate builds it -- with the DEPENDENCY'S OWN build task, never a recipe copied
% here (R1), which also smoke-tests the gateway it just built. The owner's host usually has it
% already and this is then a no-op.
if ~gjToolbox.nativeStepper("available")
    odeRoot = char(fullfile(options.DependencyRoot, "odeSuite"));
    fprintf("ciGate: building gjToolbox's native stepping kernel in %s\n", odeRoot);
    back = pwd;
    try
        cd(odeRoot);
        buildtool("mex");
        cd(back);
    catch err
        cd(back);
        rethrow(err);
    end
    if ~gjToolbox.nativeStepper("available")
        error("NonGravForces:ciGate:depMissing", ...
            "gjToolbox's native stepping kernel is still not available after buildtool mex in %s -> the native binding cannot be gated without it", odeRoot);
    end
end
% The model sources are fetched when missing (buildMex cannot build without them); the data are NOT:
% every test reads the committed excerpts, and a 147 MB archive per CI run would test the network.
info = setupNonGravForces('Verbose', true, 'FetchData', "never");
mexReport = NonGravForces.buildMex();          % O-14: a runner without a compiler turns the gate red here, by design
for m = 1:numel(mexReport)                     % one element per target (msis, msis21, hwm14, rhs, kernel)
    fprintf("ciGate: MEX %s (%s) present=%d, built with %s\n", mexReport(m).file, mexReport(m).target, ...
        mexReport(m).upToDate, mexReport(m).compiler);
end
% buildMex REPORTS an unbuildable rhs instead of raising, so that a user without gjToolbox still
% gets the other four. The gate is the place where that report becomes a verdict.
notBuilt = [mexReport(~[mexReport.upToDate]).target];
if ~isempty(notBuilt)
    error("NonGravForces:ciGate:depMissing", "the MEX target(s) %s did not build -> read the table above", strjoin(notBuilt, ", "));
end
% A FETCH THAT FAILED ITS CONTENT CHECK IS A RED GATE. setupNonGravForces reports per item and
% nothing here read that report until 2026-09-16, when a kernel whose SHA-256 did not match was
% printed as "failed" three lines above "GATE TRUE". A "missing (... never)" item is not a failure:
% the gate deliberately does not fetch observations.
failed = info.data.items(startsWith([info.data.items.action], "failed"));
if ~isempty(failed)
    error("NonGravForces:ciGate:dataMissing", "%d data fetch(es) failed their content check -> %s", ...
        numel(failed), strjoin(compose("%s: %s", [failed.name]', [failed.action]'), "; "));
end
req = info.dependencies([info.dependencies.required]);
if ~all([req.present])
    error("NonGravForces:ciGate:depMissing", "a required dependency did not answer after its setup -> read the table above");
end
% AND THE MINIMUM IS A VERDICT HERE. dependencyReport measures .belowMinimum and never
% refuses (R5); the gate is where a report becomes a judgement, as it is for buildMex and for
% the fetches above. The column was printed and read by nothing from 2026-08-27 to
% 2026-09-17, which is exactly as much protection as not having it (O-25).
old = info.dependencies([info.dependencies.belowMinimum]);
if ~isempty(old)
    error("NonGravForces:ciGate:depMissing", "%s -> install at least the declared minimum", ...
        strjoin(compose("%s %s is below the declared minimum %s", [old.name], [old.version], [old.minimum]), "; "));
end
r = runNgfTests('ExcludeSlow', options.ExcludeSlow, 'Strict', false);
reasons = strings(0, 1);
if r.failed > 0,     reasons(end+1) = sprintf("%d failed", r.failed); end
if r.incomplete > 0, reasons(end+1) = sprintf("%d incomplete", r.incomplete); end
short = r.audit(~r.audit.complete, :);          % an exemption for one category never excuses another
if ~isempty(short), reasons(end+1) = "category audit shortfall: " + strjoin(short.unit, ", "); end
if ~isfile(r.logFile), reasons(end+1) = "tests/last_run.txt was not written"; end
if r.passed == 0,    reasons(end+1) = "no test passed (suite not discovered?)"; end
fprintf("ciGate: %d passed, %d failed, %d incomplete, %d exempt units, GATE %s\n", ...
    r.passed, r.failed, r.incomplete, nnz(r.audit.exempt), string(isempty(reasons)));
if ~isempty(reasons)
    error("NonGravForces:ciGate:red", "gate is RED -> %s", strjoin(reasons, "; "));
end
end
