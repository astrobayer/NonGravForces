function ciGate(options)
%CIGATE The ONE CI entry point for NonGravForces: install dependencies, run the suite, judge.
%
%   ciGate                       % what .github/workflows/ci.yml runs
%   ciGate('DependencyRoot', p)  % folder holding referenceSystem/ and gracefoReading/
%
%   Usable locally as a postflight: matlab -batch ciGate. It is never a bare
%   runNgfTests call: runNgfTests RETURNS a result and errors only with
%   'Strict'; the verdict is made here, in one place, and a red run raises so
%   that matlab -batch exits non-zero.
%
% NAME-VALUE ARGUMENTS
%   DependencyRoot  1x1 string  Folder that contains referenceSystem/ and
%                               gracefoReading/ as siblings. Default: the
%                               parent of this repository, which is the CI
%                               checkout layout and the owner's local layout.
%   ExcludeSlow     1x1 logical Skip the 'speed' tag (default false in CI:
%                               the speed tests are ratios, not host times).
%
% ERRORS
%   NonGravForces:ciGate:depMissing   a required dependency setup was not found
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
deps = ["referenceSystem/setupReferenceSystem.m", "gracefoReading/setupGraceFOReading.m"];
for k = 1:numel(deps)
    f = fullfile(options.DependencyRoot, deps(k));
    if ~isfile(f)
        error("NonGravForces:ciGate:depMissing", "%s not found -> check out the dependency next to this repository", f);
    end
    run(f);
end
cd(here);
info = setupNonGravForces('Verbose', true);
mexReport = NonGravForces.buildMex();          % O-14: a runner without a compiler turns the gate red here, by design
fprintf("ciGate: MEX %s built with %s\n", mexReport.file, mexReport.compiler);
req = info.dependencies([info.dependencies.required]);
if ~all([req.present])
    error("NonGravForces:ciGate:depMissing", "a required dependency did not answer after its setup -> read the table above");
end
r = runNgfTests('ExcludeSlow', options.ExcludeSlow, 'Strict', false);
reasons = strings(0, 1);
if r.failed > 0,     reasons(end+1) = sprintf("%d failed", r.failed); end
if r.incomplete > 0, reasons(end+1) = sprintf("%d incomplete", r.incomplete); end
short = r.audit(~r.audit.complete & ~r.audit.exempt, :);
if ~isempty(short), reasons(end+1) = "category audit shortfall: " + strjoin(short.unit, ", "); end
if ~isfile(r.logFile), reasons(end+1) = "tests/last_run.txt was not written"; end
if r.passed == 0,    reasons(end+1) = "no test passed (suite not discovered?)"; end
fprintf("ciGate: %d passed, %d failed, %d incomplete, %d exempt units, GATE %s\n", ...
    r.passed, r.failed, r.incomplete, nnz(r.audit.exempt), string(isempty(reasons)));
if ~isempty(reasons)
    error("NonGravForces:ciGate:red", "gate is RED -> %s", strjoin(reasons, "; "));
end
end
