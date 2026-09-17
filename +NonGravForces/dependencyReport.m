function deps = dependencyReport(options)
%DEPENDENCYREPORT Presence and version of the four ecosystem toolboxes.
%
%   deps = NonGravForces.dependencyReport
%   deps = NonGravForces.dependencyReport('Verbose', false)
%
%   Measured, never declared (rule R3): each version is obtained by calling
%   the dependency's own version function. Absence WARNS and is reported; it
%   never errors, because a library that refuses to load is not protecting
%   the user (rule R5). Consumers that need a dependency call
%   NonGravForces.fail("<unit>", "depMissing", ...) at the point of use.
%
%   THE MINIMUM IS COMPARED, AND IT WAS NOT UNTIL 2026-09-17. For three weeks this
%   column was printed and never read, which made it a guarantee the toolbox did not
%   give (O-25, rule R4). It is compared over the components the MINIMUM states -- "2.1"
%   asks nothing about a patch, "1.0.1" does -- and the outcome is reported in
%   .belowMinimum. ciGate refuses on it; this unit still only reports. A version that is
%   not a dotted number cannot be ordered and is never called too old: MICE answers
%   CSPICE_N0067, so its .comparable is false and its minimum is documentation.
%
% NAME-VALUE ARGUMENTS
%   Verbose   1x1  logical  -  Print a table (default true).
%
% OUTPUT ARGUMENTS
%   Name     Size   Type   Units  Description
%   ----     ----   ----   -----  -----------
%   deps     1x4    struct -      .name .required (logical) .present (logical)
%                                 .version (string, "" if absent) .minimum
%                                 (string) .comparable (logical: the two are dotted
%                                 numbers and can be ordered) .belowMinimum (logical:
%                                 present, comparable and older than .minimum)
%                                 .setup (string: setup function) .role (string)
%
% ERRORS
%   None. Warning NonGravForces:dependencyReport:depMissing per absent
%   REQUIRED dependency.
%
% EXAMPLE
%   d = NonGravForces.dependencyReport('Verbose', false); [d.present]
%
% See also setupNonGravForces.
%
% KEYWORDS dependencies versions ecosystem
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    options.Verbose (1,1) logical = true
end
spec = { ...
% Minimums are MAJOR.MINOR per astrobayer/ci-actions/VERSIONING.md. gracefoReading is
% the one exception and carries a PATCH: io/private/l1bLoad routes every read through
% readL1bMulti, which raised MATLAB:badsubscript on an all-empty TimeRange until
% gracefoReading 1.0.1 (its issue #6, this repository's O-12). Declaring 1.0 would state
% a minimum this toolbox does not actually run on. This column was REPORTED and never COMPARED
% from 2026-08-27 to 2026-09-17; .belowMinimum is measured here now and ciGate refuses on it (O-25).
% lsqvce was listed here as an optional dependency for app.calibrate until 2026-09-09, when the
% owner ruled it overkill for this estimation problem and WP-19 got the toolbox's own solver. A
% dependency table that names a toolbox nothing calls is a table that will be believed, so the row
% is gone rather than kept as a courtesy. If a future unit needs variance components, it comes back.
  "referenceSystem", true,  "2.1", "setupReferenceSystem", "time, frames, ephemerides, attitude kinematics", @() string(referenceSystem.toolboxVersion().version); ...
  "gracefoReading",  true,  "1.0.1", "setupGraceFOReading",  "GRACE-FO Level-1B reading",                        @() string(gracefoReading.version()); ...
% gjToolbox stays OPTIONAL: the toolbox computes forces without an integrator. It is required by
% exactly one unit, app.DynamicsSource, which subclasses gjToolbox.DynamicsSource and therefore
% cannot even load without it, and by the buildMex target "rhs", which compiles against its
% csrc/gj_native.h. ciGate requires it anyway, because a gate that skips those is not a gate.
% versionInfo("string") returns a SENTENCE -- "Gauss-Jackson Toolbox 0.1.0 (2026-09-06)" -- and this
% column had been printing it as a version since 2026-09-08. Nothing noticed until the minimum became
% comparable on 2026-09-17 and the sentence could not be ranked. .Version is the version.
  "gjToolbox",       false, "0.1", "gjToolboxSetup",       "orbit integration; app.DynamicsSource and the rhs MEX", @() string(gjToolbox.versionInfo().Version); ...
% MICE is NAIF's compiled MATLAB toolkit and is NOT a repository anyone can check out: it ships as
% platform-specific MEX and sites manage it centrally, which is why referenceSystem DETECTS it and
% never adds it to the path. It is listed here because since 2026-09-16 the suite evaluates
% sun.source rs.de440, which cannot work without it, and because the ONE thing worse than a missing
% dependency is one that is present only because a particular machine happens to carry it on a saved
% pathdef. There is no setup function to name; "on the MATLAB path" is the instruction.
% cspice_tkvrsn is a MEX and takes a CHAR: called with a string it throws, which this probe would
% have reported as "mice absent" on a machine that has it. Measured 2026-09-16, when the table said
% false while ciGate had just printed CSPICE_N0067 two lines above.
  "mice",            false, "N0067", "(on the MATLAB path)", "SPICE kernels for sun.source rs.de440", @() string(cspice_tkvrsn('TOOLKIT'))};
n = size(spec,1);
deps = repmat(struct('name',"",'required',false,'present',false,'version',"",'minimum',"", ...
    'comparable',false,'belowMinimum',false,'setup',"",'role',""), 1, n);
for k = 1:n
    deps(k).name = spec{k,1}; deps(k).required = spec{k,2}; deps(k).minimum = spec{k,3};
    deps(k).setup = spec{k,4}; deps(k).role = spec{k,5};
    % exist() does not resolve package functions (finding V-9); probe by calling.
    try
        v = spec{k,6}();
        deps(k).present = true;
        deps(k).version = string(v);
    catch
        deps(k).present = false;
        if deps(k).required
            warning("NonGravForces:dependencyReport:depMissing", ...
                "%s not reachable -> run %s from its own folder first", deps(k).name, deps(k).setup);
        end
    end
    % THE MINIMUM IS COMPARED NOW, not merely printed. It was decorative from 2026-08-27 to
    % 2026-09-17 and read as a guarantee it did not give (O-25, rule R4); ciGate turns the
    % comparison into a verdict, because this unit reports and never refuses (R5).
    % NOT EVERY VERSION IS A NUMBER. MICE answers CSPICE_N0067, which no dotted comparison can
    % order, so `comparable` is reported and a version nobody can order is never called too old.
    if deps(k).present
        [deps(k).comparable, deps(k).belowMinimum] = localBelowMinimum(deps(k).version, deps(k).minimum);
    end
end
if options.Verbose
    fprintf('%-18s %-9s %-9s %-12s %-9s %s\n','dependency','required','present','version','minimum','role');
    for k = 1:n
        minShown = deps(k).minimum;
        if deps(k).present && ~deps(k).comparable
            minShown = minShown + "*";     % not orderable; the star says the column is documentation
        elseif deps(k).belowMinimum
            minShown = minShown + " BELOW";
        end
        fprintf('%-18s %-9s %-9s %-12s %-9s %s\n', deps(k).name, string(deps(k).required), ...
            string(deps(k).present), deps(k).version, minShown, deps(k).role);
    end
end
end

function [comparable, below] = localBelowMinimum(version, minimum)
%LOCALBELOWMINIMUM Order two dotted versions over the components the MINIMUM states.
%   The minimum is the shorter statement and decides the precision of the question: "2.1" asks
%   about major and minor and says nothing about the patch, so 2.1.216 satisfies it; gracefoReading
%   states "1.0.1" and therefore does ask about the patch (O-12). A version that is not a dotted
%   number at all -- MICE's CSPICE_N0067 -- is not comparable, and this says so rather than
%   inventing an order: R3, and the honest answer to a question nobody can ask.
comparable = false;
below = false;
vParts = str2double(split(string(version), "."));
mParts = str2double(split(string(minimum), "."));
if any(isnan(vParts)) || any(isnan(mParts)) || isempty(vParts) || isempty(mParts)
    return                                      % not orderable; comparable stays false
end
comparable = true;
nCompare = numel(mParts);
vParts(end+1:nCompare) = 0;                     % "2" against a minimum of "2.1" is 2.0
for i = 1:nCompare
    if vParts(i) < mParts(i)
        below = true;
        return
    elseif vParts(i) > mParts(i)
        return
    end
end
end

% localContentsVersion lived here until 2026-09-09. It read a version out of a package's
% Contents.m and lsqvce was its only caller, so dropping that dependency left it dead. Code
% Analyzer caught it on the first run (DEFNU) and it is removed rather than kept "in case":
% every other dependency publishes a version function, which is the better probe anyway.
