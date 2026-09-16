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
% NAME-VALUE ARGUMENTS
%   Verbose   1x1  logical  -  Print a table (default true).
%
% OUTPUT ARGUMENTS
%   Name     Size   Type   Units  Description
%   ----     ----   ----   -----  -----------
%   deps     1x4    struct -      .name .required (logical) .present (logical)
%                                 .version (string, "" if absent) .minimum
%                                 (string) .setup (string: setup function)
%                                 .role (string)
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
% a minimum this toolbox does not actually run on. NOTE (R3): this column is REPORTED,
% never COMPARED - nothing here or in ciGate refuses an older dependency. See O-25.
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
  "gjToolbox",       false, "0.1", "gjToolboxSetup",       "orbit integration; app.DynamicsSource and the rhs MEX", @() string(gjToolbox.versionInfo("string")); ...
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
deps = repmat(struct('name',"",'required',false,'present',false,'version',"",'minimum',"",'setup',"",'role',""), 1, n);
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
end
if options.Verbose
    fprintf('%-18s %-9s %-9s %-10s %s\n','dependency','required','present','version','role');
    for k = 1:n
        fprintf('%-18s %-9s %-9s %-10s %s\n', deps(k).name, string(deps(k).required), ...
            string(deps(k).present), deps(k).version, deps(k).role);
    end
end
end

% localContentsVersion lived here until 2026-09-09. It read a version out of a package's
% Contents.m and lsqvce was its only caller, so dropping that dependency left it dead. Code
% Analyzer caught it on the first run (DEFNU) and it is removed rather than kept "in case":
% every other dependency publishes a version function, which is the better probe anyway.
