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
  "referenceSystem", true,  "2.1", "setupReferenceSystem", "time, frames, ephemerides, attitude kinematics", @() string(referenceSystem.toolboxVersion().version); ...
  "gracefoReading",  true,  "1.0", "setupGraceFOReading",  "GRACE-FO Level-1B reading",                        @() string(gracefoReading.version()); ...
  "gjToolbox",       false, "0.1", "gjToolboxSetup",       "orbit integration (Gauss-Jackson)",                @() string(gjToolbox.versionInfo("string")); ...
  "lsqvce",          false, "1.0", "lsqvce_setup",         "least squares, VCE, whitening, snooping",          @() localContentsVersion("lsqvce")};
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

function v = localContentsVersion(pkg)
% Read "Version x.y" from the package Contents.m - measured, not declared (R3).
f = which(pkg + ".Contents");
if isempty(f), f = which("+" + pkg + "/Contents.m"); end
if isempty(f), error("NonGravForces:dependencyReport:depMissing", "no Contents.m for %s", pkg); end
tok = regexp(fileread(f), 'Version\s+([0-9][0-9.]*)', 'tokens', 'once');
if isempty(tok), error("NonGravForces:dependencyReport:depMissing", "no Version line in %s", f); end
v = string(tok{1});
end
