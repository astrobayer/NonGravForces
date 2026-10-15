function p = provenance(unit, options)
%PROVENANCE Metadata record that travels with every NonGravForces result (rule R2).
%
%   p = NonGravForces.provenance("Model.accel")
%   p = NonGravForces.provenance("io.readACC1B", 'Files', files, 'Options', opts)
%
%   A result without provenance cannot be audited. The record answers: which
%   toolbox version, which dependency versions, which option hash, which data
%   files, which warnings, when. Two results may be combined only if their
%   provenance is compatible (see NonGravForces.app.forceBudget for an example).
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   unit      1x1   string  -      Producing unit, qualified inside the package
%
% NAME-VALUE ARGUMENTS
%   Files     1xN   string  -      Data files consumed (default empty)
%   Options   1x1   struct  -      Resolved options; hashed (default empty)
%   Warnings  1xN   string  -      Warning identifiers issued (default empty)
%   Note      1x1   string  -      Free text (default "")
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   p         1x1   struct  -      .unit .toolboxVersion .dependencies (name,
%                                  version) .optionHash (string, 8 hex) .files
%                                  .warnings .note .created (datetime UTC)
%
% ERRORS
%   None.
%
% EXAMPLE
%   p = NonGravForces.provenance("env.eclipseFactor");
%
% See also NonGravForces.dependencyReport, NonGravForces.Options.
%
% KEYWORDS provenance metadata audit
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    unit (1,1) string
    options.Files (1,:) string = string.empty(1,0)
    options.Options struct = struct()
    options.Warnings (1,:) string = string.empty(1,0)
    options.Note (1,1) string = ""
end
% The dependency versions are read once per MATLAB session and cached: dependencyReport
% shells out to each toolbox and is far too slow to run per provenance record. A session
% that swaps a dependency underneath itself must `clear functions` to see the new version.
persistent depCache
if isempty(depCache)
    d = NonGravForces.dependencyReport('Verbose', false);
    depCache = struct('name', {d.name}, 'version', {d.version});
end
p.unit           = unit;
p.toolboxVersion = string(NonGravForces.toolboxVersion('string'));
p.dependencies   = depCache;
p.optionHash     = localHash(options.Options);
p.files          = options.Files;
p.warnings       = options.Warnings;
p.note           = options.Note;
p.created        = datetime('now', 'TimeZone', 'UTC');
end

function h = localHash(s)
% Deterministic 32-bit FNV-1a over the canonical text of the struct.
if isempty(fieldnames(s)), h = "00000000"; return; end
txt = uint8(char(jsonencode(s)));
x = uint32(2166136261);
for b = txt
    x = bitxor(x, uint32(b));
    x = uint32(mod(double(x) * 16777619, 2^32));
end
h = string(lower(dec2hex(x, 8)));
end
