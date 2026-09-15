function v = toolboxVersion(format)
%TOOLBOXVERSION Semantic version of NonGravForces - the single authority.
%
%   v = NonGravForces.toolboxVersion            % struct
%   s = NonGravForces.toolboxVersion('string')  % '0.10.0'
%
%   THE NUMERIC FIELDS ARE THE AUTHORITY AND THE STRING IS COMPOSED FROM THEM
%   (rule R1, decision of 2026-09-06): an edit changes major/minor/patch below
%   and nothing else. Version rule (P-7 as re-based): minor = number of
%   levels delivered (LP-2 -> 0.2.0), patch = HANDOVER amendment without a code
%   change. Every copy of the version elsewhere (Contents.m, README.md,
%   doc/html footers) is compared against this function by auditPackage.py;
%   tests compare, never pin a literal.
%
% INPUT ARGUMENTS
%   Name     Size  Type   Units  Description
%   ----     ----  ----   -----  -----------
%   format   1x1   string -      "struct" (default) or "string".
%
% OUTPUT ARGUMENTS
%   Name     Size  Type   Units  Description
%   ----     ----  ----   -----  -----------
%   v        1x1   struct -      .major .minor .patch (uint16), .semantic
%                                (string, composed), .date (string, ISO),
%                                .concept (uint8, 3 = this redesign), .model
%                                (string: authoring model of this revision).
%                                With "string": char '0.10.0'.
%
% ERRORS
%   MATLAB:validators:mustBeMember   Unknown format keyword.
%
% EXAMPLE
%   NonGravForces.toolboxVersion('string')
%
% See also NonGravForces.dependencyReport.
%
% KEYWORDS version semver
% Created 2026-08-27 by Claude Fable 5 (Anthropic); re-based 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    format (1,1) string {mustBeMember(format, ["struct","string"])} = "struct"
end
major = uint16(0); minor = uint16(10); patch = uint16(0);
v = struct('major', major, 'minor', minor, 'patch', patch, ...
           'semantic', string(sprintf('%d.%d.%d', major, minor, patch)), ...
           'date', "2026-09-16", 'concept', uint8(3), ...
           'model', "Claude Opus 5 (Anthropic)");
if format == "string"
    v = char(v.semantic);
end
end
