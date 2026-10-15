function v = toolboxVersion(format)
%TOOLBOXVERSION Semantic version of NonGravForces - the single authority.
%
%   v = NonGravForces.toolboxVersion            % struct
%   s = NonGravForces.toolboxVersion('string')  % '0.15.0'
%
%   THE NUMERIC FIELDS ARE THE AUTHORITY AND THE STRING IS COMPOSED FROM THEM
%   (rule R1, decision of 2026-09-06): an edit changes major/minor/patch below
%   and nothing else. Every copy of the version elsewhere (Contents.m,
%   README.md, doc/html footers) is compared against this function by
%   auditPackage.py; tests compare, never pin a literal.
%
%   WHAT THE THREE COMPONENTS MEAN (owner's ruling, 2026-09-21). One scheme
%   governs every repository of this organisation: MAJOR is a breaking
%   interface, MINOR is added functionality, PATCH is a fix. See ci-actions
%   VERSIONING.md, which is the authority; this help states the rule, it does
%   not own it.
%
%   THE LOCAL RULE THIS REPLACED, named because a dead rule that is not named
%   stays alive: minor used to be the number of LEVELS DELIVERED (LP-2 gave
%   0.2.0) and patch a HANDOVER amendment without a code change. Under that
%   rule a document edit moved the version and a delivered level moved the
%   minor whether or not a caller gained anything. 0.12.0 is NOT renumbered -
%   the ruling says FROM NOW ON, and at 0.y.z there is no stability promise to
%   revise, which makes this the cheapest possible moment to change meaning.
%   The levels delivered remain a fact worth recording; they are recorded in
%   the HANDOVER, not in the version field.
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
%                                With "string": char '0.15.0'.
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
major = uint16(0); minor = uint16(15); patch = uint16(0);
v = struct('major', major, 'minor', minor, 'patch', patch, ...
           'semantic', string(sprintf('%d.%d.%d', major, minor, patch)), ...
           'date', "2026-10-15", 'concept', uint8(3), ...
           'model', "Claude Opus 5.5 (Anthropic)");
if format == "string"
    v = char(v.semantic);
end
end
