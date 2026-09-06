function idx = l1bColumns(unit, raw, names, required)
%L1BCOLUMNS Column indices of named observation variables, measured from the header (rule R3).
%
%   idx = l1bColumns(unit, raw, names)            % every name required
%   idx = l1bColumns(unit, raw, names, false)     % missing names give 0
%
%   Columns are located by the variable names the file's own YAML header
%   declares (raw.obsNames), never by position: a release that reorders or
%   adds a column cannot silently shift a coordinate.
%
% INPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   unit      1x1   string   -      Calling unit, for the error identifier
%   raw       1x1   struct   -      readL1bMulti aux (field obsNames)
%   names     1xN   string   -      Variable names to locate
%   required  1x1   logical  -      Error when a name is absent (default true)
%
% OUTPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   idx       1xN   double   -      Column index into obs, 0 where absent and not required
%
% ERRORS
%   NonGravForces:<unit>:dataColumn   a required variable is not in the header
%
% KEYWORDS private columns header
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
if nargin < 4, required = true; end
obsNames = string(raw.obsNames(:)');
[found, idx] = ismember(names, obsNames);
if required && ~all(found)
    NonGravForces.fail(unit, "dataColumn", "variable(s) %s not declared in the file header (%s) -> check the product and release", ...
        strjoin(names(~found), ","), strjoin(obsNames, " "));
end
idx(~found) = 0;
end
