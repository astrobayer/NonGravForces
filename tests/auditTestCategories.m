function T = auditTestCategories(options)
%AUDITTESTCATEGORIES Count, per public unit, which of the seven categories exist.
%
%   T = auditTestCategories
%
%   A unit is a "done" or "partial" row of doc/functionlist.csv. A test method
%   covers a unit when its class carries the tag "unit:<name>" or its method
%   name contains the unit name (case-insensitive) and the method carries a
%   category tag. Counting is mechanical (rule T1); exemptions come from
%   tests/testCategoryExemptions.m with one reason each (rule T2).
%
% OUTPUT ARGUMENTS
%   T   table   unit, categories (1x7 logical), complete, exempt, missing (string)
%
% KEYWORDS audit test categories
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    options.Verbose (1,1) logical = true
end
cats = ["contract","reference","precision","vectorization","robustness","metamorphic","speed"];
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
fl = readtable(fullfile(root, 'doc', 'functionlist.csv'), 'TextType', 'string');
units = fl.name(fl.status ~= "stub" & fl.status ~= "interface");
suite = matlab.unittest.TestSuite.fromFolder(here);
names = string({suite.Name}); tags = {suite.Tags};
ex = testCategoryExemptions();
n = numel(units);
have = false(n, numel(cats)); exempt = false(n, 1); missing = strings(n, 1);
for u = 1:n
    hit = contains(lower(names), lower(units(u)));
    for c = 1:numel(cats)
        have(u, c) = any(hit & cellfun(@(t) any(string(t) == cats(c)), tags));
    end
    if isfield(ex, units(u))
        exCats = string(fieldnames(ex.(units(u))));
        have(u, ismember(cats, exCats)) = true;
        exempt(u) = true;
    end
    missing(u) = strjoin(cats(~have(u, :)), ",");
end
T = table(units, have, all(have, 2), exempt, missing, 'VariableNames', {'unit','categories','complete','exempt','missing'});
if options.Verbose
    disp(T(:, {'unit','complete','exempt','missing'}));
end
end
