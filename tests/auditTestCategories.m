function T = auditTestCategories(options)
%AUDITTESTCATEGORIES Count, per public unit, which of the seven categories exist.
%
%   T = auditTestCategories
%
%   A unit is a "done" or "partial" row of doc/functionlist.csv. A test method
%   covers a unit when the suite name "<Class>/<method>" CONTAINS the unit name,
%   case-insensitively, and the method carries a category tag. That is the whole
%   rule: there is no "unit:<name>" class tag, and this help claimed one until
%   2026-09-07. Two consequences the reader must know, because the instrument
%   cannot see them: a substring match credits a unit whose name is contained in
%   another (no such pair exists today, checked 2026-09-07), and a method testing
%   one unit credits any other unit whose name appears in its method name.
%   Counting is mechanical (rule T1); exemptions come from
%   tests/testCategoryExemptions.m with one reason each (rule T2). An exemption
%   credits ONLY the categories it names: the caller must filter on `complete`,
%   never on `~complete & ~exempt`, or every unit holding one exemption becomes
%   invisible for all seven categories.
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
