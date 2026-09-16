function T = auditTestCategories(options)
%AUDITTESTCATEGORIES Count, per public unit, which of the seven categories exist.
%
%   T = auditTestCategories
%
%   A unit is a "done" or "partial" row of doc/functionlist.csv. A test method
%   covers a unit when that unit's name is the LONGEST catalogued name that the
%   METHOD name begins with (case-insensitively), and the method carries a
%   category tag. The class name plays no part.
%
%   Until 2026-09-15 the rule was "the suite name <Class>/<method> contains the
%   unit name", and this help asserted that no two unit names collide. Four
%   pairs did (HANDOVER LP-C, finding F8): every densityInversion* method credited
%   env.density, every macroModel* method credited Model, every method of the class
%   tNgfQcScreen credited qc.screen, and any *Window* method credited env.wind, so a
%   unit could lose all its own tests and stay "complete". Under the longest-prefix
%   rule a method is credited to exactly one unit. Fault-injected: renaming
%   densityContract drops env.density's contract category.
%
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
% the method part of "<Class>/<method>" (a parameterised method carries "(...)" after it), and its owning unit:
% the longest catalogued name it begins with
method = lower(extractAfter(names, "/"));
method = extractBefore(method + "(", "(");
owner = strings(size(method));
bestLen = zeros(size(method));
lu = lower(units);
for u = 1:n
    starts = startsWith(method, lu(u)) & strlength(lu(u)) > bestLen;
    owner(starts) = units(u);
    bestLen(starts) = strlength(lu(u));
end
have = false(n, numel(cats)); exempt = false(n, 1); missing = strings(n, 1);
for u = 1:n
    hit = owner == units(u);
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
