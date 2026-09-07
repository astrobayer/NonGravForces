function result = runNgfTests(options)
%RUNNGFTESTS The one runner: discover, run, report, and write tests/last_run.txt.
%
%   result = runNgfTests
%   result = runNgfTests('ExcludeSlow', false, 'Name', "tNgfQc*")
%
%   GREEN means: every discovered class loaded, nothing failed, nothing
%   incomplete, the category audit found no public unit short of a category
%   (or exempted with a reason), and last_run.txt was written. The pass
%   count alone is not the gate.
%
% NAME-VALUE ARGUMENTS
%   ExcludeSlow  1x1 logical  skip tag "speed" (default true)
%   Name         1x1 string   class-name pattern (default "tNgf*")
%   Strict       1x1 logical  error when not GREEN (default false)
%
% OUTPUT ARGUMENTS
%   result       1x1 struct   .passed .failed .incomplete .duration_s .green
%                             .records (table) .audit (table) .logFile
%
% KEYWORDS test runner gate report
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    options.ExcludeSlow (1,1) logical = true
    options.Name (1,1) string = "tNgf*"
    options.Strict (1,1) logical = false
end
import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag
here = fileparts(mfilename('fullpath'));
NgfTestCase.recordStore('reset');
suite = TestSuite.fromFolder(here, 'Name', options.Name + "/*");
if isempty(suite)
    % Never silently widen to the whole folder: a mistyped pattern would then report a green
    % run of tests the caller did not ask for.
    NonGravForces.fail("runNgfTests", "inputName", "'Name' pattern '%s' selects no test class in %s -> check the pattern (the default is tNgf*)", options.Name, here);
end
if options.ExcludeSlow
    suite = suite.selectIf(~HasTag("speed"));
end
runner = matlab.unittest.TestRunner.withTextOutput('OutputDetail', 2);
t0 = tic;
res = runner.run(suite);
result.duration_s = toc(t0);
result.passed = nnz([res.Passed]); result.failed = nnz([res.Failed]); result.incomplete = nnz([res.Incomplete]);
result.records = NgfTestCase.recordTable();
result.audit = auditTestCategories('Verbose', false);
% A unit that holds an exemption for ONE category is not thereby excused for the others:
% auditTestCategories already credits the exempted categories into `complete`, so the only
% honest filter is `~complete`. Until 2026-09-07 the extra `& ~exempt` term hid three units.
shortfall = result.audit(~result.audit.complete, :);
result.green = result.failed == 0 && result.incomplete == 0 && isempty(shortfall);
% report
result.logFile = fullfile(here, 'last_run.txt');
fid = fopen(result.logFile, 'w');
fprintf(fid, 'NonGravForces %s  run %s  MATLAB %s\n', NonGravForces.toolboxVersion('string'), char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm z')), version);
fprintf(fid, 'tests %d  passed %d  failed %d  incomplete %d  duration %.1f s  GATE %s\n', numel(res), result.passed, result.failed, result.incomplete, result.duration_s, string(result.green).upper());
fprintf(fid, '\n== category audit (units short of a category) ==\n');
if isempty(shortfall), fprintf(fid, 'none\n'); else
    for k = 1:height(shortfall), fprintf(fid, '%s: missing %s\n', shortfall.unit(k), shortfall.missing(k)); end
end
fprintf(fid, '\n== recorded bounds ==\n');
for k = 1:height(result.records)
    fprintf(fid, '%-40s %12.4e <= %12.4e %-6s %s\n', result.records.name(k), result.records.measured(k), result.records.bound(k), result.records.units(k), result.records.note(k));
end
fprintf(fid, '\n== failures ==\n');
for k = find([res.Failed])
    fprintf(fid, '%s\n', res(k).Name);
end
fclose(fid);
fprintf('\nGATE %s  (%d passed, %d failed, %d incomplete, %d audit shortfalls) -> %s\n', ...
    string(result.green).upper(), result.passed, result.failed, result.incomplete, height(shortfall), result.logFile);
if options.Strict && ~result.green
    error('NonGravForces:runNgfTests:gate', 'gate is not GREEN -> read %s', result.logFile);
end
end
