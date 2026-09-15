function T = benchmark(options)
%BENCHMARK Time the paths a user waits on, and write the numbers down.
%
%   T = benchmark
%   T = benchmark('Epochs', 20, 'Write', false)
%
%   IT LIVES IN tools/ AND NOT IN THE PACKAGE, on purpose. It is a development
%   instrument, not part of the toolbox's scientific surface: it is not
%   catalogued, it carries no seven test categories (a metamorphic test of a
%   stopwatch would be theatre), and it is not shipped in the .mltbx. Putting it
%   in +NonGravForces would have bought a catalogue row and six meaningless
%   tests in exchange for nothing.
%
%   WHY THIS EXISTS AND WHAT IT IS NOT. Every unit already carries a `speed`
%   test whose measurement lands in tests/last_run.txt, and this does not repeat
%   any of them: those are RATIOS, chosen because a ratio survives a change of
%   host and an absolute does not (rule T5). What no test produces is the number
%   a user actually asks for -- how long does an epoch cost -- because that
%   number is host-dependent by nature and has no business gating anything.
%
%   It exists because the figure this project quotes most often, the 88 ms per
%   epoch of the standard preset in HANDOVER O-29, was measured once by hand and
%   has been repeated ever since. A number that cannot be regenerated is a number
%   nobody can check. This regenerates it.
%
%   IT IS NOT RUN BY THE GATE. It is slow, it measures the machine as much as
%   the code, and a gate that fails because a laptop was busy teaches people to
%   ignore gates. Run it when a performance claim needs refreshing, and paste
%   the table into the record with the host beside it.
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type     Units  Description
%   ----     ----  ----     -----  -----------
%   Epochs   1x1   double   -      Arc length for the per-epoch timings, >= 1; default 10.
%                                  Ten is deliberate: at 88 ms the standard preset alone is a
%                                  second, and a benchmark nobody waits for is a benchmark
%                                  nobody runs
%   Repeats  1x1   double   -      Timed repetitions, >= 1; default 3. The median is reported,
%                                  because one scheduling hiccup should not become a published
%                                  figure
%   Write    1x1   logical  -      Append the table to tools/benchmark_last.txt (default true)
%
% OUTPUT ARGUMENTS
%   Name  Size  Type   Units  Description
%   ----  ----  -----  -----  -----------
%   T     Nx4   table  -      .path, .seconds (median wall time), .perEpoch_ms, .note
%
% ERRORS
%   NonGravForces:benchmark:inputRange   Epochs or Repeats below 1
%   NonGravForces:benchmark:dataMissing  a sample file this needs is not in data/
%
% EXAMPLE
%   T = benchmark('Epochs', 20);
%   T(T.perEpoch_ms > 10, :)      % what is worth profiling
%
% See also NonGravForces.Model, NonGravForces.kernel.singleEpoch, NonGravForces.app.forceBudget,
% NonGravForces.app.DynamicsSource.
%
% KEYWORDS benchmark timing performance per-epoch cost profiling
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    options.Epochs (1,1) double = 10
    options.Repeats (1,1) double = 3
    options.Write (1,1) logical = true
end
unit = "benchmark";
if options.Epochs < 1 || options.Repeats < 1 || ~isfinite(options.Epochs) || ~isfinite(options.Repeats)
    NonGravForces.fail(unit, "inputRange", "'Epochs' and 'Repeats' must be finite and at least 1 -> got %g and %g", options.Epochs, options.Repeats);
end
n = round(options.Epochs);
reps = round(options.Repeats);

root = fileparts(fileparts(mfilename('fullpath')));   % tools/ -> repository root
need = {fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'), ...
        fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'), ...
        fullfile(root, 'data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), ...
        fullfile(root, 'data', 'samples', 'SCA1B_2020-01-01_C_04.txt')};
absent = need(~cellfun(@isfile, need));
if ~isempty(absent)
    NonGravForces.fail(unit, "dataMissing", "sample data missing: %s -> the benchmark uses the same 2020-01-01 GRACE-C arc the tests use", strjoin(string(absent), ", "));
end

eop = referenceSystem.readEop(need{1});
sw = NonGravForces.io.spaceWeather(need{2});
[ep, r, v] = NonGravForces.io.readGNV1B(need{3}, 'Eop', eop);
[epAtt, qq] = NonGravForces.io.readSCA1B(need{4});
att = NonGravForces.geom.AttitudeProvider(epAtt, qq);
macro = NonGravForces.geom.graceFO("C");
expo = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
fastModel = NonGravForces.Model(macro, NonGravForces.Options("fast"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'Exponential', expo);
stdModel = NonGravForces.Model(macro, NonGravForces.Options("standard"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw);
epArc = referenceSystem.daySecondToEpoch([ep.day(1:n) ep.sec(1:n)], ep.scale);
rArc = r(:, 1:n); vArc = v(:, 1:n); qArc = att.at(epArc);

% The compiled single-epoch path needs its inputs composed once, outside the timing: what is
% being measured is the kernel, not env.radiationSources.
o = fastModel.options;
vrel = NonGravForces.env.wind(o.drag.wind, epArc, rArc, vArc, 'Eop', eop);
R = referenceSystem.rotationFromQuaternion(qArc);
vrelBody = R(:, :, 1)' * vrel(:, 1);
Ritrs = referenceSystem.rotationItrsToGcrs(epArc, 'Eop', eop);
lla = referenceSystem.cartesianToGeodetic((Ritrs(:, :, 1)' * rArc(:, 1))', 'Ellipsoid', "WGS84")';
[rho, sp, Tk] = NonGravForces.env.density("exponential", ...
    referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale), lla, [], 'Exponential', expo);
src1 = NonGravForces.env.radiationSources( ...
    referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale), rArc(:, 1), qArc(:, 1), o);
gas1 = struct('rho', rho, 'species', sp, 'T_K', Tk(2), 'vrel_body', vrelBody, 'shadow', []);

cases = { ...
    "Model.accel, fast preset",      @() fastModel.accel(epArc, rArc, vArc, qArc),                       n, "exponential density, constant Cd, conical shadow, no ERP or thermal"; ...
    "Model.accel, standard preset",  @() stdModel.accel(epArc, rArc, vArc, qArc),                        n, "the figure HANDOVER O-29 quotes; dominated by geom.selfShadow"; ...
    "Model.budget, fast preset",     @() fastModel.budget(epArc, rArc, vArc, qArc),                      n, "one accel plus a reduction"; ...
    "app.forceBudget, fast, 1 chunk",@() NonGravForces.app.forceBudget(fastModel, epArc, rArc, vArc, []), n, "the chunked reduction over the same accel"; ...
    "kernel.singleEpoch (compiled)", @() NonGravForces.kernel.singleEpoch(macro, gas1, src1, [], 'Model', o.drag.gsi, 'Cd', o.drag.cd), 1, "one epoch in C, environment composed outside the timing"};

% The integrator's path (LP-C WP-C9): the WHOLE standard-preset chain for one state in C, over
% tables built once in prepare(). It is the third evaluation path (P11) and the only one a
% Gauss-Jackson step actually calls, so a per-epoch table without it is a table about
% post-processing only. It needs gjToolbox and the rhs MEX; both are optional, so its absence is
% REPORTED as a row and never degraded into silence.
ngfSrc = [];
tai = referenceSystem.epochTo(referenceSystem.convertTimeScale( ...
    referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale), "TAI"), "datetime");
tSpan = [60, double(ep.sec(min(n, numel(ep.sec))) - ep.sec(1)) + 60];
rhsReport = NonGravForces.buildMex('Check', true, 'Target', "rhs");
if isempty(which('gjToolbox.DynamicsSource')) || ~rhsReport.upToDate
    cases(end+1, :) = {"app.DynamicsSource.accel (native)", [], 1, ...
        "NOT MEASURED: gjToolbox absent or the rhs MEX not built -> run gjToolboxSetup and NonGravForces.buildMex"};
else
    ngfSrc = NonGravForces.app.DynamicsSource(stdModel, 'Epoch', tai, 'Margin_s', 120, 'TableStep_s', 60);
    ngfSrc = ngfSrc.prepare(struct('Span', tSpan, 'StencilSpan', tSpan, 'MaxColumns', 1));
    rNat = rArc(:, 1); vNat = vArc(:, 1); tNat = mean(tSpan);
    cases(end+1, :) = {"app.DynamicsSource.accel (native)", @() ngfSrc.accel(tNat, rNat, vNat), 1, ...
        "the whole standard-preset chain for one state in C, over tables built once in prepare()"};
end
cleanupSrc = onCleanup(@() localRelease(ngfSrc));

paths = strings(size(cases, 1), 1);
secs = zeros(size(cases, 1), 1);
per = zeros(size(cases, 1), 1);
notes = strings(size(cases, 1), 1);
for c = 1:size(cases, 1)
    fcn = cases{c, 2};
    if isempty(fcn)                          % a path this host cannot measure: reported, never silent
        paths(c) = cases{c, 1}; secs(c) = NaN; per(c) = NaN; notes(c) = cases{c, 4};
        continue
    end
    fcn();                                   % warm: the first call pays for lookup and JIT
    t = zeros(1, reps);
    for k = 1:reps
        t0 = tic; fcn(); t(k) = toc(t0);
    end
    paths(c) = cases{c, 1};
    secs(c) = median(t);
    per(c) = 1e3 * median(t) / cases{c, 3};
    notes(c) = cases{c, 4};
end
T = table(paths, secs, per, notes, 'VariableNames', {'path', 'seconds', 'perEpoch_ms', 'note'});
T.Properties.VariableUnits = ["", "s", "ms", ""];
T.Properties.Description = sprintf('NonGravForces %s benchmark, %d epochs, median of %d, %s', ...
    NonGravForces.toolboxVersion('string'), n, reps, version('-release'));

if options.Write
    out = fullfile(root, 'tools', 'benchmark_last.txt');
    fid = fopen(out, 'a');
    if fid < 0
        NonGravForces.fail(unit, "dataMissing", "cannot append to %s -> check the folder is writable", out);
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, '\n%s  NonGravForces %s  MATLAB %s  %d epochs, median of %d\n', ...
        string(datetime('now', 'Format', 'uuuu-MM-dd HH:mm')), ...
        NonGravForces.toolboxVersion('string'), version('-release'), n, reps);
    for c = 1:height(T)
        fprintf(fid, '  %-34s %9.4f s  %8.3f ms/epoch   %s\n', T.path(c), T.seconds(c), T.perEpoch_ms(c), T.note(c));
    end
end
end

function localRelease(src)
%LOCALRELEASE Free the native context of the integrator source, if one was built.
if ~isempty(src) && isa(src, 'NonGravForces.app.DynamicsSource')
    src.release();
end
end
