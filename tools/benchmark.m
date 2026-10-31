function T = benchmark(options)
%BENCHMARK Time the paths a user waits on, and write the numbers down.
%
%   T = benchmark
%   T = benchmark('Epochs', [1 10 100])
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
%   IT WAS FIRST RUN ON 2026-09-18, nine days after it was written, AND IT DID NOT
%   WORK -- see the comment on the native row below. In those nine days the 88 ms
%   it was built to check was quoted into four more files and one shipped PDF,
%   while the correct value sat in tests/tNgfModel.m. An instrument nobody runs is
%   the same disease as a number nobody can check.
%
%   IT WRITES TWO FILES, and the difference between them matters. The human log
%   tools/benchmark_last.txt is APPENDED and GITIGNORED: it accumulates every run
%   on this machine, which is what a developer reads and what nobody wants in a
%   diff. tools/benchmark_reference.csv is OVERWRITTEN and COMMITTED, and is the
%   authority auditPackage.py check 13 compares every per-epoch figure in the
%   repository against. The authority has to be the committed one: the ignored log
%   does not exist in a fresh clone or on the CI audit job, and a check whose
%   authority is absent everywhere but one laptop protects nothing.
%
%   IT IS NOT RUN BY THE GATE. It is slow, it measures the machine as much as
%   the code, and a gate that fails because a laptop was busy teaches people to
%   ignore gates. Run it when a performance claim needs refreshing, and paste
%   the table into the record with the host beside it.
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type     Units  Description
%   ----     ----  ----     -----  -----------
%   Epochs   1xM   double   -      Batch sizes to time, each >= 1; default [1 10 100]. A VECTOR,
%                                  because the per-epoch cost is mostly a property of the batch
%                                  and barely one of the toolbox: the standard preset reads
%                                  15.9 ms at K=1 and 0.798 at K=100, one call paying about
%                                  15 ms of setup before any per-epoch work. The scalar this
%                                  took until 2026-09-18 is what allowed a single figure to be
%                                  quoted as the cost of the toolbox (HANDOVER O-29)
%   Repeats  1x1   double   -      Timed repetitions, >= 1; default 3. The median is reported,
%                                  because one scheduling hiccup should not become a published
%                                  figure. Five warm-up calls precede them: ONE IS NOT WARM, and
%                                  on 2026-09-18 the difference was a factor of 5.4 on
%                                  kernel.singleEpoch and 7.8 on the fast preset at K=1
%   Write    1x1   logical  -      Default true. Write both output files: tools/benchmark_last.txt
%                                  (appended, gitignored) and tools/benchmark_reference.csv
%                                  (overwritten, committed). Pass false to measure without
%                                  touching the committed authority
%
% OUTPUT ARGUMENTS
%   Name  Size  Type   Units  Description
%   ----  ----  -----  -----  -----------
%   T     Nx5   table  -      .path, .epochs (the batch size timed), .seconds (median wall
%                            time), .perEpoch_ms, .note
%
% ERRORS
%   NonGravForces:benchmark:inputRange   Epochs or Repeats below 1, or an Epochs larger than the
%                                        committed arc has states to time (the message counts them;
%                                        writing the count here would be a second authority)
%   NonGravForces:benchmark:dataMissing  a sample file this needs is not in data/
%
% EXAMPLE
%   T = benchmark('Epochs', [1 10 100]);
%   T(T.path == "Model.accel, standard preset", ["epochs" "perEpoch_ms"])   % the curve
%
% See also NonGravForces.Model, NonGravForces.kernel.singleEpoch, NonGravForces.app.forceBudget,
% NonGravForces.app.DynamicsSource.
%
% KEYWORDS benchmark timing performance per-epoch cost profiling
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    options.Epochs (1,:) double = [1 10 100]
    options.Repeats (1,1) double = 3
    options.Write (1,1) logical = true
end
unit = "benchmark";
if isempty(options.Epochs) || any(options.Epochs < 1) || any(~isfinite(options.Epochs)) || ...
        options.Repeats < 1 || ~isfinite(options.Repeats)
    NonGravForces.fail(unit, "inputRange", "'Epochs' and 'Repeats' must be finite and at least 1 -> got [%s] and %g", num2str(options.Epochs), options.Repeats);
end
% A VECTOR, and the default is the three batch sizes that matter. The single value this took until
% 2026-09-18 is what let a per-epoch cost be stated as a property of the toolbox when it is mostly
% a property of K: the standard preset reads 15.9 ms at K=1 and 0.798 at K=100, because one call
% pays about 15 ms of setup before it does any per-epoch work. One number cannot say that; a curve
% says it by itself, which is why O-29 was able to survive ten days of quotation unchallenged.
nList = unique(round(options.Epochs), 'stable');
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

% The arc is a committed 600-record excerpt, so a batch size larger than it has no states to
% time. Refused HERE, by name, rather than left to ep.day(1:n) -- which raises MATLAB's
% "Index exceeds the number of array elements", a message that names neither the fixture nor the
% limit and comes from no unit the caller can look up (P10, and the mnemonic taxonomy of
% CLAUDE.md section 3). The pre-2026-09-18 code knew about this bound and merely clamped with
% min(n, numel(ep.sec)) in a line that has since been deleted, which is a silent answer to a
% question the caller got wrong.
if any(nList > numel(ep.sec))
    NonGravForces.fail(unit, "inputRange", "'Epochs' asks for up to %d epochs and the committed arc holds %d -> lower 'Epochs', or point the benchmark at a longer GNV1B product", max(nList), numel(ep.sec));
end
att = NonGravForces.geom.AttitudeProvider(epAtt, qq);
macro = NonGravForces.geom.graceFO("C");
expo = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
fastModel = NonGravForces.Model(macro, NonGravForces.Options("fast"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'Exponential', expo);
% The standard preset takes the combined observed density scale (since 0.17.0): Swarm's daily factor
% where it has one, the F30 driver fit elsewhere. Both come from the committed excerpts of 2020-01.
f = @(varargin) fullfile(root, varargin{:});
swarm = string(f('data', 'samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
cls = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
scale = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
stdModel = NonGravForces.Model(macro, NonGravForces.Options("standard"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', scale);
% The compiled single-epoch path needs its inputs composed once, outside the timing: what is
% being measured is the kernel, not env.radiationSources. Built from epoch 1 and independent of
% 'Epochs' -- this path and the native one below take ONE state by construction.
ep1 = referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale);
q1 = att.at(ep1);
o = fastModel.options;
vrel = NonGravForces.env.wind(o.drag.wind, ep1, r(:, 1), v(:, 1), 'Eop', eop);
R = referenceSystem.rotationFromQuaternion(q1);
vrelBody = R(:, :, 1)' * vrel(:, 1);
Ritrs = referenceSystem.rotationItrsToGcrs(ep1, 'Eop', eop);
lla = referenceSystem.cartesianToGeodetic((Ritrs(:, :, 1)' * r(:, 1))', 'Ellipsoid', "WGS84")';
[rho, sp, Tk] = NonGravForces.env.density("exponential", ep1, lla, [], 'Exponential', expo);
src1 = NonGravForces.env.radiationSources(ep1, r(:, 1), q1(:, 1), o);
gas1 = struct('rho', rho, 'species', sp, 'T_K', Tk(2), 'vrel_body', vrelBody, 'shadow', []);

% The integrator's path (LP-C WP-C9): the WHOLE standard-preset chain for one state in C, over
% tables built once in prepare(). It is the third evaluation path (P11) and the only one a
% Gauss-Jackson step actually calls, so a per-epoch table without it is a table about
% post-processing only. It needs gjToolbox and the rhs MEX; both are optional, so its absence is
% REPORTED as a row and never degraded into silence.
ngfSrc = [];
tai = referenceSystem.epochTo(referenceSystem.convertTimeScale( ...
    referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale), "TAI"), "datetime");

% THE SPAN IS CENTRED IN THE ATTITUDE THE FIXTURE ACTUALLY CARRIES, measured from the SCA1B
% record rather than assumed, and it does NOT track 'Epochs': this row times ONE state, so tying
% its span to the arc length was never anything but a way to go out of range.
%
% Until 2026-09-18 the span read [60, 60 + dt]. With 'Margin_s' 120 the attitude stencil then
% reaches back to -60 s while the fixture's attitude starts at 0, so prepare() refused and THIS
% BENCHMARK ERRORED FOR EVERY VALUE OF 'Epochs', including its own default -- from the day it was
% written (2026-09-09) until the day it was first run (2026-09-18). Nine days in which the
% instrument built to regenerate the figure of O-29 could regenerate nothing, and nothing said so,
% because nothing ever ran it.
%
% The refusal is caught rather than pre-empted ON PURPOSE. buildTables snaps its grid outward
% (floor to the step, one node below, at least four nodes) before the margin is applied, so the
% attitude a span needs is wider than span +- Margin_s by an amount only that method knows.
% Re-deriving it here would put the same rule in two places and compare them never, which is the
% defect this benchmark is being repaired for in the first place. Ask, and report the answer.
tAtt = 86400 * (double(epAtt.day) - double(ep.day(1))) + (double(epAtt.sec) - double(ep.sec(1)));
tMid = 0.5 * (tAtt(1) + tAtt(end));
tSpan = [tMid, tMid];
rhsReport = NonGravForces.buildMex('Check', true, 'Target', "rhs");

% The two POINTWISE paths, timed once: they take one state whatever 'Epochs' says, so repeating
% them per batch size would present run-to-run noise as if it were batch scaling.
pointwise = { ...
    "kernel.singleEpoch (compiled)", @() NonGravForces.kernel.singleEpoch(macro, gas1, src1, [], 'Model', o.drag.gsi, 'Cd', o.drag.cd), ...
        "one epoch in C, environment composed outside the timing"};
if isempty(which('gjToolbox.DynamicsSource')) || ~rhsReport.upToDate
    pointwise(end+1, :) = {"app.DynamicsSource.accel (native)", [], ...
        "NOT MEASURED: gjToolbox absent or the rhs MEX not built -> run gjToolboxSetup and NonGravForces.buildMex"};
else
    try
        % the standard preset's transient walls (since 2026-10-31) are integrated along the fixture's own orbit
        ngfSrc = NonGravForces.app.DynamicsSource(stdModel, 'Epoch', tai, 'Margin_s', 120, 'TableStep_s', 60, ...
            'ThermalOrbit', struct('epoch', ep, 'r', r, 'v', v));
        ngfSrc = ngfSrc.prepare(struct('Span', tSpan, 'StencilSpan', tSpan, 'MaxColumns', 1));
        rNat = r(:, 1); vNat = v(:, 1); tNat = tMid;
        pointwise(end+1, :) = {"app.DynamicsSource.accel (native)", @() ngfSrc.accel(tNat, rNat, vNat), ...
            "the whole standard-preset chain for one state in C, over tables built once in prepare()"};
    catch err
        ngfSrc = [];
        pointwise(end+1, :) = {"app.DynamicsSource.accel (native)", [], ...
            sprintf("NOT MEASURED: prepare() refused the %.0f s attitude of this fixture -> %s", tAtt(end) - tAtt(1), err.message)};
    end
end
cleanupSrc = onCleanup(@() localRelease(ngfSrc));

nBatched = 4;                                   % the rows of `batched` below, per batch size
nRow = numel(nList) * nBatched + size(pointwise, 1);
paths = strings(nRow, 1); nEp = zeros(nRow, 1);
secs = zeros(nRow, 1); per = zeros(nRow, 1); notes = strings(nRow, 1);
row = 0;
for n = nList
    epArc = referenceSystem.daySecondToEpoch([ep.day(1:n) ep.sec(1:n)], ep.scale);
    rArc = r(:, 1:n); vArc = v(:, 1:n); qArc = att.at(epArc);
    batched = { ...
        "Model.accel, fast preset",      @() fastModel.accel(epArc, rArc, vArc, qArc),                       "exponential density, constant Cd, conical shadow, no ERP or thermal"; ...
        "Model.accel, standard preset",  @() stdModel.accel(epArc, rArc, vArc, qArc),                        "the figure HANDOVER O-29 quoted; selfShadow dominates the MARGINAL cost, not the K=1 cost"; ...
        "Model.budget, fast preset",     @() fastModel.budget(epArc, rArc, vArc, qArc),                      "one accel plus a reduction"; ...
        "app.forceBudget, fast, 1 chunk",@() NonGravForces.app.forceBudget(fastModel, epArc, rArc, vArc, []), "the chunked reduction over the same accel"};
    assert(size(batched, 1) == nBatched);       % the preallocation above counts these rows
    for c = 1:nBatched
        row = row + 1;
        [secs(row), per(row)] = localTime(batched{c, 2}, n, reps);
        paths(row) = batched{c, 1}; nEp(row) = n; notes(row) = batched{c, 3};
    end
end
for c = 1:size(pointwise, 1)
    row = row + 1;
    [secs(row), per(row)] = localTime(pointwise{c, 2}, 1, reps);
    paths(row) = pointwise{c, 1}; nEp(row) = 1; notes(row) = pointwise{c, 3};
end
T = table(paths, nEp, secs, per, notes, 'VariableNames', {'path', 'epochs', 'seconds', 'perEpoch_ms', 'note'});
T.Properties.VariableUnits = ["", "", "s", "ms", ""];
T.Properties.Description = sprintf('NonGravForces %s benchmark, epochs [%s], median of %d, %s', ...
    NonGravForces.toolboxVersion('string'), num2str(nList), reps, version('-release'));

if options.Write
    % THE COMMITTED AUTHORITY IS WRITTEN FIRST. The other file is the gitignored host log, and
    % writing it first meant an editor holding it open could abort the run before the file the
    % audit actually depends on was refreshed -- the artefact nothing depends on deciding the fate
    % of the one everything does.
    localWriteReference(fullfile(root, 'tools', 'benchmark_reference.csv'), T, reps, unit);
    localAppendLog(fullfile(root, 'tools', 'benchmark_last.txt'), T, reps, unit);   % ignored
end
end

function [sec, perEpoch_ms] = localTime(fcn, n, reps)
%LOCALTIME Median wall time of reps calls, after five warm-ups, and the per-epoch cost it implies.
%   A path this host cannot measure is reported as NaN, never silently omitted.
if isempty(fcn)
    sec = NaN; perEpoch_ms = NaN;
    return
end
% FIVE warm-up calls, not one. One is not warm, and the difference is not small: measured
% 2026-09-18, kernel.singleEpoch read 2.44 ms after a single warm-up against the 0.456 ms in
% tools/benchmark_reference.csv, and app.DynamicsSource.accel 0.72 against 0.544 -- factors of
% 5.4 and 1.3 of pure measurement artefact. (Take the settled figures from that file and not
% from here; this comment is the explanation, the CSV is the authority.) The evidence was
% sitting in three consecutive invocations inside one MATLAB process: the first reported the
% high value and the second and third agreed on the low one, which is what a JIT reaching
% steady state looks like and not what the code costs.
for w = 1:5
    fcn();
end
t = zeros(1, reps);
for k = 1:reps
    t0 = tic; fcn(); t(k) = toc(t0);
end
sec = median(t);
perEpoch_ms = 1e3 * sec / n;
end

function localAppendLog(out, T, reps, unit)
%LOCALAPPENDLOG The human log: append-only, per host, and GITIGNORED on purpose.
%   It accumulates every run on this machine, which is what a developer wants to read and exactly
%   what nobody wants in a diff.
fid = fopen(out, 'a');
if fid < 0
    NonGravForces.fail(unit, "dataMissing", "cannot append to %s -> check the folder is writable", out);
end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '\n%s  NonGravForces %s  MATLAB %s  median of %d\n', ...
    string(datetime('now', 'Format', 'uuuu-MM-dd HH:mm')), ...
    NonGravForces.toolboxVersion('string'), version('-release'), reps);
for c = 1:height(T)
    fprintf(fid, '  %-34s K=%-5d %9.4f s  %8.3f ms/epoch   %s\n', ...
        T.path(c), T.epochs(c), T.seconds(c), T.perEpoch_ms(c), T.note(c));
end
end

function localWriteReference(out, T, reps, unit)
%LOCALWRITEREFERENCE The machine-readable authority: OVERWRITTEN, and COMMITTED.
%   auditPackage.py check 13 compares every per-epoch figure quoted anywhere in the repository
%   against this file, so it has to travel with the source. The append-only log cannot serve:
%   it is gitignored, so on a fresh checkout -- the CI audit job, or any new clone -- it does not
%   exist, and a check whose authority is absent everywhere except the machine that wrote it
%   protects nothing. That is the same shape as the 2026-09-16 finding about a digest check that
%   only ran on download.
%
%   Overwritten rather than appended, so re-running the benchmark shows up as a readable diff:
%   this file is meant to be committed alongside the prose it licenses. A run on a slower host
%   will move the numbers, which is correct -- the prose should then be re-checked, and check 13
%   allows 25 % before it complains.
fid = fopen(out, 'w');
if fid < 0
    NonGravForces.fail(unit, "dataMissing", "cannot write %s -> check the folder is writable", out);
end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '# NonGravForces per-epoch costs. WRITTEN BY tools/benchmark.m -- do not hand-edit.\n');
fprintf(fid, '# The authority for every "<x> ms per epoch" figure quoted in this repository;\n');
fprintf(fid, '# auditPackage.py check 13 compares each of them against the perEpoch_ms column.\n');
fprintf(fid, '# generated %s  toolbox %s  MATLAB %s  median of %d\n', ...
    string(datetime('now', 'Format', 'uuuu-MM-dd HH:mm')), ...
    NonGravForces.toolboxVersion('string'), version('-release'), reps);
fprintf(fid, 'path,epochs,seconds,perEpoch_ms\n');
for c = 1:height(T)
    fprintf(fid, '"%s",%d,%.6f,%.4f\n', T.path(c), T.epochs(c), T.seconds(c), T.perEpoch_ms(c));
end
end

function localRelease(src)
%LOCALRELEASE Free the native context of the integrator source, if one was built.
if ~isempty(src) && isa(src, 'NonGravForces.app.DynamicsSource')
    src.release();
end
end
