%% Tutorial 7 - The force model: presets, budgets, series, torques and partials
% |NonGravForces.Model| is the single point of access of the toolbox: it composes the
% environment (wind, density, radiation sources, occlusion, panel temperatures) and the kernel of
% Tutorial 6, once per call, for every epoch of an arc. This tutorial uses all of its methods on
% the ten committed minutes of GRACE-C, 2020-01-01: |accel| per effect for several presets,
% |budget| and |app.forceBudget| (which effect carries the arc), |series| (long arcs in chunks,
% with the transient wall temperatures carried across them, and the |'T0_K'| option behind
% that), |torque|, and |partials| - and what the partials are for: estimating scale factors.
%
% *What you will learn*
%
% * how the fast and standard presets differ on a real arc, and what one changed option
%   (|thermal.model|) does;
% * the per-effect table of |Model.budget| and the multi-model table of |app.forceBudget|;
% * how |Model.series| chunks an arc and why a transient thermal model must carry its wall
%   temperatures across chunks (|'T0_K'|);
% * the surface torques of |Model.torque| and the opt-in gravity-gradient torque;
% * the parameter partials of |Model.partials| (cdScale, crScale, alpha), checked by finite
%   differences, and a scale estimated with them.
%
% *Requirements:* none - every file is committed in |data/|. The precise preset is described but
% not evaluated: it needs model sources and data the data pool provides. Run time: under a
% minute.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Three models on one arc
% The data are those of Tutorial 1. Three models are built on the same spacecraft:
%
% * |fast|: exponential atmosphere (it needs an |'Exponential'| profile: rho0 at h0 and a scale
%   height H, here 3e-12 kg/m^3 at 450 km and 60 km, as in |demos/demoForceBudget.m|), constant
%   drag coefficient, no Earth radiation, no thermal recoil, no self-shadowing;
% * |standard| as shipped: NRLMSISE-00 times the observed daily density scale (|'DensityScale'|),
%   Sentman, Earth radiation from cap tiles, self-shadowing, TRANSIENT wall temperatures;
% * |standard| with ONE option changed, |thermal.model = "steadyState"|: the walls follow the
%   heating instantly, as the preset did until 0.18.0.
%
% |precise| (NRLMSIS 2.1, HWM14 winds, CERES fluxes, the DE440 Sun) is not evaluated here: its
% model sources are fetched by |setupNonGravForces| and its data live in the data pool.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
[ep, r, vel] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
macro = NonGravForces.geom.graceFO("C");
expo = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);

optF = NonGravForces.Options("fast");
optS = NonGravForces.Options("standard");
optT = optS.set("thermal.model", "steadyState");
Mf = NonGravForces.Model(macro, optF, 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'Exponential', expo);
Ms = NonGravForces.Model(macro, optS, 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
Mt = NonGravForces.Model(macro, optT, 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);

fields = ["drag.density" "drag.densityScale" "drag.gsi" "erp.backend" "geom.selfShadow" "thermal.model"];
optP = NonGravForces.Options("precise");
opts = {optF, optS, optT, optP};
fprintf('%-18s %-12s %-12s %-14s %-12s\n', 'option', 'fast', 'standard', 'std, steady', 'precise');
for k = 1:numel(fields)
    part = split(fields(k), ".");
    val = strings(1, 4);
    for i = 1:4, val(i) = string(opts{i}.(part(1)).(part(2))); end
    fprintf('%-18s %-12s %-12s %-14s %-12s\n', fields(k), val);
end

K = 600; idx = 1:K;
e = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
rr = r(:, idx); vv = vel(:, idx);
t = (ep.sec(idx) - ep.sec(1))';

%% Step 2 - accel, per effect, for the three models
% |[a, info] = M.accel(epoch, r, v)| returns the acceleration (m/s^2, body frame, 3 x K) and, in
% |info.drag|, |info.radiation| and |info.thermal|, the three effects that add up to it exactly.
% Without a quaternion argument the Model asks its attitude provider. The table is the RMS per
% axis (x along track, y cross track, z radial) in nm/s^2.
models = {Mf, Ms, Mt}; labels = ["fast" "standard" "std, steady"];
infos = cell(1, 3); accs = cell(1, 3);
rms3 = @(x) sqrt(mean(x .^ 2, 2))' * 1e9;
fprintf('%-12s %-10s %8s %8s %8s\n', 'model', 'effect', 'x', 'y', 'z');
for i = 1:3
    [accs{i}, infos{i}] = models{i}.accel(e, rr, vv);
    for eff = ["drag" "radiation" "thermal"]
        fprintf('%-12s %-10s %8.3f %8.3f %8.3f\n', labels(i), eff, rms3(infos{i}.(eff)));
    end
end
fprintf('standard: accel minus the sum of its effects, max |.| = %.1e m/s^2\n', ...
    max(abs(accs{2} - (infos{2}.drag + infos{2}.radiation + infos{2}.thermal)), [], 'all'));
figure
subplot(2, 1, 1)
plot(t, 1e9 * [infos{1}.drag(1, :); infos{2}.drag(1, :)]), grid on
ylabel('nm/s^2'), title('drag, along track'), legend("fast", "standard")
subplot(2, 1, 2)
plot(t, 1e9 * [infos{2}.thermal(3, :); infos{3}.thermal(3, :)]), grid on
ylabel('nm/s^2'), title('thermal recoil, radial'), legend("transient walls", "steady state")
xlabel('seconds since 2020-01-01 00:00:00 GPS')

%% Step 3 - Model.budget: which effect carries the arc
% |M.budget(epoch, r, v)| reduces the same per-effect arrays to a table: one row per effect plus
% the total, the RMS and the maximum magnitude per axis and the norm of the RMS (all m/s^2). It is
% evaluated from ONE call to |accel|, so the total row is bitwise the RMS of what |accel| returns.
B = Ms.budget(e, rr, vv);
disp(B(:, {'effect', 'rms_x', 'rms_y', 'rms_z', 'rms_norm'}))
fprintf('total row vs RMS of accel: max difference %.1e m/s^2\n', ...
    max(abs(B{4, {'rms_x', 'rms_y', 'rms_z'}} - sqrt(mean(accs{2} .^ 2, 2))')));

%% Step 4 - app.forceBudget: several models side by side
% |app.forceBudget(models, epoch, r, v, q)| evaluates in chunks of |'Chunk'| epochs and keeps only
% sums of squares and running maxima, so its memory does not grow with the arc - a day at 1 Hz is
% 86400 epochs. With a cell array of models it stacks their tables; |'Labels'| names them. The
% |share| column is each effect's |rms_norm| over the total's; it can sum to more than one when
% effects partly cancel. q = [] takes each model's own attitude provider.
[Bm, bi] = NonGravForces.app.forceBudget({Mf, Ms, Mt}, e, rr, vv, [], 'Chunk', 200, 'Labels', labels);
disp(Bm(:, {'model', 'effect', 'rms_norm', 'share'}))
fprintf('%d chunks of %d epochs\n', bi.chunks, bi.chunkSize);

%% Step 5 - Model.series: a long arc in chunks
% |M.series(epoch, r, v)| returns a timetable (|ax|, |ay|, |az| in m/s^2; with
% |'Quantities', ["accel" "torque"]| also |tx|, |ty|, |tz|) and evaluates in chunks of
% |series.chunk| epochs (8640 by default; 100 here to make the point). The result does not depend on
% the chunking beyond rounding, transient walls included: each chunk after the first starts one
% epoch early, at the previous chunk's last epoch and temperatures, and drops that epoch. Its
% Description records how it was computed.
Mc = NonGravForces.Model(macro, optS.set("series.chunk", 100), 'Attitude', att, 'Eop', eop, ...
    'SpaceWeather', sw, 'DensityScale', S);
Ts = Mc.series(e, rr, vv);
fprintf('%s\n', Ts.Properties.Description);
fprintf('series in chunks of 100 vs one accel call: max difference %.1e m/s^2 (signal %.1e)\n', ...
    max(abs([Ts.ax Ts.ay Ts.az]' - accs{2}), [], 'all'), max(abs(accs{2}), [], 'all'));

%% Step 6 - Carrying the walls yourself: 'T0_K'
% A TRANSIENT thermal model integrates the panel temperatures over the epochs of a call, starting
% from the steady state of the call's FIRST epoch. Cutting an arc into separate calls therefore
% restarts the walls at every cut. |accel(..., 'T0_K', T0)| starts them from given temperatures
% (one per panel, K) instead; |info.temperature_K(:, end)| of one call is the start of the next.
% Here the arc is cut into three calls twice: once naively, and once as |series| does it (the
% next call starts at the previous call's last epoch, with its temperatures, and drops it).
cuts = [1 200; 200 400; 400 600];
sub = @(i) referenceSystem.daySecondToEpoch([ep.day(i) ep.sec(i)], ep.scale);
naive = zeros(3, K); carried = zeros(3, K); T0 = [];
for c = 1:3
    i0 = cuts(c, 1) + (c > 1); i1 = cuts(c, 2);           % 1:200, 201:400, 401:600
    [~, in] = Ms.accel(sub(i0:i1), rr(:, i0:i1), vv(:, i0:i1));
    naive(:, i0:i1) = in.thermal;
    j0 = cuts(c, 1);                                      % 1:200, 200:400, 400:600
    [~, ic] = Ms.accel(sub(j0:i1), rr(:, j0:i1), vv(:, j0:i1), [], 'T0_K', T0);
    th = ic.thermal;
    if c > 1, th = th(:, 2:end); end
    carried(:, i0:i1) = th;
    T0 = ic.temperature_K(:, end);
end
fprintf('thermal recoil vs one call: restarted walls %.2e, carried walls %.2e m/s^2 (max difference)\n', ...
    max(abs(naive - infos{2}.thermal), [], 'all'), max(abs(carried - infos{2}.thermal), [], 'all'));
figure
plot(t, 1e9 * [infos{2}.thermal(3, :); naive(3, :); carried(3, :)]), grid on
xline(t(cuts(2:3, 1)), ':')
ylabel('nm/s^2'), xlabel('seconds since 2020-01-01 00:00:00 GPS')
title('radial thermal recoil, transient walls'), legend("one call", "three calls, restarted", "three calls, T0\_K carried")

%% Step 7 - Torques
% |M.torque(epoch, r, v)| returns the torque about the centre of mass in N m, body frame: the
% surface torques of drag, radiation and thermal recoil (|info.drag|, |info.radiation|,
% |info.thermal|, reductions of the same per-panel forces as |accel|) plus the BODY torques,
% which no preset enables. |torques.gravityGradient = true| adds 3 GM/|r|^5 r x (I r) with the
% handbook inertia, in |info.bodyTorque|; the magnetic torque would also need |'B_body'| and
% |'Dipole_Am2'| at the call. The gravity gradient acts on the attitude, never on the
% accelerometer.
Mg = NonGravForces.Model(macro, optS.set("torques.gravityGradient", true), 'Attitude', att, 'Eop', eop, ...
    'SpaceWeather', sw, 'DensityScale', S);
[tau, tq] = Mg.torque(e, rr, vv);
fprintf('%-18s %10s %10s %10s  (RMS, N m)\n', 'torque', 'x', 'y', 'z');
for nm = ["drag" "radiation" "thermal" "bodyTorque"]
    fprintf('%-18s %10.3e %10.3e %10.3e\n', nm, sqrt(mean(tq.(nm) .^ 2, 2)));
end
fprintf('%-18s %10.3e %10.3e %10.3e\n', 'total', sqrt(mean(tau .^ 2, 2)));

%% Step 8 - Model.partials: derivatives with respect to the model parameters
% |[J, names] = M.partials(epoch, r, v, 'WithRespectTo', ...)| returns d a / d p as 3 x P x K
% (m/s^2 per unit of the parameter, body frame). Three scalar parameters have closed forms:
% |cdScale| and |crScale| scale the whole drag and the whole radiation force, so their partials
% ARE those forces divided by the current scale; |alpha| (only for |drag.accommodation = "fixed"|)
% goes through |kernel.partials|. |'r'| and |'v'| are central differences of |accel|. The
% |'Method', "finiteDiff"| route differentiates the scalar options numerically and exists to
% check the analytic one.
wrt = ["cdScale" "crScale" "alpha"];
[J, names] = Ms.partials(e, rr, vv, 'WithRespectTo', wrt);
Jn = Ms.partials(e, rr, vv, 'WithRespectTo', wrt, 'Method', "finiteDiff");
fprintf('%-9s %14s %12s\n', 'parameter', 'RMS |da/dp|', 'analytic-FD');
for i = 1:3
    Ji = squeeze(J(:, i, :)); Jni = squeeze(Jn(:, i, :));
    fprintf('%-9s %14.4e %12.1e\n', names(i), sqrt(mean(vecnorm(Ji) .^ 2)), max(abs(Ji - Jni), [], 'all') / max(abs(Jni), [], 'all'));
end
fprintf('d a / d cdScale vs info.drag: max difference %.1e m/s^2\n', max(abs(squeeze(J(:, 1, :)) - infos{2}.drag), [], 'all'));

%% Step 9 - What the partials are for: estimating scale factors
% Calibration and drag studies ask "which drag scale does the measurement want?" Here the
% answer is planted: the "measurement" is the standard model with |drag.cdScale| 0.9 and
% |srp.crScale| 1.05. The difference to the nominal model is LINEAR in both scales, so one
% least-squares solve with the two partial columns, stacked over all epochs and axes, recovers
% the planted changes to rounding. On real data the same columns enter the calibration of
% Tutorial 9 beside a bias and an accelerometer scale, and |tools/dragScaleMonth.m| estimates
% k_drag over whole months this way.
Mp = NonGravForces.Model(macro, optS.set("drag.cdScale", 0.9).set("srp.crScale", 1.05), 'Attitude', att, ...
    'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
aMeas = Mp.accel(e, rr, vv);
A = [reshape(J(:, 1, :), [], 1), reshape(J(:, 2, :), [], 1)];
dp = A \ reshape(aMeas - accs{2}, [], 1);
fprintf('recovered cdScale %.12f (planted 0.9), crScale %.12f (planted 1.05)\n', 1 + dp(1), 1 + dp(2));
