%% Tutorial 13 - Beyond the electrostatic accelerometer: atom interferometers and drag-free control
% The physics core of NonGravForces returns TRUE dynamics: the non-gravitational acceleration in the
% body frame, the torque, the angular rate. What an instrument reports is something else, and the
% |sensor| layer maps the one to the other. GRACE-FO's electrostatic accelerometer is one such map.
% This tutorial takes the modelled acceleration of the reference day through two others: a
% cold-atom interferometer, which measures a PHASE, once per cycle, modulo 2 pi, and reads the
% rotation of the satellite as an acceleration; and a drag-free control loop, whose output is not
% a reading at all but the residual the thrusters FAIL to remove.
%
% *What you will learn*
%
% * the common interface |sensor.Sensor|: |observe|, |partials|, |noiseModel|, |numParameters|;
% * |sensor.AtomInterferometer|: the phase k_eff a T^2, the triangular sensitivity (a phase is an
%   average, not a sample), the dead time between cycles, the 2 pi ambiguity, and the Coriolis
%   phase driven by the ATOM velocity relative to the instrument;
% * |sensor.DragFree|: the loop's high-pass sensitivity, its start-up transient, the thruster
%   noise it injects, and saturation when the drag exceeds the actuators' authority.
%
% *Requirements:* none - the acceleration is the standard Model's along the committed ten minutes
% of 2020-01-01. Run time: about ten seconds.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - The true acceleration on a 10 Hz grid
% The input every sensor takes is a struct |dyn| with the true acceleration |dyn.a| (3xK, body
% frame, m/s^2) and, for the instruments that need them, the angular rate |dyn.omega| and its
% derivative |dyn.omegaDot| (3xK, rad/s and rad/s^2, body frame). Here |a| is the standard Model
% (Tutorial 1) along the committed ten minutes of GRACE-C at 1 Hz, interpolated by pchip onto a
% 10 Hz grid, because the atom interferometer's sensitivity window and the drag-free loop are
% resolved in tenths of a second. The rate is held constant at the Earth-pointing value (0,
% 1.10e-3, 0) rad/s that |geom.attitudeQuality| measures on these days (CLAUDE.md, N-12).
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
[ep, r, v] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("standard"), 'Attitude', att, ...
    'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
idx = 1:600;
e = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
a1 = M.accel(e, r(:, idx), v(:, idx));
t1 = double(ep.sec(idx) - ep.sec(1));
fs = 10;
t10 = (t1(1):1 / fs:t1(end)).';
a = interp1(t1(:), a1.', t10, 'pchip').';
K = numel(t10);
omega = repmat([0; 1.10e-3; 0], 1, K);
dyn = struct('a', a, 'omega', omega, 'omegaDot', zeros(3, K));
fprintf('%d samples at %d Hz over %.0f s; |a| from %.3e to %.3e m/s^2, mean a_x %.3e\n', K, fs, t10(end) - t10(1), ...
    min(vecnorm(a)), max(vecnorm(a)), mean(a(1, :)));

%% Step 2 - One interface, three instruments
% Every instrument subclasses |NonGravForces.sensor.Sensor|, so code that simulates a reading is
% written once against |observe(dyn)|, which returns |.value|, a |.mask| of the samples that are
% valid, the |.frame| of the value and an |.info| struct; |partials(dyn)| gives the derivatives of
% the reading with respect to the instrument's own calibration parameters, |noiseModel(f)| its
% amplitude spectral density. The three built here: an electrostatic accelerometer without its
% CRN filter ([] - Tutorial 8 uses the filtered one), a rubidium atom interferometer with its
% Raman beam along x, and a drag-free loop of 0.1 Hz bandwidth. Extra fields of |dyn| that an
% instrument does not read are ignored.
sensors = {NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), zeros(3, 1), []), ...
    NonGravForces.sensor.AtomInterferometer(1.610e7, 0.5, 2, [1; 0; 0], 'SampleRate', fs), ...
    NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', fs)};
fprintf('%-20s %-20s %6s %-14s %10s %s\n', 'class', 'name', 'params', 'frame', 'size', 'parameter names');
for k = 1:numel(sensors)
    s = sensors{k};
    y = s.observe(dyn);
    [~, pn] = s.partials(dyn);
    fprintf('%-20s %-20s %6d %-14s %10s %s\n', extractAfter(class(s), "sensor."), s.name, s.numParameters(), ...
        y.frame, mat2str(size(y.value)), strjoin(pn(1:min(4, end)), ' '));
end
fprintf('all are a NonGravForces.sensor.Sensor: %d\n', all(cellfun(@(s) isa(s, 'NonGravForces.sensor.Sensor'), sensors)));

%% Step 3 - The atom interferometer: a phase, once per cycle
% A Mach-Zehnder atom interferometer reads dphi = k_eff <n, a> T^2, with k_eff = 1.610e7 rad/m the
% two-photon Raman wave vector on the rubidium D2 line, T = 0.5 s the interrogation time and n the
% beam direction. Three things make it unlike an accelerometer, and the class keeps all three
% explicit: the phase is a TRIANGULAR AVERAGE of a(t) over 2T (|info.weighted|), not a sample;
% the instrument is BLIND between cycles - with a 2 s cycle one reading in twenty samples, the
% dead time |info.deadTime_s| - and |y.mask| is true only there; and the readout gives the phase
% MODULO 2 pi (|info.wrapped|), so the integer fringe number (|info.fringe|) must be resolved by
% other means, usually a classical accelerometer beside it.
ai = sensors{2};
y = ai.observe(dyn);
gain = 1.610e7 * 0.5 ^ 2;
m = y.mask;
fprintf('%d readings in %d samples, dead time %.1f s per cycle\n', nnz(m), K, y.info.deadTime_s);
fprintf('mean phase %.4f rad = %.4e m/s^2; fringe numbers %d .. %d; wrapped phase in [%.3f, %.3f] rad\n', ...
    mean(y.value(m)), mean(y.value(m)) / gain, min(y.info.fringe(m)), max(y.info.fringe(m)), ...
    min(y.info.wrapped(m)), max(y.info.wrapped(m)));
fprintf('weighted vs sampled a_x at the readings: largest difference %.3e m/s^2\n', max(abs(y.info.weighted(m) - a(1, m))));
asd = ai.noiseModel([1e-3 0.1 1]);
fprintf('phase ASD at 1e-3 / 0.1 / 1 Hz: %s rad/rtHz (Inf above half the cycle rate, 0.25 Hz)\n', mat2str(asd, 3));
figure
plot(t10, 1e9 * a(1, :), t10(m), 1e9 * y.value(m) / gain, 'o'), grid on
xlabel('s since 2020-01-01 00:00 GPS'), ylabel('nm/s^2')
legend(["a_x, true" "phase / (k_{eff} T^2), one per cycle"], 'Location', 'best')
title('An atom interferometer samples once per cycle')

%% Step 4 - The Coriolis phase
% A rotating instrument reads the transverse velocity of the ATOMS as an acceleration:
% -2 k_eff T^2 <n, omega x v_at>. v_at is the velocity of the atom cloud RELATIVE TO THE
% INSTRUMENT - its launch and release residual, millimetres per second - and NOT the orbital
% velocity, which atoms and satellite share. With 'Rotation' true the class needs |dyn.omega| and
% |dyn.atomVelocity| (3xK, body frame, m/s); here 1 mm/s along z, perpendicular to both the beam
% (x) and the rotation (y), the worst case. The printed equivalent acceleration is to be read
% against the drag of Step 1: the rotation term is not a small correction.
sr = NonGravForces.sensor.AtomInterferometer(1.610e7, 0.5, 2, [1; 0; 0], 'SampleRate', fs, 'Rotation', true);
dynR = dyn;
dynR.atomVelocity = repmat([0; 0; 1e-3], 1, K);
yr = sr.observe(dynR);
cor = mean(yr.value(yr.mask) - y.value(y.mask));
fprintf('Coriolis phase %.3f rad = %.2f fringes = %.3e m/s^2 equivalent (mean |a_x| %.3e)\n', cor, cor / (2 * pi), ...
    abs(cor) / gain, mean(abs(a(1, :))));
dynR.atomVelocity = repmat([0; 1e-3; 0], 1, K);
yp = sr.observe(dynR);
fprintf('the same velocity ALONG the rotation axis: Coriolis phase %.3g rad\n', mean(yp.value(yp.mask) - y.value(y.mask)));

%% Step 5 - Drag-free control: the residual a loop fails to remove
% A drag-free satellite fires its thrusters against the drag measured by its own sensor. Inside
% the control bandwidth the loop cancels the acceleration; outside it cannot act. |DragFree|
% models what reaches the payload as a cascade of first-order high-pass sections (corner 0.1 Hz
% here, one section), so the residual is the drag's content above the corner. The GRACE-like drag
% of Step 1 varies on the orbital time scale, far below 0.1 Hz, so almost nothing is left - which
% is the point of drag-free flight. The filter starts from rest; the first five time constants
% are masked as a start-up transient, not trimmed. |magnitudeResponse(f)| is |S(f)| of the loop
% and |noiseModel(f)| the thruster noise it injects - inside the bandwidth, where the loop acts.
df = sensors{3};
yd = df.observe(dyn);
ok = yd.mask(1, :);
fprintf('settling: %d samples masked; |S| at f_c = %.4f (one section: 1/sqrt(2) = %.4f)\n', yd.info.settleSamples, ...
    yd.info.gainAtBandwidth, 1 / sqrt(2));
fprintf('RMS of the true a: [%.3e %.3e %.3e], of the residual: [%.3e %.3e %.3e] m/s^2\n', ...
    sqrt(mean(a(:, ok) .^ 2, 2)), sqrt(mean(yd.value(:, ok) .^ 2, 2)));
fq = [1 / 5640, 1e-2, 0.1, 1];
fprintf('|S(f)| at 1/rev, 0.01, 0.1, 1 Hz: %s\n', mat2str(df.magnitudeResponse(fq), 3));
nz = df.noiseModel(fq);
fprintf('thruster noise ASD at the same f (x axis): %s m/s^2/rtHz\n', mat2str(nz(1, :), 3));

%% Step 6 - Saturation
% A thruster set has finite AUTHORITY, the largest acceleration it can command. When the drag
% exceeds it the loop commands everything it has and the excess reaches the payload, a hard
% nonlinearity that is exactly what a geomagnetic storm does to a drag-free mission. Here the
% authority is set to half the largest |a| of Step 1, so the loop saturates over part of the arc;
% |info.saturated| marks those samples.
authority = 0.5 * max(vecnorm(a));
dfSat = NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', fs, 'Authority', authority);
ys = dfSat.observe(dyn);
oks = ys.mask(1, :);
fprintf('authority %.3e m/s^2: %d of %d samples saturated; residual |a| up to %.3e m/s^2\n', authority, ...
    nnz(ys.info.saturated), K, max(vecnorm(ys.value(:, oks))));
figure
plot(t10(ok), 1e9 * vecnorm(yd.value(:, ok)), t10(oks), 1e9 * vecnorm(ys.value(:, oks)), ...
    t10, 1e9 * vecnorm(a), ':'), grid on
xlabel('s since 2020-01-01 00:00 GPS'), ylabel('|a| (nm/s^2)')
legend(["residual, ideal loop" "residual, saturated loop" "true drag + radiation + thermal"], 'Location', 'best')
title('What a drag-free loop leaves behind')

%% Step 7 - Calibration parameters of the other sensors
% Each instrument names its own calibration parameters and returns the partial derivatives of
% its reading with respect to them, the design matrix a calibration of that instrument would
% use. For the atom interferometer: a phase bias and a relative gain error of k_eff T^2; for the
% drag-free loop: a residual bias per axis and a gain error. Both enter linearly, so the partials
% are exact. Printed: the shapes, and the gain column against the phase it multiplies.
[Jai, nai] = ai.partials(dyn);
[Jdf, ndf] = df.partials(dyn);
fprintf('AtomInterferometer: J %s, parameters %s\n', mat2str(size(Jai)), strjoin(nai, ', '));
fprintf('DragFree:           J %s, parameters %s\n', mat2str(size(Jdf)), strjoin(ndf, ', '));
fprintf('AI gain column equals the bias-free phase at the readings: largest difference %.2e rad\n', ...
    max(abs(Jai(m, 2) - y.value(m).')));
