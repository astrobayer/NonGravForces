%% Tutorial 8 - From true dynamics to what the accelerometer reads
% |Model.accel| returns the TRUE non-gravitational acceleration at the centre of mass. An
% accelerometer does not read that. Its proof mass sits a fraction of a millimetre away from the
% centre of mass of a rotating spacecraft, it has a bias and a scale of its own, and the Level-1B
% product ACT1B is the 10 Hz instrument stream passed through a low-pass filter (the CRN filter of
% the GRACE L1B processing) and sampled at 1 Hz. |sensor.Electrostatic| is the adapter between the
% two (principle P8: sensors are adapters). This tutorial builds it step by step, puts the model
% through it, compares the result with ACT1B on the committed ten minutes of GRACE-C, and looks at
% what the thruster firings in THR1B do to the measurement.
%
% *What you will learn*
%
% * the sensor equation y = bias + S R (a + omegaDot x dr + omega x (omega x dr)) and how
%   |sensor.Electrostatic| takes bias, the 3x3 scale S and the proof-mass offset dr;
% * the sign convention of omega, and how to get it from the attitude;
% * the CRN filter (|timeSeriesPreparation.crnKernel("ACC")|): 10 Hz grid, edge samples, DC gain;
% * how to put a 1 Hz model through the filter and sample it like ACT1B (|'OutputRate', 1|);
% * what THR1B says about thruster firings and how they show in ACT1B;
% * why a bias and a scale must be estimated - which Tutorial 9 does properly.
%
% *Requirements:* none - every file is committed in |data/|; timeSeriesPreparation provides the CRN
% kernel. Run time: about twenty seconds.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - The true dynamics: the model along ten minutes
% As in Tutorial 1: the standard preset with the observed density scale, evaluated at the 600
% GNV1B epochs (1 Hz). |a| is the non-gravitational acceleration at the centre of mass, body frame
% (the science reference frame, SRF), m/s^2. ACT1B and THR1B are read beside it; THR1B is the
% WHOLE day, ACT1B the first ten minutes.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
[ep, r, vel] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
[epa, act] = NonGravForces.io.readACC1B(f('samples', 'ACT1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
macro = NonGravForces.geom.graceFO("C");
M = NonGravForces.Model(macro, NonGravForces.Options("standard"), 'Attitude', att, 'Eop', eop, ...
    'SpaceWeather', sw, 'DensityScale', S);
K = 600; idx = 1:K;
e = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
a = M.accel(e, r(:, idx), vel(:, idx));
t1 = (ep.sec(idx) - ep.sec(1))';                       % 1xK, seconds since the first epoch
fprintf('model: %d epochs over %.0f s, |a| RMS %.3e m/s^2\n', K, t1(end), sqrt(mean(vecnorm(a) .^ 2)));
fprintf('ACT1B: %d records; mean per axis [%.3e %.3e %.3e] m/s^2 (the instrument bias dominates)\n', ...
    numel(epa.sec), mean(act, 2));

%% Step 2 - The sensor equation: bias, scale, frame
% |sensor.Electrostatic(bias, S, dr, kernel)| models
%
%   y = bias + S R a_sensed,   a_sensed = a + omegaDot x dr + omega x (omega x dr)
%
% with |bias| (3x1, m/s^2, accelerometer frame), |S| (3x3, dimensionless; off-diagonal entries are
% axis couplings, not a rotation), |dr| (3x1, m, the proof-mass offset from the centre of mass in
% the body frame) and R the body-to-accelerometer rotation (|'Rotation'|, default identity: for
% GRACE-FO ACT1B is given in the SRF). |kernel| [] means no filter; then any sample rate is
% accepted (|'SampleRate'| 1 here). |observe(dyn)| takes a struct |dyn| with |a|, |omega| and
% |omegaDot| (3 x K each, body frame) and returns |y.value| (3 x K) and a validity |y.mask|. With
% dr = 0 and no rotation the reading is exactly bias + S a.
b = [1e-6; -3e-6; 2e-7]; Sc = diag([1.05 0.96 1.01]);
s0 = NonGravForces.sensor.Electrostatic(b, Sc, zeros(3, 1), [], 'SampleRate', 1);
dyn1 = struct('a', a, 'omega', zeros(3, K), 'omegaDot', zeros(3, K));
y0 = s0.observe(dyn1);
fprintf('%s sensor, %d parameters, frame "%s", filter "%s"\n', s0.name, s0.numParameters(), y0.frame, y0.info.kernel);
fprintf('max |y - (bias + S a)| = %.1e m/s^2\n', max(abs(y0.value - (b + Sc * a)), [], 'all'));

%% Step 3 - The rotation: omega in the sensor's convention, and the proof-mass offset
% The spacecraft rotates once per revolution to keep pointing at its partner, and a proof mass off
% the centre of mass feels the centrifugal and the Euler acceleration of that rotation. The sign of
% omega is part of the sensor's contract: for R mapping BODY -> INERTIAL, skew(omega) = R' Rdot
% (the |Electrostatic| help explains why only that sign makes the equation kinematics). It is
% formed here from the SCA1B attitude at the model epochs by central differences, with the
% symmetric part of R' Rdot - the differencing error - projected out, as |tests/tNgfSensor.m|
% does. |referenceSystem.angularVelocityFromRotationRate(R, Rdot, 'Resolution', "reference")|
% returns the NEGATIVE of this omega (its rotation runs the other way), which is checked too. The
% offset is the handbook's centre of gravity -> accelerometer vector, |macro.proofMassOffset|
% (JPL D-56935 Table 4).
qk = att.at(e);
R = referenceSystem.rotationFromQuaternion(qk);       % 3x3xK, body -> GCRS
Rdot = zeros(3, 3, K);
Rdot(:, :, 2:K-1) = (R(:, :, 3:K) - R(:, :, 1:K-2)) / 2;
Rdot(:, :, 1) = Rdot(:, :, 2); Rdot(:, :, K) = Rdot(:, :, K-1);
Mw = pagemtimes(pagetranspose(R), Rdot);
Om = (Mw - pagetranspose(Mw)) / 2;
w = [squeeze(Om(3, 2, :)), squeeze(Om(1, 3, :)), squeeze(Om(2, 1, :))]';
wd = zeros(3, K);
wd(:, 2:K-1) = (w(:, 3:K) - w(:, 1:K-2)) / 2;
wd(:, 1) = wd(:, 2); wd(:, K) = wd(:, K-1);
wRef = referenceSystem.angularVelocityFromRotationRate(R, Rdot, 'Resolution', "reference");
fprintf('mean omega [%.3e %.3e %.3e] rad/s (once per revolution of ~94 min: %.3e rad/s)\n', mean(w, 2), 2 * pi / (94 * 60));
fprintf('max |omega + referenceSystem''s| / max |omega| = %.1e\n', max(abs(w + wRef), [], 'all') / max(abs(w), [], 'all'));
dr = macro.proofMassOffset;
sC = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), dr, [], 'SampleRate', 1);
yC = sC.observe(struct('a', zeros(3, K), 'omega', w, 'omegaDot', wd));
fprintf('proof-mass offset [%.2f %.2f %.2f] mm -> coupling RMS [%.2e %.2e %.2e] m/s^2\n', ...
    1e3 * dr, sqrt(mean(yC.value(:, 2:end-1) .^ 2, 2)));

%% Step 4 - The CRN filter
% ACT1B is not a sample of the acceleration; it is the 10 Hz instrument stream convolved with the
% CRN low-pass of the GRACE L1B ATBD (Wu et al. 2006, App. B). Its one authority is
% |timeSeriesPreparation.crnKernel("ACC")|: 1407 symmetric taps on the 10 Hz grid (140.7 s), a
% pass band to about 35 mHz, normalised to unit gain at f0 = 0.37 mHz with its DC gain of
% 1 + 1.1e-7 KEPT (never renormalise it: every calibrated scale would move by it). The response is
% evaluated here from the taps directly. Because the taps span 140.7 s, (1407 - 1)/2 = 703
% samples - 70.3 s - at each end of an arc are invalid.
kc = timeSeriesPreparation.crnKernel("ACC");
fprintf('CRN "%s": %d taps at %g Hz, bandwidth %.3f Hz, f0 %.2e Hz, DC gain 1 + %.3e\n', ...
    kc.preset, numel(kc.F), kc.fs, kc.B, kc.f0, kc.dcGain - 1);
fr = logspace(-4, log10(0.2), 200);
H = abs(exp(-2i * pi * fr(:) * (kc.lags(:)' / kc.fs)) * kc.F(:));
figure
semilogx(fr, 20 * log10(H)), grid on
xlabel('frequency (Hz)'), ylabel('gain (dB)'), title('CRN filter of ACT1B (timeSeriesPreparation.crnKernel("ACC"))')
xline(kc.B, ':', 'B'), ylim([-120 5])

%% Step 5 - The model through the sensor, sampled like ACT1B
% The CRN default needs 10 Hz dynamics, and |observe| checks that from |dyn.epoch| (here seconds,
% 1 x N): a 1 Hz arc convolved with 0.1 s taps would be a filter of ten times the span. The model
% is therefore interpolated from its 1 Hz epochs onto the 10 Hz grid by pchip, as |app.calibrate|
% does (its help measures the interpolation error after the filter at 9e-13 m/s^2), and so are
% omega and omegaDot. |'OutputRate', 1| then SAMPLES every tenth filtered value - L1B forms ACT1B
% the same way; it never averages. dr is left at zero for the comparison with ACT1B; Step 3 showed
% its size.
t10 = (0:10 * t1(end)) / 10;                            % 10 Hz grid over the same span
a10 = interp1(t1, a', t10, 'pchip')';
w10 = interp1(t1, w', t10, 'pchip')';
wd10 = interp1(t1, wd', t10, 'pchip')';
sF = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), zeros(3, 1), "crn", 'OutputRate', 1);
dyn10 = struct('a', a10, 'omega', w10, 'omegaDot', wd10, 'epoch', t10);
yF = sF.observe(dyn10);
ok = yF.mask(1, :);
tY = round(yF.epoch);
fprintf('%d samples at 10 Hz -> %d at 1 Hz, %d valid (%.1f s lost at each end)\n', numel(t10), numel(tY), nnz(ok), (numel(kc.F) - 1) / 2 / kc.fs);
[~, iM] = ismember(tY, t1);
dF = yF.value(:, ok) - a(:, iM(ok));
fprintf('filtered minus unfiltered model, RMS per axis [%.2e %.2e %.2e] m/s^2\n', sqrt(mean(dF .^ 2, 2)));

%% Step 6 - The model against ACT1B: one bias and one scale per axis
% ACT1B is matched to the filtered model by epoch. Its bias is unknown and large (Step 1), and its
% scale is not exactly 1, so per axis the simplest calibration is fitted by least squares:
% ACT1B = b + s y_model on the valid samples. Over ten minutes this is only illustrative - the
% along-track signal hardly varies, so its scale is poorly determined - but the cross-track and
% radial shapes already agree. The proper estimate (one scale over many days, one bias per day,
% uncertainties by bootstrap) is |app.calibrate|, Tutorial 9; the month baseline of 0.19.0 found
% S_xx 0.999 and S_zz 1.092 for GRACE-C (CLAUDE.md section 5).
tA = round(epa.sec - ep.sec(1))';
[inA, iA] = ismember(tY, tA);
use = ok & inA;
yA = act(:, iA(use)); yM = yF.value(:, use); tU = tY(use);
axn = ["x (along track)" "y (cross track)" "z (radial)"];
fit = zeros(2, 3); res = zeros(size(yA));
for j = 1:3
    X = [ones(nnz(use), 1), yM(j, :)'];
    fit(:, j) = X \ yA(j, :)';
    res(j, :) = yA(j, :) - (X * fit(:, j))';
    fprintf('%-16s bias %10.3e m/s^2, scale %7.3f, residual RMS %.2f nm/s^2\n', axn(j), fit(1, j), fit(2, j), ...
        1e9 * sqrt(mean(res(j, :) .^ 2)));
end

%% Step 7 - Thruster firings: THR1B
% The attitude is held by cold-gas thrusters. |io.readTHR1B| returns an event table (one row per
% record with a firing; timeSeriesPreparation's |qc.eventTable| schema) and, in its second output,
% every record: |epoch|, the thruster names, |onTime_s| (thrusters x records, seconds) and
% |fired|. A firing is a short impulse that the accelerometer sees, and through the CRN filter it
% is smeared over +-70 s. The test |tNgfModel/modelAct1bReference| therefore cuts the window
% [t_fire - 0.5 - 70.3, t_fire + on-time + 7 + 70.3] s around every firing; the same is done here,
% and the fit of Step 6 is repeated on what is left.
[thrEvents, thr] = NonGravForces.io.readTHR1B(f('samples', 'THR1B_2020-01-01_C_04.txt'));
kf = any(thr.fired, 1);
tFire = reshape(thr.epoch.sec(kf), 1, []) - ep.sec(1);
onT = max(thr.onTime_s(:, kf), [], 1);
inArc = tFire >= 0 & tFire <= t1(end);
fprintf('THR1B: %d records, %d with a firing over the day (%d thrusters); total on-time %.1f s\n', ...
    thr.numRecords, height(thrEvents), numel(thr.thruster), sum(thr.onTime_s(:)));
fprintf('firings inside the ACT1B excerpt at t = %s s, on-time %s s\n', mat2str(tFire(inArc), 5), mat2str(onT(inArc), 3));
halfCrn = (numel(kc.F) - 1) / 2 / kc.fs;
cut = any(tU(:) >= tFire - 0.5 - halfCrn & tU(:) <= tFire + onT + 7 + halfCrn, 2)';
fprintf('%d of %d compared samples lie within a firing window\n', nnz(cut), numel(tU));
if nnz(~cut) > 30
    for j = 1:3
        X = [ones(nnz(~cut), 1), yM(j, ~cut)'];
        c2 = X \ yA(j, ~cut)';
        r2 = yA(j, ~cut) - (X * c2)';
        fprintf('%-16s without the firing windows: residual RMS %.2f nm/s^2 (all samples %.2f)\n', axn(j), ...
            1e9 * sqrt(mean(r2 .^ 2)), 1e9 * sqrt(mean(res(j, :) .^ 2)));
    end
end
figure
for j = 2:3
    subplot(2, 1, j - 1)
    plot(tU, 1e9 * res(j, :)), grid on, hold on
    xline(tFire(inArc), 'r:')
    ylabel('nm/s^2'), title(axn(j) + ": ACT1B minus fitted model (red: thruster firings)")
end
xlabel('seconds since 2020-01-01 00:00:00 GPS')

%% Step 8 - The sensor's partials: why a bias and a scale are estimated
% The reading is linear in the bias and in S (and in dr for a fixed S: the coupling term is
% S times C dr, a product), so |[J, names] = s.partials(dyn)| - (3M) x 15 for
% theta = [bias; S(:); dr], axis fastest, through the same filter and sampling - is the design
% matrix of a calibration. With dr = 0, as here, J theta reproduces the reading exactly. That is
% the structure |app.fitCalibration| exploits (with a bias that may vary in time, many days and a
% day bootstrap). Step 6 showed what is at stake: a bias of the order of micro-m/s^2 against a
% signal of tens of nano-m/s^2, and scales a few per cent from 1. Next: Tutorial 9.
th = [b; Sc(:); zeros(3, 1)];
sP = NonGravForces.sensor.Electrostatic(b, Sc, zeros(3, 1), "crn", 'OutputRate', 1);
[Jp, pn] = sP.partials(dyn10);
yP = sP.observe(dyn10);
yJ = reshape(Jp * th, 3, []);
fprintf('partials: %d x %d, columns %s ... %s\n', size(Jp, 1), size(Jp, 2), pn(1), pn(end));
fprintf('max |J theta - y| on valid samples = %.1e m/s^2 (signal %.1e)\n', ...
    max(abs(yJ(:, yP.mask(1, :)) - yP.value(:, yP.mask(1, :))), [], 'all'), max(abs(yP.value(:, yP.mask(1, :))), [], 'all'));
