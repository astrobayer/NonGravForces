%% Tutorial 9 - Calibrating the accelerometer: bias and scale
% An electrostatic accelerometer does not report the non-gravitational acceleration a; it reports
% y = b(t) + S a, with a bias b(t) of order 1e-7 m/s^2 that drifts and a scale S near one, and
% Level-1B (ACT1B) has passed the reading through the CRN low-pass filter. Calibrating means
% estimating b(t) and S against an a-priori model of a. This tutorial does it twice: first on a
% SYNTHETIC pair whose bias and scale are known because they were planted, so that every number
% the estimator returns can be checked against the truth; then on the real GRACE-C excerpt of
% 2020-01-01, against the standard Model, with the thruster firings of THR1B masked.
%
% *What you will learn*
%
% * what |app.syntheticPair| plants and how it passes the truth through the real Level-1B
%   (CRN) chain;
% * how |app.calibrate| composes the instrument stages (|timeSeriesPreparation.run|) and the
%   estimator (|app.fitCalibration| over |lsqvce.fit|);
% * why the bootstrap sigma, not the formal one, is the uncertainty to quote, and what the
%   inflation and the errors-in-variables (EIV) bound say;
% * how to use an |app.Calibration|: scale, bias curve, |apply|, |invert|, |save|/|load|,
%   |summaryTable|, and its figure |app.plotCalibration|;
% * how to calibrate real ACT1B against a Model along GNV1B with SCA1B attitude, with THR1B
%   firings masked, and how to declare a noise law with |lsqvce.Noise|.
%
% *Requirements:* none - the real data are the committed excerpts in |data/|; the synthetic pair
% is generated. timeSeriesPreparation and lsqvce must be installed (Tutorial 1, Step 1). Run time:
% one to two minutes.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - A planted truth: app.syntheticPair
% |app.syntheticPair| builds, on the 10 Hz raw grid of the instrument, an orbital-band model
% acceleration (an along-track drag-like signal, a small cross-track and a radial part), forms
% the RAW reading |S model + b(t) + noise + spikes|, and passes it through the real Level-1B CRN
% filter (|timeSeriesPreparation.crnKernel("ACC")|), sampled at the whole seconds as ACT1B is.
% It returns three things:
%
% * |obs| - the 1 Hz observation, a |timeSeriesPreparation.Arc| (epochs |obs.t| in int64 GPS
%   nanoseconds since 2000-01-01 12:00, values |obs.obs| K x 3 in m/s^2, body frame), NaN where
%   the 140.7 s filter window met an edge or the planted gap;
% * |ref| - the model on the 10 Hz grid, the reference the estimator filters exactly as the
%   observation was filtered;
% * |truth| - every planted number: the scale matrix, the bias at the start and its drift per
%   day, the noise, the spikes, the gap.
%
% Everything random comes from |'Seed'|, so the same seed gives the same arcs bit for bit.
[obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 7);
fprintf('observation: %d epochs at 1 Hz (%d with a NaN), product %s\n', numel(obs.t), ...
    nnz(any(isnan(obs.obs), 2)), obs.product);
fprintf('reference:   %d epochs at 10 Hz\n', numel(ref.t));
fprintf('planted scale diag  [%.4f %.4f %.4f]\n', diag(truth.scale));
fprintf('planted bias at t0  [%.3e %.3e %.3e] m/s^2, drift [%.2e %.2e %.2e] m/s^2/day\n', ...
    truth.bias0, truth.biasDrift_per_day);
fprintf('planted spikes %d, gap %s\n', height(truth.spikes), mat2str(~isempty(truth.gap)));
hObs = double(obs.t - obs.t(1)) / 3600e9;
hRef = double(ref.t - obs.t(1)) / 3600e9;
figure
for j = 1:3
    subplot(3, 1, j)
    plot(hRef(1:10:end), 1e9 * ref.obs(1:10:end, j), hObs, 1e9 * obs.obs(:, j))
    ylabel('nm/s^2'), grid on
    title("axis " + j)
end
xlabel('hours since the first observation epoch')
legend(["model (reference)" "observation y = b + S a, CRN-filtered"], 'Location', 'southoutside')

%% Step 2 - Calibrate the synthetic pair
% |app.calibrate(acc, ref)| takes the observation Arc and a reference - here the 10 Hz model Arc
% of Step 1, later a |NonGravForces.Model|. It runs |timeSeriesPreparation.run| (mask, the CRN
% filter of the reference, alignment, quality control; what QC removes becomes a zero weight) and
% then |app.fitCalibration|, which solves, per axis,
%
% $$y_i = H b_i + \sum_j S_{ij} (H a_j) + e_i$$
%
% with H the CRN filter applied to every deterministic column. 'Fit' passes options to the
% estimator: the bias basis is |"linear"| by default; 'FixCrossTrack' false estimates S_yy too
% (by default it is held at one, because the real cross-track signal is too small to calibrate a
% scale against); 'NBoot' sets the number of bootstrap replicates and 'Seed' their seed.
res = NonGravForces.app.calibrate(obs, ref, 'Fit', {'FixCrossTrack', false, 'NBoot', 50, 'Seed', 2});
cal = res.calibration;
fprintf('stages run: %s\n', strjoin(res.diagnostics.stages, ', '));
fprintf('rows fitted per axis: %s of %d\n', mat2str(sum(res.diagnostics.fit.rows, 1)), numel(obs.t));
fprintf('%-5s %10s %10s %12s %10s\n', 'axis', 'truth', 'estimate', 'boot sigma', 'error/sig');
for i = 1:3
    sb = cal.uncertainty.sigmaScaleBootstrap(i, i);
    fprintf('S_%d%d  %10.5f %10.5f %12.2e %10.2f\n', i, i, truth.scale(i, i), cal.scale(i, i), sb, ...
        abs(cal.scale(i, i) - truth.scale(i, i)) / sb);
end

%% Step 3 - Bootstrap against formal sigma, inflation and the EIV bound
% The least-squares solve gives FORMAL sigmas, which assume white noise. A CRN-filtered residual is
% anything but white - it is band-limited at 0.035 Hz - so the formal sigmas are far too small.
% With 'NBoot' > 0 the estimator therefore refits moving-block bootstrap replicates of the
% residual, and the covariance the Calibration carries is the BOOTSTRAP one;
% |uncertainty.inflation| is the ratio bootstrap / formal. The second check is the
% errors-in-variables bound: the model a is itself uncertain, which attenuates a regression scale
% towards zero by a factor lambda. |uncertainty.eiv.lambda| is a lower bound on that factor
% (computed with the residual variance as an upper bound on the model error), and
% |eiv.scaleBiasBound| = 1 - lambda the worst-case relative bias of the scale. It is a bound,
% not a correction. The gates report the observability checks the fit passed.
u = cal.uncertainty;
fprintf('%-5s %12s %12s %10s %12s %12s\n', 'axis', 'formal', 'bootstrap', 'inflation', 'EIV lambda', 'resid RMS');
for i = 1:3
    fprintf('%-5d %12.2e %12.2e %10.1f %12.6f %12.3e\n', i, u.sigmaScaleFormal(i, i), ...
        u.sigmaScaleBootstrap(i, i), u.inflation(i, i), u.eiv.lambda(i), u.residualRms(i));
end
fprintf('noise model %s, %d replicates; gates: %s (condition numbers %s)\n', u.noise, u.nBoot, ...
    strjoin(cal.gates.verdict, ' / '), mat2str(cal.gates.cond, 3));

%% Step 4 - The Calibration object: scale and bias curve
% |res.calibration| is an |app.Calibration|. It holds the OBSERVATION direction y = b(t) + S a
% (|cal.convention|), the 3x3 |scale| with its 9x9 covariance |covScale|, and the bias as a CURVE:
% |bias.basis| names the basis, |bias.coef| (3xP) its coefficients, and |biasAt(t)| evaluates it at
% int64 GPS-ns epochs. For the linear basis the coefficients are the value at |bias.t0| and the
% slope per day, which can be put beside the planted bias directly. |summaryTable| lists every
% parameter with its sigma.
fprintf('convention %s, bias basis %s, %d coefficients per axis\n', cal.convention, cal.bias.basis, size(cal.bias.coef, 2));
b0 = cal.biasAt(obs.t(1));
bTrue0 = truth.bias(:, 1);
fprintf('bias at the first epoch: estimate [%.4e %.4e %.4e], planted [%.4e %.4e %.4e] m/s^2\n', b0, bTrue0);
disp(cal.summaryTable())

%% Step 5 - apply and invert
% |cal.apply(arc)| inverts the observation equation epoch by epoch, a = S \ (y - b(t)), and returns
% a new Arc of calibrated accelerations. For the synthetic pair the right answer is known: the
% CRN-filtered model (|timeSeriesPreparation.crnApply| of the reference at the observation
% epochs); what is left is the planted noise and spikes. |cal.invert()| returns the same
% calibration in the CORRECTION direction a = b'(t) + S' y (the one GRACE TN-02 uses), with the
% covariance propagated to first order; applying either gives the same accelerations. A number
% quoted without its direction is the commonest way two calibrations disagree.
aCal = cal.apply(obs);
K = timeSeriesPreparation.crnKernel("ACC");
Href = timeSeriesPreparation.crnApply(ref, K, 'OutputEpochs', obs.t);
d = aCal.obs - Href.obs;
ok = all(isfinite(d), 2);
fprintf('calibrated minus filtered model, RMS over %d epochs: [%.3e %.3e %.3e] m/s^2\n', nnz(ok), ...
    sqrt(mean(d(ok, :) .^ 2, 1)));
calInv = cal.invert();
aInv = calInv.apply(obs);
fprintf('inverted calibration: convention %s, S''_xx = %.6f (1/S_xx = %.6f)\n', calInv.convention, ...
    calInv.scale(1, 1), 1 / cal.scale(1, 1));
fprintf('apply(invert) vs apply, largest difference: %.3e m/s^2\n', max(abs(aInv.obs(ok, :) - aCal.obs(ok, :)), [], 'all'));
figure
plot(hObs, 1e9 * d(:, 1), hObs, 1e9 * d(:, 3)), grid on
xlabel('hours'), ylabel('nm/s^2'), legend(["x" "z"])
title('calibrated observation minus the filtered model')

%% Step 6 - save and load
% A Calibration is saved as a versioned MAT file (|schemaVersion|), so a file written today is
% read by a later release - with a warning if it was written by a newer schema. Write only to
% |tempdir| in a tutorial.
file = fullfile(tempdir, 'ngfTutorial09_calibration.mat');
cal.save(file);
cal2 = NonGravForces.app.Calibration.load(file);
fprintf('saved to %s and loaded: scale identical %d, bias identical %d\n', file, ...
    isequal(cal2.scale, cal.scale), isequal(cal2.bias.coef, cal.bias.coef));
delete(file);

%% Step 7 - The calibration figure
% |app.plotCalibration(cal, t)| draws, per axis, the scale with its 1-sigma band, the bias curve
% with its band and - with 'Residual' - the residual's autocorrelation before and after a
% low-order AR whitening, with the 95 % white-noise band. The first ACF is why the formal sigmas
% lie; the second whether a simple AR law would have explained the colour. It opens its own
% figure and draws only numbers the Calibration's methods return.
fig = NonGravForces.app.plotCalibration(cal, obs.t, 'Residual', res.diagnostics.fit.residual);
fprintf('figure "%s" with %d axes\n', fig.Name, numel(findobj(fig, 'Type', 'axes')));

%% Step 8 - The real day: ACT1B, GNV1B, SCA1B and THR1B
% Now the same on real data: the first ten minutes (600 records at 1 Hz) of GRACE-C's ACT1B of
% 2020-01-01. The observation must be a |timeSeriesPreparation.Arc| with three columns in m/s2;
% its product name ("ACT1B") selects the CRN preset. The reference comes from a Model evaluated
% along the GNV1B orbit (rotated into the GCRS with the EOP; |'Orbit'| is a struct with .epoch,
% .r and .v) with the SCA1B attitude (BODY -> GCRS). THR1B, the thruster firings, is read with
% |gracefoReading.readL1b| into an Arc whose |on_time_...| columns (ms) the mask stage reads: every
% firing, widened by the preset's margin, is cut before the fit.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
[ep, r, v] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
orb = struct('epoch', ep, 'r', r, 'v', v);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
[epC, aMeas] = NonGravForces.io.readACC1B(f('samples', 'ACT1B_2020-01-01_C_04.txt'));
names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", units, "fileName", "ACT1B_2020-01-01_C_04.txt", "product", "ACT1B");
acc = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), aMeas.', aux);
[tt, ot, at] = gracefoReading.readL1b(f('samples', 'THR1B_2020-01-01_C_04.txt'), "TimeFormat", "gpsns");
thr = timeSeriesPreparation.Arc(tt, ot, at);
on = startsWith(thr.names, "on_time_");
fired = any(thr.obs(:, on) > 0, 2);
inArc = thr.t >= acc.t(1) & thr.t <= acc.t(end);
fprintf('ACT1B: %d epochs; GNV1B: %d epochs; SCA1B: %d quaternions\n', numel(acc.t), numel(ep.sec), size(q, 2));
fprintf('THR1B: %d records over the day, %d with a firing, %d of those inside the ACT1B excerpt\n', ...
    numel(thr.t), nnz(fired), nnz(fired & inArc));

%% Step 9 - Calibrate the real day against the standard Model
% The a priori is the standard preset as shipped, so it needs its observed density scale
% (Tutorial 1, Step 6). |app.calibrate(acc, M, 'Orbit', orb, 'Attitude', att, 'Thrusters', thr)|
% evaluates the Model per effect at the orbit epochs over the arc plus half a CRN window,
% interpolates it onto the 10 Hz raw grid (pchip), and fits as in Step 2. S_yy is held at one
% (the default). With a Model reference there is also an ATTRIBUTION: the residual of each axis
% is regressed on the filtered drag, radiation and thermal series, and the coefficients say
% which a-priori effect would need correcting - or "unresolved" when the effects are collinear,
% which on ten minutes of data they will usually be.
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("standard"), 'Attitude', att, ...
    'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
resR = NonGravForces.app.calibrate(acc, M, 'Orbit', orb, 'Attitude', att, 'Thrusters', thr, 'Fit', {'NBoot', 50});
calR = resR.calibration;
fprintf('samples with zero weight after the stages (filter edges, firings, QC): %d of %d\n', ...
    nnz(resR.arc.weights == 0), numel(resR.arc.weights));
fprintf('S_xx %.4f +- %.4f, S_zz %.4f +- %.4f (bootstrap), S_yy held at %.1f\n', calR.scale(1, 1), ...
    calR.uncertainty.sigmaScaleBootstrap(1, 1), calR.scale(3, 3), calR.uncertainty.sigmaScaleBootstrap(3, 3), calR.scale(2, 2));
fprintf('bias at the first epoch [%.3e %.3e %.3e] m/s^2\n', calR.biasAt(acc.t(1)));
fprintf('residual RMS [%.3e %.3e %.3e] m/s^2\n', calR.uncertainty.residualRms);
A = resR.attribution;
fprintf('attribution per axis: %s\n', strjoin(A.verdict, ' / '));
disp(array2table(A.relativeCorrection, 'RowNames', cellstr(A.effects), 'VariableNames', {'x', 'y', 'z'}))
pred = resR.diagnostics.fit.prediction;
tm = double(acc.t - acc.t(1)) / 1e9;
figure
ax = ["x (along track)" "y (cross track)" "z (radial)"];
for j = 1:3
    subplot(3, 1, j)
    plot(tm, 1e9 * acc.obs(:, j), tm, 1e9 * pred.obs(:, j)), grid on
    ylabel('nm/s^2'), title(ax(j))
end
xlabel('seconds since 2020-01-01 00:00 GPS')
legend(["ACT1B" "b + S H a (fitted)"], 'Location', 'southoutside')

%% Step 10 - A declared noise law
% The default 'Noise' is "white", made honest by the bootstrap. Alternatively a noise LAW can be
% declared for every axis with an |lsqvce.Noise| object (lsqvce >= 7.0; the struct of lsqvce 6 is
% refused): |lsqvce.Noise.arma(ar, sigma)| is a known whitening filter with AR coefficients
% [1 a1 ... ap], |lsqvce.Noise.white(sigma)| white noise, |lsqvce.Noise.psd(f, S)| a measured
% spectrum. Here an AR(1) law, without a bootstrap, so the sigmas are the formal ones under that
% law. The reference is handed back as the Arc the previous call fitted
% (|diagnostics.rawReference|), which gives the same calibration as the Model itself would and
% skips a second Model evaluation.
law = lsqvce.Noise.arma([1 -0.5], 1);
resN = NonGravForces.app.calibrate(acc, resR.diagnostics.rawReference, 'Thrusters', thr, ...
    'Fit', {'NBoot', 0, 'Noise', law});
calN = resN.calibration;
fprintf('noise %s: S_xx %.4f (formal sigma %.2e), S_zz %.4f (formal sigma %.2e)\n', calN.uncertainty.noise, ...
    calN.scale(1, 1), calN.uncertainty.sigmaScaleFormal(1, 1), calN.scale(3, 3), calN.uncertainty.sigmaScaleFormal(3, 3));
fprintf('white + bootstrap (Step 9): S_xx %.4f, S_zz %.4f\n', calR.scale(1, 1), calR.scale(3, 3));
