function res = calibrate(acc, reference, options)
%CALIBRATE Calibrate an accelerometer arc against a NonGravForces Model or a reference arc, through timeSeriesPreparation.run and lsqvce.fit.
%
%   res = NonGravForces.app.calibrate(acc, M, 'Orbit', orb, 'Attitude', att)   % a priori from the Model
%   res = NonGravForces.app.calibrate(acc, refArc)                              % an external reference
%   res = NonGravForces.app.calibrate(acc, M, 'Orbit', orb, 'Fit', {'BiasBasis', "pspline", 'Scale', "full"})
%
%   THE APPLICATION OF CONCEPT 4 AND 6.2: one call from an observation arc to
%   a calibration, with every instrument stage on the one code path of
%   timeSeriesPreparation.run (v3 F4) and every number estimated by
%   lsqvce.fit (through app.fitCalibration). It replaces the matrix-in,
%   struct-out app.calibrate of 0.12.0 outright (owner ruling 2026-10-14).
%
%   1. REFERENCE. With a Model, the a priori acceleration is evaluated at the
%      ORBIT's epochs -- body frame, per effect (drag, radiation, thermal: the
%      arrays Model.budget sums) -- over the observation span widened by half
%      a CRN window, and interpolated onto the kernel's raw grid (10 Hz for
%      "ACC") by pchip, which is shape-preserving and so does not overshoot an
%      eclipse transition. The Model is not evaluated at 10 Hz: that would cost
%      ten times as much for a signal the CRN band (0.035 Hz) cannot see the
%      difference in. With an Arc, the reference is taken as given and must
%      already lie on the raw grid.
%   2. INSTRUMENT STAGES. timeSeriesPreparation.run(acc, Reference = ref):
%      mask, crn (the reference filtered like the L1B data), align, qc. Its
%      layers become zero weights of the arc the fit sees. Its "jump" events
%      are converted from UTC to GPS ns through referenceSystem and reported
%      in diagnostics.jumps; with 'Jumps' true they become step parameters of
%      the bias (CONCEPT 6.1 row 2: a jump is a Heaviside column, never a
%      mask). THAT IS OPT-IN, measured: qc screens the UNCALIBRATED residual,
%      obs minus model, where a scale error and a bias drift look like level
%      changes to the change-point detector. On a 0.25-day synthetic arc the
%      spurious jumps it found, fitted as steps, put the scale 22 bootstrap
%      sigma off the truth and the fit took 17 minutes (CI run of 033ac76).
%      Declare the jumps you trust with 'Fit', {'Jumps', t}, or pass 'Jumps'
%      true after reading diagnostics.jumps.
%   3. ESTIMATION. app.fitCalibration(run's arc, raw reference) with the
%      kernel of the run's crn preset and the jumps; 'Fit' passes any of its
%      name-value arguments on.
%   4. ATTRIBUTION (per arc; the per-revolution regression on drivers is
%      N-03). With a Model, the fit residual of each axis is regressed on the
%      CRN-filtered per-effect series plus a constant, through lsqvce.fit: its
%      coefficient divided by S_ii is the RELATIVE correction each effect's
%      a priori would need. When the effect columns are collinear (VIF above
%      'VifMax' from timeSeriesPreparation.observability) the axis's verdict is
%      "unresolved" and no correction is claimed (CONCEPT 6.2). With an Arc
%      reference there are no effects and the verdict is "unavailable".
%
% INPUT ARGUMENTS
%   Name       Size  Type                          Units  Description
%   ----       ----  ----                          -----  -----------
%   acc        1x1   timeSeriesPreparation.Arc     m/s^2  Observation, three axes in m/s2 (e.g. ACT1B,
%                                                         lin_accl_x/y/z); its product selects the preset
%   reference  1x1   NonGravForces.Model or Arc    m/s^2  The a priori: a Model (needs 'Orbit'), or an
%                                                         Arc with acc's three columns on the raw grid
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type    Units  Description
%   ----       ----  ----    -----  -----------
%   Orbit      1x1   struct  -      .epoch (referenceSystem epoch, GPS), .r and .v (3xK, GCRS, m and
%                                   m/s), e.g. from io.readGNV1B; required with a Model, covering the
%                                   arc plus half a CRN window
%   Attitude   1x1   struct  -      geom.AttitudeProvider for the orbit epochs; [] (default) uses the
%                                   Model's own provider
%   Config     1x1   struct  -      timeSeriesPreparation.Config; default Config(acc), with
%                                   run.filterReference set true
%   QualFlags  1x1   struct  -      Passed to run's mask stage; default none
%   Thrusters  1x1   struct  -      THR1B Arc, passed to run's mask stage; default none
%   Jumps      1x1   logical -      Turn run's "jump" events into bias steps; default false
%                                   (they are reported in diagnostics.jumps either way)
%   VifMax     1x1   double  -      Attribution's VIF gate; default 10
%   Fit        1xN   cell    -      Name-value pairs for app.fitCalibration (BiasBasis, Scale,
%                                   ScaleBasis for a time-varying scale, Noise, NBoot, ...); default {}
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   res   1x1   struct  -      timeSeriesPreparation.run's fields (arc, reference, diagnostics,
%                              correlation, flags, config, provenance), diagnostics.fit
%                              (fitCalibration's info), diagnostics.rawReference (the reference
%                              on the raw grid, as fitted) and diagnostics.jumps (1xJ int64, qc's
%                              jump epochs in GPS ns, fitted only with 'Jumps' true), then .calibration (app.Calibration) and
%                              .attribution (.effects 1xE string, .rms ExE3 a priori RMS of the
%                              filtered effects over the fitted rows, .share, .relativeCorrection
%                              Ex3, .sigma Ex3, .vif Ex3, .verdict 1x3 "resolved" | "unresolved" |
%                              "unavailable", .note)
%
% ACCURACY
%   The interpolation of the Model from the orbit epochs (1 s for GNV1B) to
%   the 10 Hz raw grid is pchip, third order in the step on smooth stretches.
%   MEASURED across a shadow entry (standard preset, 490 km, against the Model
%   evaluated at 10 Hz; tNgfCalibrate/calibrateReferenceInterpolationPrecision):
%   1.1e-10 m/s^2 on the raw grid at the shadow boundary, 9.3e-13 m/s^2 after
%   the CRN filter, 3e-5 of the filtered signal -- far below the Model's own
%   error (density, 10-30 %). The estimation accuracy is
%   app.fitCalibration's; the stages' are timeSeriesPreparation's.
%
% ERRORS
%   NonGravForces:calibrate:inputArc         acc is not a three-axis Arc in m/s2
%   NonGravForces:calibrate:inputReference   reference is neither a Model nor such an Arc
%   NonGravForces:calibrate:inputOrbit       a Model without 'Orbit', a malformed Orbit, or an orbit
%                                            that does not cover the arc
%   NonGravForces:calibrate:inputConfig      the Config names no CRN preset
%   NonGravForces:fitCalibration:*           the estimate
%   NonGravForces:Model:*                    the Model evaluation
%   timeSeriesPreparation:*, lsqvce:*, referenceSystem:*   propagated unwrapped
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md 4 (the entry and the result fields),
%   6.1 (the stages), 6.2 (the observation equation and the attribution gate);
%   Fritsch F., Carlson R. (1980) SIAM J. Numer. Anal. 17:238 (pchip);
%   Klinger B., Mayer-Guerr T. (2016) Adv. Space Res. 58:1597.
%
% EXAMPLE
%   eop = referenceSystem.readEop(eopFile);  sw = NonGravForces.io.spaceWeather(swFile);
%   [ep, r, v] = NonGravForces.io.readGNV1B(gnvFile, 'Eop', eop);
%   [epA, q] = NonGravForces.io.readSCA1B(scaFile);  att = NonGravForces.geom.AttitudeProvider(epA, q);
%   M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("standard"), 'Eop', eop, 'SpaceWeather', sw);
%   res = NonGravForces.app.calibrate(actArc, M, 'Orbit', struct('epoch', ep, 'r', r, 'v', v), 'Attitude', att);
%   res.calibration.summaryTable()
%
% See also NonGravForces.app.fitCalibration, NonGravForces.app.Calibration, timeSeriesPreparation.run, NonGravForces.Model.
%
% KEYWORDS calibration accelerometer application run lsqvce Model reference attribution jumps
% Created 2026-09-09 by Claude Opus 5 (Anthropic); rewritten 2026-10-15 by Claude (Anthropic) for Matthias Weigelt (N-02c).
arguments
    acc (1,1)
    reference (1,1)
    options.Orbit = []
    options.Attitude = []
    options.Config = []
    options.QualFlags = []
    options.Thrusters = []
    options.Jumps (1,1) logical = false
    options.VifMax (1,1) double {mustBePositive} = 10
    options.Fit cell = {}
end
unit = "calibrate";
if ~isa(acc, 'timeSeriesPreparation.Arc') || size(acc.obs, 2) ~= 3 || any(acc.units ~= "m/s2")
    NonGravForces.fail(unit, "inputArc", "acc must be one timeSeriesPreparation.Arc with three acceleration columns in m/s2 -> got a %s", class(acc));
end
cfg = options.Config;
if isempty(cfg), cfg = timeSeriesPreparation.Config(acc); end
c = toStruct(cfg);
if ~isfield(c, 'crn') || ~isfield(c.crn, 'preset') || isempty(c.crn.preset) || strlength(string(c.crn.preset)) == 0
    NonGravForces.fail(unit, "inputConfig", "the Config names no CRN preset, so the reference cannot be filtered like the L1B data -> set(cfg, ""crn"", ""preset"", ""ACC"")");
end
K = timeSeriesPreparation.crnKernel(string(c.crn.preset));
cfg = set(cfg, "run", "filterReference", true);

% --- 1. the reference on the raw grid ------------------------------------------------------
effects = struct('name', {}, 'arc', {});
if isa(reference, 'NonGravForces.Model')
    [ref, effects] = localModelReference(reference, acc, K, options.Orbit, options.Attitude, unit);
elseif isa(reference, 'timeSeriesPreparation.Arc') && size(reference.obs, 2) == 3 && all(reference.units == "m/s2")
    ref = reference;
else
    NonGravForces.fail(unit, "inputReference", "the reference must be a NonGravForces.Model or a three-axis timeSeriesPreparation.Arc in m/s2 -> got a %s", class(reference));
end

% --- 2. the instrument stages --------------------------------------------------------------
runArgs = {'Reference', ref, 'Config', cfg};
if ~isempty(options.QualFlags), runArgs = [runArgs, {'QualFlags', options.QualFlags}]; end
if ~isempty(options.Thrusters), runArgs = [runArgs, {'Thrusters', options.Thrusters}]; end
res = timeSeriesPreparation.run(acc, runArgs{:});

jumps = zeros(1, 0, 'int64');
ev = res.flags.events;
if ~isempty(ev) && height(ev) > 0
    isJump = ev.class == "jump";
    if any(isJump)
        e = referenceSystem.convertTimeScale(referenceSystem.datetimeToEpoch(ev.tStart(isJump), "UTC"), "GPS");
        jumps = unique(referenceSystem.epochTo(e, "gpsNanoseconds")).';
    end
end

% --- 3. the estimate -----------------------------------------------------------------------
fitJumps = zeros(1, 0, 'int64');
if options.Jumps, fitJumps = jumps; end
fitArgs = localMerge({'Kernel', K, 'Jumps', fitJumps}, options.Fit);
[cal, info] = NonGravForces.app.fitCalibration(res.arc, ref, fitArgs{:});
res.diagnostics.fit = info;
res.diagnostics.jumps = jumps;
res.diagnostics.rawReference = ref;
res.calibration = cal;

% --- 4. attribution ------------------------------------------------------------------------
res.attribution = localAttribution(effects, cal, info, res.arc, K, options.VifMax);
end

% =========================================================================================
function [ref, effects] = localModelReference(M, acc, K, orb, att, unit)
%LOCALMODELREFERENCE The Model at the orbit epochs, per effect, pchip onto the kernel's raw grid.
if ~isstruct(orb) || ~all(isfield(orb, {'epoch', 'r', 'v'}))
    NonGravForces.fail(unit, "inputOrbit", "a Model reference needs 'Orbit' with .epoch, .r and .v (e.g. io.readGNV1B) -> the a priori is evaluated along the orbit");
end
ep = orb.epoch;
if string(ep.scale) ~= "GPS", ep = referenceSystem.convertTimeScale(ep, "GPS"); end
tOrb = referenceSystem.epochTo(ep, "gpsNanoseconds");
n = numel(tOrb);
if ~isequal(size(orb.r), [3 n]) || ~isequal(size(orb.v), [3 n]) || any(diff(tOrb) <= 0)
    NonGravForces.fail(unit, "inputOrbit", "'Orbit' needs .r and .v of 3x%d and strictly increasing epochs -> check the orbit product", n);
end
step = int64(round(1e9 / K.fs));
half = int64(K.Nh) * step;
t0 = acc.t(1);
lo = max(acc.t(1) - half, tOrb(1)); hi = min(acc.t(end) + half, tOrb(end));
k0 = localCeilDiv(lo - t0, step); k1 = localFloorDiv(hi - t0, step);
if k1 <= k0
    NonGravForces.fail(unit, "inputOrbit", "the orbit (%d epochs) does not overlap the arc -> supply the orbit of the same span", n);
end
tRaw = t0 + (k0:k1).' * step;
use = find(tOrb >= tRaw(1), 1) - 1:find(tOrb <= tRaw(end), 1, 'last') + 1;
use = use(use >= 1 & use <= n);
epUse = referenceSystem.daySecondToEpoch([ep.day(use), ep.sec(use)], ep.scale);
q = [];
if ~isempty(att), q = att.at(epUse); end
[a, ai] = M.accel(epUse, orb.r(:, use), orb.v(:, use), q);
ts = double(tOrb(use) - t0) / 1e9; tq = double(tRaw - t0) / 1e9;
names = cellstr(acc.names);
ref = localArc(tRaw, pchip(ts, a, tq).', names, "MODEL");
eff = ["drag" "radiation" "thermal"];
effects = struct('name', {}, 'arc', {});
for k = 1:numel(eff)
    effects(k).name = eff(k);
    effects(k).arc = localArc(tRaw, pchip(ts, ai.(eff(k)), tq).', names, "MODEL_" + upper(eff(k)));
end
end

function att = localAttribution(effects, cal, info, arc, K, vifMax)
%LOCALATTRIBUTION Per axis, the fit residual on the filtered per-effect a priori, through lsqvce.fit.
if isempty(effects)
    att = struct('effects', strings(1, 0), 'rms', zeros(0, 3), 'share', zeros(0, 3), ...
        'relativeCorrection', zeros(0, 3), 'sigma', zeros(0, 3), 'vif', zeros(0, 3), ...
        'verdict', repmat("unavailable", 1, 3), ...
        'note', "the reference is an Arc: per-effect series need a NonGravForces.Model");
    return
end
E = numel(effects);
H = cell(1, E);
for k = 1:E
    H{k} = timeSeriesPreparation.crnApply(effects(k).arc, K, OutputEpochs = arc.t).obs;
end
rmsE = nan(E, 3); relCorr = nan(E, 3); sig = nan(E, 3); vif = nan(E, 3); verdict = strings(1, 3);
for i = 1:3
    rows = info.rows(:, i);
    X = zeros(nnz(rows), E);
    for k = 1:E, X(:, k) = H{k}(rows, i); end
    rmsE(:, i) = sqrt(mean(X.^2, 1)).';
    live = rmsE(:, i) > 0;                          % an effect the preset switched off has no column
    r = info.residual.obs(rows, i);
    D = [ones(nnz(rows), 1), X(:, live)];
    full = nan(numel(arc.t), size(D, 2)); full(rows, :) = D;
    rep = timeSeriesPreparation.observability(localOneAxis(arc.t, arc.obs(:, i)), full, MinWeightFraction = eps);
    v = rep.vif(2:end, 1);
    vif(live, i) = v;
    s = sqrt(mean(D.^2, 1)); s(s == 0) = 1; sl = sqrt(mean(r.^2)); if sl == 0, sl = 1; end
    R = lsqvce.fit(struct('A', D ./ s, 'l', r / sl), lsqvce.Options(Verbose = false));
    x = R.x(:) * sl ./ s(:);
    u = sl^2 ./ (s(:) * s(:).');
    C = R.sol.Cxx .* u;
    relCorr(live, i) = x(2:end) / cal.scale(i, i);
    sig(live, i) = sqrt(diag(C(2:end, 2:end))) / abs(cal.scale(i, i));
    if any(~(v <= vifMax)), verdict(i) = "unresolved"; else, verdict(i) = "resolved"; end
end
total = sqrt(sum(rmsE.^2, 1));
att = struct('effects', [effects.name], 'rms', rmsE, 'share', rmsE ./ total, ...
    'relativeCorrection', relCorr, 'sigma', sig, 'vif', vif, 'verdict', verdict, ...
    'note', "per arc: residual = sum_k c_k H e_k + const; relativeCorrection = c_k / S_ii, formal sigma under white noise; unresolved where an effect column's VIF exceeds the gate");
end

function a = localOneAxis(t, y)
%LOCALONEAXIS One axis as a one-column arc, for observability.
names = {'y'};
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", struct('y', 'm/s2'), "fileName", "ATTRIBUTION", "product", "ATTRIBUTION");
a = timeSeriesPreparation.Arc(t, y, aux);
end

function arc = localArc(t, X, names, product)
%LOCALARC A reference arc with the observation's column names, as timeSeriesPreparation.Arc builds one.
u = struct();
for k = 1:numel(names), u.(names{k}) = 'm/s2'; end
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", u, "fileName", product, "product", product);
arc = timeSeriesPreparation.Arc(t, X, aux);
end

function args = localMerge(defaults, given)
%LOCALMERGE Name-value pairs: the defaults, each replaced by a later pair of the same name.
names = string(defaults(1:2:end));
args = defaults;
for k = 1:2:numel(given)
    j = find(strcmpi(names, string(given{k})), 1);
    if isempty(j)
        args = [args, given(k:k + 1)]; %#ok<AGROW>
    else
        args{2 * j} = given{k + 1};
    end
end
end

function q = localCeilDiv(a, b)
q = idivide(a, b, 'ceil');
end

function q = localFloorDiv(a, b)
q = idivide(a, b, 'floor');
end
