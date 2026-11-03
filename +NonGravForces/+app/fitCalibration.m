function [cal, info] = fitCalibration(obs, ref, options)
%FITCALIBRATION Estimate an accelerometer calibration from a filtered observation arc and a raw-rate model arc, through lsqvce.fit.
%
%   cal = NonGravForces.app.fitCalibration(obs, ref)
%   [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'Scale', "full", 'BiasBasis', "pspline", 'NBoot', 0)
%   cal = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear")   % S(t), N-03a
%   cal = NonGravForces.app.fitCalibration(obsDays, refDays)                  % 1xN: one scale, N biases (N-09)
%
%   The estimator of the calibration (N-02b; CONCEPT 6.2). Per axis i, with H
%   the CRN kernel of timeSeriesPreparation.crnKernel:
%
%       y_i(t) = H*b_i(t) + sum_j S_ij (H*a_j)(t) [+ K2_i (H*a_i^2)(t), demeaned] + e_i(t)
%
%   y is the observation `obs` (e.g. ACT1B at 1 Hz), a the model `ref` on the
%   kernel's raw grid (10 Hz), b_i(t) = B(t) c_i the bias curve on the basis
%   NonGravForces.app.Calibration.basis evaluates. EVERY deterministic column
%   is filtered by H exactly as the observation was -- the bias basis, the
%   model and the K2 regressor -- by timeSeriesPreparation.crnApply at the
%   observation epochs, so the fitted curve is the one the instrument carried
%   before the L1B filter (CONCEPT 6.2). The K2 column is H applied to a_i^2,
%   because the nonlinearity acts on the instrument's raw signal, before the
%   filter; (H a_i)^2 differs from it by the filter's smoothing of a square.
%   Its mean over the fitted rows goes into the bias, so K2 is reported with a
%   bias that carries -K2 * mean(H a_i^2).
%
%   THE THREE AXES ARE THREE PROBLEMS (CONCEPT 6.6), each one lsqvce.fit call
%   on declaration kind A: one system, one day, one arc per contiguous run of
%   usable rows (a gap or a NaN window of the filter ends an arc). An entry of
%   S that is not estimated is held at its a-priori value and moved to the
%   observation side. The columns are equilibrated to unit RMS before the
%   solve and the estimates scaled back, so a bias of 1e-7 m/s^2 and a scale
%   of 1 are solved at the same precision.
%
%   SCALE. "diag" estimates S_ii, "full" all nine, a 3x3 logical mask any
%   pattern. The cross-track scale S_yy is held at its a-priori value by
%   default ('FixCrossTrack', true; CONCEPT 6.2: the cross-track signal is too
%   small to calibrate a scale against). A row that estimates more than one
%   model column is GATED: when the variance inflation factor of any of its
%   model columns (timeSeriesPreparation.observability) exceeds 'VifMax', the
%   call refuses rather than report a shear the data cannot separate (v3 I.5-3).
%
%   TIME-VARYING SCALE (N-03a; CONCEPT 6.2, the v3 plan's WP-14). With
%   'ScaleBasis' other than "constant" every estimated entry becomes a curve,
%   S_ij(t) = B_s(t) g_ij, on the basis NonGravForces.app.Calibration.basis
%   evaluates, and its Q columns are H(B_s,q .* a_j): the scale multiplies the
%   RAW signal, before the filter, exactly as the bias curve does. The
%   Calibration then carries the curve (scaleCurve), its pointwise sigma band
%   and a Wald test of constancy (drift, isDriftSignificant); `scale` is the
%   curve at the middle of the fitted span. Identifiability against the bias is
%   by BASIS SIZE, not by a penalty: the scale must vary no faster than the
%   bias, so a hat or B-spline scale with a knot spacing below the bias's is
%   refused (v3 WP-14 asked the same asymmetry of a smoothness prior with
%   lambdas by blocked CV; a coarse basis states it without a weight to tune,
%   and lsqvce priors remain the route if a penalty is ever needed). Entries
%   held fixed stay constant.
%
%   MEASURED REGRESSORS (N-11a; CONCEPT 6.2). 'Regressors' is an Arc of K
%   measured series -- accelerometer temperatures from AHK1B, magnetorquer
%   currents from MAG1B, any housekeeping -- on its own grid. Each is
%   interpolated LINEARLY onto the kernel's raw grid (no overshoot between
%   housekeeping samples; NaN outside its span, so those rows are not fitted),
%   demeaned over that span, filtered by H like the bias, and enters every axis
%   as y_i += c_ik (z_k - m_k). With one arc its columns join the bias block;
%   with several in common mode its coefficients are COMMON, a physical coupling
%   not changing from day to day, and each day keeps its own mean m_k. A series
%   with no finite value or no variation on the raw span is refused. That check
%   does NOT catch every dead channel: TICUR reads -273.205 to -273.202 on
%   2020-01-01 -- absolute zero plus a 3 mK quantisation wobble -- so it varies
%   and would be fitted; choosing live channels is the caller's. MEASURED: the GRACE-FO sensor-unit
%   temperatures span 0.022 C over that whole day (thermal control), so within
%   one day a temperature column is nearly the constant bias column and its
%   coefficient is weak; it is a campaign quantity.
%
%   MANY ARCS, ONE SCALE (N-09; CONCEPT 6.6). With obs a 1xN Arc array (one
%   arc per day, in time order, not overlapping) and ref 1xN or one arc
%   covering them all, 'ScaleMode' "common" (the default for N > 1) solves,
%   per axis, ONE lsqvce kind-A system whose days are the arcs and whose arcs
%   are each day's gap-free runs: the scale (and K2) columns are COMMON, each
%   day's bias basis is a 'day'-scope LOCAL block, pre-eliminated by lsqvce
%   and recovered by back-substitution. The result is a 1xN Calibration
%   array: the common scale in every element, each with its own day's bias.
%   "perArc" is N independent single-arc fits (seed Seed + k - 1 for arc k).
%   Common mode supports a "constant" or "linear" bias and a constant scale
%   only, and refuses 'Jumps' (unsupportedMultiArc): a day's local block must
%   have one size, and a common curve over days is not defined here. Its
%   uncertainty is a bootstrap over DAYS: timeSeriesPreparation.bootstrap on
%   the day index with block length 1 -- days are independent where the
%   samples of one day are not. A replicate is the weighted sum of the drawn
%   days' reduced normal equations (lsqvce >= 6.35: R.contributions and
%   lsqvce.solveContributions, a day drawn twice counted twice), so the data
%   pass is paid once and a replicate costs O(N p^2). Under "white" noise and
%   no 'Robust' that IS the joint refit of the drawn days; an estimated noise
%   model or robust weights are held at the full fit's values in every
%   replicate (a bootstrap under fixed weights). Below 10 days it warns (fewArcs). The per-day biases carry
%   lsqvce's FORMAL local covariance (s0^2 Nll^-1 plus the propagated common
%   part), not a bootstrap one. The gates are computed on the model columns
%   with each day's bias basis projected out -- the columns the pre-eliminated
%   solve actually sees.
%
%   NOISE. "white" (default) solves under unit weights; the formal sigmas
%   then assume white noise, which the CRN-filtered residual is not. That is
%   why 'NBoot' (default 200) adds a moving-block bootstrap of every parameter
%   through timeSeriesPreparation.bootstrap, refitting each replicate with
%   lsqvce.fit: the covariance the Calibration carries is the BOOTSTRAP one
%   when NBoot > 0, and uncertainty.inflation reports bootstrap/formal sigma.
%   "estimate" hands the noise to lsqvce's joint variance-component estimator
%   (Estimate "noise", the law of each axis estimated, AR whitener); an
%   lsqvce.Noise law (lsqvce >= 7.0) is declared for every axis as given (e.g.
%   lsqvce.Noise.arma(ar, sigma) from timeSeriesPreparation.noise's AR fit).
%   Since 2026-11-03 the fit is declared as a lsqvce.TaskDefinition (one
%   system, the axis's commons, day-scoped bias locals, one stream "acc");
%   the v6 covariance struct is refused (inputNoise). The white default is a
%   decision, not a convenience: the filtered residual is band-limited at
%   0.035 Hz, so a whitening law estimated from it is ill-conditioned above
%   the band, while the block bootstrap needs no law at all.
%
%   TEMPERATURE-DEPENDENT SCALE (N-11b; CONCEPT 6.2). With 'ScaleTemperature'
%   naming a column T of 'Regressors', every ESTIMATED scale entry S_ij gains a
%   coefficient S1_ij = dS_ij/dT: the column H((T - T0) a_j), the measured series
%   (linear onto the raw grid, as for the regressors) times the model axis before
%   the filter. T0 is the mean of T over the raw grid of the arc -- over ALL arcs
%   in common mode, where S1 is common: a per-day reference would remove exactly
%   the day-to-day temperature differences that make S1 observable on a
%   thermally controlled instrument (GRACE-FO's TSU span 0.022 C in a day). The
%   calibration's `scale` is S(T0). The columns enter the gate; a VIF above
%   'VifMax' gives the verdict "warn", not a refusal, since S1 then carries the
%   large sigma it deserves. It is refused beside a time-varying scale
%   (unsupportedScaleTemperature), with which it is not separable on one arc.
%
%   EIV BOUND (v3 D4). The model carries the dominant error, so a regression
%   scale is attenuated toward zero by lambda = var(H a_j) / (var(H a_j) +
%   sigma_e^2) (Fuller 1987, sect. 1.1). With the residual variance of axis j,
%   divided by S_jj^2, as an UPPER bound on sigma_e^2, uncertainty.eiv.lambda
%   is a lower bound on the attenuation factor and scaleBiasBound = 1 - lambda
%   the worst-case relative bias of every scale on model axis j. It is a
%   bound, not a correction.
%
% INPUT ARGUMENTS
%   Name  Size  Type                        Units  Description
%   ----  ----  ----                        -----  -----------
%   obs   1xN   timeSeriesPreparation.Arc   m/s^2  Observation, three axes, units m/s2, already
%                                                  through the L1B chain (e.g. ACT1B); its
%                                                  epochs lie on ref's grid
%   ref   1xN   timeSeriesPreparation.Arc   m/s^2  Model accelerations, three axes, units m/s2, on
%                                                  the kernel's raw grid around obs (10 Hz for
%                                                  "ACC"); with 'Kernel' [] on obs's epochs
%
% NAME-VALUE ARGUMENTS
%   Name               Size  Type    Units  Description
%   ----               ----  ----    -----  -----------
%   Kernel             1x1   string  -      "crn" (default): crnKernel("ACC"); a crnKernel struct
%                                           (fields F, fs, Nf, Nh, Tf, preset); or [] for an
%                                           unfiltered observation
%   BiasBasis          1x1   string  -      "constant" | "linear" (default) | "hat" | "pspline"
%   BiasKnotSpacing_s  1x1   double  s      Knot spacing of "hat" and "pspline"; default 10800.
%                                           Below the kernel span Tf it warns
%   Scale              3x3   string  -      "diag" (default) | "full" | a 3x3 logical mask
%   ScaleBasis         1x1   string  -      "constant" (default) | "linear" | "hat" | "pspline":
%                                           the time dependence of every estimated entry
%   ScaleKnotSpacing_s 1x1   double  s      Knot spacing of a "hat" or "pspline" scale; default
%                                           43200, at least the bias's when both are knotted
%   ScaleMode          1x1   string  -      With several arcs: "common" (default) or "perArc"
%   Regressors         1xN   Arc     -      Measured series, one Arc of K columns (one per obs arc,
%                                           or one covering all), interpolated linearly but never
%                                           across a gap in a series (Calibration.interpolateSeries);
%                                           default none
%   ScaleTemperature   1x1   string  -      Name of one 'Regressors' column the scale depends on
%                                           (N-11b): S(T) = S + S1 (T - T0); default "" (none)
%   Span               1x2   int64   ns     [first last] epoch the bias and scale bases must cover;
%                                           default [] = the first and last usable observation epoch.
%                                           Cross-validation passes the WHOLE arc's span, so that a
%                                           B-spline is defined over a held-out block (N-03b)
%   FixCrossTrack      1x1   logical -      Hold S_yy at its a-priori value; default true
%   ScaleApriori       3x3   double  -      Values of the entries not estimated; default eye(3)
%   Quadratic          1x1   logical -      Estimate K2 per axis; default false
%   Noise              1x1   string  -      "white" (default) | "estimate" | one
%                                           lsqvce.Noise law (lsqvce >= 7.0)
%   NBoot              1x1   double  -      Bootstrap replicates, 0 for none; default 200
%   Seed               1x1   double  -      Bootstrap seed; default 0
%   VifMax             1x1   double  -      VIF gate for a row with several model columns;
%                                           default 10
%   MinWeightFraction  1x1   double  -      Minimum-data gate of observability; default 0.5
%   Jumps              1xJ   int64   ns     Epochs (GPS ns) of steps in the bias, one Heaviside
%                                           column each, filtered like the basis; default none
%   Robust             1x1   string  -      lsqvce M-estimation: "none" (default) | "huber" |
%                                           "igg3" | "bisquare"; a single arc's bootstrap refits
%                                           with it, the day bootstrap keeps its final weights
%
% OUTPUT ARGUMENTS
%   Name  Size  Type                            Units  Description
%   ----  ----  ----                            -----  -----------
%   cal   1xN   NonGravForces.app.Calibration   -      One per arc. Convention "observation"; scaleCurve when
%                                                      'ScaleBasis' is not "constant"; covScale and bias
%                                                      cov from the bootstrap when NBoot > 0, else
%                                                      formal; .uncertainty (.sigmaScaleFormal,
%                                                      .sigmaScaleBootstrap, .inflation 3x3,
%                                                      .residualRms 1x3 m/s^2, .eiv, .noise, .nBoot,
%                                                      .stopReason); .gates (.vif, .cond,
%                                                      .maxAbsCorr, .weightFraction, .vifMax,
%                                                      .verdict)
%   info  1x1   struct                          -      .residual (Arc at obs epochs, NaN where not
%                                                      fitted) .names (1x3 cell of column names)
%                                                      .rows (Nx3 logical) .Qxx (1x3 cell, the
%                                                      unit-weight cofactor of each axis, in the
%                                                      unscaled parameters) .results (1x3 cell of
%                                                      lsqvce.Results) .observability (1x3 cell)
%                                                      .prediction (Arc at obs epochs: the fitted
%                                                      observation, fixed scale entries included, at
%                                                      EVERY row whose design is finite, fitted or not;
%                                                      NaN elsewhere). With several arcs in common
%                                                      mode: .mode "common", .residual (1xN cell of
%                                                      Arc), .rows (1xN cell), .names, .Qxx and
%                                                      .results per axis; in perArc mode: .mode,
%                                                      .arcs (1xN cell of single-arc infos)
%
% ACCURACY
%   The estimate is lsqvce's QR solve on equilibrated columns: rounding of
%   order cond * eps with cond the equilibrated condition number the gates
%   report. Each filtered column carries the rounding of one Nf-term dot
%   product (crnApply). The formal sigmas assume the declared noise; the
%   bootstrap sigmas are Monte-Carlo estimates (timeSeriesPreparation.
%   bootstrap's accuracy) with the block length of the Politis-White rule on
%   the residual of the axis. The
%   bias is a curve before the filter; Calibration.apply evaluates it
%   unfiltered (its help names that approximation).
%
% ERRORS
%   NonGravForces:fitCalibration:depMissing        lsqvce >= 7.0 (lsqvce.TaskDefinition) or
%                                                  timeSeriesPreparation is not on the path
%   NonGravForces:fitCalibration:inputNoise        'Noise' neither "white", "estimate" nor one
%                                                  lsqvce.Noise (a v6 covariance struct)
%   NonGravForces:fitCalibration:inputArc          obs or ref not a three-axis Arc in m/s2; one obs
%                                                  arc with several refs; several obs arcs with
%                                                  neither one ref nor one per arc, or not in time
%                                                  order, or overlapping
%   NonGravForces:fitCalibration:inputScale        a mask that is neither "diag", "full" nor 3x3,
%                                                  or a non-finite 'ScaleApriori'
%   NonGravForces:fitCalibration:inputGrid         'Kernel' [] and an obs epoch ref does not carry
%   NonGravForces:fitCalibration:inputKernel       a kernel that is neither "crn", a struct nor []
%   NonGravForces:fitCalibration:inputJumps        'Jumps' not int64
%   NonGravForces:fitCalibration:inputScaleBasis   a knotted scale basis finer than the knotted
%                                                  bias basis
%   NonGravForces:fitCalibration:inputSpan         'Span' not an increasing int64 pair
%   NonGravForces:fitCalibration:inputRegressors   'Regressors' not an Arc (or one per arc), or a
%                                                  series with no finite value or no variation on
%                                                  the arc
%   NonGravForces:fitCalibration:unsupportedMultiArc  common mode with a knotted bias, 'Jumps' or a
%                                                  time-varying scale
%   NonGravForces:fitCalibration:inputScaleTemperature  'ScaleTemperature' names no column of
%                                                  'Regressors'
%   NonGravForces:fitCalibration:unsupportedScaleTemperature  'ScaleTemperature' beside a
%                                                  time-varying scale
%   NonGravForces:fitCalibration:numericCollinear  the VIF gate refused a row with several model
%                                                  columns
%   NonGravForces:Calibration:inputBasis           knots the basis refuses
%   timeSeriesPreparation:*                        crnApply, observability, bootstrap
%   lsqvce:*                                       the solve
%
% WARNINGS
%   NonGravForces:fitCalibration:knotsBelowFilterSpan  knot spacing below the kernel span Tf: the
%                                                      filtered basis columns are nearly dependent
%   NonGravForces:fitCalibration:fewArcs               a bootstrap over fewer than 10 days
%   NonGravForces:fitCalibration:jumpOutsideFit        a jump whose step is constant over the fitted
%                                                      rows; it is dropped from the bias
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md 6.2 (the observation equation, the
%   gates, the EIV bound) and 6.6 (three axes, three problems); the v3 plan
%   (docs/inputs/HANDOVER_AccelCal_v3.md) D4, I.5-3, WP-08, WP-09, WP-13, WP-14;
%   Fuller W. (1987) Measurement Error Models, Wiley, sect. 1.1 (attenuation);
%   Kuensch H. (1989) Ann. Statist. 17:1217 (moving-block bootstrap);
%   Kutner M. et al. (2005) Applied Linear Statistical Models, ch. 10 (VIF);
%   Klinger B., Mayer-Guerr T. (2016) Adv. Space Res. 58:1597 (K2 term).
%
% EXAMPLE
%   [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 3);
%   cal = NonGravForces.app.fitCalibration(obs, ref, 'FixCrossTrack', false, 'NBoot', 50);
%   [diag(cal.scale), diag(truth.scale)]
%
% See also NonGravForces.app.Calibration, NonGravForces.app.syntheticPair, lsqvce.fit, timeSeriesPreparation.crnApply.
%
% KEYWORDS calibration accelerometer scale bias estimator lsqvce CRN bootstrap EIV VIF
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    obs (1,:)
    ref (1,:)
    options.Kernel = "crn"
    options.BiasBasis (1,1) string {mustBeMember(options.BiasBasis, ["constant" "linear" "hat" "pspline"])} = "linear"
    options.BiasKnotSpacing_s (1,1) double {mustBePositive, mustBeFinite} = 10800
    options.Scale = "diag"
    options.ScaleBasis (1,1) string {mustBeMember(options.ScaleBasis, ["constant" "linear" "hat" "pspline"])} = "constant"
    options.ScaleKnotSpacing_s (1,1) double {mustBePositive, mustBeFinite} = 43200
    options.Span = int64([])
    options.ScaleMode (1,1) string {mustBeMember(options.ScaleMode, ["common" "perArc"])} = "common"
    options.Regressors = []
    options.ScaleTemperature (1,1) string = ""
    options.FixCrossTrack (1,1) logical = true
    options.ScaleApriori (3,3) double = eye(3)
    options.Quadratic (1,1) logical = false
    options.Noise = "white"
    options.NBoot (1,1) double {mustBeInteger, mustBeNonnegative} = 200
    options.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.VifMax (1,1) double {mustBePositive} = 10
    options.MinWeightFraction (1,1) double {mustBePositive, mustBeLessThanOrEqual(options.MinWeightFraction, 1)} = 0.5
    options.Jumps = zeros(1, 0, 'int64')
    options.Robust (1,1) string {mustBeMember(options.Robust, ["none" "huber" "igg3" "bisquare"])} = "none"
end
unit = "fitCalibration";
if isempty(which('lsqvce.TaskDefinition')) || isempty(which('timeSeriesPreparation.crnApply'))
    NonGravForces.fail(unit, "depMissing", "the estimate goes through lsqvce.fit on a lsqvce.TaskDefinition (lsqvce >= 7.0) and the filter through timeSeriesPreparation.crnApply -> run lsqvce_setup (7.0 or later) and setupTimeSeriesPreparation");
end
if ~(isa(options.Noise, 'lsqvce.Noise') && isscalar(options.Noise)) && ...
        ~((isstring(options.Noise) || ischar(options.Noise)) && any(string(options.Noise) == ["white" "estimate"]))
    NonGravForces.fail(unit, "inputNoise", "'Noise' must be ""white"", ""estimate"" or one lsqvce.Noise law -> got a %s (a v6 covariance struct is now lsqvce.Noise.<law>(...), e.g. lsqvce.Noise.arma([1 -0.5], 1))", class(options.Noise));
end
if ~isa(options.Noise, 'lsqvce.Noise'), options.Noise = string(options.Noise); end
if ~isa(obs, 'timeSeriesPreparation.Arc') || ~isa(ref, 'timeSeriesPreparation.Arc')
    NonGravForces.fail(unit, "inputArc", "obs and ref must be timeSeriesPreparation.Arc objects with three acceleration columns in m/s2 -> got a %s and a %s", class(obs), class(ref));
end
for arc = [obs, ref]
    if size(arc.obs, 2) ~= 3 || any(arc.units ~= "m/s2")
        NonGravForces.fail(unit, "inputArc", "obs and ref must each have three acceleration columns in m/s2 -> got %d columns in %s", size(arc.obs, 2), strjoin(arc.units, ", "));
    end
end
if numel(obs) > 1
    [cal, info] = localMultiArc(obs, ref, options, unit);
    return
end
if numel(ref) ~= 1
    NonGravForces.fail(unit, "inputArc", "one observation arc takes one reference arc -> got %d", numel(ref));
end
mask = localMask(options.Scale, unit);
if options.FixCrossTrack, mask(2, 2) = false; end
S0 = options.ScaleApriori;
if ~all(isfinite(S0(:)))
    NonGravForces.fail(unit, "inputScale", "'ScaleApriori' must be finite -> it supplies every entry not estimated");
end
K = localKernel(options.Kernel, unit);

% --- the bias basis, on the grid the filter reads -----------------------------------------
% The knots span the USABLE epochs, the first and the last observation epoch that carries a finite
% value (B-splines: widened by half a kernel, below). Spanning the whole arc instead put a hat
% or a B-spline whose support lies outside every fitted window, a zero column the solve rightly
% refused (lsqvce:qr:singularR, CI run 38690853810); spanning less would leave the curve
% undefined where a window looks.
t1 = obs.t(1);
live = find(obs.weights > 0 & ~obs.filled & any(isfinite(obs.obs), 2));
if isempty(live), live = 1:numel(obs.t); end
tA = obs.t(live(1)); tB = obs.t(live(end));
if ~isempty(options.Span)
    if ~isa(options.Span, 'int64') || numel(options.Span) ~= 2 || options.Span(2) <= options.Span(1)
        NonGravForces.fail(unit, "inputSpan", "'Span' must be an increasing int64 pair [first last] of GPS ns -> got %s %s", class(options.Span), mat2str(size(options.Span)));
    end
    tA = options.Span(1); tB = options.Span(2);
end
tA0 = tA; tB0 = tB;                                % the span itself, before any widening
knotted = ["hat" "pspline"];
if ismember(options.ScaleBasis, knotted) && ismember(options.BiasBasis, knotted) && options.ScaleKnotSpacing_s < options.BiasKnotSpacing_s
    NonGravForces.fail(unit, "inputScaleBasis", "a scale knot spacing of %g s is finer than the bias's %g s: a scale that varies faster than the bias is not separable from it -> raise 'ScaleKnotSpacing_s' to at least 'BiasKnotSpacing_s'", ...
        options.ScaleKnotSpacing_s, options.BiasKnotSpacing_s);
end
if ~isempty(K) && options.BiasBasis == "pspline"
    % B-splines are a partition of unity only between knot 4 and knot end-3, so their span is
    % widened to every raw sample a window reads. Hats are not widened: their linear extension
    % beyond the end knots is still a partition of unity, and a hat centred outside the usable
    % epochs would be seen by a window's edge taps alone.
    half = int64(K.Nh) * int64(round(1e9 / K.fs));
    tA = tA - half; tB = tB + half;
end
if ~isa(options.Jumps, 'int64')
    NonGravForces.fail(unit, "inputJumps", "'Jumps' must be int64 GPS ns, as timeSeriesPreparation.Arc stores time -> got %s", class(options.Jumps));
end
jumps = unique(options.Jumps(:)');                 % sorted, each once
bias = struct('basis', options.BiasBasis, 't0', t1, 'knots', int64([]), 'jumps', jumps, 'coef', [], 'cov', {{[], [], []}});
h = int64(round(options.BiasKnotSpacing_s * 1e9));
switch options.BiasBasis
    case "hat"
        bias.knots = tA + int64(0:ceil(double(tB - tA) / double(h))) * h;
    case "pspline"
        bias.knots = tA + int64(-3:ceil(double(tB - tA) / double(h)) + 3) * h;
end
if ismember(options.BiasBasis, ["hat" "pspline"]) && ~isempty(K) && options.BiasKnotSpacing_s < K.Tf
    warning("NonGravForces:fitCalibration:knotsBelowFilterSpan", ...
        "knot spacing %g s is below the kernel span %g s: the filtered basis columns are nearly dependent -> raise 'BiasKnotSpacing_s'", options.BiasKnotSpacing_s, K.Tf);
end

% --- the scale basis (N-03a): the same construction on the same span -----------------------
sc = struct('basis', options.ScaleBasis, 't0', t1, 'knots', int64([]), 'tRef', tA0 + (tB0 - tA0) / 2);
hs = int64(round(options.ScaleKnotSpacing_s * 1e9));
switch options.ScaleBasis
    case "hat"
        sc.knots = tA0 + int64(0:ceil(double(tB0 - tA0) / double(hs))) * hs;
    case "pspline"
        sA = tA0; sB = tB0;
        if ~isempty(K)
            sA = sA - int64(K.Nh) * int64(round(1e9 / K.fs)); sB = sB + int64(K.Nh) * int64(round(1e9 / K.fs));
        end
        sc.knots = sA + int64(-3:ceil(double(sB - sA) / double(hs)) + 3) * hs;
end
Q = size(NonGravForces.app.Calibration.basis(sc, t1), 2);    % 1 for a constant scale
varying = options.ScaleBasis ~= "constant";

% --- the filtered columns at the observation epochs ----------------------------------------
if isempty(K)
    [found, at] = ismember(obs.t, ref.t);
    if ~all(found)
        NonGravForces.fail(unit, "inputGrid", "with 'Kernel' [] every observation epoch must be a model epoch -> %d of %d are not", nnz(~found), numel(found));
    end
    Ha = ref.obs(at, :);
    HB = NonGravForces.app.Calibration.basis(bias, obs.t);
    Hq = ref.obs(at, :).^2;
else
    Ha = timeSeriesPreparation.crnApply(ref, K, OutputEpochs = obs.t).obs;
    Braw = NonGravForces.app.Calibration.basis(bias, ref.t);
    HB = timeSeriesPreparation.crnApply(localArc(ref.t, Braw, "b"), K, OutputEpochs = obs.t).obs;
    Hq = [];
    if options.Quadratic
        Hq = timeSeriesPreparation.crnApply(localArc(ref.t, ref.obs.^2, "q"), K, OutputEpochs = obs.t).obs;
    end
end
% The scale-curve columns H(B_s,q .* a_j), one block of Q per model axis an entry needs.
HaS = cell(1, 3);
if varying
    for j = find(any(mask, 1))
        if isempty(K)
            HaS{j} = NonGravForces.app.Calibration.basis(sc, obs.t) .* ref.obs(at, j);
        else
            Braw = NonGravForces.app.Calibration.basis(sc, ref.t) .* ref.obs(:, j);
            HaS{j} = timeSeriesPreparation.crnApply(localArc(ref.t, Braw, "s"), K, OutputEpochs = obs.t).obs;
        end
    end
end
% A jump whose filtered step is the same on every usable row -- all of them before it (0) or all
% after it (the DC gain, bitwise: one dot product of the same taps with ones) -- is the constant
% column again, or nothing. It cannot be estimated and would make the design singular, so it is
% dropped and said so, never fitted as noise.
if ~isempty(jumps)
    anyRow = obs.weights > 0 & ~obs.filled & any(isfinite(obs.obs), 2) & all(isfinite(HB), 2);
    nB = size(HB, 2) - numel(jumps);
    keep = true(1, numel(jumps));
    for k = 1:numel(jumps)
        v = HB(anyRow, nB + k);
        keep(k) = ~isempty(v) && any(v ~= v(1));
    end
    if ~all(keep)
        warning("NonGravForces:fitCalibration:jumpOutsideFit", ...
            "%d of %d jumps lie outside the fitted rows (their step is constant there) and are not estimated -> nothing to do unless a jump inside the arc was expected", ...
            nnz(~keep), numel(keep));
        jumps = jumps(keep); bias.jumps = jumps;
        HB = HB(:, [true(1, nB), keep]);
    end
end
P = size(HB, 2);
nJ = numel(jumps);
bNames = ["b_" + string(1:P - nJ), "jump_" + string(1:nJ)];
if isempty(K), gridT = obs.t; else, gridT = ref.t; end
[HR, reg] = localRegressors(options.Regressors, gridT, obs.t, K, unit);
Kr = size(HR, 2); rNames = "r_" + reg.names;
axisNames = ["x" "y" "z"];
% The temperature-dependent scale (N-11b): one column H((T - T0) a_j) per model axis an entry needs.
[tCol, T0] = localTemperatureColumn(options, reg, unit);
HT = cell(1, 3);
if tCol > 0
    if varying
        NonGravForces.fail(unit, "unsupportedScaleTemperature", "a scale that depends on time and on temperature is not separable on one arc -> use 'ScaleBasis' ""constant"" with 'ScaleTemperature'");
    end
    Zt = NonGravForces.app.Calibration.interpolateSeries(options.Regressors, tCol, gridT) - T0;
    for j = find(any(mask, 1))
        if isempty(K)
            HT{j} = Zt .* ref.obs(at, j);
        else
            HT{j} = timeSeriesPreparation.crnApply(localArc(ref.t, Zt .* ref.obs(:, j), "t"), K, OutputEpochs = obs.t).obs;
        end
    end
end
coefT = zeros(3); covT = nan(9);

% --- three problems ----------------------------------------------------------------------
N = numel(obs.t); N0 = N;
usable = obs.weights > 0 & ~obs.filled;
scale = S0; covScale = zeros(9); coefB = zeros(3, P); covB = cell(1, 3);
constRow = ones(1, Q); if options.ScaleBasis == "linear", constRow = [1 0]; end
coefS = S0(:) * constRow; covS = cell(1, 9);            % a fixed entry is its constant a priori
Bref = NonGravForces.app.Calibration.basis(sc, sc.tRef);
quad = nan(3, 1);
sigF = nan(3); sigB = nan(3); resRms = nan(1, 3); stopReason = strings(1, 3);
res = nan(N, 3); rowsUsed = false(N, 3); pred = nan(N, 3);
info = struct('residual', [], 'names', {cell(1, 3)}, 'rows', [], 'Qxx', {cell(1, 3)}, ...
    'results', {cell(1, 3)}, 'observability', {cell(1, 3)});
gates = struct('vif', {cell(1, 3)}, 'cond', nan(1, 3), 'maxAbsCorr', nan(1, 3), ...
    'weightFraction', nan(1, 3), 'vifMax', options.VifMax, 'verdict', strings(1, 3));
blockLen = nan(1, 3);
for i = 1:3
    est = find(mask(i, :)); fix = find(~mask(i, :));
    if varying
        Xs = [HaS{est}];
        sNames = reshape(("S_" + axisNames(i) + axisNames(est)) + "_" + string((1:Q)'), 1, []);
    else
        Xs = Ha(:, est);
        sNames = "S_" + axisNames(i) + axisNames(est);
    end
    Xt = zeros(N0, 0); tNames = strings(1, 0);
    if tCol > 0
        Xt = [HT{est}]; tNames = "dS_" + axisNames(i) + axisNames(est) + "_dT";
    end
    if options.Quadratic
        X = [HB, HR, Xs, Xt, Hq(:, i)];
        names = [bNames, rNames, sNames, tNames, "K2_" + axisNames(i)];
    else
        X = [HB, HR, Xs, Xt];
        names = [bNames, rNames, sNames, tNames];
    end
    l = obs.obs(:, i) - Ha(:, fix) * S0(i, fix).';
    ok = usable & isfinite(l) & all(isfinite(X), 2);
    if options.Quadratic
        X(:, end) = X(:, end) - mean(X(ok, end));
    end
    rowsUsed(:, i) = ok;

    % gate: identifiability before the fit (timeSeriesPreparation's metrics)
    yArc = localArc(obs.t, obs.obs(:, i), "y");
    Xg = X; Xg(~ok, :) = NaN;
    rep = timeSeriesPreparation.observability(yArc, Xg, Names = names, MinWeightFraction = options.MinWeightFraction);
    info.observability{i} = rep;
    gates.vif{i} = rep.vif(:, 1).'; gates.cond(i) = rep.cond(1);
    gates.maxAbsCorr(i) = rep.maxAbsCorr(1); gates.weightFraction(i) = rep.weightFraction(1);
    modelVif = rep.vif(P + Kr + (1:numel(est) * Q), 1);
    if numel(est) > 1 && any(~(modelVif <= options.VifMax))
        NonGravForces.fail(unit, "numericCollinear", "axis %s estimates %d model columns and their VIF reaches %.3g > %g: the shear is not separable on this arc -> use 'Scale' ""diag"" or a longer arc", ...
            axisNames(i), numel(est), max(modelVif), options.VifMax);
    end
    gates.verdict(i) = "pass";
    if any(~(modelVif <= options.VifMax)), gates.verdict(i) = "warn"; end
    tVif = rep.vif(P + Kr + numel(est) * Q + (1:numel(tNames)), 1);
    if any(~(tVif <= options.VifMax)), gates.verdict(i) = "warn"; end

    % the solve: lsqvce.fit, kind A, one arc per contiguous run
    [x, Qxx, Cxx, R] = localSolve(X, l, ok, obs.t, options.Noise, options.Robust);
    info.results{i} = R; info.Qxx{i} = Qxx; info.names{i} = names;
    stopReason(i) = string(R.stopReason);
    v = nan(N, 1); v(ok) = l(ok) - X(ok, :) * x;
    % the fitted observation wherever the design is defined, held-out rows included (N-03b)
    rp = all(isfinite(X), 2);
    pred(rp, i) = X(rp, :) * x + Ha(rp, fix) * S0(i, fix).';
    res(:, i) = v;
    resRms(i) = sqrt(mean(v(ok).^2));
    if options.NBoot > 0
        % The block length follows the dependence of the ERRORS, so it is the Politis-White rule on
        % the residual (bootstrap's own automatic rule, read off a two-replicate call), not on the
        % regressors, whose constant bias column has no autocorrelation to measure.
        [~, bl] = timeSeriesPreparation.bootstrap(v, @(r) mean(r, 1, "omitnan"), NBoot = 2, Stage = "fitCalibration.block");
        z = [l, X]; z(~ok, :) = NaN;
        [~, bi] = timeSeriesPreparation.bootstrap(z, @(zz) localBootStat(zz, options.Robust), NBoot = options.NBoot, BlockLength = bl.blockLength, ...
            Seed = options.Seed, Stage = "fitCalibration." + axisNames(i));
        rows = all(isfinite(bi.replicates), 2);
        Cuse = cov(bi.replicates(rows, :));
        blockLen(i) = bi.blockLength;
    else
        Cuse = Cxx;
    end

    coefB(i, :) = x(1:P).';
    covB{i} = Cuse(1:P, 1:P);
    reg.coef(i, :) = reshape(x(P + (1:Kr)), 1, []); reg.cov{i} = Cuse(P + (1:Kr), P + (1:Kr));
    sIdx = P + Kr + (1:numel(est) * Q);
    lin = (est - 1) * 3 + i;                         % positions of S(i, est) in S(:)
    % J maps the curve coefficients of row i to S(i, est) at tRef; for a constant scale it is the
    % identity and every line below reduces to what it was before N-03a.
    J = kron(eye(numel(est)), Bref);
    % x(sIdx) is forced to a column: x is a SCALAR when the bias is "constant" and the axis estimates
    % no scale (S_yy held), and a scalar indexed by an empty index is 1x0, not 0x1 (found by N-03b).
    scale(i, est) = (J * reshape(x(sIdx), [], 1)).';
    covScale(lin, lin) = J * Cuse(sIdx, sIdx) * J.';
    sigF(i, est) = sqrt(diag(J * Cxx(sIdx, sIdx) * J.')).';
    if options.NBoot > 0, sigB(i, est) = sqrt(diag(J * Cuse(sIdx, sIdx) * J.')).'; end
    if varying
        for e = 1:numel(est)
            q = sIdx((e - 1) * Q + (1:Q));
            coefS(lin(e), :) = x(q).';
            covS{lin(e)} = Cuse(q, q);
        end
    end
    if tCol > 0
        tIdx = P + Kr + numel(est) * Q + (1:numel(est));
        coefT(i, est) = reshape(x(tIdx), 1, []);
        covT(lin, lin) = Cuse(tIdx, tIdx);
    end
    if options.Quadratic, quad(i) = x(end); end
end

% --- the EIV bound per model axis --------------------------------------------------------
lambda = nan(1, 3);
for j = 1:3
    ok = rowsUsed(:, j);
    vs = var(Ha(ok, j));
    lambda(j) = vs / (vs + var(res(ok, j)) / scale(j, j)^2);
end

bias.coef = coefB; bias.cov = covB;
uncertainty = struct('sigmaScaleFormal', sigF, 'sigmaScaleBootstrap', sigB, 'inflation', sigB ./ sigF, ...
    'residualRms', resRms, 'eiv', struct('lambda', lambda, 'scaleBiasBound', 1 - lambda, ...
    'note', "lower bound on the attenuation factor, from the residual variance as an upper bound on the model error; a bound, not a correction"), ...
    'noise', localNoiseName(options.Noise), 'nBoot', options.NBoot, 'blockLength', blockLen, 'stopReason', stopReason);
prov = NonGravForces.provenance("app.fitCalibration", 'Options', localHashable(options), ...
    'Note', "obs " + obs.product + ", ref " + ref.product + ", kernel " + localKernelName(K));
sc.coef = coefS; sc.cov = covS;
if ~varying, sc = struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'tRef', int64(0), 'coef', zeros(9, 0), 'cov', {cell(1, 9)}); end
if Kr == 0, reg = []; end
st = localScaleTemperature(options, tCol, T0, coefT, covT);
cal = NonGravForces.app.Calibration(struct('scale', scale, 'scaleCurve', sc, 'covScale', covScale, 'bias', bias, 'regressors', reg, ...
    'scaleTemperature', st, 'quadratic', quad, 'uncertainty', uncertainty, 'gates', gates, 'convention', "observation", 'provenance', prov));
info.residual = obs.derive(obs.t, res, "NonGravForces.app.fitCalibration", struct('basis', options.BiasBasis));
info.rows = rowsUsed;
info.prediction = obs.derive(obs.t, pred, "NonGravForces.app.fitCalibration.prediction", struct('basis', options.BiasBasis));
end

% =========================================================================================
function [tCol, T0] = localTemperatureColumn(options, reg, unit)
%LOCALTEMPERATURECOLUMN The 'Regressors' column 'ScaleTemperature' names (0: none) and its raw-grid mean.
tCol = 0; T0 = NaN;
if options.ScaleTemperature == "", return, end
[found, tCol] = ismember(options.ScaleTemperature, reg.names);
if ~found
    NonGravForces.fail(unit, "inputScaleTemperature", "'ScaleTemperature' names %s, which is not a column of 'Regressors' (%s) -> pass the temperature series in 'Regressors' too", ...
        options.ScaleTemperature, strjoin(reg.names, ", "));
end
T0 = reg.mean(tCol);
end

function st = localScaleTemperature(options, tCol, T0, coefT, covT)
%LOCALSCALETEMPERATURE The Calibration's scaleTemperature struct ("" name when not estimated).
st = struct('name', "", 'units', "", 'reference', NaN, 'coef', zeros(3), 'cov', nan(9));
if tCol == 0, return, end
units = "-";
R = options.Regressors(1);
if numel(R.units) >= tCol, units = string(R.units(tCol)); end
st = struct('name', options.ScaleTemperature, 'units', units, 'reference', T0, 'coef', coefT, 'cov', covT);
end

% =========================================================================================
function [HR, reg] = localRegressors(R, gridT, outT, K, unit)
%LOCALREGRESSORS The measured series on the raw grid, demeaned, filtered like the bias (N-11a).
reg = struct('names', strings(1, 0), 'units', strings(1, 0), 'mean', zeros(1, 0), 'coef', zeros(3, 0), 'cov', {cell(1, 3)});
HR = zeros(numel(outT), 0);
if isempty(R), return, end
if ~isa(R, 'timeSeriesPreparation.Arc') || ~isscalar(R)
    NonGravForces.fail(unit, "inputRegressors", "'Regressors' must be a timeSeriesPreparation.Arc of measured series -> got a %s", class(R));
end
Kr = numel(R.names);
Z = NonGravForces.app.Calibration.interpolateSeries(R, 1:Kr, gridT);
for k = 1:Kr
    z = Z(isfinite(Z(:, k)), k);
    if isempty(z) || all(z == z(1))
        NonGravForces.fail(unit, "inputRegressors", "regressor %s has no finite value or no variation over the arc -> remove it from 'Regressors'", R.names(k));
    end
end
m = mean(Z, 1, "omitnan");
Zc = Z - m;
if isempty(K)
    HR = Zc;
else
    HR = timeSeriesPreparation.crnApply(localArc(gridT, Zc, "r"), K, OutputEpochs = outT).obs;
end
reg.names = reshape(string(R.names), 1, []); reg.units = reshape(string(R.units), 1, []); reg.mean = m;
reg.coef = zeros(3, Kr);
end

% =========================================================================================
function [cal, info] = localMultiArc(obs, ref, options, unit)
%LOCALMULTIARC Several arcs: N independent fits, or one common scale with per-day local biases (N-09).
N = numel(obs);
if isscalar(ref), refs = repmat(ref, 1, N);
elseif numel(ref) == N, refs = ref;
else
    NonGravForces.fail(unit, "inputArc", "ref must be one arc covering every observation arc, or one per arc -> got %d for %d", numel(ref), N);
end
for k = 2:N
    if obs(k).t(1) <= obs(k - 1).t(end)
        NonGravForces.fail(unit, "inputArc", "the observation arcs must be in time order and must not overlap -> arc %d starts before arc %d ends", k, k - 1);
    end
end
if options.ScaleMode == "perArc"
    c = cell(1, N); infos = cell(1, N);
    for k = 1:N
        o = rmfield(options, 'ScaleMode'); o.Seed = options.Seed + k - 1;
        if numel(options.Regressors) == N, o.Regressors = options.Regressors(k); end
        nv = namedargs2cell(o);
        [c{k}, infos{k}] = NonGravForces.app.fitCalibration(obs(k), refs(k), nv{:});
    end
    cal = [c{:}];
    info = struct('mode', "perArc", 'arcs', {infos});
    return
end
if ~ismember(options.BiasBasis, ["constant" "linear"]) || ~isempty(options.Jumps) || options.ScaleBasis ~= "constant"
    NonGravForces.fail(unit, "unsupportedMultiArc", "a common scale over several arcs takes a ""constant"" or ""linear"" bias per day, a constant scale and no 'Jumps' -> use 'ScaleMode' ""perArc"" for anything else");
end
mask = localMask(options.Scale, unit);
if options.FixCrossTrack, mask(2, 2) = false; end
S0 = options.ScaleApriori;
if ~all(isfinite(S0(:)))
    NonGravForces.fail(unit, "inputScale", "'ScaleApriori' must be finite -> it supplies every entry not estimated");
end
K = localKernel(options.Kernel, unit);
if options.NBoot > 0 && N < 10
    warning("NonGravForces:fitCalibration:fewArcs", "a bootstrap over %d days estimates the common scale's spread from %d independent draws -> read its sigma as indicative, or use 10 days or more", N, N);
end

% --- per arc: the filtered columns at the observation epochs (as the single-arc path) ----------
if ~(isempty(options.Regressors) || isscalar(options.Regressors) || numel(options.Regressors) == N)
    NonGravForces.fail(unit, "inputRegressors", "'Regressors' must be one Arc covering every arc or one per arc -> got %d for %d", numel(options.Regressors), N);
end
arcs = struct('t', cell(1, N), 'y', [], 'usable', [], 'Ha', [], 'HB', [], 'Hq', [], 'bias', [], 'HR', [], 'reg', []);
for k = 1:N
    o = obs(k); r = refs(k);
    b = struct('basis', options.BiasBasis, 't0', o.t(1), 'knots', int64([]), 'jumps', zeros(1, 0, 'int64'), 'coef', [], 'cov', {{[], [], []}});
    if isempty(K)
        [found, at] = ismember(o.t, r.t);
        if ~all(found)
            NonGravForces.fail(unit, "inputGrid", "with 'Kernel' [] every observation epoch must be a model epoch -> arc %d: %d of %d are not", k, nnz(~found), numel(found));
        end
        Ha = r.obs(at, :); HB = NonGravForces.app.Calibration.basis(b, o.t); Hq = r.obs(at, :).^2;
    else
        Ha = timeSeriesPreparation.crnApply(r, K, OutputEpochs = o.t).obs;
        HB = timeSeriesPreparation.crnApply(localArc(r.t, NonGravForces.app.Calibration.basis(b, r.t), "b"), K, OutputEpochs = o.t).obs;
        Hq = [];
        if options.Quadratic
            Hq = timeSeriesPreparation.crnApply(localArc(r.t, r.obs.^2, "q"), K, OutputEpochs = o.t).obs;
        end
    end
    if numel(options.Regressors) > 1, Rk = options.Regressors(k); else, Rk = options.Regressors; end
    if isempty(K), gridT = o.t; else, gridT = r.t; end
    [HR, regk] = localRegressors(Rk, gridT, o.t, K, unit);
    arcs(k) = struct('t', o.t, 'y', o.obs, 'usable', o.weights > 0 & ~o.filled, 'Ha', Ha, 'HB', HB, 'Hq', Hq, 'bias', b, 'HR', HR, 'reg', regk);
end
% The temperature-dependent scale (N-11b): ONE reference T0 over every arc's raw grid, so that S1 sees
% the day-to-day temperature differences, and S1 common like the scale.
[tCol, ~] = localTemperatureColumn(options, arcs(1).reg, unit);
HT = cell(1, N); T0 = NaN;
if tCol > 0
    Zt = cell(1, N); grids = cell(1, N);
    for k = 1:N
        if numel(options.Regressors) > 1, Rk = options.Regressors(k); else, Rk = options.Regressors; end
        if isempty(K), grids{k} = obs(k).t; else, grids{k} = refs(k).t; end
        Zt{k} = NonGravForces.app.Calibration.interpolateSeries(Rk, tCol, grids{k});
    end
    zAll = cell2mat(Zt(:));
    T0 = mean(zAll(isfinite(zAll)));
    for k = 1:N
        HT{k} = nan(numel(obs(k).t), 3);
        for j = find(any(mask, 1))
            if isempty(K)
                [~, at] = ismember(obs(k).t, refs(k).t);
                HT{k}(:, j) = (Zt{k} - T0) .* refs(k).obs(at, j);
            else
                HT{k}(:, j) = timeSeriesPreparation.crnApply(localArc(refs(k).t, (Zt{k} - T0) .* refs(k).obs(:, j), "t"), K, OutputEpochs = obs(k).t).obs;
            end
        end
    end
end
coefT = zeros(3); covT = nan(9);
Kr = size(arcs(1).HR, 2); rNames = "r_" + arcs(1).reg.names;
regCoef = zeros(3, Kr); regCov = cell(1, 3);
P = size(arcs(1).HB, 2);
axisNames = ["x" "y" "z"];
bNames = "b_" + string(1:P);

% --- three problems, each one joint system -----------------------------------------------------
scale = S0; covScale = zeros(9); coefB = cell(1, N); covB = cell(1, N);
for k = 1:N, coefB{k} = zeros(3, P); covB{k} = {[], [], []}; end
quad = nan(3, 1); sigF = nan(3); sigB = nan(3); resRms = nan(1, 3); stopReason = strings(1, 3);
gates = struct('vif', {cell(1, 3)}, 'cond', nan(1, 3), 'maxAbsCorr', nan(1, 3), ...
    'weightFraction', nan(1, 3), 'vifMax', options.VifMax, 'verdict', strings(1, 3));
res = cell(1, N); rowsUsed = cell(1, N);
for k = 1:N, res{k} = nan(numel(arcs(k).t), 3); rowsUsed{k} = false(numel(arcs(k).t), 3); end
names = cell(1, 3); results = cell(1, 3); Qxx = cell(1, 3);
lambda = nan(1, 3);
for i = 1:3
    est = find(mask(i, :)); fix = find(~mask(i, :));
    sNames = "S_" + axisNames(i) + axisNames(est);
    if options.Quadratic, sNames = [sNames, "K2_" + axisNames(i)]; end %#ok<AGROW>
    sNames = [sNames, rNames]; %#ok<AGROW>
    if tCol > 0, sNames = [sNames, "dS_" + axisNames(i) + axisNames(est) + "_dT"]; end %#ok<AGROW>
    names{i} = [bNames, sNames];
    % the per-arc rows of this axis: common columns A, local columns B, observation l
    A = cell(1, N); B = cell(1, N); L = cell(1, N); ok = cell(1, N); tS = cell(1, N);
    for k = 1:N
        a = arcs(k);
        A{k} = a.Ha(:, est);
        if options.Quadratic, A{k} = [A{k}, a.Hq(:, i)]; end
        A{k} = [A{k}, a.HR];
        if tCol > 0, A{k} = [A{k}, HT{k}(:, est)]; end
        B{k} = a.HB;
        L{k} = a.y(:, i) - a.Ha(:, fix) * S0(i, fix).';
        ok{k} = a.usable & isfinite(L{k}) & all(isfinite(A{k}), 2) & all(isfinite(B{k}), 2);
        tS{k} = double(a.t - arcs(1).t(1)) / 1e9;
        rowsUsed{k}(:, i) = ok{k};
    end
    if options.Quadratic                                   % K2 demeaned over every fitted row
        qc = numel(est) + 1;
        q = cell2mat(cellfun(@(Ak, okk) Ak(okk, qc), A, ok, 'UniformOutput', false)');
        for k = 1:N, A{k}(:, qc) = A{k}(:, qc) - mean(q); end
    end

    % gate: the common columns with each day's bias basis projected out
    Xp = []; lp = []; tp = int64([]);
    for k = 1:N
        r = ok{k};
        Bk = B{k}(r, :);
        Xp = [Xp; A{k}(r, :) - Bk * (Bk \ A{k}(r, :))]; %#ok<AGROW>
        lp = [lp; L{k}(r) - Bk * (Bk \ L{k}(r))]; %#ok<AGROW>
        tp = [tp; arcs(k).t(r)]; %#ok<AGROW>
    end
    rep = timeSeriesPreparation.observability(localArc(tp, lp, "y"), Xp, Names = sNames, MinWeightFraction = options.MinWeightFraction);
    gates.vif{i} = rep.vif(:, 1).'; gates.cond(i) = rep.cond(1);
    gates.maxAbsCorr(i) = rep.maxAbsCorr(1); gates.weightFraction(i) = rep.weightFraction(1);
    modelVif = rep.vif(1:numel(est), 1);
    if numel(est) > 1 && any(~(modelVif <= options.VifMax))
        NonGravForces.fail(unit, "numericCollinear", "axis %s estimates %d model columns and their VIF reaches %.3g > %g once each day's bias is eliminated -> use 'Scale' ""diag"" or more days", ...
            axisNames(i), numel(est), max(modelVif), options.VifMax);
    end
    gates.verdict(i) = "pass";
    if any(~(modelVif <= options.VifMax)), gates.verdict(i) = "warn"; end
    if tCol > 0 && any(~(rep.vif(end - numel(est) + 1:end, 1) <= options.VifMax)), gates.verdict(i) = "warn"; end

    % the joint solve, and the bootstrap over days on its per-day reduced contributions
    [x, y, Cxx, Cyy, Rr, Qc, sc] = localJoint(A, B, L, ok, tS, 1:N, options.Noise, options.Robust);
    results{i} = Rr; Qxx{i} = Qc; stopReason(i) = string(Rr.stopReason);
    for k = 1:N
        v = nan(numel(L{k}), 1);
        v(ok{k}) = L{k}(ok{k}) - A{k}(ok{k}, :) * x - B{k}(ok{k}, :) * y(:, k);
        res{k}(:, i) = v;
        coefB{k}(i, :) = y(:, k).';
        covB{k}{i} = Cyy{k};
    end
    allv = cell2mat(cellfun(@(r) r(:, i), res, 'UniformOutput', false)');
    resRms(i) = sqrt(mean(allv(isfinite(allv)).^2));
    if options.NBoot > 0 && (~isempty(est) || options.Quadratic || Kr > 0)
        % A replicate is the weighted sum of the days' reduced normal equations (lsqvce >= 6.35,
        % R.contributions / solveContributions): one data pass, then O(N p^2) per replicate. Under the
        % default white noise and no robust weights it IS the refit of the drawn days; an estimated
        % noise model or robust weights stay at the full fit's values (lsqvce: fixed weights).
        Cn = Rr.contributions();
        [~, bi] = timeSeriesPreparation.bootstrap((1:N)', ...
            @(ii) localReplicate(Cn, ii, N, sc).', ...
            NBoot = options.NBoot, BlockLength = 1, Seed = options.Seed, Stage = "fitCalibration.days." + axisNames(i));
        rows = all(isfinite(bi.replicates), 2);
        Cuse = cov(bi.replicates(rows, :));
    else
        Cuse = Cxx;
    end
    sIdx = 1:numel(est);
    scale(i, est) = x(sIdx).';
    lin = (est - 1) * 3 + i;
    covScale(lin, lin) = Cuse(sIdx, sIdx);
    sigF(i, est) = sqrt(diag(Cxx(sIdx, sIdx))).';
    if options.NBoot > 0, sigB(i, est) = sqrt(diag(Cuse(sIdx, sIdx))).'; end
    if options.Quadratic, quad(i) = x(numel(est) + 1); end
    rIdx = numel(est) + options.Quadratic + (1:Kr);
    regCoef(i, :) = reshape(x(rIdx), 1, []); regCov{i} = Cuse(rIdx, rIdx);
    if tCol > 0
        tIdx = numel(est) + options.Quadratic + Kr + (1:numel(est));
        coefT(i, est) = reshape(x(tIdx), 1, []);
        covT(lin, lin) = Cuse(tIdx, tIdx);
    end
end

% --- the EIV bound per model axis, over every fitted row ----------------------------------------
for j = 1:3
    h = []; v = [];
    for k = 1:N
        h = [h; arcs(k).Ha(rowsUsed{k}(:, j), j)]; %#ok<AGROW>
        v = [v; res{k}(rowsUsed{k}(:, j), j)]; %#ok<AGROW>
    end
    lambda(j) = var(h) / (var(h) + var(v) / scale(j, j)^2);
end
uncertainty = struct('sigmaScaleFormal', sigF, 'sigmaScaleBootstrap', sigB, 'inflation', sigB ./ sigF, ...
    'residualRms', resRms, 'eiv', struct('lambda', lambda, 'scaleBiasBound', 1 - lambda, ...
    'note', "lower bound on the attenuation factor, from the residual variance as an upper bound on the model error; a bound, not a correction"), ...
    'noise', localNoiseName(options.Noise), 'nBoot', options.NBoot, 'blockLength', nan(1, 3), 'stopReason', stopReason, ...
    'mode', "common", 'nArcs', N, 'bootstrap', "days", ...
    'biasCovariance', "formal (lsqvce local covariance), not bootstrap");
c = cell(1, N); infoRes = cell(1, N);
for k = 1:N
    b = arcs(k).bias; b.coef = coefB{k}; b.cov = covB{k};
    regk = [];
    if Kr > 0, regk = arcs(k).reg; regk.coef = regCoef; regk.cov = regCov; end
    prov = NonGravForces.provenance("app.fitCalibration", 'Options', localHashable(options), ...
        'Note', sprintf("arc %d of %d, common scale; obs %s, ref %s, kernel %s", k, N, obs(k).product, refs(k).product, localKernelName(K)));
    c{k} = NonGravForces.app.Calibration(struct('scale', scale, 'covScale', covScale, 'bias', b, 'regressors', regk, ...
        'scaleTemperature', localScaleTemperature(options, tCol, T0, coefT, covT), 'quadratic', quad, 'uncertainty', uncertainty, 'gates', gates, 'convention', "observation", 'provenance', prov));
    infoRes{k} = obs(k).derive(obs(k).t, res{k}, "NonGravForces.app.fitCalibration", struct('basis', options.BiasBasis, 'arc', k));
end
cal = [c{:}];
info = struct('mode', "common", 'residual', {infoRes}, 'rows', {rowsUsed}, 'names', {names}, ...
    'Qxx', {Qxx}, 'results', {results});
end

function x = localReplicate(Cn, ii, N, sc)
%LOCALREPLICATE One bootstrap replicate: the drawn days' multiplicities on the full fit's contributions.
c = accumarray(ii(:), 1, [N 1]).';
out = lsqvce.solveContributions(Cn, c);
x = out.x(:) * sc.sl ./ sc.sA(:);
end

function [x, y, Cxx, Cyy, R, Qc, sc] = localJoint(A, B, L, ok, tS, days, noise, robust)
%LOCALJOINT The common columns with each day's bias as lsqvce 'day' locals: one kind-A system.
%   days lists the arcs to use, in order; they are relabelled 1..numel(days). The day bootstrap
%   no longer refits through here (localReplicate), so every caller passes 1:N.
nd = numel(days);
if size(A{1}, 2) == 0
    % No common column (a held scale and no K2): the days decouple into independent bias fits,
    % each through the single-arc solver; lsqvce takes no system without a common parameter.
    x = zeros(0, 1);
    if nargout < 2, return, end
    P = size(B{1}, 2); y = zeros(P, nd); Cyy = cell(1, nd); Rs = cell(1, nd);
    for j = 1:nd
        k = days(j);
        [xk, ~, Ck, Rs{j}] = localSolve(B{k}, L{k}, ok{k}, int64(round(tS{k} * 1e9)), noise, robust);
        y(:, j) = xk; Cyy{j} = Ck;
    end
    Cxx = zeros(0); Qc = zeros(0); R = Rs{end}; sc = struct('sA', zeros(1, 0), 'sl', 1);
    return
end
Ac = cell2mat(cellfun(@(Ak, okk) Ak(okk, :), A(days), ok(days), 'UniformOutput', false)');
Bc = cell2mat(cellfun(@(Bk, okk) Bk(okk, :), B(days), ok(days), 'UniformOutput', false)');
Lc = cell2mat(cellfun(@(Lk, okk) Lk(okk), L(days), ok(days), 'UniformOutput', false)');
sA = sqrt(mean(Ac.^2, 1)); sA(sA == 0) = 1;
sB = sqrt(mean(Bc.^2, 1)); sB(sB == 0) = 1;
sl = sqrt(mean(Lc.^2)); if sl == 0, sl = 1; end
P = size(Bc, 2);
D = cell(1, nd);
for j = 1:nd
    k = days(j);
    idx = find(ok{k});
    cut = [0; find(diff(idx) > 1); numel(idx)];
    nArc = numel(cut) - 1;
    Dj = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, 'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
    for a = 1:nArc
        r = idx(cut(a) + 1:cut(a + 1));
        % this day's local block only: lsqvce places it in the system's B by the locals' names
        Dj.A{a} = A{k}(r, :) ./ sA; Dj.B{a} = B{k}(r, :) ./ sB; Dj.l{a} = L{k}(r) / sl;
        Dj.t{a} = tS{k}(r);                                % lags matter inside an arc only
    end
    D{j} = Dj;
end
p = size(Ac, 2);
R = lsqvce.fit(localTask(p, P, nd, @(day, ctx) localBuildDay(ctx, D{day}, P), noise), localFitOptions(noise, robust));
x = R.x(:) * sl ./ sA(:);
sc = struct('sA', sA, 'sl', sl);
if nargout < 2, return, end
Cxx = R.sol.Cxx .* (sl^2 ./ (sA(:) * sA(:).'));
Qc = R.acc.commonCofactor() ./ (sA(:) * sA(:).');
rr = R.residuals('Verbose', false);
yl = rr.sys(1).y(:); Cl = rr.sys(1).Cyy;
y = zeros(P, nd); Cyy = cell(1, nd);
for j = 1:nd
    q = (j - 1) * P + (1:P);
    y(:, j) = yl(q) * sl ./ sB(:);
    Cyy{j} = Cl(q, q) .* (sl^2 ./ (sB(:) * sB(:).'));
end
end

% =========================================================================================
function [x, Qxx, Cxx, R] = localSolve(X, l, ok, t, noise, robust)
%LOCALSOLVE One axis through lsqvce.fit, kind A: one system, one day, one arc per run of rows.
s = sqrt(mean(X(ok, :).^2, 1)); s(s == 0) = 1;
sl = sqrt(mean(l(ok).^2)); if sl == 0, sl = 1; end
Xe = X ./ s; le = l / sl;
idx = find(ok);
cut = [0; find(diff(idx) > 1); numel(idx)];
nArc = numel(cut) - 1;
D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, 'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
tSec = double(t - t(1)) / 1e9;
for k = 1:nArc
    r = idx(cut(k) + 1:cut(k + 1));
    D.A{k} = Xe(r, :); D.B{k} = zeros(numel(r), 0); D.l{k} = le(r); D.t{k} = tSec(r);
end
p = size(X, 2);
R = lsqvce.fit(localTask(p, 0, 1, @(day, ctx) localBuildDay(ctx, D, 0), noise), localFitOptions(noise, robust));
% Back to the unscaled parameters. Q is the solve's own cofactor (sol.Cxx = sol.s0sq * Q), read
% from the accumulator rather than divided out of Cxx, which a noise-free arc leaves 0 or NaN.
% The cofactor (X'X)^-1 scales with the COLUMNS only; the covariance also carries the observation
% scale sl^2, because its s0^2 was estimated in units of l / sl. Scaling the cofactor by sl^2 too
% made it 4e-14 too small and a bound built on it 5e6 too tight (CI run 38690853810).
x = R.x(:) * sl ./ s(:);
Qxx = R.acc.commonCofactor() ./ (s(:) * s(:).');
Cxx = R.sol.Cxx .* (sl^2 ./ (s(:) * s(:).'));
end

function task = localTask(p, P, nd, build, noise)
%LOCALTASK The lsqvce 7 task: p commons, P day-scoped locals of one system, one stream "acc" over
% days 1..nd. It compiles to the v6 declaration this unit built by hand until 2026-11-03, which
% lsqvce fits bitwise alike (its test_taskDefinition): commons and locals in the same order, the
% stream's law the same spec, its Estimate set only for the estimated law (v6 set 'estimate' there).
task = lsqvce.TaskDefinition(Days = 1:nd);
task.addCommon("c" + (1:p));
sys = task.addSystem("acc");
if P > 0, sys.addLocal("b" + (1:P), Scope = "day"); end
est = string.empty;
if isa(noise, 'lsqvce.Noise')
    law = noise;
else
    law = lsqvce.Noise.white(1);
    if noise == "estimate", est = "law"; end
end
sys.addStream("acc", Build = build, Noise = law, Estimate = est);
end

function out = localBuildDay(ctx, D, P)
%LOCALBUILDDAY One day's arcs; the day's local block by name (lsqvce places the columns).
out = ctx.newDay();
names = "b" + (1:P);
for a = 1:numel(D.A)
    if P > 0
        out.addArc(D.A{a}, D.l{a}, D.t{a}, Local = {names, D.B{a}});
    else
        out.addArc(D.A{a}, D.l{a}, D.t{a});
    end
end
end

function o = localFitOptions(noise, robust)
%LOCALFITOPTIONS One GLS solve under the declared law, or lsqvce's joint VCE with an AR whitener.
if ~isa(noise, 'lsqvce.Noise') && noise == "estimate"
    o = lsqvce.Options(Estimate = "noise", Pool = "stream", ArcFactors = false, Robust = robust, Verbose = false);
    o.Noise = struct('Whitener', "ar");
else
    o = lsqvce.Options(Robust = robust, Verbose = false);
end
end

function x = localBootStat(z, robust)
%LOCALBOOTSTAT One bootstrap replicate: the same equilibrated solve, through lsqvce.fit (kind B).
z = z(all(isfinite(z), 2), :);
l = z(:, 1); X = z(:, 2:end);
s = sqrt(mean(X.^2, 1)); s(s == 0) = 1;
sl = sqrt(mean(l.^2)); if sl == 0, sl = 1; end
R = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(X ./ s, l / sl), lsqvce.Options(Robust = robust, Verbose = false));
x = (R.x(:) * sl ./ s(:)).';
end

function mask = localMask(scale, unit)
%LOCALMASK The estimated entries of S as a 3x3 logical.
if (ischar(scale) || isstring(scale)) && isscalar(string(scale))
    switch string(scale)
        case "diag", mask = logical(eye(3));
        case "full", mask = true(3);
        otherwise
            NonGravForces.fail(unit, "inputScale", "'Scale' must be ""diag"", ""full"" or a 3x3 logical mask -> got ""%s""", string(scale));
    end
elseif (islogical(scale) || isnumeric(scale)) && isequal(size(scale), [3 3])
    mask = logical(scale);
else
    NonGravForces.fail(unit, "inputScale", "'Scale' must be ""diag"", ""full"" or a 3x3 logical mask -> got a %s %s", mat2str(size(scale)), class(scale));
end
end

function K = localKernel(kernel, unit)
%LOCALKERNEL Resolve 'Kernel': the L1B ACC kernel, a given crnKernel struct, or none.
if isempty(kernel)
    K = [];
elseif (ischar(kernel) || isstring(kernel)) && string(kernel) == "crn"
    K = timeSeriesPreparation.crnKernel("ACC");
elseif isstruct(kernel) && isscalar(kernel) && all(isfield(kernel, {'F', 'fs', 'Nf', 'Nh', 'Tf', 'preset'}))
    K = kernel;
else
    NonGravForces.fail(unit, "inputKernel", "'Kernel' must be ""crn"", a timeSeriesPreparation.crnKernel struct or [] -> got a %s", class(kernel));
end
end

function arc = localArc(t, X, prefix)
%LOCALARC A design arc of named columns, the shape timeSeriesPreparation.Arc builds from.
names = cellstr(prefix + "_" + string(1:size(X, 2)));
u = struct();
for k = 1:numel(names), u.(names{k}) = '-'; end
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", u, "fileName", "DESIGN", "product", "DESIGN");
arc = timeSeriesPreparation.Arc(t, X, aux);
end

function s = localNoiseName(noise)
if isa(noise, 'lsqvce.Noise'), s = "declared " + noise.Law; else, s = string(noise); end
end

function s = localKernelName(K)
if isempty(K), s = "none"; else, s = string(K.preset); end
end

function o = localHashable(options)
%LOCALHASHABLE The options as provenance hashes them: a kernel struct by its preset only.
o = options;
if isstruct(o.Kernel), o.Kernel = "struct " + string(o.Kernel.preset); end
if isa(o.Noise, 'lsqvce.Noise'), o.Noise = "declared " + o.Noise.Law; end
if ~isempty(o.Regressors), o.Regressors = "regressors " + strjoin(unique([o.Regressors.names]), " "); end
end
