classdef Calibration
%CALIBRATION An accelerometer calibration: scale, time-varying bias, uncertainty, gates and provenance, applicable to an arc.
%
%   cal = NonGravForces.app.Calibration(fit)
%   b   = cal.biasAt(t)                 % 3xK bias at int64 GPS-ns epochs
%   S   = cal.scaleAt(t)                % 3x3xK scale at int64 GPS-ns epochs (N-03a)
%   d   = cal.drift()                   % Wald test of a time-varying scale against a constant
%   a   = cal.apply(obsArc)             % calibrated accelerations, a timeSeriesPreparation.Arc
%   a   = cal.apply(obsArc, 'Regressors', R)   % with measured regressor terms (N-11a), and the
%                                              % temperature of a temperature-dependent scale (N-11b)
%   S   = cal.scaleAtTemperature(T)     % 3x3xK scale at temperatures T (N-11b)
%   inv = cal.invert()                  % the same calibration in the correction direction
%   T   = cal.summaryTable()
%   cal.save(file);  cal = NonGravForces.app.Calibration.load(file)
%   B   = NonGravForces.app.Calibration.basis(bias, t)   % the bias basis, shared with the estimator
%
%   WHAT IT HOLDS. The parameters of the OBSERVATION equation the calibration
%   estimates (CONCEPT 6.2), in the direction this toolbox uses everywhere
%   (sensor.Electrostatic, qc.residuals):
%
%       y(t) = b(t) + S a(t)            convention "observation"
%
%   y is what the accelerometer reports, a the true (modelled) acceleration.
%   Published calibrations come in other directions -- GRACE TN-02 writes
%   a = b' + S' y, ZARM writes a = s (y + b) -- and a number quoted without its
%   direction is the commonest way two calibrations disagree by a factor near
%   two in the bias. So the direction is a property, `invert` converts between
%   the two this class knows, and `apply` always returns calibrated a.
%
%   THE BIAS IS A CURVE, not a vector: b(t) = B(t) c with the basis named in
%   bias.basis and its coefficients in bias.coef (3xP). The basis is evaluated
%   by the static method `basis`, and the estimator builds its design with the
%   SAME method, so the curve a calibration reports is the one that was fitted
%   (one authority). Bases:
%     "constant"  P = 1
%     "linear"    P = 2: c(:,1) + c(:,2) * (t - t0) / 1 day
%     "hat"       P = numel(knots): piecewise linear through the knot values
%     "pspline"   P = numel(knots) - 4: cubic B-splines on the FULL knot vector
%                 (Cox-de Boor), so a uniform interior spacing h comes with
%                 three extra knots beyond each end
%   followed by one unit step per entry of bias.jumps (int64, optional): a
%   jump in the bias is a Heaviside column, never a mask (CONCEPT 6.1, row 2).
%   Knots and t0 are int64 GPS ns, exactly as timeSeriesPreparation.Arc stores
%   time, so no epoch is rounded through a double before it is differenced.
%
%   THE SCALE MAY BE A CURVE TOO (N-03a; CONCEPT 6.2 "time-varying s(t) only
%   with its pointwise sigma band and an isDriftSignificant flag"). Then
%   scaleCurve holds S_ij(t) = B_s(t) g_ij on a basis of the same four kinds,
%   evaluated by the same `basis` method, with g_ij in row (j-1)*3+i of
%   scaleCurve.coef (the order of S(:)); an entry the fit held fixed carries
%   its constant value and an empty covariance. `scale` is then S at
%   scaleCurve.tRef, the middle of the fitted span, and must agree with the
%   curve there: one authority, checked at construction. `scaleSigmaAt`
%   gives the pointwise band sqrt(B_s C B_s') and `drift` a Wald test of
%   "the curve is constant": W = (D g)' (D C D')^+ (D g) with D the first
%   differences of g (partition-of-unity bases, where a constant curve is
%   equal coefficients) or the slope ("linear"), chi-square with rank(D C D')
%   degrees of freedom. The test is the band made into one number: it asks
%   whether the band excludes every constant, jointly over the curve.
%
%   MEASURED REGRESSORS (N-11a; CONCEPT 6.2). A calibration may carry terms in
%   measured series z_k (accelerometer temperatures, magnetorquer currents):
%   y = b(t) + sum_k c_k (z_k(t) - m_k) + S a, with the coefficients 3xK in
%   regressors.coef, their means m_k in regressors.mean and their names and
%   units beside them. apply then needs the series ('Regressors', an Arc with
%   those column names) and REFUSES without them: dropping a term the fit
%   estimated would apply a different calibration than the one fitted. The
%   series are interpolated linearly onto the arc's epochs, as the estimator
%   interpolated them onto the raw grid, and never across a gap in the series
%   (interpolateSeries: NaN where the samples around an epoch are more than 10x
%   the 90th-percentile spacing apart).
%
%   A TEMPERATURE-DEPENDENT SCALE (N-11b; CONCEPT 6.2 "an optional
%   temperature-dependent, fully-populated scale matrix"): S(T) = S + S1 (T - T0),
%   with S1 = scaleTemperature.coef (per unit of the series, the estimated entries
%   of S only), T the measured series named scaleTemperature.name and T0 =
%   scaleTemperature.reference, the temperature at which `scale` holds. apply then
%   needs that series in 'Regressors' and refuses without it, for the reason the
%   regressors do. It cannot be combined with a time-varying scale.
%
%   APPLY inverts the observation equation per epoch: a = S \ (y - b(t)), with
%   S = S(t) when the scale is a curve and S = S(T(t)) when it depends on a
%   temperature. The
%   bias is evaluated, not filtered: the estimator fits H*b, and for a bias
%   that varies slower than the 140.7 s CRN window H*b = dcGain * b to within
%   the curvature of b over that window (exactly so for "constant" and
%   "linear", whose curves the symmetric kernel passes unchanged up to the DC
%   gain 1 + 1.1e-7). A K2 (quadratic) term is not inverted: `apply` refuses a
%   calibration that carries one, rather than invert it approximately and
%   say nothing.
%
%   INVERT returns the correction direction, a = b'(t) + S' y with S' = S^-1
%   and b' = -S^-1 b, and propagates the covariance to first order. It is NOT
%   the regression of a on y: inverting an estimate is not estimating the
%   inverse, and with noise in the model (errors in variables) the two differ
%   by the attenuation the EIV bound reports. A time-varying scale is NOT
%   inverted: S(t)^-1 is no curve on the same basis, so `invert` refuses
%   rather than return an approximation it does not name; nor is a temperature-
%   dependent one (S(T)^-1 is not linear in T). The bias covariance of the
%   inverted calibration ignores the scale's uncertainty and the bias-scale
%   cross-covariance, which this class does not hold; the help says so
%   instead of returning a number that looks complete.
%
% CONSTRUCTOR ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   fit   1x1   struct  -      .scale 3x3; .bias struct (.basis, .t0 int64, .knots int64 1xQ,
%                              optional .jumps int64 1xJ, .coef 3xP with P including the J
%                              steps, optional .cov 1x3 cell of PxP); optional .covScale
%                              9x9 (cov of scale(:)), .quadratic 3x1, .uncertainty struct,
%                              .gates struct, .convention "observation"|"correction",
%                              .provenance struct; optional .scaleCurve struct (.basis,
%                              .t0, .knots, .tRef int64, .coef 9xQ, .cov 1x9 cell of QxQ,
%                              [] where not estimated) for a time-varying scale; optional
%                              .regressors struct (.names 1xK, .units 1xK, .mean 1xK, .coef
%                              3xK, .cov 1x3 cell of KxK); optional .scaleTemperature struct
%                              (.name, .units, .reference, .coef 3x3, .cov 9x9)
%
% NAME-VALUE ARGUMENTS
%   None.
%
% PROPERTIES (read-only)
%   scale        3x3  double  -      S; for a time-varying scale, S(scaleCurve.tRef)
%   scaleCurve   1x1  struct  -      .basis ("constant" when the scale is not a curve) .t0
%                                    .knots .tRef .coef (9xQ) .cov (1x9 cell)
%   covScale     9x9  double  -      covariance of scale(:), NaN where not estimated
%   bias         1x1  struct  m/s^2  .basis .t0 .knots .jumps .coef (3xP) .cov (1x3 cell)
%   quadratic    3x1  double  s^2/m  K2 per axis; NaN = not estimated
%   regressors   1x1  struct  -      .names .units .mean (1xK) .coef (3xK, m/s^2 per unit) .cov
%                                    (1x3 cell, KxK); K = 0 without measured regressors
%   scaleTemperature 1x1 struct -    .name ("" without) .units .reference T0 .coef (3x3, dS/dT, 1 per
%                                    unit, 0 where not estimated) .cov (9x9, cov of coef(:), NaN
%                                    where not estimated)
%   uncertainty  1x1  struct  -      as the estimator reports it (sigma, inflation, EIV bound, ...)
%   gates        1x1  struct  -      observability and preflight verdicts
%   convention   1x1  string  -      "observation" (y = b + S a) or "correction" (a = b + S y)
%   provenance   1x1  struct  -      NonGravForces.provenance of the fit, plus the estimator's
%
% METHODS
%   b = biasAt(obj, t)          t Kx1 int64 GPS ns; b 3xK in m/s^2
%   S = scaleAt(obj, t)         3x3xK; the constant S repeated when not time-varying
%   s = scaleSigmaAt(obj, t)    3x3xK pointwise sigma of S(t); NaN where not estimated
%   tf = isTimeVarying(obj)     scaleCurve.basis is not "constant"
%   tf = isTemperatureDependent(obj)   scaleTemperature.name is not ""
%   S = scaleAtTemperature(obj, T)   3x3xK: scale + coef (T - reference), T 1xK in the series' unit
%   r = regressorTermAt(obj, t, R)   3xK: sum_k c_k (z_k(t) - m_k) from the measured series R
%   T = temperatureAt(obj, t, R)     1xK: the series the scale depends on, at int64 epochs t
%   Z = Calibration.interpolateSeries(R, col, t)   columns col of Arc R, linearly at int64 t (NaN
%                               outside); the one rule of the estimator and of apply
%   d = drift(obj, alpha)       .statistic .dof .pValue (3x3, NaN where not tested)
%                               .significant (3x3 logical) .alpha; alpha default 0.05
%   tf = isDriftSignificant(obj, alpha)   drift(obj, alpha).significant
%   out = apply(obj, arc, 'Regressors', R)   R needed when regressors.names is not empty;
%   out = apply(obj, arc)       arc: timeSeriesPreparation.Arc with 3 columns in m/s2; out: the
%                               calibrated arc (derive: same epochs, provenance appended)
%   inv = invert(obj)           the other convention, first-order covariance
%   T = summaryTable(obj)       one row per scale entry and bias coefficient, plus one per
%                               scale-curve coefficient when time-varying and one per estimated
%                               temperature coefficient (dS_ij_dT)
%   save(obj, file)             versioned MAT file (schemaVersion "4"; "1" to "3" files load)
%   cal = Calibration.load(file)
%   B = Calibration.basis(bias, t)   KxP basis matrix at int64 epochs
%
% OUTPUT ARGUMENTS
%   cal   1x1  NonGravForces.app.Calibration
%
% ACCURACY
%   apply/invert are exact algebra to rounding (cond(S) * eps of the reading);
%   biasAt is the basis evaluation, exact to rounding (partition of unity of the
%   B-splines holds to ~1e-15). The approximation of `apply` is the one named
%   above: the bias is evaluated where the fit used it filtered. A scale that
%   varies within the 140.7 s CRN window -- S(t), or S(T(t)) -- is applied to the
%   filtered reading where the fit filtered S a: the difference is at most about
%   |S1| times the temperature change inside one window, relative to the signal
%   (2.2e-4 measured for |S1| = 2e-3 per C and a 0.5 C orbital swing; GRACE-FO's
%   sensor-unit temperatures move ~0.02 C in a DAY, which makes it negligible). The p-value of
%   `drift` is gammainc(W/2, dof/2, 'upper'), exact for the chi-square law;
%   the law itself is the large-sample one of the estimate, and holds for the
%   bootstrap covariance only as well as that covariance is estimated.
%   covScale and the bias covariances are PRECISIONS GIVEN THE FORCE MODEL: they do
%   not contain the model's own error, which is correlated with the signal and is
%   absorbed by the scale. Measured on January 2020: the steady-state to transient
%   thermal model moved S_zz from 0.820 to 1.094 (about 100 bootstrap sigmas), the
%   observed density scale moved S_xx from 0.762 to 0.966, cutting the thruster
%   firings moved GRACE-D's S_zz by 3.6 sigmas (HANDOVER, PR #126, #122, #134). The
%   per-day bias covariance is the formal (white-noise) one of a CRN-filtered, red
%   residual and is optimistic by the same kind of factor as the scale's
%   uncertainty.inflation.
%
% ERRORS
%   NonGravForces:Calibration:inputFit       the fit struct lacks a field, or one has the wrong size
%   NonGravForces:Calibration:inputBasis     unknown basis, knots or jumps not increasing, too few knots,
%                                            or epochs that are not int64 GPS ns
%   NonGravForces:Calibration:inputArc       apply: not a 3-column arc in m/s2
%   NonGravForces:Calibration:unsupportedQuadratic  apply with a non-zero K2
%   NonGravForces:Calibration:dataSchema     load: not a calibration file, or no schemaVersion
%   NonGravForces:Calibration:numericSingular invert/apply: S is singular to working precision
%                                            (at some epoch, for a time-varying scale)
%   NonGravForces:Calibration:unsupportedTimeVarying  invert of a time-varying scale
%   NonGravForces:Calibration:inputRegressors  apply without the regressor series the calibration
%                                            carries (or the temperature its scale depends on), or
%                                            an Arc that lacks one of their columns
%   NonGravForces:Calibration:unsupportedTemperatureScale  invert of a temperature-dependent scale,
%                                            or one combined with a time-varying scale
%
% WARNINGS
%   NonGravForces:Calibration:newerSchema    load: written by a newer schema; read field by field
%   NonGravForces:Calibration:formalDrift    drift: the covariance is the formal white-noise one
%                                            (NBoot 0, Noise "white"), which a CRN-filtered
%                                            residual violates; the test is then far too liberal
%
% REFERENCES
%   Bettadpur S. (2009) GRACE TN-02 (correction direction, a = b + S f);
%   Woeske F. et al. (2024) Adv. Space Res. 74:4517 (ZARM direction);
%   Fuller W. (1987) Measurement Error Models, sect. 1.1 (why an inverted
%   regression is not the reverse regression); de Boor C. (1978) A Practical
%   Guide to Splines, ch. IX (Cox-de Boor recursion); Wald A. (1943) Trans.
%   Amer. Math. Soc. 54:426 (the test of a linear hypothesis on an estimate).
%
% EXAMPLE
%   [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Spikes', 0, 'Gap_s', 0);
%   bias = struct('basis', "linear", 't0', obs.t(1), 'knots', int64([]), ...
%                 'coef', [truth.bias0, truth.biasDrift_per_day]);
%   cal = NonGravForces.app.Calibration(struct('scale', truth.scale, 'bias', bias));
%   a = cal.apply(obs);                  % the CRN-filtered model, to the planted noise
%
% See also NonGravForces.app.calibrate, NonGravForces.app.syntheticPair, NonGravForces.sensor.Electrostatic.
%
% KEYWORDS calibration accelerometer bias scale apply invert spline provenance
% Created 2026-10-14 by Claude (Anthropic) for Matthias Weigelt; temperature-dependent scale 2026-10-20.
    properties (Constant)
        SchemaVersion = "4"
    end
    properties (SetAccess = private)
        scale (3,3) double = eye(3)
        scaleCurve (1,1) struct = struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'tRef', int64(0), 'coef', zeros(9, 0), 'cov', {cell(1, 9)})
        covScale (9,9) double = nan(9)
        bias (1,1) struct = struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'jumps', zeros(1, 0, 'int64'), 'coef', zeros(3, 1), 'cov', {{[], [], []}})
        quadratic (3,1) double = nan(3, 1)
        regressors (1,1) struct = struct('names', strings(1, 0), 'units', strings(1, 0), 'mean', zeros(1, 0), 'coef', zeros(3, 0), 'cov', {cell(1, 3)})
        scaleTemperature (1,1) struct = struct('name', "", 'units', "", 'reference', NaN, 'coef', zeros(3), 'cov', nan(9))
        uncertainty (1,1) struct = struct()
        gates (1,1) struct = struct()
        convention (1,1) string = "observation"
        provenance (1,1) struct = struct()
    end
    methods
        function obj = Calibration(fit)
            arguments
                fit (1,1) struct
            end
            unit = "Calibration";
            if ~all(isfield(fit, {'scale', 'bias'}))
                NonGravForces.fail(unit, "inputFit", "the fit must carry .scale and .bias -> got fields %s", strjoin(string(fieldnames(fit))', ", "));
            end
            if ~isequal(size(fit.scale), [3 3]) || ~all(isfinite(fit.scale(:)))
                NonGravForces.fail(unit, "inputFit", ".scale must be a finite 3x3 -> got %s", mat2str(size(fit.scale)));
            end
            b = fit.bias;
            if ~isstruct(b) || ~all(isfield(b, {'basis', 't0', 'knots', 'coef'}))
                NonGravForces.fail(unit, "inputFit", ".bias must be a struct with basis, t0, knots and coef -> see the help");
            end
            jumps = NonGravForces.app.Calibration.jumpsOf(b);
            P = NonGravForces.app.Calibration.numBasis(b.basis, b.knots) + numel(jumps);
            if ~isequal(size(b.coef), [3 P]) || ~all(isfinite(b.coef(:)))
                NonGravForces.fail(unit, "inputFit", ".bias.coef must be a finite 3x%d for basis %s -> got %s", P, b.basis, mat2str(size(b.coef)));
            end
            if ~isfield(b, 'cov') || isempty(b.cov), b.cov = {[], [], []}; end
            obj.bias = struct('basis', string(b.basis), 't0', int64(b.t0), 'knots', int64(b.knots(:)'), ...
                'jumps', jumps, 'coef', b.coef, 'cov', {b.cov});
            obj.scale = fit.scale;
            if isfield(fit, 'scaleCurve') && ~isempty(fit.scaleCurve)
                obj.scaleCurve = NonGravForces.app.Calibration.checkScaleCurve(fit.scaleCurve, fit.scale);
            end
            if isfield(fit, 'covScale') && ~isempty(fit.covScale)
                if ~isequal(size(fit.covScale), [9 9])
                    NonGravForces.fail(unit, "inputFit", ".covScale must be 9x9 (the covariance of scale(:)) -> got %s", mat2str(size(fit.covScale)));
                end
                obj.covScale = fit.covScale;
            end
            if isfield(fit, 'quadratic') && ~isempty(fit.quadratic)
                if numel(fit.quadratic) ~= 3
                    NonGravForces.fail(unit, "inputFit", ".quadratic must be 3x1 -> got %s", mat2str(size(fit.quadratic)));
                end
                obj.quadratic = fit.quadratic(:);
            end
            if isfield(fit, 'regressors') && ~isempty(fit.regressors)
                obj.regressors = NonGravForces.app.Calibration.checkRegressors(fit.regressors);
            end
            if isfield(fit, 'scaleTemperature') && ~isempty(fit.scaleTemperature)
                obj.scaleTemperature = NonGravForces.app.Calibration.checkScaleTemperature(fit.scaleTemperature);
                if obj.isTemperatureDependent() && obj.isTimeVarying()
                    NonGravForces.fail(unit, "unsupportedTemperatureScale", "a scale that depends on both time and temperature is not delivered -> fit one or the other");
                end
            end
            if isfield(fit, 'uncertainty'), obj.uncertainty = fit.uncertainty; end
            if isfield(fit, 'gates'), obj.gates = fit.gates; end
            if isfield(fit, 'convention')
                if ~ismember(string(fit.convention), ["observation", "correction"])
                    NonGravForces.fail(unit, "inputFit", ".convention must be ""observation"" or ""correction"" -> got ""%s""", string(fit.convention));
                end
                obj.convention = string(fit.convention);
            end
            if isfield(fit, 'provenance') && ~isempty(fit.provenance)
                obj.provenance = fit.provenance;
            else
                obj.provenance = NonGravForces.provenance("app.Calibration");
            end
        end

        function b = biasAt(obj, t)
            %BIASAT The bias curve at int64 GPS-ns epochs, 3xK.
            B = NonGravForces.app.Calibration.basis(obj.bias, t);
            b = obj.bias.coef * B.';
        end

        function out = apply(obj, arc, nv)
            %APPLY Calibrated accelerations from a measured arc, in the arc's own epochs.
            arguments
                obj
                arc
                nv.Regressors = []
            end
            unit = "Calibration";
            if ~isa(arc, 'timeSeriesPreparation.Arc') || ~isscalar(arc) || size(arc.obs, 2) ~= 3
                NonGravForces.fail(unit, "inputArc", "apply takes one timeSeriesPreparation.Arc with three acceleration columns -> got a %s", class(arc));
            end
            if any(arc.units ~= "m/s2")
                NonGravForces.fail(unit, "inputArc", "the arc's units must be m/s2 -> got %s", strjoin(arc.units, ", "));
            end
            if any(obj.quadratic ~= 0 & ~isnan(obj.quadratic))
                NonGravForces.fail(unit, "unsupportedQuadratic", "a calibration with a K2 term is not inverted approximately -> apply the scale and bias yourself, or fit without 'Quadratic'");
            end
            y = arc.obs.';                                     % 3xK
            b = obj.biasAt(arc.t);
            if ~isempty(obj.regressors.names)
                b = b + obj.regressorTermAt(arc.t, nv.Regressors);
            end
            if obj.isTimeVarying() || obj.isTemperatureDependent()
                % One 3x3 solve per epoch, as pages: the curve, or the temperature, gives every epoch
                % its own S.
                if obj.isTemperatureDependent()
                    S3 = obj.scaleAtTemperature(obj.temperatureAt(arc.t, nv.Regressors));
                else
                    S3 = obj.scaleAt(arc.t);
                end
                K = size(y, 2);
                if obj.convention == "observation"
                    obj.requireRegularPages(S3);
                    a = reshape(pagemldivide(S3, reshape(y - b, 3, 1, K)), 3, K);
                else
                    a = b + reshape(pagemtimes(S3, reshape(y, 3, 1, K)), 3, K);
                end
            elseif obj.convention == "observation"
                obj.requireRegular();
                a = obj.scale \ (y - b);
            else
                a = b + obj.scale * y;
            end
            out = arc.derive(arc.t, a.', "NonGravForces.app.Calibration.apply", ...
                struct('convention', obj.convention, 'basis', obj.bias.basis));
        end

        function out = invert(obj)
            %INVERT The same calibration in the other convention; covariance to first order.
            if obj.isTimeVarying()
                NonGravForces.fail("Calibration", "unsupportedTimeVarying", "S(t)^-1 is not a curve on the scale basis, so a time-varying scale has no exact inverse here -> apply it to the arc instead, or refit with 'ScaleBasis' ""constant""");
            end
            if obj.isTemperatureDependent()
                NonGravForces.fail("Calibration", "unsupportedTemperatureScale", "S(T)^-1 is not linear in T, so a temperature-dependent scale has no exact inverse here -> apply it to the arc instead");
            end
            obj.requireRegular();
            Si = inv3(obj.scale);
            fit = struct('scale', Si, 'quadratic', obj.quadratic, 'uncertainty', obj.uncertainty, ...
                'gates', obj.gates, 'provenance', obj.provenance);
            % vec(S^-1) = -(S^-T kron S^-1) vec(dS) to first order.
            J = -kron(Si.', Si);
            fit.covScale = J * obj.covScale * J.';
            nb = obj.bias;
            nb.coef = -Si * obj.bias.coef;
            % Bias covariance: -S^-1 mixes the axes; with only per-axis covariances held, the
            % propagated one is diagonal by axis and ignores the uncertainty of S (see the help).
            if all(cellfun(@(c) ~isempty(c), obj.bias.cov))
                P = size(obj.bias.coef, 2);
                newCov = cell(1, 3);
                for i = 1:3
                    C = zeros(P);
                    for j = 1:3, C = C + Si(i, j)^2 * obj.bias.cov{j}; end
                    newCov{i} = C;
                end
                nb.cov = newCov;
            end
            fit.bias = nb;
            if ~isempty(obj.regressors.names)
                g = obj.regressors;
                g.coef = -Si * g.coef;
                if all(cellfun(@(c) ~isempty(c), g.cov))
                    Kr = numel(g.names); nc = cell(1, 3);
                    for i = 1:3
                        C = zeros(Kr);
                        for j = 1:3, C = C + Si(i, j)^2 * obj.regressors.cov{j}; end
                        nc{i} = C;
                    end
                    g.cov = nc;
                end
                fit.regressors = g;
            end
            if obj.convention == "observation", fit.convention = "correction"; else, fit.convention = "observation"; end
            out = NonGravForces.app.Calibration(fit);
        end

        function T = summaryTable(obj)
            %SUMMARYTABLE One row per scale entry and per bias coefficient: value and sigma.
            [ii, jj] = ndgrid(1:3, 1:3);
            names = compose("S_%d%d", ii(:), jj(:));
            vals = obj.scale(:);
            sig = sqrt(diag(obj.covScale));
            P = size(obj.bias.coef, 2);
            ax = ["x", "y", "z"];
            for i = 1:3
                names = [names; compose("b_%s_%d", ax(i), (1:P)')]; %#ok<AGROW>
                vals = [vals; obj.bias.coef(i, :).']; %#ok<AGROW>
                if isempty(obj.bias.cov{i}), s = nan(P, 1); else, s = sqrt(diag(obj.bias.cov{i})); end
                sig = [sig; s]; %#ok<AGROW>
            end
            for k = 1:numel(obj.regressors.names)
                for i = 1:3
                    names = [names; "r_" + obj.regressors.names(k) + "_" + ax(i)]; %#ok<AGROW>
                    vals = [vals; obj.regressors.coef(i, k)]; %#ok<AGROW>
                    if isempty(obj.regressors.cov{i}), sg = NaN; else, sg = sqrt(obj.regressors.cov{i}(k, k)); end
                    sig = [sig; sg]; %#ok<AGROW>
                end
            end
            if obj.isTimeVarying()
                Q = size(obj.scaleCurve.coef, 2);
                for k = 1:9
                    if isempty(obj.scaleCurve.cov{k}), continue; end
                    names = [names; compose("s_%d%d_%d", ii(k), jj(k), (1:Q)')]; %#ok<AGROW>
                    vals = [vals; obj.scaleCurve.coef(k, :).']; %#ok<AGROW>
                    sig = [sig; sqrt(diag(obj.scaleCurve.cov{k}))]; %#ok<AGROW>
                end
            end
            if obj.isTemperatureDependent()
                st = obj.scaleTemperature;
                cs = sqrt(diag(st.cov));
                for k = find(isfinite(cs))'
                    names = [names; compose("dS_%d%d_dT", ii(k), jj(k))]; %#ok<AGROW>
                    vals = [vals; st.coef(k)]; %#ok<AGROW>
                    sig = [sig; cs(k)]; %#ok<AGROW>
                end
            end
            T = table(names, vals, sig, 'VariableNames', {'parameter', 'value', 'sigma'});
        end

        function r = regressorTermAt(obj, t, R)
            %REGRESSORTERMAT sum_k c_k (z_k(t) - m_k), 3xK, from the measured series R.
            g = obj.regressors;
            if isempty(R) || ~isa(R, 'timeSeriesPreparation.Arc')
                NonGravForces.fail("Calibration", "inputRegressors", "this calibration carries regressor terms in %s -> pass 'Regressors', an Arc with those columns", strjoin(g.names, ", "));
            end
            [found, col] = ismember(g.names, R.names);
            if ~all(found)
                NonGravForces.fail("Calibration", "inputRegressors", "'Regressors' lacks %s -> pass the series the calibration was fitted with", strjoin(g.names(~found), ", "));
            end
            Z = NonGravForces.app.Calibration.interpolateSeries(R, col, t) - g.mean;
            r = g.coef * Z.';
        end

        function tf = isTemperatureDependent(obj)
            %ISTEMPERATUREDEPENDENT Whether the scale depends on a measured temperature (N-11b).
            tf = obj.scaleTemperature.name ~= "";
        end

        function S = scaleAtTemperature(obj, T)
            %SCALEATTEMPERATURE The scale at temperatures T (1xK, the series' unit), 3x3xK.
            T = reshape(T, 1, 1, []);
            S = obj.scale + obj.scaleTemperature.coef .* (T - obj.scaleTemperature.reference);
        end

        function T = temperatureAt(obj, t, R)
            %TEMPERATUREAT The series the scale depends on, interpolated at int64 epochs t (1xK).
            name = obj.scaleTemperature.name;
            if isempty(R) || ~isa(R, 'timeSeriesPreparation.Arc')
                NonGravForces.fail("Calibration", "inputRegressors", "this calibration's scale depends on %s -> pass 'Regressors', an Arc with that column", name);
            end
            [found, col] = ismember(name, R.names);
            if ~found
                NonGravForces.fail("Calibration", "inputRegressors", "'Regressors' lacks %s, the temperature this scale depends on -> pass the series the calibration was fitted with", name);
            end
            T = reshape(NonGravForces.app.Calibration.interpolateSeries(R, col, t), 1, []);
        end

        function tf = isTimeVarying(obj)
            %ISTIMEVARYING Whether the scale is a curve rather than a constant matrix.
            tf = obj.scaleCurve.basis ~= "constant";
        end

        function S = scaleAt(obj, t)
            %SCALEAT The scale at int64 GPS-ns epochs, 3x3xK.
            K = numel(t);
            if ~obj.isTimeVarying()
                S = repmat(obj.scale, 1, 1, K);
                return
            end
            B = NonGravForces.app.Calibration.basis(obj.scaleCurve, t);
            S = reshape(obj.scaleCurve.coef * B.', 3, 3, K);
        end

        function s = scaleSigmaAt(obj, t)
            %SCALESIGMAAT The pointwise sigma of S(t), 3x3xK; NaN where the entry was not estimated.
            K = numel(t);
            if ~obj.isTimeVarying()
                s = repmat(reshape(sqrt(diag(obj.covScale)), 3, 3), 1, 1, K);
                return
            end
            B = NonGravForces.app.Calibration.basis(obj.scaleCurve, t);
            s = nan(9, K);
            for k = 1:9
                C = obj.scaleCurve.cov{k};
                if isempty(C), continue; end
                s(k, :) = sqrt(max(sum((B * C) .* B, 2), 0)).';
            end
            s = reshape(s, 3, 3, K);
        end

        function d = drift(obj, alpha)
            %DRIFT Wald test, per entry, of the hypothesis that the scale curve is constant.
            arguments
                obj (1,1) NonGravForces.app.Calibration
                alpha (1,1) double {mustBePositive, mustBeLessThan(alpha, 1)} = 0.05
            end
            d = struct('statistic', nan(3), 'dof', nan(3), 'pValue', nan(3), 'significant', false(3), 'alpha', alpha);
            if ~obj.isTimeVarying(), return; end
            u = obj.uncertainty;
            if isfield(u, 'nBoot') && isfield(u, 'noise') && u.nBoot == 0 && string(u.noise) == "white"
                % Measured (N-03a): on a CRN-filtered residual the bootstrap sigma is ~26x the formal
                % white one, and the formal test flagged a planted-zero drift at p = 1e-18.
                warning("NonGravForces:Calibration:formalDrift", "the scale covariance is the formal white-noise one, which the filtered residual violates by an order of magnitude -> refit with 'NBoot' > 0 before reading the drift test");
            end
            Q = size(obj.scaleCurve.coef, 2);
            if obj.scaleCurve.basis == "linear"
                D = [0 1];                          % c1 + c2 * days is constant iff c2 = 0
            else
                D = diff(eye(Q));                   % partition of unity: constant iff all equal
            end
            for k = 1:9
                C = obj.scaleCurve.cov{k};
                if isempty(C), continue; end
                g = obj.scaleCurve.coef(k, :).';
                M = D * C * D.';
                r = rank(M);
                if r == 0, continue; end
                W = (D * g).' * pinv(M) * (D * g);
                d.statistic(k) = W; d.dof(k) = r;
                d.pValue(k) = gammainc(W / 2, r / 2, 'upper');
                d.significant(k) = d.pValue(k) < alpha;
            end
        end

        function tf = isDriftSignificant(obj, alpha)
            %ISDRIFTSIGNIFICANT 3x3 logical: entries whose curve the Wald test says is not constant.
            arguments
                obj (1,1) NonGravForces.app.Calibration
                alpha (1,1) double {mustBePositive, mustBeLessThan(alpha, 1)} = 0.05
            end
            tf = obj.drift(alpha).significant;
        end

        function save(obj, file)
            %SAVE A versioned MAT file: every property plus schemaVersion, under one variable.
            arguments
                obj (1,1) NonGravForces.app.Calibration
                file (1,1) string
            end
            calibration = obj.toStruct(); %#ok<NASGU>
            builtin('save', file, 'calibration', '-v7');
        end
    end

    methods (Static)
        function obj = load(file)
            %LOAD Read a calibration saved by save; a newer schema is read field by field and warns.
            arguments
                file (1,1) string
            end
            s = builtin('load', file);
            if ~isfield(s, 'calibration') || ~isfield(s.calibration, 'schemaVersion')
                NonGravForces.fail("Calibration", "dataSchema", "%s holds no calibration with a schemaVersion -> save it with NonGravForces.app.Calibration.save", file);
            end
            c = s.calibration;
            if str2double(c.schemaVersion) > str2double(NonGravForces.app.Calibration.SchemaVersion)
                warning("NonGravForces:Calibration:newerSchema", "%s was written by schema %s, newer than this toolbox's %s; the fields this version knows are read and the rest ignored -> update NonGravForces", ...
                    file, c.schemaVersion, NonGravForces.app.Calibration.SchemaVersion);
            end
            obj = NonGravForces.app.Calibration(rmfield(c, 'schemaVersion'));
        end

        function Z = interpolateSeries(R, col, t)
            %INTERPOLATESERIES Columns col of the Arc R, linearly at int64 epochs t (NaN outside).
            %   One rule for the estimator and for apply: linear, so a housekeeping series sampled every
            %   few seconds is never overshot between its samples -- and never bridged across a GAP: where
            %   the two samples around an epoch lie more than 10 times the 90th percentile of the sample
            %   spacing apart, the value is NaN, so the rows there are not fitted (a 20 min housekeeping
            %   gap, interpolated, moved a temperature-dependent scale by 3.7 %; found by N-11b,
            %   2026-10-20). A percentile, not the median: real series come in BURSTS -- MAG1B in pairs
            %   0.5 ms apart once a second, AHK1B's TSU in 0.1 s bursts 4.7 s apart -- whose median
            %   spacing is the burst's, which made every epoch between bursts a "gap".
            tr = double(R.t - R.t(1)) / 1e9; tq = double(t(:) - R.t(1)) / 1e9;
            Z = zeros(numel(tq), numel(col));
            for k = 1:numel(col)
                z = R.obs(:, col(k)); ok = isfinite(z);
                Z(:, k) = interp1(tr(ok), z(ok), tq, 'linear', NaN);
                if nnz(ok) > 1
                    tk = tr(ok);
                    span = interp1(tk, tk, tq, 'next', NaN) - interp1(tk, tk, tq, 'previous', NaN);
                    d = sort(diff(tk));
                    Z(span > 10 * d(ceil(0.9 * numel(d))), k) = NaN;
                end
            end
        end

        function B = basis(bias, t)
            %BASIS The KxP bias basis at int64 GPS-ns epochs t (the estimator uses this too).
            t = t(:);
            if ~isa(t, 'int64')
                NonGravForces.fail("Calibration", "inputBasis", "epochs must be int64 GPS ns, as timeSeriesPreparation.Arc stores them -> got %s", class(t));
            end
            NonGravForces.app.Calibration.numBasis(bias.basis, bias.knots);   % validates
            days = double(t - int64(bias.t0)) / 86400e9;
            switch string(bias.basis)
                case "constant"
                    B = ones(numel(t), 1);
                case "linear"
                    B = [ones(numel(t), 1), days];
                case "hat"
                    kn = double(int64(bias.knots(:)') - int64(bias.t0)) / 86400e9;
                    B = zeros(numel(t), numel(kn));
                    for p = 1:numel(kn)
                        e = zeros(1, numel(kn)); e(p) = 1;
                        B(:, p) = interp1(kn, e, days, 'linear', 'extrap');
                    end
                case "pspline"
                    kn = double(int64(bias.knots(:)') - int64(bias.t0)) / 86400e9;
                    B = NonGravForces.app.Calibration.coxDeBoor(kn, days, 4);
            end
            % Jumps: one unit step per epoch, 0 before it and 1 from it on (the sample AT the jump
            % already carries the new level), appended after the smooth basis.
            jumps = NonGravForces.app.Calibration.jumpsOf(bias);
            if ~isempty(jumps)
                B = [B, double(t >= jumps)];
            end
        end
    end

    methods (Access = private)
        function requireRegular(obj)
            if rcond(obj.scale) < eps
                NonGravForces.fail("Calibration", "numericSingular", "the scale matrix is singular to working precision (rcond %g) -> the calibration cannot be applied or inverted", rcond(obj.scale));
            end
        end

        function requireRegularPages(~, S3)
            % Every page at once, closed form: |det| over the product of the row norms is 0 for a
            % singular page and 1 for an orthogonal one (Hadamard), with no per-epoch call.
            c = @(i, j) reshape(S3(i, j, :), 1, []);
            dt = c(1,1) .* (c(2,2) .* c(3,3) - c(2,3) .* c(3,2)) - c(1,2) .* (c(2,1) .* c(3,3) - c(2,3) .* c(3,1)) ...
                + c(1,3) .* (c(2,1) .* c(3,2) - c(2,2) .* c(3,1));
            rn = sqrt(sum(S3.^2, 2));
            q = abs(dt) ./ reshape(prod(rn, 1), 1, []);
            if any(~(q >= eps))
                NonGravForces.fail("Calibration", "numericSingular", "the scale curve is singular to working precision at %d of %d epochs -> the calibration cannot be applied there", nnz(~(q >= eps)), numel(q));
            end
        end

        function s = toStruct(obj)
            s = struct('schemaVersion', NonGravForces.app.Calibration.SchemaVersion, 'scale', obj.scale, 'scaleCurve', obj.scaleCurve, 'regressors', obj.regressors, 'scaleTemperature', obj.scaleTemperature, ...
                'covScale', obj.covScale, 'bias', obj.bias, 'quadratic', obj.quadratic, ...
                'uncertainty', obj.uncertainty, 'gates', obj.gates, 'convention', obj.convention, ...
                'provenance', obj.provenance);
        end
    end

    methods (Static, Access = private)
        function g = checkRegressors(g)
            %CHECKREGRESSORS Validate a regressor struct: names, units, means, coefficients, covariances.
            unit = "Calibration";
            if ~isstruct(g) || ~isscalar(g) || ~all(isfield(g, {'names', 'mean', 'coef'}))
                NonGravForces.fail(unit, "inputFit", ".regressors must be a struct with names, mean and coef -> see the help");
            end
            Kr = numel(g.names);
            if ~isequal(size(g.coef), [3 Kr]) || numel(g.mean) ~= Kr || ~all(isfinite([g.coef(:); g.mean(:)]))
                NonGravForces.fail(unit, "inputFit", ".regressors.coef must be a finite 3x%d and .mean 1x%d -> got %s and %s", Kr, Kr, mat2str(size(g.coef)), mat2str(size(g.mean)));
            end
            if ~isfield(g, 'units') || isempty(g.units), g.units = repmat("-", 1, Kr); end
            if ~isfield(g, 'cov') || isempty(g.cov), g.cov = cell(1, 3); end
            g = struct('names', reshape(string(g.names), 1, []), 'units', reshape(string(g.units), 1, []), ...
                'mean', reshape(g.mean, 1, []), 'coef', g.coef, 'cov', {reshape(g.cov, 1, 3)});
        end

        function st = checkScaleTemperature(st)
            %CHECKSCALETEMPERATURE Validate a temperature-dependent scale: name, reference, coef, cov.
            unit = "Calibration";
            if ~isstruct(st) || ~isscalar(st) || ~all(isfield(st, {'name', 'reference', 'coef'}))
                NonGravForces.fail(unit, "inputFit", ".scaleTemperature must be a struct with name, reference and coef -> see the help");
            end
            if string(st.name) == ""
                st = struct('name', "", 'units', "", 'reference', NaN, 'coef', zeros(3), 'cov', nan(9));
                return
            end
            if ~isequal(size(st.coef), [3 3]) || ~all(isfinite(st.coef(:))) || ~isscalar(st.reference) || ~isfinite(st.reference)
                NonGravForces.fail(unit, "inputFit", ".scaleTemperature needs a finite 3x3 coef and a finite scalar reference -> got %s and %s", mat2str(size(st.coef)), mat2str(size(st.reference)));
            end
            if ~isfield(st, 'units') || isempty(st.units), st.units = "-"; end
            if ~isfield(st, 'cov') || isempty(st.cov), st.cov = nan(9); end
            if ~isequal(size(st.cov), [9 9])
                NonGravForces.fail(unit, "inputFit", ".scaleTemperature.cov must be 9x9 (the covariance of coef(:)) -> got %s", mat2str(size(st.cov)));
            end
            st = struct('name', string(st.name), 'units', string(st.units), 'reference', double(st.reference), 'coef', st.coef, 'cov', st.cov);
        end

        function sc = checkScaleCurve(sc, scale)
            %CHECKSCALECURVE Validate a scale curve and its agreement with S at tRef.
            unit = "Calibration";
            if ~isstruct(sc) || ~isscalar(sc) || ~all(isfield(sc, {'basis', 't0', 'knots', 'coef'}))
                NonGravForces.fail(unit, "inputFit", ".scaleCurve must be a struct with basis, t0, knots and coef -> see the help");
            end
            if string(sc.basis) == "constant" && isempty(sc.coef)
                sc = struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'tRef', int64(0), 'coef', zeros(9, 0), 'cov', {cell(1, 9)});
                return
            end
            Q = NonGravForces.app.Calibration.numBasis(sc.basis, sc.knots);
            if ~isequal(size(sc.coef), [9 Q]) || ~all(isfinite(sc.coef(:)))
                NonGravForces.fail(unit, "inputFit", ".scaleCurve.coef must be a finite 9x%d for basis %s -> got %s", Q, string(sc.basis), mat2str(size(sc.coef)));
            end
            if ~isfield(sc, 'cov') || isempty(sc.cov), sc.cov = cell(1, 9); end
            if ~iscell(sc.cov) || numel(sc.cov) ~= 9 || any(cellfun(@(c) ~isempty(c) && ~isequal(size(c), [Q Q]), sc.cov))
                NonGravForces.fail(unit, "inputFit", ".scaleCurve.cov must be a 1x9 cell of %dx%d or [] -> got another shape", Q, Q);
            end
            if ~isfield(sc, 'tRef') || isempty(sc.tRef), sc.tRef = sc.t0; end
            sc = struct('basis', string(sc.basis), 't0', int64(sc.t0), 'knots', int64(sc.knots(:)'), ...
                'tRef', int64(sc.tRef), 'coef', sc.coef, 'cov', {reshape(sc.cov, 1, 9)});
            Sref = reshape(sc.coef * NonGravForces.app.Calibration.basis(sc, sc.tRef).', 3, 3);
            if any(abs(Sref(:) - scale(:)) > 1e-12 * max(1, abs(scale(:))))
                NonGravForces.fail(unit, "inputFit", ".scale must be the curve at .scaleCurve.tRef (one authority) -> they differ by %g", max(abs(Sref(:) - scale(:))));
            end
        end

        function jumps = jumpsOf(bias)
            %JUMPSOF The jump epochs of a bias struct as a strictly increasing int64 row (none: 1x0).
            jumps = zeros(1, 0, 'int64');
            if isfield(bias, 'jumps') && ~isempty(bias.jumps)
                if ~isa(bias.jumps, 'int64') || any(diff(bias.jumps(:)') <= 0)
                    NonGravForces.fail("Calibration", "inputBasis", "bias.jumps must be strictly increasing int64 GPS ns -> sort them and drop duplicates");
                end
                jumps = bias.jumps(:)';
            end
        end

        function P = numBasis(basis, knots)
            %NUMBASIS The number of basis functions, validating the basis and its knots.
            unit = "Calibration";
            kn = int64(knots(:)');
            if any(diff(kn) <= 0)
                NonGravForces.fail(unit, "inputBasis", "knots must be strictly increasing int64 GPS ns -> sort them");
            end
            switch string(basis)
                case "constant", P = 1;
                case "linear",   P = 2;
                case "hat"
                    P = numel(kn);
                    if P < 2, NonGravForces.fail(unit, "inputBasis", "a hat basis needs at least 2 knots -> got %d", P); end
                case "pspline"
                    P = numel(kn) - 4;
                    if P < 4, NonGravForces.fail(unit, "inputBasis", "a cubic B-spline basis needs the full knot vector, at least 8 knots -> got %d", numel(kn)); end
                otherwise
                    NonGravForces.fail(unit, "inputBasis", "unknown bias basis ""%s"" -> use constant, linear, hat or pspline", string(basis));
            end
        end

        function B = coxDeBoor(kn, x, order)
            %COXDEBOOR B-spline basis of the given order on the full knot vector kn, at x.
            x = x(:); m = numel(kn);
            B = double(x >= kn(1:m - 1) & x < kn(2:m));                % order 1
            last = kn(m - order + 1);                                    % right end of the support
            B(x == last, :) = 0; B(x == last, find(kn(1:m - 1) < last, 1, 'last')) = 1;
            for k = 2:order
                Bn = zeros(numel(x), m - k);
                for i = 1:m - k
                    d1 = kn(i + k - 1) - kn(i); d2 = kn(i + k) - kn(i + 1);
                    a = 0; c = 0;
                    if d1 > 0, a = (x - kn(i)) / d1 .* B(:, i); end
                    if d2 > 0, c = (kn(i + k) - x) / d2 .* B(:, i + 1); end
                    Bn(:, i) = a + c;
                end
                B = Bn;
            end
        end
    end
end

function Si = inv3(S)
%INV3 The inverse of a regular 3x3 by a solve, not inv(): one LU, the same rounding as S \ I.
Si = S \ eye(3);
end
