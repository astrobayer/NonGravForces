classdef Calibration
%CALIBRATION An accelerometer calibration: scale, time-varying bias, uncertainty, gates and provenance, applicable to an arc.
%
%   cal = NonGravForces.app.Calibration(fit)
%   b   = cal.biasAt(t)                 % 3xK bias at int64 GPS-ns epochs
%   a   = cal.apply(obsArc)             % calibrated accelerations, a timeSeriesPreparation.Arc
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
%   Knots and t0 are int64 GPS ns, exactly as timeSeriesPreparation.Arc stores
%   time, so no epoch is rounded through a double before it is differenced.
%
%   APPLY inverts the observation equation per epoch: a = S \ (y - b(t)). The
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
%   by the attenuation the EIV bound reports. The bias covariance of the
%   inverted calibration ignores the scale's uncertainty and the bias-scale
%   cross-covariance, which this class does not hold; the help says so
%   instead of returning a number that looks complete.
%
% CONSTRUCTOR ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   fit   1x1   struct  -      .scale 3x3; .bias struct (.basis, .t0 int64, .knots int64 1xQ,
%                              .coef 3xP, optional .cov 1x3 cell of PxP); optional .covScale
%                              9x9 (cov of scale(:)), .quadratic 3x1, .uncertainty struct,
%                              .gates struct, .convention "observation"|"correction",
%                              .provenance struct
%
% NAME-VALUE ARGUMENTS
%   None.
%
% PROPERTIES (read-only)
%   scale        3x3  double  -      S
%   covScale     9x9  double  -      covariance of scale(:), NaN where not estimated
%   bias         1x1  struct  m/s^2  .basis .t0 .knots .coef (3xP) .cov (1x3 cell)
%   quadratic    3x1  double  s^2/m  K2 per axis; NaN = not estimated
%   uncertainty  1x1  struct  -      as the estimator reports it (sigma, inflation, EIV bound, ...)
%   gates        1x1  struct  -      observability and preflight verdicts
%   convention   1x1  string  -      "observation" (y = b + S a) or "correction" (a = b + S y)
%   provenance   1x1  struct  -      NonGravForces.provenance of the fit, plus the estimator's
%
% METHODS
%   b = biasAt(obj, t)          t Kx1 int64 GPS ns; b 3xK in m/s^2
%   out = apply(obj, arc)       arc: timeSeriesPreparation.Arc with 3 columns in m/s2; out: the
%                               calibrated arc (derive: same epochs, provenance appended)
%   inv = invert(obj)           the other convention, first-order covariance
%   T = summaryTable(obj)       one row per scale entry and bias coefficient
%   save(obj, file)             versioned MAT file (schemaVersion "1")
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
%   above: the bias is evaluated where the fit used it filtered.
%
% ERRORS
%   NonGravForces:Calibration:inputFit       the fit struct lacks a field, or one has the wrong size
%   NonGravForces:Calibration:inputBasis     unknown basis, knots not increasing or too few
%   NonGravForces:Calibration:inputArc       apply: not a 3-column arc in m/s2
%   NonGravForces:Calibration:unsupportedQuadratic  apply with a non-zero K2
%   NonGravForces:Calibration:dataSchema     load: not a calibration file, or no schemaVersion
%   NonGravForces:Calibration:numericSingular invert/apply: S is singular to working precision
%
% WARNINGS
%   NonGravForces:Calibration:newerSchema    load: written by a newer schema; read field by field
%
% REFERENCES
%   Bettadpur S. (2009) GRACE TN-02 (correction direction, a = b + S f);
%   Woeske F. et al. (2024) Adv. Space Res. 74:4517 (ZARM direction);
%   Fuller W. (1987) Measurement Error Models, sect. 1.1 (why an inverted
%   regression is not the reverse regression); de Boor C. (1978) A Practical
%   Guide to Splines, ch. IX (Cox-de Boor recursion).
%
% EXAMPLE
%   [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Spikes', 0, 'Gap_s', 0);
%   bias = struct('basis', "linear", 't0', obs.t(1), 'knots', int64([]), ...
%                 'coef', [truth.bias0, truth.biasDrift_per_day]);
%   cal = NonGravForces.app.Calibration(struct('scale', truth.scale, 'bias', bias));
%   a = cal.apply(obs);                  % the CRN-filtered model, to rounding
%
% See also NonGravForces.app.calibrate, NonGravForces.app.syntheticPair, NonGravForces.sensor.Electrostatic.
%
% KEYWORDS calibration accelerometer bias scale apply invert spline provenance
% Created 2026-10-14 by Claude (Anthropic) for Matthias Weigelt.
    properties (Constant)
        SchemaVersion = "1"
    end
    properties (SetAccess = private)
        scale (3,3) double = eye(3)
        covScale (9,9) double = nan(9)
        bias (1,1) struct = struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'coef', zeros(3, 1), 'cov', {{[], [], []}})
        quadratic (3,1) double = nan(3, 1)
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
            P = NonGravForces.app.Calibration.numBasis(b.basis, b.knots);
            if ~isequal(size(b.coef), [3 P]) || ~all(isfinite(b.coef(:)))
                NonGravForces.fail(unit, "inputFit", ".bias.coef must be a finite 3x%d for basis %s -> got %s", P, b.basis, mat2str(size(b.coef)));
            end
            if ~isfield(b, 'cov') || isempty(b.cov), b.cov = {[], [], []}; end
            obj.bias = struct('basis', string(b.basis), 't0', int64(b.t0), 'knots', int64(b.knots(:)'), 'coef', b.coef, 'cov', {b.cov});
            obj.scale = fit.scale;
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

        function out = apply(obj, arc)
            %APPLY Calibrated accelerations from a measured arc, in the arc's own epochs.
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
            if obj.convention == "observation"
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
            T = table(names, vals, sig, 'VariableNames', {'parameter', 'value', 'sigma'});
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
        end
    end

    methods (Access = private)
        function requireRegular(obj)
            if rcond(obj.scale) < eps
                NonGravForces.fail("Calibration", "numericSingular", "the scale matrix is singular to working precision (rcond %g) -> the calibration cannot be applied or inverted", rcond(obj.scale));
            end
        end

        function s = toStruct(obj)
            s = struct('schemaVersion', NonGravForces.app.Calibration.SchemaVersion, 'scale', obj.scale, ...
                'covScale', obj.covScale, 'bias', obj.bias, 'quadratic', obj.quadratic, ...
                'uncertainty', obj.uncertainty, 'gates', obj.gates, 'convention', obj.convention, ...
                'provenance', obj.provenance);
        end
    end

    methods (Static, Access = private)
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
