classdef Electrostatic < NonGravForces.sensor.Sensor
%ELECTROSTATIC Electrostatic accelerometer: bias, 3x3 scale, proof-mass offset coupling, FIR filter.
%
%   s = NonGravForces.sensor.Electrostatic(bias, S, dr, kernel)
%   s = NonGravForces.sensor.Electrostatic(bias, S, dr, kernel, 'Rotation', R, 'SampleRate', 10)
%   y = s.observe(dyn);  [J, names] = s.partials(dyn);  asd = s.noiseModel(f)
%
%   WHAT AN ELECTROSTATIC ACCELEROMETER READS. Its proof mass sits at a point
%   that is never exactly the satellite's centre of mass, and the satellite
%   rotates, so the specific force at that point is not the one at the CoM:
%
%       a_sensed = a_ng + omegaDot x dr + omega x (omega x dr) - T dr
%       y        = bias + S R a_sensed                         (then filtered)
%
%   dr is the offset of the proof mass from the CoM in the BODY frame, R takes
%   the body frame to the accelerometer frame, S is the 3x3 scale matrix (its
%   off-diagonal entries are the axis couplings, not a rotation) and bias is in
%   the accelerometer frame. The convention is the one qc.residuals already
%   applies, `y = b + S a`, so a calibration estimated once fits both.
%
%   THE SIGN OF omega IS PART OF THIS CONTRACT, because it is the one place a
%   sign error hides: the centrifugal term contains omega twice and survives any
%   convention, so only the Euler term would be wrong, by a factor of exactly
%   -1. dyn.omega is the angular velocity of the body frame with respect to
%   inertial, resolved in the BODY frame, with the sign fixed by
%
%       skew(omega) = R' Rdot        for R mapping BODY -> INERTIAL
%
%   which is what falls out of differentiating v_inertial = R v_body twice, and
%   is therefore the only sign for which the equation above is the kinematics of
%   the proof mass.
%
%   MEASURED 2026-09-08, so a caller does not have to derive it again:
%
%       omega = -referenceSystem.angularVelocityFromRotationRate( ...
%                    referenceSystem.rotationFromQuaternion(q), Rdot, ...
%                    'Resolution', "reference")
%
%   The MINUS SIGN IS NOT OPTIONAL. That unit's rotation runs navigation -> body
%   where this one runs body -> inertial, and its sense is the frame (passive)
%   one, so it returns exactly the negative of the omega defined above (ratio
%   -1.0 measured, tests/tNgfSensor/ElectrostaticReference). Dropping it leaves
%   the centrifugal term right -- it contains omega twice -- and the Euler term
%   wrong by a factor of -1, which is why this is written down instead of left
%   to be rediscovered.
%
%   THE GRAVITY GRADIENT TERM IS NOT COMPUTED HERE, and its absence is not an
%   approximation this class makes quietly. T = d(g)/d(r) needs a gravity field,
%   and this toolbox models no gravity at all (HANDOVER D.7 C4): one authority
%   per fact, and that authority is elsewhere. A caller who has T passes it in
%   `dyn.gravityGradient` and it enters exactly where the equation above puts
%   it; a caller who does not gets the other three terms and a `info.gradient`
%   of "none" saying so. For GRACE-FO's few-hundred-micrometre offset the term
%   is of order 1e-12 m/s^2 and is usually left out, but that is the caller's
%   judgement to make with a number, not this class's to make by omission.
%
%   THE FILTER IS A KERNEL YOU SUPPLY, NOT ONE THIS CLASS INVENTS. L1B
%   accelerations are the ACT1A stream through the CRN low-pass of Wu et al.
%   (2006), and a model compared against them must go through the same filter or
%   the comparison is between two different bandwidths. The coefficients of that
%   filter are defined in that paper, which has NOT been read for this
%   implementation, so no kernel is constructed here and none is guessed: pass
%   the coefficients as `kernel` and they are applied by convolution, with the
%   (L-1)/2 samples at each end marked invalid in `y.mask` because the kernel
%   there overlaps data that does not exist. Pass [] and nothing is filtered.
%   The partials go through the same kernel, or they would describe a signal
%   the observation is not.
%
% CONSTRUCTOR ARGUMENTS
%   Name     Size  Type    Units    Description
%   ----     ----  ----    -----    -----------
%   bias     3x1   double  m/s^2    Additive bias in the accelerometer frame; default zeros
%   S        3x3   double  -        Scale matrix; a 3x1 is taken as its diagonal; default eye(3)
%   dr       3x1   double  m        Proof-mass offset from the CoM, body frame; default zeros
%   kernel   1xL   double  -        FIR coefficients, L odd, sum ~= 0; default [] (no filter)
%
% NAME-VALUE ARGUMENTS
%   Name         Size  Type    Units       Description
%   ----         ----  ----    -----       -----------
%   Rotation     3x3   double  -           Body -> accelerometer frame; default eye(3)
%   SampleRate   1x1   double  Hz          Nominal rate, carried for noiseModel; default 1
%   NoiseAsd     3x1   double  m/s^2/rtHz  White floor per axis; default [1e-10; 1e-9; 1e-10],
%                                          the GRACE-FO order with the cross-track axis ten
%                                          times worse (Kornfeld et al. 2019)
%   NoiseCorner  1x1   double  Hz          Below it the ASD rises as 1/f; default 5e-3
%   Name         1x1   string  -           Instrument label; default "electrostatic"
%
% METHODS
%   y          = observe(obj, dyn)     dyn: .a 3xK body, .omega 3xK body, .omegaDot 3xK body,
%                                      optional .epoch, optional .gravityGradient 3x3xK.
%                                      y: .value 3xK (accelerometer frame), .mask 3xK logical,
%                                      .frame, .epoch, .info
%   [J, names] = partials(obj, dyn)    d(y)/d(theta), (3K)x15, theta = [bias; S(:); dr].
%                                      Row (k-1)*3+i is axis i of epoch k: axis fastest
%   asd        = noiseModel(obj, f_Hz) 3xN one-sided amplitude spectral density
%   n          = numParameters(obj)    15
%
% ERRORS
%   NonGravForces:Electrostatic:inputShape   an argument or a dyn field has the wrong size
%   NonGravForces:Electrostatic:inputKernel  the kernel is empty-summed, of even length, or not a vector
%   NonGravForces:Electrostatic:inputDyn     dyn lacks a, omega or omegaDot, or their epoch counts differ
%   NonGravForces:Electrostatic:numericNaN   an input contains NaN or Inf
%
% REFERENCES
%   Montenbruck O., Gill E. (2000) Satellite Orbits, sect. 3.7 (specific force at
%   a point offset from the centre of mass); Frommknecht B. (2007) PhD, TU
%   Munich (GRACE accelerometer calibration model); Klinger B., Mayer-Guerr T.
%   (2016) Adv. Space Res. 58:1597 (bias and scale parametrisation);
%   Wu S.-C. et al. (2006) GRACE ACC filter (the CRN kernel, NOT implemented
%   here); JPL D-56935 GRACE-FO L1 Handbook (frames, proof-mass geometry);
%   Kornfeld R. et al. (2019) J. Spacecraft Rockets 56:931 (noise floor).
%
% EXAMPLE
%   s = NonGravForces.sensor.Electrostatic([1e-6; -3e-6; 2e-7], diag([1.05 0.96 1.01]), [2e-4; 0; -1e-4]);
%   dyn = struct('a', aBody, 'omega', w, 'omegaDot', wdot);
%   y = s.observe(dyn);
%   r = NonGravForces.qc.residuals(measured', y.value', 'Whitening', "highpass", 'Window', 601);
%
% See also NonGravForces.sensor.Sensor, NonGravForces.qc.residuals, NonGravForces.app.calibrate, referenceSystem.angularVelocityFromRotationRate.
%
% KEYWORDS accelerometer electrostatic bias scale factor centre of mass offset calibration
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
    properties (SetAccess = private)
        bias (3,1) double = zeros(3, 1)
        scale (3,3) double = eye(3)
        dr_m (3,1) double = zeros(3, 1)
        rotation (3,3) double = eye(3)
        kernel (1,:) double = []
        noiseAsd (3,1) double = [1e-10; 1e-9; 1e-10]
        noiseCorner (1,1) double = 5e-3
    end
    methods
        function obj = Electrostatic(bias, S, dr, kernel, options)
            arguments
                bias double = zeros(3, 1)
                S double = eye(3)
                dr double = zeros(3, 1)
                kernel double = []
                options.Rotation double = eye(3)
                options.SampleRate (1,1) double {mustBePositive} = 1
                options.NoiseAsd double = [1e-10; 1e-9; 1e-10]
                options.NoiseCorner (1,1) double {mustBePositive} = 5e-3
                options.Name (1,1) string = "electrostatic"
            end
            unit = "Electrostatic";
            b = bias(:);
            if ~isequal(size(b), [3 1])
                NonGravForces.fail(unit, "inputShape", "bias must be 3x1 -> got %s", mat2str(size(bias)));
            end
            if isvector(S) && numel(S) == 3
                S = diag(S(:));                 % a 3x1 is the diagonal, as qc.residuals also accepts
            end
            if ~isequal(size(S), [3 3])
                NonGravForces.fail(unit, "inputShape", "S must be 3x3 or 3x1 (its diagonal) -> got %s", mat2str(size(S)));
            end
            d = dr(:);
            if ~isequal(size(d), [3 1])
                NonGravForces.fail(unit, "inputShape", "dr must be 3x1 -> got %s", mat2str(size(dr)));
            end
            if ~isequal(size(options.Rotation), [3 3])
                NonGravForces.fail(unit, "inputShape", "'Rotation' must be 3x3 -> got %s", mat2str(size(options.Rotation)));
            end
            k = kernel(:)';
            if ~isempty(k)
                if ~isvector(kernel) || mod(numel(k), 2) ~= 1
                    NonGravForces.fail(unit, "inputKernel", "the filter kernel must be a vector of ODD length so it has a centre sample -> got %s", mat2str(size(kernel)));
                end
                if ~all(isfinite(k)) || sum(k) == 0
                    NonGravForces.fail(unit, "inputKernel", "the filter kernel must be finite and must not sum to zero -> a zero-sum kernel removes the signal this sensor measures");
                end
            end
            asd = options.NoiseAsd(:);
            if ~isequal(size(asd), [3 1]) || any(asd < 0)
                NonGravForces.fail(unit, "inputShape", "'NoiseAsd' must be 3x1 and non-negative -> got %s", mat2str(size(options.NoiseAsd)));
            end
            if any(~isfinite([b; S(:); d; options.Rotation(:); asd]))
                NonGravForces.fail(unit, "numericNaN", "a constructor argument contains NaN/Inf -> screen the calibration before building the sensor");
            end
            obj.bias = b;
            obj.scale = S;
            obj.dr_m = d;
            obj.rotation = options.Rotation;
            obj.kernel = k;
            obj.noiseAsd = asd;
            obj.noiseCorner = options.NoiseCorner;
            obj.name = options.Name;
            obj.axes = ["x", "y", "z"];
            obj.sampleRate = options.SampleRate;
        end

        function y = observe(obj, dyn)
            %OBSERVE The accelerometer reading for the given true dynamics.
            [a, w, wd, K, T] = obj.checkDyn(dyn);
            aSensed = obj.sensedAcceleration(a, w, wd, T);
            raw = obj.bias + obj.scale * (obj.rotation * aSensed);
            [val, mask] = obj.applyFilter(raw);
            edge = 0; if ~isempty(obj.kernel), edge = (numel(obj.kernel) - 1) / 2; end
            grad = "none"; if ~isempty(T), grad = "supplied"; end
            info = struct('K', K, 'filtered', ~isempty(obj.kernel), 'edgeSamples', edge, ...
                'gradient', grad, 'provenance', NonGravForces.provenance("sensor.Electrostatic"));
            ep = []; if isfield(dyn, 'epoch'), ep = dyn.epoch; end
            y = struct('value', val, 'mask', mask, 'frame', "accelerometer", 'epoch', ep, 'info', info);
        end

        function [J, names] = partials(obj, dyn)
            %PARTIALS d(reading)/d(bias, scale, proof-mass offset), (3K)x15, axis fastest.
            [a, w, wd, K, T] = obj.checkDyn(dyn);
            aSensed = obj.sensedAcceleration(a, w, wd, T);
            aAf = obj.rotation * aSensed;                      % 3xK, what the scale multiplies
            P = 15;
            J = zeros(3 * K, P);
            for k = 1:K
                rows = (k - 1) * 3 + (1:3);
                J(rows, 1:3) = eye(3);                         % d/d bias
                % d/d S(i,j) = a_af(j) in row i. MATLAB's S(:) runs down the columns, so
                % parameter 3 + (j-1)*3 + i is S(i,j) and its column of J is e_i a_af(j).
                for j = 1:3
                    for i = 1:3
                        J(rows(i), 3 + (j - 1) * 3 + i) = aAf(j, k);
                    end
                end
                % d/d dr: the coupling is linear in dr, so its Jacobian is the operator itself.
                Om = NonGravForces.sensor.Electrostatic.skew(w(:, k));
                Od = NonGravForces.sensor.Electrostatic.skew(wd(:, k));
                Ck = Od + Om * Om;
                if ~isempty(T), Ck = Ck - T(:, :, k); end
                J(rows, 13:15) = obj.scale * obj.rotation * Ck;
            end
            % The observation is filtered, so its derivatives are filtered too: a Jacobian of the
            % unfiltered signal would describe a quantity nobody observes.
            if ~isempty(obj.kernel)
                for p = 1:P
                    col = reshape(J(:, p), 3, K);
                    J(:, p) = reshape(obj.applyFilter(col), 3 * K, 1);
                end
            end
            names = ["bias_x", "bias_y", "bias_z", ...
                "scale_11", "scale_21", "scale_31", "scale_12", "scale_22", "scale_32", ...
                "scale_13", "scale_23", "scale_33", "dr_x", "dr_y", "dr_z"];
        end

        function S = noiseModel(obj, f_Hz)
            %NOISEMODEL One-sided amplitude spectral density, 3xN: white above the corner, 1/f below.
            f = f_Hz(:)';
            if isempty(f) || any(f <= 0) || any(~isfinite(f))
                NonGravForces.fail("Electrostatic", "inputShape", "f_Hz must be positive and finite -> the ASD is not defined at or below zero frequency");
            end
            S = obj.noiseAsd .* sqrt(1 + (obj.noiseCorner ./ f).^2);
        end

        function n = numParameters(~)
            n = 15;
        end
    end

    methods (Access = private)
        function [a, w, wd, K, T] = checkDyn(obj, dyn) %#ok<INUSL>
            %CHECKDYN The three required fields, their shapes, and the optional gradient.
            unit = "Electrostatic";
            if ~isstruct(dyn) || ~isscalar(dyn) || ~all(isfield(dyn, {'a', 'omega', 'omegaDot'}))
                NonGravForces.fail(unit, "inputDyn", "dyn must be a 1x1 struct with fields a, omega and omegaDot (3xK each) -> omega and omegaDot come from referenceSystem.angularVelocityFromRotationRate");
            end
            a = dyn.a; w = dyn.omega; wd = dyn.omegaDot;
            if ~ismatrix(a) || size(a, 1) ~= 3
                NonGravForces.fail(unit, "inputShape", "dyn.a must be 3xK -> got %s", mat2str(size(a)));
            end
            K = size(a, 2);
            if ~isequal(size(w), [3 K]) || ~isequal(size(wd), [3 K])
                NonGravForces.fail(unit, "inputDyn", "dyn.omega and dyn.omegaDot must be 3x%d to match dyn.a -> got %s and %s", K, mat2str(size(w)), mat2str(size(wd)));
            end
            T = [];
            if isfield(dyn, 'gravityGradient') && ~isempty(dyn.gravityGradient)
                T = dyn.gravityGradient;
                if ~isequal(size(T), [3 3 K])
                    NonGravForces.fail(unit, "inputShape", "dyn.gravityGradient must be 3x3x%d -> got %s", K, mat2str(size(T)));
                end
            end
            if any(~isfinite([a(:); w(:); wd(:); T(:)]))
                NonGravForces.fail(unit, "numericNaN", "dyn contains NaN/Inf -> screen the dynamics before observing them");
            end
        end

        function aSensed = sensedAcceleration(obj, a, w, wd, T)
            %SENSEDACCELERATION The specific force at the proof mass, body frame.
            K = size(a, 2);
            d = repmat(obj.dr_m, 1, K);
            coupling = cross(wd, d) + cross(w, cross(w, d));
            if ~isempty(T)
                coupling = coupling - reshape(pagemtimes(T, obj.dr_m), 3, K);
            end
            aSensed = a + coupling;
        end

        function [val, mask] = applyFilter(obj, raw)
            %APPLYFILTER Convolve each axis with the kernel; the ends are marked, never faked.
            K = size(raw, 2);
            mask = true(3, K);
            if isempty(obj.kernel)
                val = raw;
                return
            end
            L = numel(obj.kernel);
            h = (L - 1) / 2;
            val = zeros(3, K);
            for i = 1:3
                c = conv(raw(i, :), obj.kernel, 'same');
                val(i, :) = c;
            end
            if K <= 2 * h
                mask(:) = false;                 % every sample touches an edge
            else
                mask(:, [1:h, K - h + 1:K]) = false;
            end
        end
    end

    methods (Static)
        function S = skew(v)
            %SKEW The 3x3 matrix with S w = v x w.
            S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
        end
    end
end
