classdef Electrostatic < NonGravForces.sensor.Sensor
%ELECTROSTATIC Electrostatic accelerometer: bias, 3x3 scale, proof-mass offset coupling; L1B CRN filter by default.
%
%   s = NonGravForces.sensor.Electrostatic(bias, S, dr)                       % CRN-filtered, 10 Hz input
%   s = NonGravForces.sensor.Electrostatic(bias, S, dr, "crn", 'OutputRate', 1) % as ACT1B: 1 Hz samples
%   s = NonGravForces.sensor.Electrostatic(bias, S, dr, [])                   % unfiltered, any rate
%   s = NonGravForces.sensor.Electrostatic(bias, S, dr, kernel, 'Rotation', R, 'SampleRate', fs)
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
%   the accelerometer frame. The convention is the one timeSeriesPreparation.qc.residuals
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
%   THE DEFAULT FILTER IS THE L1B CRN KERNEL, AND IT IS NOT BUILT HERE (N-04,
%   owner ruling 2026-10-14). L1B accelerations are the 10 Hz ACT1A stream
%   through the CRN low-pass of the GRACE L1B ATBD (Wu et al. 2006, App. B), and
%   a model compared against them must go through the same filter or the
%   comparison is between two different bandwidths. The kernel's one authority
%   is timeSeriesPreparation.crnKernel("ACC"): 1407 taps on the 10 Hz grid,
%   normalised at f0 with its DC gain of 1 + 1.1e-7 kept. This class takes its
%   taps from there and convolves; it computes none of them.
%     Because the taps are defined on a 10 Hz grid, the CRN default needs 10 Hz
%   dynamics: 'SampleRate' defaults to 10, a kernel whose rate differs from
%   'SampleRate' is refused, and dyn.epoch is REQUIRED and checked to be on that
%   grid, because a 1 Hz arc convolved with 0.1 s taps is a filter of a
%   different bandwidth that nothing would report. 'OutputRate' then samples the
%   filtered series (never averages it), as L1B forms 1 Hz ACT1B from it.
%     Any other FIR is still accepted as a numeric `kernel` and applied at
%   whatever rate the caller supplies; [] filters nothing and accepts any rate.
%   The (L-1)/2 samples at each end are marked invalid in `y.mask`, because the
%   kernel there overlaps data that does not exist -- 703 samples, 70.3 s, for
%   the CRN kernel. The partials go through the same kernel and the same
%   sampling, or they would describe a signal the observation is not.
%
% CONSTRUCTOR ARGUMENTS
%   Name     Size  Type    Units    Description
%   ----     ----  ----    -----    -----------
%   bias     3x1   double  m/s^2    Additive bias in the accelerometer frame; default zeros
%   S        3x3   double  -        Scale matrix; a 3x1 is taken as its diagonal; default eye(3)
%   dr       3x1   double  m        Proof-mass offset from the CoM, body frame; default zeros
%   kernel   1x1   string  -        "crn" (default): timeSeriesPreparation.crnKernel("ACC")
%            1x1   struct  -        a timeSeriesPreparation.crnKernel result (e.g. other options);
%                                   its .fs must equal 'SampleRate'
%            1xL   double  -        any FIR, L odd, sum ~= 0, applied at the caller's rate;
%                                   [] = no filter
%
% NAME-VALUE ARGUMENTS
%   Name         Size  Type    Units       Description
%   ----         ----  ----    -----       -----------
%   Rotation     3x3   double  -           Body -> accelerometer frame; default eye(3)
%   SampleRate   1x1   double  Hz          Rate of dyn; must equal a CRN kernel's fs; default 10
%   OutputRate   1x1   double  Hz          Rate of y: every (SampleRate/OutputRate)-th filtered
%                                          sample from the first, an integer ratio; default
%                                          'SampleRate' (no sampling). 1 reproduces ACT1B's rate
%   NoiseAsd     3x1   double  m/s^2/rtHz  White floor per axis; default [1e-10; 1e-9; 1e-10],
%                                          the GRACE-FO order with the cross-track axis ten
%                                          times worse (Kornfeld et al. 2019)
%   NoiseCorner  1x1   double  Hz          Below it the ASD rises as 1/f; default 5e-3
%   Name         1x1   string  -           Instrument label; default "electrostatic"
%
% METHODS
%   y          = observe(obj, dyn)     dyn: .a 3xK body, .omega 3xK body, .omegaDot 3xK body,
%                                      .epoch (1xK s, datetime, or a referenceSystem epoch with
%                                      .day/.sec; REQUIRED with a CRN kernel, optional otherwise),
%                                      optional .gravityGradient 3x3xK.
%                                      y: .value 3xM (accelerometer frame, M = ceil(K/stride)),
%                                      .mask 3xM logical, .frame, .epoch (sampled like .value),
%                                      .info (.kernel "crn ACC" | "user" | "none", .outputStride)
%   [J, names] = partials(obj, dyn)    d(y)/d(theta), (3M)x15, theta = [bias; S(:); dr].
%                                      Row (k-1)*3+i is axis i of epoch k: axis fastest
%   asd        = noiseModel(obj, f_Hz) 3xN one-sided amplitude spectral density
%   n          = numParameters(obj)    15
%
% ERRORS
%   NonGravForces:Electrostatic:inputShape   an argument or a dyn field has the wrong size
%   NonGravForces:Electrostatic:inputKernel  the kernel is empty-summed, of even length, not a vector,
%                                            text other than "crn", a struct without .F/.fs, or a
%                                            CRN kernel whose .fs differs from 'SampleRate'
%   NonGravForces:Electrostatic:depMissing   "crn" asked for and timeSeriesPreparation is not on the path
%   NonGravForces:Electrostatic:inputDyn     dyn lacks a, omega or omegaDot, or their epoch counts differ;
%                                            with a CRN kernel: dyn.epoch missing, of the wrong length
%                                            or type, or not on the kernel's grid
%   NonGravForces:Electrostatic:numericNaN   an input contains NaN or Inf
%
% REFERENCES
%   Montenbruck O., Gill E. (2000) Satellite Orbits, sect. 3.7 (specific force at
%   a point offset from the centre of mass); Frommknecht B. (2007) PhD, TU
%   Munich (GRACE accelerometer calibration model); Klinger B., Mayer-Guerr T.
%   (2016) Adv. Space Res. 58:1597 (bias and scale parametrisation);
%   Wu S.-C., Kruizinga G., Bertiger W. (2006) Algorithm Theoretical Basis
%   Document for GRACE Level-1B Data Processing, JPL D-27672, App. B (the CRN
%   kernel; its taps come from timeSeriesPreparation.crnKernel); JPL D-56935
%   GRACE-FO L1 Handbook (frames, proof-mass geometry);
%   Kornfeld R. et al. (2019) J. Spacecraft Rockets 56:931 (noise floor).
%
% EXAMPLE
%   % 10 Hz dynamics through the L1B CRN kernel, sampled at 1 Hz like ACT1B:
%   s = NonGravForces.sensor.Electrostatic([1e-6; -3e-6; 2e-7], diag([1.05 0.96 1.01]), ...
%           [2e-4; 0; -1e-4], "crn", 'OutputRate', 1);
%   dyn = struct('a', aBody, 'omega', w, 'omegaDot', wdot, 'epoch', t10Hz_s);
%   y = s.observe(dyn);
%   r = timeSeriesPreparation.qc.residuals(measured', y.value', 'Whitening', "highpass", 'Window', 601);
%
% See also timeSeriesPreparation.crnKernel, NonGravForces.sensor.Sensor, timeSeriesPreparation.qc.residuals, NonGravForces.app.calibrate, referenceSystem.angularVelocityFromRotationRate.
%
% KEYWORDS accelerometer electrostatic bias scale factor centre of mass offset calibration
% Created 2026-09-08 by Claude Opus 5 (Anthropic); CRN default kernel (N-04) 2026-10-14 by Claude (Anthropic).
    properties (SetAccess = private)
        bias (3,1) double = zeros(3, 1)
        scale (3,3) double = eye(3)
        dr_m (3,1) double = zeros(3, 1)
        rotation (3,3) double = eye(3)
        kernel (1,:) double = []
        kernelSource (1,1) string = "none"
        outputRate (1,1) double = 10
        noiseAsd (3,1) double = [1e-10; 1e-9; 1e-10]
        noiseCorner (1,1) double = 5e-3
    end
    methods
        function obj = Electrostatic(bias, S, dr, kernel, options)
            arguments
                bias double = zeros(3, 1)
                S double = eye(3)
                dr double = zeros(3, 1)
                kernel = "crn"
                options.Rotation double = eye(3)
                options.SampleRate (1,1) double {mustBePositive, mustBeFinite} = 10
                options.OutputRate (1,1) double = NaN
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
                S = diag(S(:));                 % a 3x1 is the diagonal, as timeSeriesPreparation.qc.residuals also accepts
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
            [k, src] = NonGravForces.sensor.Electrostatic.resolveKernel(kernel, options.SampleRate);
            outRate = options.OutputRate;
            if isnan(outRate), outRate = options.SampleRate; end
            stride = options.SampleRate / outRate;
            if ~(isfinite(stride) && stride >= 1 && abs(stride - round(stride)) <= 1e-9 * stride)
                NonGravForces.fail(unit, "inputShape", "'OutputRate' must divide 'SampleRate' (%g Hz) by a whole number -> got %g Hz", options.SampleRate, outRate);
            end
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
            obj.kernelSource = src;
            obj.outputRate = outRate;
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
            idx = obj.outputIndex(K);
            val = val(:, idx); mask = mask(:, idx);
            edge = 0; if ~isempty(obj.kernel), edge = (numel(obj.kernel) - 1) / 2; end
            grad = "none"; if ~isempty(T), grad = "supplied"; end
            info = struct('K', K, 'filtered', ~isempty(obj.kernel), 'kernel', obj.kernelSource, ...
                'edgeSamples', edge, 'outputStride', round(obj.sampleRate / obj.outputRate), ...
                'gradient', grad, 'provenance', NonGravForces.provenance("sensor.Electrostatic"));
            ep = []; if isfield(dyn, 'epoch'), ep = NonGravForces.sensor.Electrostatic.sampleEpoch(dyn.epoch, idx, K); end
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
            % ...and sampled exactly as the observation is.
            idx = obj.outputIndex(K);
            J = J(reshape((idx - 1) * 3 + (1:3)', [], 1), :);
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
        function [a, w, wd, K, T] = checkDyn(obj, dyn)
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
            if startsWith(obj.kernelSource, "crn")
                obj.checkGrid(dyn, K);
            end
        end

        function checkGrid(obj, dyn, K)
            %CHECKGRID A CRN kernel is defined on its own grid: the dynamics must be on it, measured.
            unit = "Electrostatic";
            fixHint = "-> pass dynamics sampled at %g Hz with dyn.epoch, or build the sensor with kernel [] for unfiltered dynamics at another rate";
            if ~isfield(dyn, 'epoch') || isempty(dyn.epoch)
                NonGravForces.fail(unit, "inputDyn", "the CRN kernel is defined on a %g Hz grid and dyn.epoch is needed to check that the dynamics are on it " + fixHint, obj.sampleRate, obj.sampleRate);
            end
            t = NonGravForces.sensor.Electrostatic.epochSeconds(dyn.epoch);
            if numel(t) ~= K
                NonGravForces.fail(unit, "inputDyn", "dyn.epoch must hold one epoch per sample (%d) -> got %d", K, numel(t));
            end
            step = diff(t);
            if K > 1 && any(abs(step - 1 / obj.sampleRate) > 1e-6)
                NonGravForces.fail(unit, "inputDyn", "dyn.epoch is not on the kernel's grid: steps %g to %g s, the kernel needs %g s " + fixHint, ...
                    min(step), max(step), 1 / obj.sampleRate, obj.sampleRate);
            end
        end

        function idx = outputIndex(obj, K)
            %OUTPUTINDEX The samples y keeps: every stride-th from the first. Sampled, never averaged.
            idx = 1:round(obj.sampleRate / obj.outputRate):K;
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

    methods (Static, Access = private)
        function [k, src] = resolveKernel(kernel, sampleRate)
            %RESOLVEKERNEL The taps and where they came from: "crn <preset>", "user" or "none".
            unit = "Electrostatic";
            crn = [];
            if ischar(kernel) || isstring(kernel)
                if ~(isscalar(string(kernel)) && string(kernel) == "crn")
                    NonGravForces.fail(unit, "inputKernel", "a kernel given as text must be ""crn"" -> got ""%s""; pass a numeric FIR or [] otherwise", strjoin(string(kernel), " "));
                end
                % exist() does not resolve package functions (CLAUDE.md section 7); which() does.
                if isempty(which('timeSeriesPreparation.crnKernel'))
                    NonGravForces.fail(unit, "depMissing", "the default kernel is timeSeriesPreparation.crnKernel(""ACC"") and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation, or pass [] for an unfiltered sensor");
                end
                crn = timeSeriesPreparation.crnKernel("ACC");
            elseif isstruct(kernel)
                if ~isscalar(kernel) || ~all(isfield(kernel, {'F', 'fs'}))
                    NonGravForces.fail(unit, "inputKernel", "a kernel given as a struct must be a timeSeriesPreparation.crnKernel result with fields F and fs -> got fields %s", strjoin(string(fieldnames(kernel))', ", "));
                end
                crn = kernel;
            elseif ~isnumeric(kernel)
                NonGravForces.fail(unit, "inputKernel", "the kernel must be ""crn"", a crnKernel struct, a numeric FIR or [] -> got a %s", class(kernel));
            end
            if isempty(crn)
                k = double(kernel(:)');
                src = "user"; if isempty(k), src = "none"; end
                return
            end
            if crn.fs ~= sampleRate
                NonGravForces.fail(unit, "inputKernel", "the CRN kernel is defined at %g Hz and 'SampleRate' is %g Hz -> pass dynamics on the kernel's grid ('SampleRate', %g), or [] for no filter", crn.fs, sampleRate, crn.fs);
            end
            k = crn.F(:)';
            src = "crn";
            if isfield(crn, 'preset'), src = "crn " + string(crn.preset); end
        end

        function t = epochSeconds(e)
            %EPOCHSECONDS Seconds since the first epoch, from any of the three accepted forms.
            if isstruct(e) && isscalar(e) && all(isfield(e, {'day', 'sec'}))
                t = (double(e.day(:)) - double(e.day(1))) * 86400 + (double(e.sec(:)) - double(e.sec(1)));
            elseif isdatetime(e)
                t = seconds(e(:) - e(1));
            elseif isnumeric(e)
                t = double(e(:)) - double(e(1));
            else
                NonGravForces.fail("Electrostatic", "inputDyn", "dyn.epoch must be seconds, a datetime, or a referenceSystem epoch with .day and .sec -> got a %s", class(e));
            end
        end

        function ep = sampleEpoch(e, idx, K)
            %SAMPLEEPOCH The epochs of the samples y keeps, in the form they were given.
            if isempty(e)
                ep = e;
            elseif isstruct(e) && isscalar(e)
                ep = e;
                for f = string(fieldnames(e))'
                    if numel(e.(f)) == K, ep.(f) = e.(f)(idx); end
                end
            elseif numel(e) == K
                ep = e(idx);
            else
                ep = e;
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
