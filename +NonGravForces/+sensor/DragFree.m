classdef DragFree < NonGravForces.sensor.Sensor
%DRAGFREE Residual acceleration left by a drag-free control loop: sensitivity high-pass and actuator saturation.
%
%   s = NonGravForces.sensor.DragFree(bandwidth_Hz, order)
%   s = NonGravForces.sensor.DragFree(1e-3, 2, 'Authority', 1e-6, 'SampleRate', 10)
%   y = s.observe(dyn);  [J, names] = s.partials(dyn);  asd = s.noiseModel(f)
%
%   WHAT A DRAG-FREE SATELLITE ACTUALLY DELIVERS is not zero non-gravitational
%   acceleration but a RESIDUAL, and its shape is the loop's sensitivity
%   function. Inside the control bandwidth the thrusters cancel the drag; outside
%   it they cannot act at all, and the drag passes through untouched:
%
%       a_residual = bias + (1 + g) [ S{a_ng}  +  saturation excess ]
%
%   with S a high-pass of corner f_c and the stated order. So this class is the
%   opposite of an accelerometer: `sensor.Electrostatic` reports what a sensor
%   reads, and this reports what a controller FAILS to remove. It is in the
%   sensor layer because that is what a science instrument downstream sees.
%
%   THE FILTER IS A CASCADE, NOT A BUTTERWORTH, and the difference is stated
%   rather than left for someone to assume. `order` first-order sections are
%   applied in series, each the bilinear transform of s/(s + 2 pi f_c), so the
%   magnitude is |f/f_c|^n / (1 + (f/f_c)^2)^(n/2) and the phase is n times a
%   first-order phase. A real loop's shape is set by its controller and is the
%   caller's to model; what this offers is the right ASYMPTOTES -- n poles of
%   rejection below f_c, unity above it -- with no toolbox dependency, since
%   C2 rules out the Signal Processing Toolbox and a hand-written IIR is four
%   lines.
%
%   SATURATION IS MODELLED BECAUSE IT IS THE FAILURE MODE THAT MATTERS. A
%   thruster set has finite authority; when the drag exceeds it the loop commands
%   everything it has and the excess reaches the payload. That is a hard
%   nonlinearity, it is what a storm does to a drag-free mission, and leaving it
%   out would make the residual look linear in the drag when it stops being so
%   exactly where the mission cares. With `'Authority'` at Inf (the default) the
%   loop is ideal and the model is linear again.
%
%   THE SETTLING TRANSIENT IS MASKED, NOT HIDDEN. The filter starts from rest, so
%   the first samples carry a start-up transient that is an artefact of where the
%   arc was cut and not a property of the loop. `y.mask` is false over five time
%   constants per section; nothing is trimmed or back-filled.
%
% CONSTRUCTOR ARGUMENTS
%   Name           Size  Type    Units  Description
%   ----           ----  ----    -----  -----------
%   bandwidth_Hz   1x1   double  Hz     Control bandwidth f_c, > 0. Below it the loop rejects
%   order          1x1   double  -      Number of first-order sections, integer >= 1; default 1
%
% NAME-VALUE ARGUMENTS
%   Name          Size  Type    Units       Description
%   ----          ----  ----    -----       -----------
%   Authority     1x1   double  m/s^2       Maximum acceleration the actuators can command,
%                                           > 0; default Inf (an ideal loop)
%   SampleRate    1x1   double  Hz          Sample rate of dyn, > 0; default 1. Must exceed
%                                           2 f_c or the bandwidth is not on the grid
%   Bias          3x1   double  m/s^2       Additive residual offset; default zeros
%   GainError     1x1   double  -           Relative error of the residual; default 0
%   ThrusterAsd   3x1   double  m/s^2/rtHz  Actuation noise floor injected inside the
%                                           bandwidth; default 1e-11 per axis
%   Name          1x1   string  -           Label; default "dragFree"
%
% METHODS
%   y          = observe(obj, dyn)     dyn: .a 3xK non-gravitational acceleration, body frame.
%                                      y: .value 3xK residual, .mask 3xK, .frame "body",
%                                      .info (.commanded .saturated .settleSamples .gainAtBandwidth)
%   [J, names] = partials(obj, dyn)    d(residual)/d([bias; gainError]), (3K)x4, axis fastest
%   asd        = noiseModel(obj, f_Hz) 3xN residual ASD: thruster noise shaped by the
%                                      complementary sensitivity, i.e. injected where the loop
%                                      is active and absent where it is not
%   n          = numParameters(obj)    4
%   m          = magnitudeResponse(obj, f_Hz)   |S(f)| of the cascade, evaluated on the unit
%                                      circle from the difference-equation coefficients, so a
%                                      caller (and a test) can ask what this loop does to a
%                                      frequency without running a series through it
%
% ERRORS
%   NonGravForces:DragFree:inputShape   an argument or dyn.a has the wrong size
%   NonGravForces:DragFree:inputRange   bandwidth <= 0, order < 1 or not an integer,
%                                       Authority <= 0, SampleRate <= 2 f_c
%   NonGravForces:DragFree:inputDyn     dyn lacks a
%   NonGravForces:DragFree:numericNaN   an input contains NaN or Inf
%
% REFERENCES
%   Lange B. (1964) AIAA J. 2:1590 (the drag-free concept); Canuto E. (2008)
%   Automatica 44:1766 (drag-free control design and the sensitivity function);
%   Tapley B. et al. (2004) GRL 31:L09607 for the drag magnitudes this must
%   survive; Armano M. et al. (2018) PRL 120:061101 (LISA Pathfinder: what a
%   real residual looks like). Oppenheim A., Schafer R. (2010) Discrete-Time
%   Signal Processing, sect. 7.1 for the bilinear transform used here.
%   HANDOVER.md Part G LP-7 (WP-29).
%
% EXAMPLE
%   [~, ai] = M.accel(ep, r, v, q);
%   s = NonGravForces.sensor.DragFree(1e-3, 2, 'Authority', 5e-7, 'SampleRate', 1);
%   y = s.observe(struct('a', ai.drag + ai.radiation + ai.thermal));
%   rmsResidual = sqrt(mean(y.value(:, y.mask(1, :)).^2, 2));   % rms() is Signal Processing Toolbox
%
% See also NonGravForces.sensor.Sensor, NonGravForces.sensor.Electrostatic, NonGravForces.Model, NonGravForces.app.forceBudget.
%
% KEYWORDS drag-free control residual acceleration sensitivity function thruster authority saturation
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
    properties (SetAccess = private)
        bandwidth_Hz (1,1) double = 1e-3
        order (1,1) double = 1
        authority (1,1) double = Inf
        bias (3,1) double = zeros(3, 1)
        gainError (1,1) double = 0
        thrusterAsd (3,1) double = 1e-11 * ones(3, 1)
    end
    methods
        function obj = DragFree(bandwidth_Hz, order, options)
            arguments
                bandwidth_Hz (1,1) double = 1e-3
                order (1,1) double = 1
                options.Authority (1,1) double = Inf
                options.SampleRate (1,1) double = 1
                options.Bias double = zeros(3, 1)
                options.GainError (1,1) double = 0
                options.ThrusterAsd double = 1e-11 * ones(3, 1)
                options.Name (1,1) string = "dragFree"
            end
            unit = "DragFree";
            b = options.Bias(:);
            asd = options.ThrusterAsd(:);
            if isscalar(asd), asd = repmat(asd, 3, 1); end
            if ~isequal(size(b), [3 1])
                NonGravForces.fail(unit, "inputShape", "'Bias' must be 3x1 -> got %s", mat2str(size(options.Bias)));
            end
            if ~isequal(size(asd), [3 1])
                NonGravForces.fail(unit, "inputShape", "'ThrusterAsd' must be 3x1 or scalar -> got %s", mat2str(size(options.ThrusterAsd)));
            end
            if any(~isfinite([bandwidth_Hz; order; b; asd; options.SampleRate; options.GainError]))
                NonGravForces.fail(unit, "numericNaN", "a constructor argument contains NaN/Inf -> Authority may be Inf, nothing else may");
            end
            if ~(bandwidth_Hz > 0)
                NonGravForces.fail(unit, "inputRange", "the control bandwidth must be positive -> got %g Hz", bandwidth_Hz);
            end
            if order < 1 || mod(order, 1) ~= 0
                NonGravForces.fail(unit, "inputRange", "'order' must be an integer >= 1: it is the number of first-order sections in the cascade -> got %g", order);
            end
            if ~(options.Authority > 0)
                NonGravForces.fail(unit, "inputRange", "'Authority' must be positive (Inf for an ideal loop) -> a zero-authority loop is no loop, and passing the drag through unchanged is what Authority = eps already does");
            end
            if any(asd < 0)
                NonGravForces.fail(unit, "inputRange", "'ThrusterAsd' must be non-negative -> got %s", mat2str(asd'));
            end
            if options.SampleRate <= 2 * bandwidth_Hz
                NonGravForces.fail(unit, "inputRange", "'SampleRate' %g Hz must exceed twice the bandwidth %g Hz -> below that the corner is not on the grid and the bilinear transform warps it past Nyquist", options.SampleRate, bandwidth_Hz);
            end
            obj.bandwidth_Hz = bandwidth_Hz;
            obj.order = order;
            obj.authority = options.Authority;
            obj.bias = b;
            obj.gainError = options.GainError;
            obj.thrusterAsd = asd;
            obj.name = options.Name;
            obj.axes = ["x", "y", "z"];
            obj.sampleRate = options.SampleRate;
        end

        function y = observe(obj, dyn)
            %OBSERVE The residual acceleration the payload sees.
            [a, K] = obj.checkDyn(dyn);
            lin = obj.highpass(a);                        % what the loop cannot reject
            cmd = a - lin;                                % what it therefore has to command
            excess = zeros(3, K);
            nrm = vecnorm(cmd);
            over = nrm > obj.authority;
            if any(over)
                % The command is clipped in its own direction; the difference reaches the payload.
                excess(:, over) = cmd(:, over) .* (1 - obj.authority ./ nrm(over));
            end
            res = obj.bias + (1 + obj.gainError) * (lin + excess);
            settle = obj.settleSamples();
            mask = true(3, K);
            mask(:, 1:min(settle, K)) = false;
            info = struct('commanded', cmd - excess, 'saturated', over, 'settleSamples', settle, ...
                'gainAtBandwidth', obj.magnitudeResponse(obj.bandwidth_Hz), ...
                'provenance', NonGravForces.provenance("sensor.DragFree"));
            ep = []; if isfield(dyn, 'epoch'), ep = dyn.epoch; end
            y = struct('value', res, 'mask', mask, 'frame', "body", 'epoch', ep, 'info', info);
        end

        function [J, names] = partials(obj, dyn)
            %PARTIALS d(residual)/d([bias; gainError]), (3K)x4, row (k-1)*3+i is axis i of epoch k.
            [a, K] = obj.checkDyn(dyn);
            lin = obj.highpass(a);
            cmd = a - lin;
            excess = zeros(3, K);
            nrm = vecnorm(cmd);
            over = nrm > obj.authority;
            if any(over)
                excess(:, over) = cmd(:, over) .* (1 - obj.authority ./ nrm(over));
            end
            core = lin + excess;                          % what the gain multiplies
            J = zeros(3 * K, 4);
            for k = 1:K
                rows = (k - 1) * 3 + (1:3);
                J(rows, 1:3) = eye(3);
                J(rows, 4) = core(:, k);
            end
            names = ["bias_x", "bias_y", "bias_z", "gainError"];
        end

        function S = noiseModel(obj, f_Hz)
            %NOISEMODEL Residual ASD: thruster noise shaped by the COMPLEMENTARY sensitivity.
            % The loop injects its actuation noise where it is active and injects nothing where it
            % is not, so the noise follows the low-pass complement of the rejection, not the
            % rejection itself. That is the opposite shape from the signal path and is the reason
            % a wider bandwidth is not free.
            f = f_Hz(:)';
            if isempty(f) || any(f <= 0) || any(~isfinite(f))
                NonGravForces.fail("DragFree", "inputShape", "f_Hz must be positive and finite -> the ASD is not defined at or below zero frequency");
            end
            ratio = f / obj.bandwidth_Hz;
            S = obj.thrusterAsd .* (1 ./ sqrt(1 + ratio.^2)).^obj.order;
        end

        function n = numParameters(~)
            n = 4;
        end

        function m = magnitudeResponse(obj, f_Hz)
            %MAGNITUDERESPONSE |S(f)| of the cascade, from the DIFFERENCE EQUATION coefficients.
            % Published so a caller can check what this loop does to a given frequency without
            % running a series through it, and so a test has something to compare against that is
            % not the filter itself.
            f = f_Hz(:)';
            dt = 1 / obj.sampleRate;
            c = 2 * pi * obj.bandwidth_Hz * dt;
            alpha = 2 / (2 + c);
            beta = (2 - c) / (2 + c);
            z = exp(-2i * pi * f * dt);                   % z^-1 on the unit circle
            h = alpha * (1 - z) ./ (1 - beta * z);
            m = abs(h).^obj.order;
        end
    end

    methods (Access = private)
        function [a, K] = checkDyn(~, dyn)
            unit = "DragFree";
            if ~isstruct(dyn) || ~isscalar(dyn) || ~isfield(dyn, 'a')
                NonGravForces.fail(unit, "inputDyn", "dyn must be a 1x1 struct with field a (3xK), the non-gravitational acceleration the loop is fighting -> Model.accel returns it");
            end
            a = dyn.a;
            if ~ismatrix(a) || size(a, 1) ~= 3
                NonGravForces.fail(unit, "inputShape", "dyn.a must be 3xK -> got %s", mat2str(size(a)));
            end
            K = size(a, 2);
            if any(~isfinite(a(:)))
                NonGravForces.fail(unit, "numericNaN", "dyn.a contains NaN/Inf -> screen the accelerations before filtering them");
            end
        end

        function n = settleSamples(obj)
            %SETTLESAMPLES Five time constants per section, the transient of starting from rest.
            tau = 1 / (2 * pi * obj.bandwidth_Hz);
            n = ceil(5 * tau * obj.sampleRate) * obj.order;
        end

        function out = highpass(obj, x)
            %HIGHPASS The cascade, bilinear, written out rather than borrowed from a toolbox.
            % One section: y(k) = alpha (x(k) - x(k-1)) + beta y(k-1), which is s/(s + w_c) under
            % the bilinear transform with w_c = 2 pi f_c.
            dt = 1 / obj.sampleRate;
            c = 2 * pi * obj.bandwidth_Hz * dt;
            alpha = 2 / (2 + c);
            beta = (2 - c) / (2 + c);
            out = x;
            for s = 1:obj.order
                in = out;
                out = zeros(size(in));
                for k = 2:size(in, 2)
                    out(:, k) = alpha * (in(:, k) - in(:, k - 1)) + beta * out(:, k - 1);
                end
            end
        end
    end
end
