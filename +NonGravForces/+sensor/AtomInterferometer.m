classdef AtomInterferometer < NonGravForces.sensor.Sensor
%ATOMINTERFEROMETER Cold-atom interferometer: phase = k_eff a T^2, with dead time, rotation and the 2 pi ambiguity.
%
%   s = NonGravForces.sensor.AtomInterferometer(kEff, T, cycle, axis)
%   s = NonGravForces.sensor.AtomInterferometer(kEff, T, cycle, axis, 'Rotation', true, 'SampleRate', 10)
%   y = s.observe(dyn);  [J, names] = s.partials(dyn);  asd = s.noiseModel(f)
%
%   WHAT A MACH-ZEHNDER ATOM INTERFEROMETER READS is a phase, not an
%   acceleration, and the difference matters at three points that this class
%   makes explicit rather than hides:
%
%       dphi = phaseBias + (1 + g) k_eff T^2 <n, a>  -  2 (1 + g) k_eff T^2 <n, omega x v_at>
%
%   with n the Raman beam direction in the body frame, T the interrogation time
%   and k_eff the effective wave vector. The second term is the Coriolis phase:
%   a rotating instrument reads the transverse velocity of the ATOMS as an
%   acceleration, which is why an atom interferometer on a satellite needs the
%   attitude rate as an input and not only the acceleration.
%
%   v_at IS THE VELOCITY OF THE ATOM CLOUD RELATIVE TO THE INSTRUMENT, in the body
%   frame -- its launch velocity and the residual of its release, millimetres to
%   centimetres per second -- and NOT the orbital velocity. Satellite and atoms
%   share the orbital velocity, so it cancels from the phase. Until 2026-09-15
%   the field was called `velocity`, the demo fed it 7.6 km/s, and this toolbox
%   recorded an apparent 16.7 m/s^2 that no instrument would see (HANDOVER LP-C,
%   finding F6). At 1.1 mrad/s, 1 mm/s of atom velocity is an apparent 2.2e-6
%   m/s^2, still seven times a 3e-7 drag: the term remains essential, at the
%   right size.
%
%   ONE: THE PHASE IS AN AVERAGE, NOT A SAMPLE. The pulse sequence responds to
%   a(t) through a triangular sensitivity of half-width T about the mirror pulse,
%   so a varying acceleration enters weighted, not sampled. This class applies
%   that weighting as a normalised triangular kernel on the given sample grid.
%   The weights sum to one, so a constant acceleration reproduces k_eff a T^2
%   exactly -- which is the identity the reference test asserts.
%
%   TWO: DEAD TIME. The instrument is blind between cycles. A measurement exists
%   only once per cycle time, and `y.mask` is true only on those epochs and only
%   where the interrogation window fits inside the arc. Nothing is interpolated
%   across the gaps; a caller who wants a continuous series must say how, and
%   that is a decision about their problem, not about this instrument.
%
%   THREE: THE 2 PI AMBIGUITY IS REAL AND IS NOT SWEPT UP. The fringe readout
%   gives the phase modulo 2 pi and nothing more. `y.value` carries the FULL
%   modelled phase, because a model knows the fringe it is on; `y.info.wrapped`
%   carries what an instrument would actually return, and `y.info.fringe` the
%   integer that separates them. A caller comparing against real data compares
%   the wrapped one, or resolves the ambiguity by other means first -- with a
%   classical accelerometer beside it, which is the usual answer and the reason
%   hybrid instruments exist.
%
% CONSTRUCTOR ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   kEff    1x1   double  rad/m  Effective wave vector, > 0. Two-photon Raman on the rubidium
%                                D2 line is 4 pi / 780.24 nm = 1.610e7 rad/m
%   T       1x1   double  s      Interrogation time (half the pulse separation), > 0
%   cycle   1x1   double  s      Cycle time, >= 2 T. The dead time is cycle - 2 T
%   axis    3x1   double  -      Raman beam direction in the body frame; normalised on entry
%
% NAME-VALUE ARGUMENTS
%   Name         Size  Type     Units      Description
%   ----         ----  ----     -----      -----------
%   Rotation     1x1   logical  -          Include the Coriolis phase; default false. When true,
%                                          dyn must carry .omega and .atomVelocity
%   SampleRate   1x1   double   Hz         Sample rate of dyn, > 0; default 1. Sets the grid the
%                                          triangular weighting and the cycle are counted on
%   PhaseBias    1x1   double   rad        Additive phase offset; default 0
%   GainError    1x1   double   -          Relative error of k_eff T^2; default 0
%   PhaseAsd     1x1   double   rad/rtHz   White phase noise floor; default 1e-3
%   Name         1x1   string   -          Instrument label; default "atomInterferometer"
%
% METHODS
%   y          = observe(obj, dyn)     dyn: .a 3xK body [m/s^2]; when 'Rotation', .omega 3xK
%                                      [rad/s] and .atomVelocity 3xK [m/s], the atom cloud
%                                      relative to the instrument, body frame.
%                                      y: .value 1xK full phase [rad], .mask 1xK,
%                                      .frame "beam", .info (.wrapped .fringe .weighted
%                                      .coriolis .cycleEpochs)
%   [J, names] = partials(obj, dyn)    d(phase)/d([phaseBias; gainError]), Kx2
%   asd        = noiseModel(obj, f_Hz) 1xN phase ASD [rad/sqrt(Hz)]; divide by k_eff T^2 for
%                                      the equivalent acceleration ASD
%   n          = numParameters(obj)    2
%
% ERRORS
%   NonGravForces:AtomInterferometer:inputShape   an argument or a dyn field has the wrong size
%   NonGravForces:AtomInterferometer:inputRange   kEff <= 0, T <= 0, cycle < 2T, SampleRate <= 0,
%                                                 or the axis has zero length
%   NonGravForces:AtomInterferometer:inputDyn     dyn lacks a, or lacks omega/atomVelocity when
%                                                 'Rotation' is on, or their epoch counts differ, or
%                                                 it carries the retired field velocity
%   NonGravForces:AtomInterferometer:numericNaN   an input contains NaN or Inf
%
% REFERENCES
%   Kasevich M., Chu S. (1991) Phys. Rev. Lett. 67:181 (the pulse sequence);
%   Cheinet P. et al. (2008) IEEE Trans. Instrum. Meas. 57:1141 (the sensitivity
%   function); Lan S.-Y. et al. (2012) Phys. Rev. Lett. 108:090402 (the Coriolis
%   phase and its compensation); Leveque T. et al. (2022) CEAS Space J. 14:597
%   (CARIOQA: a space-borne quantum accelerometer, its cycle time and dead time);
%   Zahzam N. et al. (2022) Remote Sens. 14:3273 (hybridisation with a classical
%   accelerometer, i.e. how the 2 pi ambiguity is resolved in practice).
%   HANDOVER.md Part G LP-7 (WP-29).
%
% EXAMPLE
%   s = NonGravForces.sensor.AtomInterferometer(1.610e7, 0.5, 2, [1; 0; 0], 'SampleRate', 10);
%   y = s.observe(struct('a', aBody));
%   phaseAtCycles = y.value(y.mask);        % one reading per cycle, nothing between
%   aEquivalent   = y.value / (1.610e7 * 0.5^2);
%
% See also NonGravForces.sensor.Sensor, NonGravForces.sensor.Electrostatic, NonGravForces.app.calibrate.
%
% KEYWORDS cold atom interferometer quantum accelerometer phase dead time Coriolis fringe
% Created 2026-09-09 by Claude Opus 5 (Anthropic); atomVelocity 2026-09-15 by Claude Opus 5 (Anthropic).
    properties (SetAccess = private)
        kEff (1,1) double = 1.610e7
        T_s (1,1) double = 0.5
        cycle_s (1,1) double = 2
        beam (3,1) double = [1; 0; 0]
        useRotation (1,1) logical = false
        phaseBias (1,1) double = 0
        gainError (1,1) double = 0
        phaseAsd (1,1) double = 1e-3
    end
    methods
        function obj = AtomInterferometer(kEff, T, cycle, axis, options)
            arguments
                kEff (1,1) double = 1.610e7
                T (1,1) double = 0.5
                cycle (1,1) double = 2
                axis double = [1; 0; 0]
                options.Rotation (1,1) logical = false
                options.SampleRate (1,1) double = 1
                options.PhaseBias (1,1) double = 0
                options.GainError (1,1) double = 0
                options.PhaseAsd (1,1) double = 1e-3
                options.Name (1,1) string = "atomInterferometer"
            end
            unit = "AtomInterferometer";
            n = axis(:);
            if ~isequal(size(n), [3 1])
                NonGravForces.fail(unit, "inputShape", "axis must be 3x1, the Raman beam direction in the body frame -> got %s", mat2str(size(axis)));
            end
            if any(~isfinite([kEff; T; cycle; n; options.SampleRate; options.PhaseBias; options.GainError; options.PhaseAsd]))
                NonGravForces.fail(unit, "numericNaN", "a constructor argument contains NaN/Inf -> screen the instrument parameters before building it");
            end
            if ~(kEff > 0) || ~(T > 0) || ~(options.SampleRate > 0) || options.PhaseAsd < 0
                NonGravForces.fail(unit, "inputRange", "kEff, T and 'SampleRate' must be positive and 'PhaseAsd' non-negative -> got %g, %g, %g, %g", kEff, T, options.SampleRate, options.PhaseAsd);
            end
            if norm(n) == 0
                NonGravForces.fail(unit, "inputRange", "the beam axis has zero length -> it has no direction to project the acceleration on");
            end
            if cycle < 2 * T
                NonGravForces.fail(unit, "inputRange", "the cycle time %g s is shorter than the interrogation 2T = %g s -> a sequence cannot finish before it ends; the dead time is cycle - 2T and cannot be negative", cycle, 2 * T);
            end
            obj.kEff = kEff;
            obj.T_s = T;
            obj.cycle_s = cycle;
            obj.beam = n / norm(n);
            obj.useRotation = options.Rotation;
            obj.phaseBias = options.PhaseBias;
            obj.gainError = options.GainError;
            obj.phaseAsd = options.PhaseAsd;
            obj.name = options.Name;
            obj.axes = "beam";
            obj.sampleRate = options.SampleRate;
        end

        function y = observe(obj, dyn)
            %OBSERVE The interferometer phase for the given true dynamics.
            [a, w, vel, K] = obj.checkDyn(dyn);
            [aw, valid] = obj.weight(obj.beam' * a);
            cor = zeros(1, K);
            if obj.useRotation
                % -2 <n, omega x v_at>: a rotating instrument reads the atoms' transverse velocity as an
                % acceleration. It goes through the SAME weighting, because it is part of the
                % acceleration the sequence integrates, not a correction applied afterwards.
                cor = obj.weight(-2 * (obj.beam' * cross(w, vel)));
            end
            gain = (1 + obj.gainError) * obj.kEff * obj.T_s^2;
            phase = obj.phaseBias + gain * (aw + cor);
            mask = valid & obj.cycleMask(K);
            wrapped = mod(phase + pi, 2 * pi) - pi;          % what an instrument returns: (-pi, pi]
            fringe = round((phase - wrapped) / (2 * pi));    % and the integer that was lost
            y = struct('value', phase, 'mask', mask, 'frame', "beam", 'epoch', [], 'info', struct());
            if isfield(dyn, 'epoch'), y.epoch = dyn.epoch; end
            y.info = struct('wrapped', wrapped, 'fringe', fringe, 'weighted', aw, 'coriolis', cor, ...
                'cycleEpochs', find(obj.cycleMask(K)), 'deadTime_s', obj.cycle_s - 2 * obj.T_s, ...
                'provenance', NonGravForces.provenance("sensor.AtomInterferometer"));
        end

        function [J, names] = partials(obj, dyn)
            %PARTIALS d(phase)/d([phaseBias; gainError]), Kx2. Both enter linearly.
            [a, w, vel, K] = obj.checkDyn(dyn);
            aw = obj.weight(obj.beam' * a);
            cor = zeros(1, K);
            if obj.useRotation
                cor = obj.weight(-2 * (obj.beam' * cross(w, vel)));
            end
            J = [ones(K, 1), (obj.kEff * obj.T_s^2 * (aw + cor))'];
            names = ["phaseBias", "gainError"];
        end

        function S = noiseModel(obj, f_Hz)
            %NOISEMODEL Phase ASD, 1xN. White to the cycle rate; above it there is no measurement.
            f = f_Hz(:)';
            if isempty(f) || any(f <= 0) || any(~isfinite(f))
                NonGravForces.fail("AtomInterferometer", "inputShape", "f_Hz must be positive and finite -> the ASD is not defined at or below zero frequency");
            end
            S = obj.phaseAsd * ones(1, numel(f));
            S(f > 1 / (2 * obj.cycle_s)) = Inf;   % beyond the Nyquist of the cycle rate nothing is measured
        end

        function n = numParameters(~)
            n = 2;
        end
    end

    methods (Access = private)
        function [a, w, vel, K] = checkDyn(obj, dyn)
            unit = "AtomInterferometer";
            if ~isstruct(dyn) || ~isscalar(dyn) || ~isfield(dyn, 'a')
                NonGravForces.fail(unit, "inputDyn", "dyn must be a 1x1 struct with field a (3xK) -> and with omega and atomVelocity when 'Rotation' is on");
            end
            a = dyn.a;
            if ~ismatrix(a) || size(a, 1) ~= 3
                NonGravForces.fail(unit, "inputShape", "dyn.a must be 3xK -> got %s", mat2str(size(a)));
            end
            K = size(a, 2);
            w = zeros(3, K); vel = zeros(3, K);
            if isfield(dyn, 'velocity')
                % Refused, not ignored: a caller still passing the orbital velocity under the old name
                % would otherwise get a phase without the Coriolis term and no sign of why.
                NonGravForces.fail(unit, "inputDyn", "dyn.velocity was retired on 2026-09-15 -> pass dyn.atomVelocity, the atom cloud's velocity RELATIVE TO THE INSTRUMENT (mm/s to cm/s), not the orbital velocity, which atoms and instrument share");
            end
            if obj.useRotation
                if ~all(isfield(dyn, {'omega', 'atomVelocity'}))
                    NonGravForces.fail(unit, "inputDyn", "'Rotation' is on, so dyn needs omega and atomVelocity (3xK each) -> the Coriolis phase is -2 k_eff T^2 <n, omega x v_at> and neither term is optional");
                end
                w = dyn.omega; vel = dyn.atomVelocity;
                if ~isequal(size(w), [3 K]) || ~isequal(size(vel), [3 K])
                    NonGravForces.fail(unit, "inputDyn", "dyn.omega and dyn.atomVelocity must be 3x%d to match dyn.a -> got %s and %s", K, mat2str(size(w)), mat2str(size(vel)));
                end
            end
            if any(~isfinite([a(:); w(:); vel(:)]))
                NonGravForces.fail(unit, "numericNaN", "dyn contains NaN/Inf -> screen the dynamics before observing them");
            end
        end

        function [out, valid] = weight(obj, x)
            %WEIGHT The triangular sensitivity of the pulse sequence, normalised to unit sum.
            % Half-width T on the sample grid. The normalisation is what makes a constant
            % acceleration give back exactly k_eff a T^2 rather than that times a kernel sum.
            K = numel(x);
            h = floor(obj.T_s * obj.sampleRate);
            valid = true(1, K);
            if h < 1
                out = x;                       % the window is shorter than one sample: it IS a sample
                return
            end
            ker = 1 - abs(-h:h) / (h + 1);
            ker = ker / sum(ker);
            out = conv(x, ker, 'same');
            if K <= 2 * h
                valid(:) = false;              % no sample has a full window
            else
                valid([1:h, K - h + 1:K]) = false;
            end
        end

        function m = cycleMask(obj, K)
            %CYCLEMASK True once per cycle: the instrument is blind between sequences.
            step = max(1, round(obj.cycle_s * obj.sampleRate));
            m = false(1, K);
            m(1:step:K) = true;
        end
    end
end
