classdef tNgfSensor < NgfTestCase
%TNGFSENSOR Tests of the sensor layer (LP-7, WP-29): sensor.Electrostatic, sensor.AtomInterferometer, sensor.DragFree.
%
%   ORACLES, none of them this class's own arithmetic:
%     - the rotational coupling against a FINITE-DIFFERENCE of the position of the
%       proof mass through a real attitude arc. omegaDot x dr + omega x (omega x dr)
%       is the second time derivative of R(t) dr resolved back in the body frame, so
%       differentiating the real SCA1B quaternions twice and comparing is a check
%       against kinematics rather than against the same formula written twice;
%     - the partials against central differences of observe itself;
%     - a rigid-body identity: a proof mass on the rotation axis feels nothing;
%     - referenceSystem.angularVelocityFromRotationRate for omega, which is where a
%       caller gets it, and which owns the sign convention this class must match.
%
%   The GRACE-FO numbers are used where a number is needed: the handbook proof-mass
%   offset is a few hundred micrometres, the yaw-steering rate is about one orbit per
%   revolution, and the resulting coupling is of order 1e-9 m/s^2 -- three orders below
%   the drag it sits on and well above the 1e-10 noise floor, which is exactly why it
%   is modelled at all.
%
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
    properties
        root, att, ep, q, w, wd
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            f = @(varargin) fullfile(tc.root, varargin{:});
            [eps2, qq] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
            tc.att = NonGravForces.geom.AttitudeProvider(eps2, qq);
            % A real attitude arc at 1 s, long enough for a centred second difference.
            n = 60;
            tc.ep = referenceSystem.daySecondToEpoch([eps2.day(1:n) eps2.sec(1:n)], eps2.scale);
            tc.q = tc.att.at(tc.ep);
            R = referenceSystem.rotationFromQuaternion(tc.q);
            [tc.w, tc.wd] = tNgfSensor.ratesFromRotation(R, 1);
        end
    end
    methods (Static)
        function [w, wd] = ratesFromRotation(R, dt)
            %RATESFROMROTATION Angular velocity and its derivative in the sensor's own convention.
            % THE CONVENTION IS WRITTEN OUT HERE RATHER THAN IMPORTED, because it is the one thing
            % in this model a sign error hides in. For R mapping BODY -> INERTIAL and a vector
            % fixed in the body, differentiating v_inertial = R v_body twice gives
            % R' Rddot = Omega_dot + Omega^2, which is exactly omegaDot x dr + omega x (omega x dr)
            % -- so the sensor's omega is the one with skew(omega) = R' Rdot, and nothing else.
            % referenceSystem.angularVelocityFromRotationRate carries a configurable sense and a
            % rotation whose direction is the transpose of this one; which of its two Resolution
            % settings lands on this definition is recorded, not assumed (HANDOVER O-32).
            K = size(R, 3);
            Rdot = zeros(3, 3, K);
            Rdot(:, :, 2:K-1) = (R(:, :, 3:K) - R(:, :, 1:K-2)) / (2 * dt);
            Rdot(:, :, 1) = Rdot(:, :, 2); Rdot(:, :, K) = Rdot(:, :, K-1);
            M = pagemtimes(pagetranspose(R), Rdot);
            Om = (M - pagetranspose(M)) / 2;      % the symmetric part is the differencing error,
            % not angular velocity, and reading the skew entries off M without projecting it out
            % leaves it in: 6.0e-7 relative on this arc, measured 2026-09-08. referenceSystem
            % projects, so a comparison against it must project too or it measures the projection.
            w = [squeeze(Om(3, 2, :)), squeeze(Om(1, 3, :)), squeeze(Om(2, 1, :))]';
            wd = zeros(3, K);
            wd(:, 2:K-1) = (w(:, 3:K) - w(:, 1:K-2)) / (2 * dt);
            wd(:, 1) = wd(:, 2); wd(:, K) = wd(:, K-1);
        end
        function R = spin(axis, rate0, accel, t)
            %SPIN A body->inertial rotation about a fixed axis whose rate ramps linearly.
            % Analytic, so the second difference below has only its own truncation error, and both
            % the centrifugal and the Euler term are non-zero.
            n = axis(:) / norm(axis);
            N = NonGravForces.sensor.Electrostatic.skew(n);
            K = numel(t);
            R = zeros(3, 3, K);
            for k = 1:K
                th = rate0 * t(k) + 0.5 * accel * t(k)^2;
                R(:, :, k) = eye(3) + sin(th) * N + (1 - cos(th)) * (N * N);   % Rodrigues
            end
        end
    end
    methods
        function dyn = scene(tc, K)
            a = 1e-7 * [sin((1:K) / 9); 0.3 * cos((1:K) / 7); 0.1 * sin((1:K) / 13)] - [3e-7; 0; 0];
            dyn = struct('a', a, 'omega', tc.w(:, 1:K), 'omegaDot', tc.wd(:, 1:K), 'epoch', []);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function ElectrostaticContract(tc)
            s = NonGravForces.sensor.Electrostatic([1e-6; -3e-6; 2e-7], diag([1.05 0.96 1.01]), [2e-4; 0; -1e-4]);
            tc.verifyEqual(s.numParameters(), 15);
            tc.verifyEqual(s.name, "electrostatic");
            tc.verifyTrue(isa(s, 'NonGravForces.sensor.Sensor'));
            dyn = tc.scene(20);
            y = s.observe(dyn);
            tc.verifySize(y.value, [3 20]); tc.verifySize(y.mask, [3 20]);
            tc.verifyEqual(y.frame, "accelerometer");
            tc.verifyTrue(all(y.mask(:)), 'an unfiltered reading has no invalid samples');
            tc.verifyEqual(y.info.gradient, "none");
            [J, names] = s.partials(dyn);
            tc.verifySize(J, [60 15]); tc.verifySize(names, [1 15]);
            tc.verifyEqual(names(1), "bias_x"); tc.verifyEqual(names(end), "dr_z");
            tc.verifySize(s.noiseModel([1e-3 1e-2 1e-1]), [3 3]);

            % A 3x1 scale is its diagonal, as qc.residuals also accepts.
            s2 = NonGravForces.sensor.Electrostatic(zeros(3, 1), [1.05; 0.96; 1.01]);
            tc.verifyEqual(s2.scale, diag([1.05 0.96 1.01]));

            E = @NonGravForces.sensor.Electrostatic;
            tc.verifyErrorId(@() E([1; 2]), 'NonGravForces:Electrostatic:inputShape');
            tc.verifyErrorId(@() E(zeros(3, 1), ones(2, 2)), 'NonGravForces:Electrostatic:inputShape');
            tc.verifyErrorId(@() E(zeros(3, 1), eye(3), [1; 2]), 'NonGravForces:Electrostatic:inputShape');
            tc.verifyErrorId(@() E(zeros(3, 1), eye(3), zeros(3, 1), [1 1 1 1]), 'NonGravForces:Electrostatic:inputKernel');
            tc.verifyErrorId(@() E(zeros(3, 1), eye(3), zeros(3, 1), [1 -2 1]), 'NonGravForces:Electrostatic:inputKernel');
            tc.verifyErrorId(@() E([NaN; 0; 0]), 'NonGravForces:Electrostatic:numericNaN');
            tc.verifyErrorId(@() E(zeros(3, 1), eye(3), zeros(3, 1), [], 'Rotation', eye(2)), 'NonGravForces:Electrostatic:inputShape');

            tc.verifyErrorId(@() s.observe(struct('a', dyn.a)), 'NonGravForces:Electrostatic:inputDyn');
            tc.verifyErrorId(@() s.observe(struct('a', dyn.a, 'omega', dyn.omega(:, 1:5), 'omegaDot', dyn.omegaDot)), 'NonGravForces:Electrostatic:inputDyn');
            bad = dyn; bad.a(1) = Inf;
            tc.verifyErrorId(@() s.observe(bad), 'NonGravForces:Electrostatic:numericNaN');
            bad = dyn; bad.gravityGradient = zeros(3, 3, 5);
            tc.verifyErrorId(@() s.observe(bad), 'NonGravForces:Electrostatic:inputShape');
            tc.verifyErrorId(@() s.noiseModel([0 1]), 'NonGravForces:Electrostatic:inputShape');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function ElectrostaticReference(tc)
            % THE COUPLING IS KINEMATICS, so it is checked against kinematics and not against the
            % same formula written a second time. A point rigidly attached at dr in the body frame
            % sits at R(t) dr in the inertial frame; its inertial acceleration is the second
            % derivative of that, and resolved back into the body frame it must be exactly what the
            % sensor computes. The rotation is ANALYTIC -- a fixed axis whose rate ramps, so both
            % the centrifugal and the Euler term are non-zero -- which leaves the second difference
            % with only its own O(dt^2) truncation instead of the attitude product's noise.
            dt = 1; t = (0:79) * dt; K = numel(t);
            dr = [2.34e-4; -1.1e-4; 3.0e-4];
            R = tNgfSensor.spin([0.3; -0.2; 0.93], 2.0e-3, 4e-6, t);
            [w, wd] = tNgfSensor.ratesFromRotation(R, dt);
            p = zeros(3, K);
            for k = 1:K, p(:, k) = R(:, :, k) * dr; end          % inertial position of the proof mass
            want = zeros(3, K);
            for k = 2:K-1
                acc = (p(:, k+1) - 2 * p(:, k) + p(:, k-1)) / dt^2;
                want(:, k) = R(:, :, k)' * acc;                  % resolved back in the body frame
            end
            s = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), dr);
            got = s.observe(struct('a', zeros(3, K), 'omega', w, 'omegaDot', wd)).value;
            idx = 4:K-3;                                         % away from the ends of both differences
            scale = max(max(abs(want(:, idx)), [], 'all'), realmin);
            tc.verifyBoundAndRecord("proof-mass coupling vs a twice-differentiated analytic rotation", ...
                max(abs(got(:, idx) - want(:, idx)), [], 'all') / scale, 1e-5, "-", ...
                "a centred second difference truncates at (omega dt)^2/12 ~ 3e-7 here; the bound is the oracle's error, not the model's");
            % The Euler term must be present, or a sign error in it would hide behind the
            % centrifugal term, which is even in omega and survives any sign convention.
            sEuler = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), dr);
            onlyEuler = sEuler.observe(struct('a', zeros(3, K), 'omega', zeros(3, K), 'omegaDot', wd)).value;
            tc.verifyBoundAndRecord("the Euler term is a real part of the reference case", ...
                max(vecnorm(want(:, idx))) / max(max(vecnorm(onlyEuler(:, idx))), realmin), 1e6, "-", ...
                "recorded so the reference case cannot silently degenerate into a pure spin, where omegaDot vanishes and half the model is untested");

            % The relation to referenceSystem's angular velocity is RECORDED, not assumed: that
            % unit carries a configurable sense and its rotation runs the other way, so which of
            % its settings lands on skew(omega) = R' Rdot is a fact about the host's convention.
            Rdot = zeros(3, 3, K);
            Rdot(:, :, 2:K-1) = (R(:, :, 3:K) - R(:, :, 1:K-2)) / (2 * dt);
            Rdot(:, :, 1) = Rdot(:, :, 2); Rdot(:, :, K) = Rdot(:, :, K-1);
            wRefBody = referenceSystem.angularVelocityFromRotationRate(R, Rdot, 'Resolution', "reference");
            rel = wRefBody(:, idx) ./ w(:, idx);
            tc.verifyBoundAndRecord("|referenceSystem omega| vs the sensor's, same magnitude", ...
                max(abs(abs(wRefBody(:, idx)) - abs(w(:, idx))), [], 'all') / max(max(abs(w(:, idx)), [], 'all'), realmin), ...
                1e-12, "-", ...
                sprintf("magnitudes must match whatever the sense; the SIGN ratio measured here is %+.1f, which is the convention fact O-32 asks for", median(rel(isfinite(rel)), 'all')));

            % Magnitude on the real arc: the reason this term is modelled at all.
            sReal = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), dr);
            yReal = sReal.observe(struct('a', zeros(3, numel(tc.ep.sec)), 'omega', tc.w, 'omegaDot', tc.wd)).value;
            tc.verifyBoundAndRecord("proof-mass coupling magnitude on the real GRACE-C arc", ...
                max(vecnorm(yReal(:, 4:end-3))), 1e-7, "m/s^2", ...
                "a few hundred micrometres at the yaw-steering rate: order 1e-9, above the 1e-10 noise floor and below the drag");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function ElectrostaticPrecision(tc)
            % The partials against central differences of observe itself. Every parameter is
            % perturbed in its own natural size: a bias in m/s^2, a scale dimensionless, an offset
            % in metres, so one step size cannot be right for all fifteen.
            b0 = [1e-6; -3e-6; 2e-7]; S0 = [1.05 0.01 -0.02; 0.005 0.96 0.008; -0.01 0.002 1.01];
            dr0 = [2.34e-4; -1.1e-4; 3.0e-4];
            dyn = tc.scene(25);
            s = NonGravForces.sensor.Electrostatic(b0, S0, dr0);
            [J, names] = s.partials(dyn);
            theta = [b0; S0(:); dr0];
            % THE MODEL IS EXACTLY LINEAR in all fifteen parameters, so a central difference is
            % exact for ANY step and the only error is cancellation -- which a LARGER step reduces,
            % the opposite of the usual rule. Steps of 1e-9/1e-6/1e-7 measured 2.2e-8 on the runner
            % (2026-09-08), pure cancellation against a 1e-6 bias; these steps put the difference
            % three orders above it.
            step = [1e-4 * ones(3, 1); 1e-2 * ones(9, 1); 1e-3 * ones(3, 1)];
            worst = 0;
            for p = 1:15
                tp = theta; tp(p) = tp(p) + step(p);
                tm = theta; tm(p) = tm(p) - step(p);
                yp = tNgfSensor.build(tp).observe(dyn).value;
                ym = tNgfSensor.build(tm).observe(dyn).value;
                num = reshape((yp - ym) / (2 * step(p)), [], 1);
                ana = J(:, p);
                sc = max(max(abs(num)), realmin);
                worst = max(worst, max(abs(ana - num)) / sc);
            end
            tc.verifyBoundAndRecord("Electrostatic partials vs central differences (15 parameters)", ...
                worst, 1e-9, "-", ...
                "exact linearity plus a step chosen to keep the difference well above eps of the reading: what is left is rounding, not truncation");
            tc.verifyEqual(numel(names), 15);
        end
    end

    % ================================================================ vectorization
    methods (Test, TestTags = {'vectorization'})
        function ElectrostaticVectorization(tc)
            % K epochs at once and K epochs one at a time are the same numbers, EXCEPT where the
            % filter is on: a filtered sample depends on its neighbours by construction, so a
            % single epoch cannot reproduce it and the mask says so. Both halves are asserted.
            s = NonGravForces.sensor.Electrostatic([1e-6; -3e-6; 2e-7], diag([1.05 0.96 1.01]), [2e-4; 0; -1e-4]);
            dyn = tc.scene(20);
            Y = s.observe(dyn).value;
            worst = 0;
            for k = 1:20
                one = struct('a', dyn.a(:, k), 'omega', dyn.omega(:, k), 'omegaDot', dyn.omegaDot(:, k));
                worst = max(worst, max(abs(s.observe(one).value - Y(:, k))));
            end
            tc.verifyBoundAndRecord("Electrostatic, one epoch at a time vs the whole arc", worst, 0, "m/s^2", ...
                "without a filter every sample is independent, so the batch may not change a single one");
            sf = NonGravForces.sensor.Electrostatic([1e-6; -3e-6; 2e-7], diag([1.05 0.96 1.01]), [2e-4; 0; -1e-4], ones(1, 5) / 5);
            yf = sf.observe(dyn);
            tc.verifyFalse(all(yf.mask(:)), 'a filtered arc must mark its edges invalid');
            tc.verifyEqual(nnz(~yf.mask(1, :)), 4, 'a 5-tap kernel invalidates two samples at each end');
            % A partial matrix must be filtered exactly as the observation is.
            Jf = sf.partials(dyn);
            tc.verifyBoundAndRecord("filtered partials vs the filtered observation, bias column", ...
                max(abs(Jf(1:3:end, 1) - conv(ones(1, 20), ones(1, 5) / 5, 'same')')), 1e-16, "-", ...
                "the bias enters every sample as 1, so its filtered column is the kernel's step response");
        end
    end

    % ==================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function ElectrostaticRobustness(tc)
            dyn = tc.scene(20);
            % A zero offset removes the coupling entirely: the reading is bias plus scaled input.
            s0 = NonGravForces.sensor.Electrostatic([1e-6; 0; 0], diag([1.05 0.96 1.01]), zeros(3, 1));
            y0 = s0.observe(dyn).value;
            tc.verifyBoundAndRecord("zero proof-mass offset removes the coupling", ...
                max(abs(y0 - ([1e-6; 0; 0] + diag([1.05 0.96 1.01]) * dyn.a)), [], 'all'), 0, "m/s^2", ...
                "with dr = 0 both cross products vanish identically, not approximately");
            % A proof mass ON the rotation axis of a pure spin feels no coupling: omega x dr = 0
            % and omegaDot is parallel to it too. This is a rigid-body identity, not a tolerance.
            K = 12; ax = [0; 0; 1];
            spin = struct('a', zeros(3, K), 'omega', repmat(2e-3 * ax, 1, K), 'omegaDot', repmat(1e-6 * ax, 1, K));
            sAx = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), 5e-4 * ax);
            tc.verifyBoundAndRecord("a proof mass on the spin axis feels nothing", ...
                max(abs(sAx.observe(spin).value), [], 'all'), 0, "m/s^2", ...
                "omega x dr and omegaDot x dr are both zero when dr is parallel to the axis");
            % One epoch is a legal arc; a filter longer than the arc leaves nothing valid.
            one = struct('a', dyn.a(:, 1), 'omega', dyn.omega(:, 1), 'omegaDot', dyn.omegaDot(:, 1));
            tc.verifySize(s0.observe(one).value, [3 1]);
            sLong = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), zeros(3, 1), ones(1, 41) / 41);
            tc.verifyFalse(any(sLong.observe(one).mask(:)), 'a 41-tap kernel on one sample leaves nothing valid');
            % The gravity gradient enters where the equation puts it and nowhere else.
            T = repmat(diag([2.7e-6, -1.35e-6, -1.35e-6]), 1, 1, 20);       % order GM/r^3 at 490 km
            dg = dyn; dg.gravityGradient = T;
            sG = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), [2e-4; 0; -1e-4]);
            yg = sG.observe(dg);
            tc.verifyEqual(yg.info.gradient, "supplied");
            tc.verifyBoundAndRecord("the gradient term equals -T dr exactly", ...
                max(abs(yg.value - sG.observe(dyn).value + reshape(pagemtimes(T, [2e-4; 0; -1e-4]), 3, 20)), [], 'all'), ...
                1e-21, "m/s^2", "one matrix-vector product, subtracted where the observation equation subtracts it; eps(1e-7) is 1.3e-23, so this is the cancellation floor with room");
        end
    end

    % =================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function ElectrostaticMetamorphic(tc)
            dyn = tc.scene(20);
            b = [1e-6; -3e-6; 2e-7]; S = diag([1.05 0.96 1.01]); dr = [2e-4; 0; -1e-4];
            s = NonGravForces.sensor.Electrostatic(b, S, dr);
            y = s.observe(dyn).value;
            % 1. The bias is additive and nothing else: changing it shifts every sample by the same
            %    vector, whatever the dynamics.
            s2 = NonGravForces.sensor.Electrostatic(b + [1e-5; 2e-5; -3e-5], S, dr);
            tc.verifyBoundAndRecord("a bias change is a rigid shift of the whole series", ...
                max(abs((s2.observe(dyn).value - y) - repmat([1e-5; 2e-5; -3e-5], 1, 20)), [], 'all'), ...
                1e-20, "m/s^2", "bias enters once, outside everything else; eps(1e-5) is 1.7e-21");
            % 2. Doubling the scale doubles the part that is not bias.
            s3 = NonGravForces.sensor.Electrostatic(b, 2 * S, dr);
            % The bound is COMPUTED, not guessed. Subtracting a 3e-6 bias from a 1e-6 reading to
            % recover a 1e-7 signal cancels, and the floor is eps of the larger term over the
            % smaller one. Writing 1e-15 here from an estimate reddened the gate at 1.4e-15
            % (2026-09-08); the floor below is 4e-15 and the factor of ten is the margin.
            floorRel = 10 * eps(max(abs([y(:); b])) ) / max(abs(y - b), [], 'all');
            tc.verifyBoundAndRecord("doubling the scale doubles the non-bias part", ...
                max(abs((s3.observe(dyn).value - b) - 2 * (y - b)), [], 'all') / max(abs(y - b), [], 'all'), ...
                floorRel, "-", "the scale multiplies the sensed acceleration and only that; the bound is the cancellation floor of removing a bias three times the signal");
            % 3. The coupling is linear in dr: two offsets added give the two readings added,
            %    minus the bias counted twice.
            d1 = [3e-4; 0; 0]; d2 = [0; -2e-4; 1e-4];
            f = @(d) NonGravForces.sensor.Electrostatic(zeros(3, 1), S, d).observe(dyn).value;
            tc.verifyBoundAndRecord("the coupling is linear in the proof-mass offset", ...
                max(abs(f(d1 + d2) - (f(d1) + f(d2) - f(zeros(3, 1)))), [], 'all') / max(abs(f(d1)), [], 'all'), ...
                1e-14, "-", "both cross-product terms are linear operators acting on dr");
            % 4. Reversing the rotation reverses the Euler term and leaves the centrifugal one:
            %    omega -> -omega flips omega x (omega x dr) not at all and omegaDot x dr not at all
            %    unless omegaDot flips too, so only flipping BOTH must leave the centrifugal part
            %    and reverse the Euler part. Checked by isolating each with a null input.
            sc = NonGravForces.sensor.Electrostatic(zeros(3, 1), eye(3), dr);
            zeroA = struct('a', zeros(3, 20), 'omega', dyn.omega, 'omegaDot', zeros(3, 20));
            flipW = zeroA; flipW.omega = -dyn.omega;
            tc.verifyBoundAndRecord("the centrifugal term is even in omega", ...
                max(abs(sc.observe(flipW).value - sc.observe(zeroA).value), [], 'all'), 1e-21, "m/s^2", ...
                "omega x (omega x dr) contains omega twice, so its sign cannot depend on omega's");
            zeroW = struct('a', zeros(3, 20), 'omega', zeros(3, 20), 'omegaDot', dyn.omegaDot);
            flipD = zeroW; flipD.omegaDot = -dyn.omegaDot;
            tc.verifyBoundAndRecord("the Euler term is odd in omegaDot", ...
                max(abs(sc.observe(flipD).value + sc.observe(zeroW).value), [], 'all'), 1e-21, "m/s^2", ...
                "omegaDot x dr contains omegaDot once");
        end
    end

    % ========================================================================= speed
    methods (Test, TestTags = {'speed'})
        function ElectrostaticSpeed(tc)
            % The reading is vectorised over epochs and the partials are a loop over them, which is
            % the honest shape for a 15-column Jacobian. What must hold is that neither is
            % quadratic: ten times the epochs must not cost a hundred times as much.
            s = NonGravForces.sensor.Electrostatic([1e-6; -3e-6; 2e-7], diag([1.05 0.96 1.01]), [2e-4; 0; -1e-4]);
            small = tc.scene(6); big = tc.scene(60);
            s.observe(big); s.partials(big);                      % warm
            tObsBig = tNgfSensor.timeIt(@() s.observe(big));
            tObsSmall = tNgfSensor.timeIt(@() s.observe(small));
            tParBig = tNgfSensor.timeIt(@() s.partials(big));
            tParSmall = tNgfSensor.timeIt(@() s.partials(small));
            rObs = tObsBig / max(tObsSmall, eps);
            rPar = tParBig / max(tParSmall, eps);
            tc.verifyBoundAndRecord("Electrostatic observe, 60 epochs / 6 epochs", rObs, 10, "-", ...
                "linear in K would be 10 and a fixed overhead pulls it below; above 10 would mean the work is not vectorised");
            tc.verifyBoundAndRecord("Electrostatic partials, 60 epochs / 6 epochs", rPar, 20, "-", ...
                "a per-epoch loop of fixed work is linear; the allowance is generous because the loop body is tiny and the timer is not");
        end
    end

    methods (Static)
        function r = rmsOf(x)
            %RMSOF Root mean square. rms() lives in the Signal Processing Toolbox, which C2 rules  % toolboxScan:allow
            % out, and a test may not depend on a product the toolbox itself refuses to need.
            r = sqrt(mean(x.^2));
        end
        function el = timeIt(f)
            %TIMEIT The median of five runs: one scheduling hiccup should not decide a ratio (T5).
            t = zeros(1, 5);
            for k = 1:5, t0 = tic; f(); t(k) = toc(t0); end
            el = median(t);
        end
        function s = build(theta)
            %BUILD A sensor from the 15-parameter vector, in the order partials documents.
            s = NonGravForces.sensor.Electrostatic(theta(1:3), reshape(theta(4:12), 3, 3), theta(13:15));
        end
    end

    % =========================================================== AtomInterferometer
    methods (Static)
        function s = ai(varargin)
            %AI A rubidium two-photon instrument on a 10 Hz grid, the CARIOQA-class configuration.
            s = NonGravForces.sensor.AtomInterferometer(1.610e7, 0.5, 2, [1; 0; 0], ...
                'SampleRate', 10, varargin{:});
        end
    end
    methods (Test, TestTags = {'contract'})
        function AtomInterferometerContract(tc)
            s = tNgfSensor.ai();
            tc.verifyEqual(s.numParameters(), 2);
            tc.verifyEqual(s.name, "atomInterferometer");
            tc.verifyTrue(isa(s, 'NonGravForces.sensor.Sensor'));
            K = 200;
            dyn = struct('a', [3e-7 * ones(1, K); zeros(2, K)]);
            y = s.observe(dyn);
            tc.verifySize(y.value, [1 K]); tc.verifySize(y.mask, [1 K]);
            tc.verifyEqual(y.frame, "beam");
            tc.verifyEqual(y.info.deadTime_s, 1, 'cycle 2 s minus 2T = 1 s of interrogation leaves 1 s blind');
            [J, names] = s.partials(dyn);
            tc.verifySize(J, [K 2]);
            tc.verifyEqual(names, ["phaseBias", "gainError"]);
            tc.verifySize(s.noiseModel([1e-3 1e-2]), [1 2]);

            A = @NonGravForces.sensor.AtomInterferometer;
            tc.verifyErrorId(@() A(1.6e7, 0.5, 0.5, [1; 0; 0]), 'NonGravForces:AtomInterferometer:inputRange');
            tc.verifyErrorId(@() A(-1, 0.5, 2, [1; 0; 0]), 'NonGravForces:AtomInterferometer:inputRange');
            tc.verifyErrorId(@() A(1.6e7, 0.5, 2, [0; 0; 0]), 'NonGravForces:AtomInterferometer:inputRange');
            tc.verifyErrorId(@() A(1.6e7, 0.5, 2, [1; 0]), 'NonGravForces:AtomInterferometer:inputShape');
            tc.verifyErrorId(@() A(NaN, 0.5, 2, [1; 0; 0]), 'NonGravForces:AtomInterferometer:numericNaN');
            tc.verifyErrorId(@() s.observe(struct('b', 1)), 'NonGravForces:AtomInterferometer:inputDyn');
            tc.verifyErrorId(@() s.observe(struct('a', ones(2, 5))), 'NonGravForces:AtomInterferometer:inputShape');
            bad = dyn; bad.a(1) = Inf;
            tc.verifyErrorId(@() s.observe(bad), 'NonGravForces:AtomInterferometer:numericNaN');
            sr = tNgfSensor.ai('Rotation', true);
            tc.verifyErrorId(@() sr.observe(dyn), 'NonGravForces:AtomInterferometer:inputDyn');
            tc.verifyErrorId(@() sr.observe(struct('a', dyn.a, 'omega', ones(3, 5), 'velocity', ones(3, K))), ...
                'NonGravForces:AtomInterferometer:inputDyn');
            tc.verifyErrorId(@() s.noiseModel(0), 'NonGravForces:AtomInterferometer:inputShape');
        end
    end
    methods (Test, TestTags = {'reference'})
        function AtomInterferometerReference(tc)
            % THE IDENTITY IS THE ORACLE: for a CONSTANT acceleration the triangular weighting has
            % unit sum, so the phase must be exactly k_eff a T^2 -- the textbook result, and the
            % thing that would break first if the kernel were normalised wrongly.
            s = tNgfSensor.ai();
            a0 = 3.2e-7; K = 400;
            y = s.observe(struct('a', [a0 * ones(1, K); zeros(2, K)]));
            want = 1.610e7 * a0 * 0.5^2;
            tc.verifyBoundAndRecord("constant acceleration gives exactly k_eff a T^2", ...
                max(abs(y.value(y.mask) - want)) / want, 1e-15, "-", ...
                "the kernel is normalised to unit sum, so a constant passes through untouched");
            % The Coriolis phase against an independent evaluation of -2 k_eff T^2 <n, omega x v>,
            % with omega and v constant so the weighting is again exact.
            sr = tNgfSensor.ai('Rotation', true);
            w0 = [0; 0; 1.1e-3]; v0 = [0; 7600; 0];
            dyn = struct('a', zeros(3, K), 'omega', repmat(w0, 1, K), 'velocity', repmat(v0, 1, K));
            yr = sr.observe(dyn);
            wantCor = -2 * 1.610e7 * 0.5^2 * dot([1; 0; 0], cross(w0, v0));
            tc.verifyBoundAndRecord("the Coriolis phase matches -2 k_eff T^2 <n, omega x v>", ...
                max(abs(yr.value(yr.mask) - wantCor)) / abs(wantCor), 1e-15, "-", ...
                "an independent evaluation of the same closed form on constant inputs");
            % THE NUMBER IS 6.73e7 RADIANS, and the first version of this test said "thousands"
            % from an estimate that was four orders out. Worth getting right, because the correct
            % figure is the whole argument: 2 |omega x v| = 16.7 m/s^2 is an apparent acceleration
            % larger than gravity, 1.07e7 fringes across the readout. An orbiting interferometer
            % that does not counter-rotate its Raman mirror does not measure a small bias on top
            % of the drag -- it has no fringe left to read.
            tc.verifyBoundAndRecord("Coriolis phase for a GRACE-like rate and speed", ...
                abs(wantCor), 1e9, "rad", ...
                sprintf("6.73e7 rad = 1.07e7 fringes = an apparent 16.7 m/s^2, against a drag of 1e-7; the bound is a smoke ceiling, the number is the point"));
            tc.verifyBoundAndRecord("Coriolis phase as an equivalent acceleration", ...
                abs(wantCor) / (1.610e7 * 0.5^2), 1e3, "m/s^2", ...
                "16.7 m/s^2: the rotation term is not a correction to the measurement, it is larger than everything the instrument is there to see");
        end
    end
    methods (Test, TestTags = {'precision'})
        function AtomInterferometerPrecision(tc)
            % The two parameters enter linearly, so central differences with a LARGE step are exact
            % and the only error is cancellation -- the lesson from the Electrostatic partials.
            K = 300;
            t = (0:K-1) / 10;
            dyn = struct('a', [3e-7 + 5e-8 * sin(2 * pi * 0.05 * t); zeros(2, K)]);
            base = [0.3; 0.02];                       % phaseBias, gainError
            s = tNgfSensor.ai('PhaseBias', base(1), 'GainError', base(2));
            J = s.partials(dyn);
            step = [1e-2; 1e-3];
            worst = 0;
            for p = 1:2
                th = base; th(p) = th(p) + step(p);
                sp = tNgfSensor.ai('PhaseBias', th(1), 'GainError', th(2));
                th = base; th(p) = th(p) - step(p);
                sm = tNgfSensor.ai('PhaseBias', th(1), 'GainError', th(2));
                num = (sp.observe(dyn).value - sm.observe(dyn).value)' / (2 * step(p));
                sc = max(max(abs(num)), realmin);
                worst = max(worst, max(abs(J(:, p) - num)) / sc);
            end
            tc.verifyBoundAndRecord("AtomInterferometer partials vs central differences", worst, 1e-10, "-", ...
                "both parameters are exactly linear, so a large step leaves only rounding");
            % The full phase, the wrapped phase and the fringe integer must reconstruct exactly.
            y = s.observe(dyn);
            tc.verifyBoundAndRecord("wrapped phase plus 2 pi times the fringe reconstructs the phase", ...
                max(abs(y.value - (y.info.wrapped + 2 * pi * y.info.fringe))), 1e-12, "rad", ...
                "the split is exact by construction; a caller resolving the ambiguity must get back what was modelled");
            tc.verifyTrue(all(abs(y.info.wrapped) <= pi + eps), 'the wrapped phase must lie in (-pi, pi]');
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function AtomInterferometerVectorization(tc)
            % THE KERNEL HAS FINITE SUPPORT, so an interior sample cannot depend on data outside
            % its own window. That is an exact statement -- perturb a distant sample and the
            % interior phase may not move at all -- and it is stronger than comparing two arcs,
            % which would only measure how the edges are handled.
            s = tNgfSensor.ai();
            K = 300; h = 5;                              % T * fs = 5 samples each side
            a = [3e-7 * ones(1, K); zeros(2, K)];
            y0 = s.observe(struct('a', a));
            aP = a; aP(1, 250) = aP(1, 250) + 1e-6;      % a spike far from sample 100
            y1 = s.observe(struct('a', aP));
            % Relative, not exact zero: conv may take an FFT route for a long input, and whether
            % it does is MATLAB's business, not this claim's. The claim is locality, and locality
            % survives either route at the rounding level.
            d = abs(y1.value - y0.value) / max(abs(y0.value));
            tc.verifyBoundAndRecord("a sample outside the window cannot move the phase", ...
                max(d(1:200)), 1e-14, "-", ...
                "the triangular kernel is supported on 2T+1 samples and nothing else enters");
            tc.verifyGreaterThan(d(250), 1e-9, 'and inside the window it must move');
            tc.verifyEqual(nnz(d > 1e-14), 2 * h + 1, ...
                'exactly the window width of samples may change');
        end
    end
    methods (Test, TestTags = {'robustness'})
        function AtomInterferometerRobustness(tc)
            K = 200;
            a = [3e-7 * ones(1, K); zeros(2, K)];
            % DEAD TIME: one reading per cycle and nothing between. 2 s at 10 Hz is every 20th.
            s = tNgfSensor.ai();
            y = s.observe(struct('a', a));
            idx = find(y.mask);
            tc.verifyEqual(unique(diff(idx)), 20, 'a 2 s cycle at 10 Hz leaves one valid sample in twenty');
            % An interrogation window shorter than one sample is not weighted at all, and says so
            % rather than silently convolving with a degenerate kernel.
            sFast = NonGravForces.sensor.AtomInterferometer(1.610e7, 0.02, 1, [1; 0; 0], 'SampleRate', 10);
            yF = sFast.observe(struct('a', a));
            tc.verifyBoundAndRecord("an unresolved window leaves the acceleration unweighted", ...
                max(abs(yF.info.weighted - a(1, :))), 0, "m/s^2", ...
                "T * fs below one sample: the window IS a sample, and pretending otherwise would invent a kernel");
            % An arc shorter than the window leaves nothing valid.
            yShort = s.observe(struct('a', a(:, 1:6)));
            tc.verifyFalse(any(yShort.mask), 'six samples cannot hold an 11-sample window');
            % Rotation with the velocity along the rotation axis produces no Coriolis phase at all.
            sr = tNgfSensor.ai('Rotation', true);
            ax = [0; 0; 1];
            yPar = sr.observe(struct('a', zeros(3, K), 'omega', repmat(2e-3 * ax, 1, K), ...
                'velocity', repmat(7600 * ax, 1, K)));
            tc.verifyBoundAndRecord("velocity along the rotation axis gives no Coriolis phase", ...
                max(abs(yPar.value)), 0, "rad", "omega x v vanishes identically when the two are parallel");
            % Above the cycle Nyquist the noise model reports Inf rather than a number.
            asd = s.noiseModel([1e-3, 1, 10]);
            tc.verifyTrue(isfinite(asd(1)) && all(isinf(asd(2:3))), ...
                'nothing is measured faster than half the cycle rate, and the ASD must say so');
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function AtomInterferometerMetamorphic(tc)
            K = 400;
            a0 = 3e-7;
            a = [a0 * ones(1, K); zeros(2, K)];
            s = tNgfSensor.ai();
            y = s.observe(struct('a', a));
            % 1. The phase is linear in the acceleration.
            y2 = s.observe(struct('a', 2 * a));
            tc.verifyBoundAndRecord("the phase doubles with the acceleration", ...
                max(abs(y2.value - 2 * y.value)) / max(abs(2 * y.value)), 1e-15, "-", ...
                "k_eff T^2 multiplies a linearly and nothing else touches it");
            % 2. It goes as T^2, which is the whole reason these instruments are sensitive.
            sT = NonGravForces.sensor.AtomInterferometer(1.610e7, 1.0, 4, [1; 0; 0], 'SampleRate', 10);
            yT = sT.observe(struct('a', a));
            tc.verifyBoundAndRecord("doubling the interrogation time quadruples the phase", ...
                abs(max(yT.value) - 4 * max(y.value)) / (4 * max(y.value)), 1e-15, "-", ...
                "T enters squared; for a constant acceleration the weighting is exact at both T");
            % 3. Reversing the beam reverses the sign and nothing else.
            sB = NonGravForces.sensor.AtomInterferometer(1.610e7, 0.5, 2, [-1; 0; 0], 'SampleRate', 10);
            yB = sB.observe(struct('a', a));
            tc.verifyBoundAndRecord("reversing the beam reverses the phase", ...
                max(abs(yB.value + y.value)) / max(abs(y.value)), 1e-15, "-", ...
                "the projection <n, a> is the only place the direction enters");
            % 4. The Coriolis phase is odd in omega: reversing the rotation reverses it.
            sr = tNgfSensor.ai('Rotation', true);
            w0 = repmat([0; 0; 1.1e-3], 1, K); v0 = repmat([0; 7600; 0], 1, K);
            yP = sr.observe(struct('a', zeros(3, K), 'omega', w0, 'velocity', v0));
            yM = sr.observe(struct('a', zeros(3, K), 'omega', -w0, 'velocity', v0));
            tc.verifyBoundAndRecord("the Coriolis phase is odd in omega", ...
                max(abs(yP.value + yM.value)) / max(abs(yP.value)), 1e-15, "-", ...
                "omega x v contains omega once");
            % 5. A phase bias is additive and reaches nothing else.
            sBias = tNgfSensor.ai('PhaseBias', 1.234);
            tc.verifyBoundAndRecord("a phase bias is a rigid shift", ...
                max(abs((sBias.observe(struct('a', a)).value - y.value) - 1.234)), 1e-12, "rad", ...
                "it enters outside the gain and outside the weighting");
        end
    end
    methods (Test, TestTags = {'speed'})
        function AtomInterferometerSpeed(tc)
            % One convolution and a handful of vector operations: linear in K.
            s = tNgfSensor.ai();
            small = struct('a', [3e-7 * ones(1, 200); zeros(2, 200)]);
            big = struct('a', [3e-7 * ones(1, 2000); zeros(2, 2000)]);
            s.observe(big);
            r = tNgfSensor.timeIt(@() s.observe(big)) / max(tNgfSensor.timeIt(@() s.observe(small)), eps);
            tc.verifyBoundAndRecord("AtomInterferometer observe, 2000 epochs / 200 epochs", r, 10, "-", ...
                "linear in K would be 10 and a fixed overhead pulls it below; the convolution is the only term that grows");
        end
    end

    % ==================================================================== DragFree
    methods (Static)
        function [s, t] = df(K, varargin)
            %DF A loop with a 0.1 Hz bandwidth on a 10 Hz grid: 1.6 s time constant, 80 to settle.
            s = NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', 10, varargin{:});
            t = (0:K-1) / 10;
        end
    end
    methods (Test, TestTags = {'contract'})
        function DragFreeContract(tc)
            [s, t] = tNgfSensor.df(500);
            tc.verifyEqual(s.numParameters(), 4);
            tc.verifyEqual(s.name, "dragFree");
            tc.verifyTrue(isa(s, 'NonGravForces.sensor.Sensor'));
            dyn = struct('a', [3e-7 * ones(size(t)); 1e-7 * sin(t); zeros(size(t))]);
            y = s.observe(dyn);
            tc.verifySize(y.value, [3 500]); tc.verifySize(y.mask, [3 500]);
            tc.verifyEqual(y.frame, "body");
            tc.verifyFalse(any(y.mask(:, 1)), 'the start-up transient must be masked, not delivered');
            tc.verifyTrue(all(y.mask(:, end)));
            tc.verifyEqual(y.info.settleSamples, ceil(5 / (2 * pi * 0.1) * 10));
            [J, names] = s.partials(dyn);
            tc.verifySize(J, [1500 4]);
            tc.verifyEqual(names, ["bias_x", "bias_y", "bias_z", "gainError"]);
            tc.verifySize(s.noiseModel([1e-3 1e-2]), [3 2]);

            D = @NonGravForces.sensor.DragFree;
            tc.verifyErrorId(@() D(0, 1, 'SampleRate', 10), 'NonGravForces:DragFree:inputRange');
            tc.verifyErrorId(@() D(0.1, 0, 'SampleRate', 10), 'NonGravForces:DragFree:inputRange');
            tc.verifyErrorId(@() D(0.1, 1.5, 'SampleRate', 10), 'NonGravForces:DragFree:inputRange');
            tc.verifyErrorId(@() D(0.1, 1, 'SampleRate', 10, 'Authority', 0), 'NonGravForces:DragFree:inputRange');
            tc.verifyErrorId(@() D(0.1, 1, 'SampleRate', 0.1), 'NonGravForces:DragFree:inputRange');
            tc.verifyErrorId(@() D(0.1, 1, 'SampleRate', 10, 'Bias', [1; 2]), 'NonGravForces:DragFree:inputShape');
            tc.verifyErrorId(@() D(NaN, 1, 'SampleRate', 10), 'NonGravForces:DragFree:numericNaN');
            tc.verifyErrorId(@() s.observe(struct('b', 1)), 'NonGravForces:DragFree:inputDyn');
            tc.verifyErrorId(@() s.observe(struct('a', ones(2, 4))), 'NonGravForces:DragFree:inputShape');
            bad = dyn; bad.a(1) = NaN;
            tc.verifyErrorId(@() s.observe(bad), 'NonGravForces:DragFree:numericNaN');
            tc.verifyErrorId(@() s.noiseModel(-1), 'NonGravForces:DragFree:inputShape');
        end
    end
    methods (Test, TestTags = {'reference'})
        function DragFreeReference(tc)
            % THE ORACLE IS THE FREQUENCY RESPONSE, derived on the unit circle from the same
            % difference-equation coefficients but evaluated analytically -- a different route
            % through the same definition, which is what a reference for a filter can be. A
            % steady-state sinusoid in must come out scaled by exactly |S(f)|.
            K = 4000;
            [s, t] = tNgfSensor.df(K);
            settled = 1001:K;                             % 60 time constants in: the transient is gone
            worst = 0;
            for f = [0.01, 0.1, 1]                        % below, at, and above the corner
                x = [sin(2 * pi * f * t); zeros(2, K)];
                y = s.observe(struct('a', x));
                got = tNgfSensor.rmsOf(y.value(1, settled)) / tNgfSensor.rmsOf(x(1, settled));
                want = s.magnitudeResponse(f);
                worst = max(worst, abs(got - want) / want);
            end
            tc.verifyBoundAndRecord("DragFree rejection vs the analytic |S(f)|, three decades", ...
                worst, 1e-6, "-", ...
                "an integer number of periods after 60 time constants; what is left is the residue of the transient, not a modelling difference");
            % The asymptotes are the whole point, so they are asserted as such: full rejection at
            % DC, none well above the corner, and 1/sqrt(2) at the corner for one section.
            tc.verifyBoundAndRecord("DragFree gain at the bandwidth is 1/sqrt(2) for one section", ...
                abs(s.magnitudeResponse(0.1) - 1 / sqrt(2)), 1e-3, "-", ...
                "the bilinear transform warps the corner slightly; 1e-3 is that warping at f_c = fs/100, not a slack bound");
            tc.verifyEqual(s.magnitudeResponse(0), 0, 'a high-pass rejects a constant completely');
            % And a genuinely constant drag is driven to zero, which is what drag-free MEANS.
            yDC = s.observe(struct('a', repmat([3e-7; 0; 0], 1, K)));
            tc.verifyBoundAndRecord("a constant drag is rejected to nothing", ...
                max(abs(yDC.value(:, settled)), [], 'all') / 3e-7, 1e-12, "-", ...
                "the loop has all the time in the world against a constant; what remains is the decayed transient");
        end
    end
    methods (Test, TestTags = {'precision'})
        function DragFreePrecision(tc)
            % Both parameters are linear, so a large step makes the central difference exact.
            K = 600;
            [~, t] = tNgfSensor.df(K);
            dyn = struct('a', [3e-7 + 1e-7 * sin(2 * pi * 0.3 * t); 5e-8 * cos(2 * pi * 0.05 * t); zeros(1, K)]);
            base = [1e-9; 2e-9; -1e-9; 0.03];
            mk = @(th) NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', 10, ...
                'Bias', th(1:3), 'GainError', th(4));
            sBase = mk(base);
            J = sBase.partials(dyn);
            step = [1e-7 * ones(3, 1); 1e-2];
            worst = 0;
            for p = 1:4
                tp = base; tp(p) = tp(p) + step(p);
                tm = base; tm(p) = tm(p) - step(p);
                sp = mk(tp); sm = mk(tm);
                yp = sp.observe(dyn); ym = sm.observe(dyn);
                num = reshape((yp.value - ym.value) / (2 * step(p)), [], 1);
                sc = max(max(abs(num)), realmin);
                worst = max(worst, max(abs(J(:, p) - num)) / sc);
            end
            tc.verifyBoundAndRecord("DragFree partials vs central differences", worst, 1e-9, "-", ...
                "bias and gain both enter linearly; the step is large enough that cancellation is not the limit");
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function DragFreeVectorization(tc)
            % The three axes go through the same filter but share no state. Perturbing one may not
            % touch another, exactly -- a shared accumulator would show here and nowhere else.
            K = 400;
            [s, t] = tNgfSensor.df(K);
            a = [3e-7 * ones(1, K); 1e-7 * sin(t); 2e-8 * cos(t)];
            y0 = s.observe(struct('a', a));
            aP = a; aP(2, :) = aP(2, :) + 4e-7 * sin(3 * t);
            y1 = s.observe(struct('a', aP));
            tc.verifyBoundAndRecord("perturbing the y axis leaves x and z untouched", ...
                max(abs(y1.value([1 3], :) - y0.value([1 3], :)), [], 'all'), 0, "m/s^2", ...
                "one filter state per axis; a shared accumulator would leak here");
            tc.verifyGreaterThan(max(abs(y1.value(2, :) - y0.value(2, :))), 0);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function DragFreeRobustness(tc)
            K = 600;
            [~, t] = tNgfSensor.df(K);
            % SATURATION is the failure mode this class exists to show. A drag far above the
            % authority passes through minus what the thrusters can cancel, and no further.
            big = [2e-5 * ones(1, K); zeros(2, K)];
            sSat = NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', 10, 'Authority', 5e-6);
            ySat = sSat.observe(struct('a', big));
            settled = 401:K;
            tc.verifyTrue(all(ySat.info.saturated(settled)), 'a 2e-5 drag against 5e-6 of authority must saturate');
            tc.verifyBoundAndRecord("the saturated residual is the drag minus the authority", ...
                max(abs(vecnorm(ySat.value(:, settled)) - (2e-5 - 5e-6))) / (2e-5 - 5e-6), 1e-9, "-", ...
                "the command is clipped in its own direction and the excess is what reaches the payload");
            % An ideal loop (Authority Inf) never saturates and is linear again.
            sIdeal = NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', 10);
            yIdeal = sIdeal.observe(struct('a', big));
            tc.verifyFalse(any(yIdeal.info.saturated));
            tc.verifyBoundAndRecord("an ideal loop rejects the same constant to nothing", ...
                max(abs(yIdeal.value(:, settled)), [], 'all') / 2e-5, 1e-12, "-", ...
                "with infinite authority the constant is entirely inside the loop's reach");
            % Zero in, zero out, exactly: a filter with zero state and zero input invents nothing.
            yZero = sIdeal.observe(struct('a', zeros(3, K)));
            tc.verifyBoundAndRecord("zero drag gives exactly zero residual", ...
                max(abs(yZero.value), [], 'all'), 0, "m/s^2", "no bias, no input, no state: nothing to produce");
            % One epoch is legal, and it is entirely transient, so nothing is valid.
            yOne = sIdeal.observe(struct('a', [3e-7; 0; 0]));
            tc.verifySize(yOne.value, [3 1]);
            tc.verifyFalse(any(yOne.mask(:)), 'a single sample is all transient and must be masked');
            % A higher order settles longer and is masked longer, which the report must say.
            s3 = NonGravForces.sensor.DragFree(0.1, 3, 'SampleRate', 10);
            y3 = s3.observe(struct('a', [3e-7 * ones(1, K); zeros(2, K)]));
            yRef = sIdeal.observe(struct('a', zeros(3, 2)));
            tc.verifyEqual(y3.info.settleSamples, 3 * yRef.info.settleSamples);
            tc.verifyTrue(all(isfinite(y3.value(:))));
            % The noise model is the COMPLEMENT of the rejection: loudest where the loop is active.
            asd = s3.noiseModel([1e-3, 0.1, 10]);
            tc.verifyTrue(asd(1, 1) > asd(1, 2) && asd(1, 2) > asd(1, 3), ...
                'thruster noise is injected inside the bandwidth and not outside it');
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function DragFreeMetamorphic(tc)
            K = 600;
            [s, t] = tNgfSensor.df(K);
            a = [3e-7 * ones(1, K); 1e-7 * sin(2 * pi * 0.3 * t); zeros(1, K)];
            y = s.observe(struct('a', a));
            settled = 201:K;
            % 1. Below saturation the loop is linear: double the drag, double the residual.
            y2 = s.observe(struct('a', 2 * a));
            tc.verifyBoundAndRecord("the residual doubles with the drag below saturation", ...
                max(abs(y2.value(:, settled) - 2 * y.value(:, settled)), [], 'all') / max(abs(2 * y.value(:, settled)), [], 'all'), ...
                1e-14, "-", "an unsaturated loop is a linear filter and nothing else");
            % 2. Reversing the drag reverses the residual.
            yM = s.observe(struct('a', -a));
            tc.verifyBoundAndRecord("reversing the drag reverses the residual", ...
                max(abs(yM.value + y.value), [], 'all') / max(abs(y.value), [], 'all'), 1e-14, "-", ...
                "linearity again, in the sign");
            % 3. A bias is additive and reaches nothing else.
            sB = NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', 10, 'Bias', [1e-9; -2e-9; 3e-9]);
            yB = sB.observe(struct('a', a));
            % The bound is COMPUTED from the cancellation floor, not written. Recovering a 1e-9
            % shift by differencing two residuals of order 1e-7 cannot beat eps of the larger, and
            % 1e-24 written from an estimate reddened the gate at 4.5e-24 -- eps(1e-7) is 1.3e-23,
            % so the guess was below the floor it was measuring against.
            biasFloor = 10 * eps(max(abs(y.value), [], 'all'));   % y is the observe struct, not the array
            tc.verifyBoundAndRecord("a bias is a rigid shift of the residual", ...
                max(abs((yB.value - y.value) - [1e-9; -2e-9; 3e-9]), [], 'all'), ...
                biasFloor, "m/s^2", "it is added outside the filter and outside the gain; the bound is ten times eps of the residual it is differenced out of");
            % 4. MORE BANDWIDTH IS MORE REJECTION at a fixed low frequency, and that is the design
            %    trade the noise model pays for: the same widening lets more thruster noise in.
            fLow = 0.02;
            xs = [sin(2 * pi * fLow * t); zeros(2, K)];
            sNarrow = NonGravForces.sensor.DragFree(0.05, 1, 'SampleRate', 10);
            sWide = NonGravForces.sensor.DragFree(0.4, 1, 'SampleRate', 10);
            yN = sNarrow.observe(struct('a', xs));
            yW = sWide.observe(struct('a', xs));
            rN = tNgfSensor.rmsOf(yN.value(1, settled));
            rW = tNgfSensor.rmsOf(yW.value(1, settled));
            tc.verifyBoundAndRecord("a wider loop rejects a low frequency harder", rW / rN, 1, "-", ...
                "eight times the bandwidth against a 0.02 Hz signal: the ratio must be below one, and how far below is the design trade");
            tc.verifyGreaterThan(sWide.noiseModel(0.1), sNarrow.noiseModel(0.1), ...
                'and it pays for that with thruster noise, which is why the bandwidth is a choice');
        end
    end
    methods (Test, TestTags = {'speed'})
        function DragFreeSpeed(tc)
            % The filter is a per-sample recursion, so it is linear in K by construction. What is
            % measured is that nothing quadratic crept into the surrounding bookkeeping.
            s = NonGravForces.sensor.DragFree(0.1, 1, 'SampleRate', 10);
            small = struct('a', ones(3, 500) * 3e-7);
            big = struct('a', ones(3, 5000) * 3e-7);
            s.observe(big);
            r = tNgfSensor.timeIt(@() s.observe(big)) / max(tNgfSensor.timeIt(@() s.observe(small)), eps);
            tc.verifyBoundAndRecord("DragFree observe, 5000 epochs / 500 epochs", r, 20, "-", ...
                "linear in K would be 10; the allowance covers an interpreted loop whose per-sample cost is small next to the timer");
        end
    end

end
