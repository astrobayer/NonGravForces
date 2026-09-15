classdef tNgfModel < NgfTestCase
%TNGFMODEL Tests of NonGravForces.Model: the single point of access (LP-6, WP-14 to WP-17).
%
%   Oracles (none built here): the real 2020-01-01 GRACE-C day — the measured ACT1B
%   accelerations for the shape of the modelled signal, the SCA1B attitude, the GNV1B
%   state; the published cross-track thermal magnitude of Siemes et al. (2023); the
%   unchunked series for the chunked one, bitwise; the per-effect budget for the total,
%   bitwise; Richardson-extrapolated central differences for the analytic partials; and,
%   underneath all of it, the independent Python kernel mirror frozen at LP-5x.
%   Every bound carries the number measured on 2026-09-08.
%
%   WHAT THE ACT1B COMPARISON CAN AND CANNOT SHOW. The fixture is a 600-record
%   excerpt, ten minutes, and ACT1B is uncalibrated: its bias is 1e-5 m/s^2 on the
%   cross-track axis against a modelled 2e-8. Only the SHAPE is comparable, and only
%   where the modelled signal rises above the instrument's own drift over ten minutes.
%   Measured: correlation 0.980 cross-track and 0.994 radial, but 0.190 along-track,
%   where the measurement is thirty times the model and dominated by bias drift. The
%   test asserts the first two and records the third rather than pretending otherwise.
%
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
    properties
        root, eop, sw, ep, r, v, q, att, macro, M, Mfast
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            f = @(varargin) fullfile(tc.root, varargin{:});
            tc.eop = referenceSystem.readEop(f('data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            tc.sw = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [tc.ep, tc.r, tc.v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', tc.eop);
            [eps2, qq] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
            tc.att = NonGravForces.geom.AttitudeProvider(eps2, qq);
            tc.macro = NonGravForces.geom.graceFO("C");
            tc.M = NonGravForces.Model(tc.macro, NonGravForces.Options("standard"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            tc.Mfast = NonGravForces.Model(tc.macro, NonGravForces.Options("fast"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, ...
                'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));
            tc.q = tc.att.at(tc.slice(1:5));
        end
    end
    methods
        function e = slice(tc, idx)
            e = referenceSystem.daySecondToEpoch([tc.ep.day(idx) tc.ep.sec(idx)], tc.ep.scale);
        end
        function [e, rr, vv] = arc(tc, n)
            idx = 1:n; e = tc.slice(idx); rr = tc.r(:, idx); vv = tc.v(:, idx);
        end
        function c = correlation(~, x, y)
            %CORRELATION Pearson correlation in base MATLAB. `corr` is Statistics Toolbox, which
            % HANDOVER D.7 C2 keeps off every code path: it exists on the owner's host and not on
            % the CI runner, so using it turned a green local run into a red CI one on 2026-09-08.
            a = x(:) - mean(x(:)); b = y(:) - mean(y(:));
            c = (a' * b) / (norm(a) * norm(b));
        end
        function q = speedRatio(~, fcn, base)
            fcn(); base();
            rr = zeros(1, 3);
            for k = 1:3
                t1 = tic; fcn(); a = toc(t1); t2 = tic; base(); b = toc(t2); rr(k) = a / b;
            end
            q = median(rr);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function modelContract(tc)
            mm = tc.macro;
            M0 = NonGravForces.Model(mm, NonGravForces.Options("fast"));
            tc.verifyEqual(M0.options.preset, "fast");
            tc.verifyEqual(M0.provenance.unit, "Model");
            tc.verifyErrorId(@() NonGravForces.Model(struct()), 'NonGravForces:Model:inputMacro');
            tc.verifyErrorId(@() NonGravForces.Model(mm, struct()), 'NonGravForces:Model:inputOptions');
            [e, rr, vv] = tc.arc(3);
            [a, info] = tc.M.accel(e, rr, vv);
            tc.verifySize(a, [3 3]); tc.verifyEqual(info.frame, "body"); tc.verifyEqual(info.K, 3);
            tc.verifySize(info.drag, [3 3]); tc.verifySize(info.radiation, [3 3]); tc.verifySize(info.thermal, [3 3]);
            tc.verifySize(info.rho, [1 3]); tc.verifySize(info.mass_kg, [1 3]); tc.verifySize(info.temperature_K, [mm.n 3]);
            tc.verifyEqual(info.backends.density, "msis00"); tc.verifyEqual(info.backends.gsi, "sentman");
            [tau, ti] = tc.M.torque(e, rr, vv);
            tc.verifySize(tau, [3 3]); tc.verifySize(ti.bodyTorque, [3 3]);
            tc.verifyEqual(ti.bodyTorque, zeros(3, 3), 'no body torque is enabled in the standard preset');
            tc.verifyClass(tc.M.rhs(), 'function_handle');
            B = tc.M.budget(e, rr, vv);
            tc.verifyEqual(B.effect', ["drag" "radiation" "thermal" "total"]);
            T = tc.M.series(e, rr, vv, 'Quantities', ["accel" "torque"]);
            tc.verifyEqual(T.Properties.VariableNames, {'ax','ay','az','tx','ty','tz'});
            % refusals, each by its documented identifier
            tc.verifyErrorId(@() tc.M.accel(e, rr(:, 1:2), vv), 'NonGravForces:Model:inputState');
            tc.verifyErrorId(@() tc.M.accel(e, rr, vv, ones(4, 2)), 'NonGravForces:Model:inputState');
            bad = rr; bad(1) = NaN;
            tc.verifyErrorId(@() tc.M.accel(e, bad, vv), 'NonGravForces:Model:inputState');
            noAtt = NonGravForces.Model(mm, NonGravForces.Options("standard"), 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            tc.verifyErrorId(@() noAtt.accel(e, rr, vv), 'NonGravForces:Model:inputAttitude');
            noEop = NonGravForces.Model(mm, NonGravForces.Options("standard"), 'Attitude', tc.att, 'SpaceWeather', tc.sw);
            tc.verifyErrorId(@() noEop.accel(e, rr, vv), 'NonGravForces:Model:dependEop');
            noSw = NonGravForces.Model(mm, NonGravForces.Options("standard"), 'Attitude', tc.att, 'Eop', tc.eop);
            tc.verifyErrorId(@() noSw.accel(e, rr, vv), 'NonGravForces:Model:dependWeather');
            tc.verifyErrorId(@() tc.M.partials(e, rr, vv, 'WithRespectTo', "mass"), 'NonGravForces:Model:unsupported');
            tc.verifyErrorId(@() tc.M.partials(e, rr, vv, 'Method', "complexStep"), 'NonGravForces:Model:unsupported');
            tc.verifyErrorId(@() tc.M.series(e, rr, vv, 'Quantities', "flux"), 'NonGravForces:Model:unsupported');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function modelReference(tc)
            % The measured ACT1B accelerations of the same real day: shape, not level.
            [epa, aMeas] = NonGravForces.io.readACC1B(fullfile(tc.root, 'data', 'samples', 'ACT1B_2020-01-01_C_04.txt'));
            n = 200; [e, rr, vv] = tc.arc(n);
            [a, info] = tc.M.accel(e, rr, vv);
            [tf, loc] = ismember(tc.ep.sec(1:n), epa.sec);
            tc.assertEqual(nnz(tf), n, 'the GNV1B and ACT1B excerpts share their second grid');
            am = aMeas(:, loc(tf));
            rho = zeros(1, 3);
            for i = 1:3
                rho(i) = tc.correlation(am(i, :), a(i, :));
            end
            tc.verifyBoundAndRecord("ACT1B cross-track correlation shortfall", 1 - rho(2), 0.10, "-", "measured 0.980; the cross-track is radiation-dominated and the model reproduces its shape");
            tc.verifyBoundAndRecord("ACT1B radial correlation shortfall", 1 - rho(3), 0.10, "-", "measured 0.994");
            % The along-track is recorded, not asserted: over a ten-minute excerpt the measured
            % along-track is 3.2e-7 RMS against a modelled 9.4e-9 and is bias-drift dominated.
            tc.verifyBoundAndRecord("ACT1B along-track correlation (recorded, not a claim)", abs(rho(1)), 1.0, "-", "measured 0.190; the excerpt is too short for the drag signal to rise above the instrument bias drift");
            % Sentman's shear produces genuine LIFT on the 50 deg flanks, so the gas force is not
            % antiparallel to the stream and must not be asserted to be. What IS a construction
            % property is that constantCd, which carries no shear, produces none: measured
            % 2026-09-08, 1e-16 for constantCd against 0.81 deg for Sentman unshadowed and 0.34 deg
            % shadowed (self-shadowing removes part of the flank that generates it).
            u = -info.vrel_body ./ vecnorm(info.vrel_body);
            cosang = sum(info.drag .* u, 1) ./ vecnorm(info.drag);
            tc.verifyBoundAndRecord("Sentman lift: drag angle off the stream", max(acosd(min(cosang, 1))), 2.0, "deg", "the flank shear of an 8-plate GRACE-FO model; measured 0.34 deg");
            oc = tc.M.options.set("drag.gsi", "constantCd");
            Mc = NonGravForces.Model(tc.macro, oc, 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            [~, ic] = Mc.accel(e, rr, vv);
            uc = -ic.vrel_body ./ vecnorm(ic.vrel_body);
            cc = sum(ic.drag .* uc, 1) ./ vecnorm(ic.drag);
            tc.verifyBoundAndRecord("constantCd drag direction vs -v_rel", max(abs(cc - 1)), 1e-14, "-", "constantCd resolves every windward plate along the stream by construction: no lift; measured 1e-16");
            % The thermal recoil against the published GRACE magnitude.
            tc.verifyBoundAndRecord("thermal recoil magnitude vs Siemes 2023 (6.4 nm/s^2)", abs(median(vecnorm(info.thermal)) - 6.4e-9), 4e-9, "m/s^2", "published cross-track RMS for GRACE; measured median 7.9e-9 over this arc");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function modelPrecision(tc)
            [e, rr, vv] = tc.arc(4);
            [a, info] = tc.M.accel(e, rr, vv);
            % The effects are the whole of the acceleration: nothing is added or lost between them.
            tc.verifyBoundAndRecord("accel minus its three effects", max(abs(a - (info.drag + info.radiation + info.thermal)), [], 'all'), 0, "m/s^2", "the same three arrays summed in the same order: bitwise");
            % budget reports exactly what accel returned.
            B = tc.M.budget(e, rr, vv);
            tc.verifyBoundAndRecord("budget total RMS vs accel", max(abs(B{4, {'rms_x','rms_y','rms_z'}}' - sqrt(mean(a.^2, 2)))), 0, "m/s^2", "bitwise: budget reduces the same array");
            % Rotating to GCRS preserves every magnitude.
            ag = tc.M.accel(e, rr, vv, [], 'Frame', "gcrs");
            tc.verifyBoundAndRecord("gcrs vs body magnitude", max(abs(vecnorm(ag) - vecnorm(a))), 1e-23, "m/s^2", "an orthogonal rotation of a 3e-8 vector; measured 7e-24");
            % rhs returns the GCRS acceleration, not a different quantity.
            f = tc.M.rhs();
            tc.verifyBoundAndRecord("rhs handle vs accel in gcrs", max(abs(f(e, rr, vv, []) - ag), [], 'all'), 0, "m/s^2", "the handle calls the same method: bitwise");
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function modelVectorization(tc)
            % A batch of epochs equals the same epochs one at a time, bitwise, and chunking
            % the series changes nothing: that is what makes a chunked or parallel run trustworthy.
            n = 6; [e, rr, vv] = tc.arc(n);
            a = tc.M.accel(e, rr, vv);
            for k = 1:n
                ak = tc.M.accel(tc.slice(k), rr(:, k), vv(:, k));
                tc.verifyEqual(a(:, k), ak, sprintf('epoch %d batched vs alone', k));
            end
            o = tc.M.options;
            Mc = NonGravForces.Model(tc.macro, o.set("series.chunk", 2), 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            M1 = NonGravForces.Model(tc.macro, o.set("series.chunk", 1000), 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            Tc = Mc.series(e, rr, vv); T1 = M1.series(e, rr, vv);
            tc.verifyEqual(Tc.Variables, T1.Variables, 'a series in chunks of 2 must equal the unchunked one bitwise');
            tc.verifyEqual(Tc.Properties.RowTimes, T1.Properties.RowTimes);
            % The descriptions must NOT match: each records the chunking it was computed with,
            % which is metadata travelling with the data (rule R2), not part of the result.
            tc.verifyTrue(contains(Tc.Properties.Description, "3 chunk(s)") && contains(T1.Properties.Description, "1 chunk(s)"));
            tc.verifyEqual([Tc.ax Tc.ay Tc.az]', a, 'the series must equal accel');
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function modelRobustness(tc)
            M0 = NonGravForces.Model(NonGravForces.geom.graceFO("D"));
            tc.verifyEqual(M0.options.preset, "standard");
            tc.verifyEmpty(M0.attitude);
            % A single epoch, and an option combination the level does not deliver.
            [e, rr, vv] = tc.arc(1);
            a = tc.M.accel(e, rr, vv);
            tc.verifySize(a, [3 1]); tc.verifyTrue(all(isfinite(a)));
            % Options that contradict each other are refused at construction (LP-C F5), and each refusal is paired with
            % the construction that must succeed (T6).
            o = tc.M.options.set("drag.gsi", "maxwell");
            mk = @(oo, varargin) NonGravForces.Model(tc.macro, oo, 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, varargin{:});
            tc.verifyErrorId(@() mk(o), 'NonGravForces:Model:inputOptions');
            am = mk(o.set("drag.specularFraction", 0.1)).accel(e, rr, vv); tc.verifyTrue(all(isfinite(am)));
            tc.verifyErrorId(@() mk(tc.M.options.set("drag.accommodation", "walker2014")), 'NonGravForces:Model:inputOptions');
            tc.verifyErrorId(@() mk(tc.M.options.set("drag.accommodation", "walker2014").set("drag.alpha", []).set("drag.surfaceMass_amu", [])), 'NonGravForces:Model:inputOptions');
            aw = mk(tc.M.options.set("drag.accommodation", "walker2014").set("drag.alpha", []).set("drag.surfaceMass_amu", 263.3223)).accel(e, rr, vv);
            ap = mk(tc.M.options.set("drag.accommodation", "pilinski2010").set("drag.alpha", [])).accel(e, rr, vv);
            tc.verifyTrue(all(isfinite([aw; ap])));
            tc.verifyErrorId(@() mk(tc.M.options.set("drag.alpha", [])), 'NonGravForces:Model:inputOptions');
            tc.verifyErrorId(@() mk(tc.M.options.set("thermal.internalPower_W", 100)), 'NonGravForces:Model:inputOptions');
            % the mass is the macro model's, or the caller's, never GRACE-C's by default (LP-C F10)
            o6 = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            bw = NonGravForces.geom.boxWing([3 1.9 0.7], o6);
            tc.verifyErrorId(@() NonGravForces.Model(bw, NonGravForces.Options("fast")), 'NonGravForces:Model:inputMass');
            tc.verifyErrorId(@() NonGravForces.Model(bw, NonGravForces.Options("fast"), 'Mass', -3), 'NonGravForces:Model:inputMass');
            tc.verifyEqual(NonGravForces.Model(bw, NonGravForces.Options("fast"), 'Mass', 150).mass.value_kg, 150);
            tc.verifyEqual(NonGravForces.Model(NonGravForces.geom.graceFO("D"), NonGravForces.Options("fast")).mass.value_kg, NonGravForces.const().graceFO.massD);
            % The fast preset must evaluate with no space weather beyond its exponential profile.
            af = tc.Mfast.accel(e, rr, vv);
            tc.verifyTrue(all(isfinite(af))); tc.verifyGreaterThan(norm(af), 0);
            % An eclipsed epoch (if the arc has one) still returns finite numbers.
            [~, info] = tc.M.accel(tc.arc(30), tc.r(:, 1:30), tc.v(:, 1:30));
            tc.verifyTrue(all(isfinite(info.nu)) && all(info.nu >= 0 & info.nu <= 1));
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function modelMetamorphic(tc)
            [e, rr, vv] = tc.arc(4);
            [a, info] = tc.M.accel(e, rr, vv);
            o = tc.M.options;
            % Doubling the drag scale doubles the drag and leaves the rest untouched.
            M2 = NonGravForces.Model(tc.macro, o.set("drag.cdScale", 2), 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            [~, i2] = M2.accel(e, rr, vv);
            tc.verifyBoundAndRecord("drag doubles with cdScale", max(abs(i2.drag - 2 * info.drag), [], 'all') / max(abs(info.drag), [], 'all'), 1e-15, "-", "an exact power of two through a linear factor");
            tc.verifyEqual(i2.radiation, info.radiation, 'the radiation must not notice the drag scale');
            % Doubling the radiation scale does the same on its side.
            M3 = NonGravForces.Model(tc.macro, o.set("srp.crScale", 2), 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            [~, i3] = M3.accel(e, rr, vv);
            tc.verifyBoundAndRecord("radiation doubles with crScale", max(abs(i3.radiation - 2 * info.radiation), [], 'all') / max(abs(info.radiation), [], 'all'), 1e-15, "-", "as above");
            tc.verifyEqual(i3.drag, info.drag, 'the drag must not notice the radiation scale');
            % The analytic parameter partials against Richardson-extrapolated central differences.
            [J, names] = tc.M.partials(e, rr, vv, 'WithRespectTo', ["cdScale", "crScale", "alpha"]);
            tc.verifyEqual(names, ["cdScale" "crScale" "alpha"]);
            [Jn, ~] = tc.M.partials(e, rr, vv, 'WithRespectTo', ["cdScale", "crScale", "alpha"], 'Method', "finiteDiff");
            for i = 1:3
                sc = max(abs(Jn(:, i, :)), [], 'all');
                tc.verifyBoundAndRecord("d(accel)/d(" + names(i) + ") analytic vs finite difference", ...
                    max(abs(J(:, i, :) - Jn(:, i, :)), [], 'all') / max(sc, realmin), 1e-6, "-", "central difference in the option, step 1e-6 relative; its own truncation dominates");
            end
            % The state partials must reproduce a small perturbation of the acceleration.
            Jr = tc.M.partials(e, rr, vv, 'WithRespectTo', "r");
            d = [30; -20; 10];
            aPred = a + reshape(pagemtimes(Jr, d), 3, []);
            aTrue = tc.M.accel(e, rr + d, vv);
            tc.verifyBoundAndRecord("accel predicted by the r-partials over a 37 m step", max(abs(aPred - aTrue), [], 'all') / max(abs(a), [], 'all'), 5e-3, "-", "second-order truncation of a linear extrapolation over 37 m, against a 500 km scale height");
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function modelSpeed(tc)
            % The fast preset exists to be fast: no self-shadowing, no Earth radiation, no thermal model. Until LP-C it was
            % 45x the cheaper (2 against 89 ms per epoch); with the O(1) cube map, flow shadowing and provenance on demand
            % the standard preset costs 1.2 ms per epoch on the host (2026-09-15, 100 epochs) and the fast one 0.3 ms, a ratio
            % of 0.25. The fast preset must still be the cheaper one by a clear margin: budget 3x the measured ratio.
            [e, rr, vv] = tc.arc(40);
            ratio = tc.speedRatio(@() tc.Mfast.accel(e, rr, vv), @() tc.M.accel(e, rr, vv));
            tc.verifyBoundAndRecord("fast preset vs standard preset, same arc (40 epochs)", ratio, 0.75, "-", "measured 0.25 on 2026-09-15 (it was 0.02 while the standard preset cost 88 ms per epoch)");
        end
    end
end
