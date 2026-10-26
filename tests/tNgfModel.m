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
            % THE THERMAL RECOIL AGAINST TWO PUBLISHED FIGURES, FOR TWO DIFFERENT SATELLITES.
            % The 7.9e-9 this comment claimed until 2026-09-17 was measured before LP-C's flowShadow
            % correction moved the thermal recoil by 1.77e-9, and nothing compared the comment with
            % the run: last_run.txt had been recording 9.25e-10 for |median - 6.4e-9|, i.e. a median
            % of 5.475e-9, for as long as the sentence said 7.9e-9. A measurement copied into prose
            % beside the instrument that records it is a second authority for one fact (R1).
            tc.verifyBoundAndRecord("thermal recoil magnitude vs Siemes 2023 (6.4 nm/s^2)", abs(median(vecnorm(info.thermal)) - 6.4e-9), 4e-9, "m/s^2", "published cross-track RMS for GRACE-A 2008; measured median 5.475e-9 over this arc");
            % [H24] Table 3 case (a) models GRACE-FO with panel geometry, Handbook coefficients and
            % NO thermal emission, and is left with 4.22 nm/s^2 of cross-track residual. The residual
            % left by omitting a term should be the SIZE of the term omitted, so that figure is a
            % magnitude check on this toolbox's thermal recoil for THIS satellite, where Siemes is
            % GRACE-A. The two bracket us: 4.22 (GRACE-FO, 2018-2023) < 5.26 (here) < 6.4 (GRACE-A,
            % 2008). The bound is deliberately loose -- different spacecraft, epochs and coefficient
            % sets -- and exists to catch a sign, a factor or a unit, not to claim agreement.
            tc.verifyBoundAndRecord("thermal recoil cross-track RMS vs Hladczuk 2024 case (a) (4.22 nm/s^2)", abs(sqrt(mean(info.thermal(2, :).^2)) - 4.22e-9), 2.5e-9, "m/s^2", "the residual left by omitting thermal emission bounds its magnitude; measured 5.263e-9, between the GRACE-FO and GRACE-A published figures");
        end
        function modelGroopsReference(tc)
            % O-19: GROOPS against this toolbox on the real day. TU Graz's operational
            % nonConservativeForces for 2020-01-01 GRACE-C (GRACEFO-1) are MODELLED by GROOPS (owner,
            % 2026-10-26), in the satellite reference frame at 10 s, and published with the attitude
            % (SRF -> CRF, scalar first) and the reduced-dynamic orbit (CRF) they were computed on.
            % This test feeds THEIR state into Model.accel, so the two differ only by force model:
            % density, gas-surface interaction, radiation and thermal. All three files are read with
            % gracefoReading.readGroopsInstrument (>= 1.1, which reads real GROOPS output), the first
            % use of that reader on TU Graz's bytes in this repository (CLAUDE.md section 7).
            %
            % Every bound is tightened from recorded runs: the per-axis ones from CI 40406425492, the
            % magnitude, residual and drag-share ones from CI 40411254932.
            S = fullfile(tc.root, 'data', 'samples', 'itsg');
            rd = @(p) gracefoReading.readGroopsInstrument(fullfile(S, "GRACEFO-1_" + p + "_2020-01-01.txt"));
            [tF, aG, iF] = rd("nonConservativeForces");
            [tA, qG] = rd("attitude");
            [tO, xO, iO] = rd("reducedDynamicOrbit");
            tc.assertSize(aG, [600 3]); tc.assertSize(qG, [600 4]); tc.assertSize(xO, [600 9]);
            tc.assertEqual(tA, tF, 'attitude and forces share their 10 s epochs');
            tc.assertEqual(tO, tF, 'orbit and forces share their 10 s epochs');
            tc.verifyEqual(string(iF.columnUnits), ["m/s^2" "m/s^2" "m/s^2"], 'the legend names the units the comparison assumes');
            tc.verifyEqual(string(iO.columnUnits(1:6)), ["m" "m" "m" "m/s" "m/s" "m/s"]);
            ep = referenceSystem.daySecondToEpoch(tF, "GPS");
            % TU Graz's quaternion is BIT-IDENTICAL to SCA1B's for the same epochs (checked 2026-10-26 on this
            % day): their "SRF -> CRF" and this toolbox's "SCA1B is ICRF -> SRF" are one rotation in two
            % quaternion conventions. It is therefore conjugated exactly once, as io.readSCA1B does. Passed as
            % is, Model received the inverse rotation (CI 40403910825: 5.3e-8 m/s^2 RMS of RADIAL drag), and
            % checkAttitudeSense is the guard that refuses that sense now.
            q = referenceSystem.invertQuaternion(qG');
            [okSense, sense] = NonGravForces.geom.checkAttitudeSense(q, xO(:, 1:3)', xO(:, 4:6)');
            tc.assertTrue(all(okSense), "TU Graz's attitude, conjugated once, must pass geom.checkAttitudeSense (see its report)");
            tc.verifyClass(sense, 'struct');
            % THE ORACLE IS CHECKED TOO: over the 30 min where the committed ACT1B excerpt overlaps, GROOPS
            % follows the measured accelerometer (an independent standard-library computation on 2026-10-26
            % gave correlations +0.990, +0.982, +0.910 for x, y, z). Correlation ignores ACT1B's bias and scale.
            [epa, am] = NonGravForces.io.readACC1B(fullfile(tc.root, 'data', 'samples', 'multiday', 'ACT1B_2020-01-01_C_04.txt'));
            [tf, loc] = ismember(round(tF(:, 2))', epa.sec);
            tc.assertEqual(nnz(tf), 180, '30 min of 10 s epochs inside the ACT1B excerpt');
            actBound = [0.05 0.05 0.2];
            for i = 1:3
                cAct = tc.correlation(aG(tf, i)', am(i, loc(tf)));
                tc.verifyBoundAndRecord("GROOPS vs ACT1B correlation shortfall, axis " + i, 1 - cAct, actBound(i), "-", ...
                    "the oracle against the measurement; independent value 2026-10-26: shortfalls 0.010, 0.018, 0.090");
            end
            [a, info] = tc.M.accel(ep, xO(:, 1:3)', xO(:, 4:6)', q);       % body frame = SRF
            g = aG';
            rho = zeros(1, 3); ratio = zeros(1, 3);
            for i = 1:3
                rho(i) = tc.correlation(a(i, :), g(i, :));
                ratio(i) = sqrt(mean(a(i, :).^2)) / sqrt(mean(g(i, :).^2));
            end
            % DECOMPOSITION, printed so the CI log carries it (a passing bound reaches only last_run.txt).
            % The first run (CI 40401231905) found the along-track ANTI-correlated (-0.686) with an RMS
            % residual of 2.5x GROOPS; this decomposition (CI 40403910825) then showed drag of 5.3e-8 RMS in
            % the RADIAL axis, which only a wrong attitude sense produces -- fixed above.
            eff = struct('total', a, 'drag', info.drag, 'radiation', info.radiation, 'thermal', info.thermal);
            fprintf('modelGroopsReference decomposition (2020-01-01 C, 600 x 10 s, m/s^2):\n');
            fprintf('  %-10s %-4s %11s %11s %8s\n', 'effect', 'axis', 'mean', 'rms', 'corr');
            for i = 1:3
                fprintf('  %-10s %-4d %+11.3e %11.3e %8s\n', 'GROOPS', i, mean(g(i, :)), sqrt(mean(g(i, :).^2)), '-');
                for fn = string(fieldnames(eff))'
                    x = eff.(fn)(i, :);
                    c = NaN; if any(x ~= x(1)), c = tc.correlation(x, g(i, :)); end
                    fprintf('  %-10s %-4d %+11.3e %11.3e %+8.3f\n', fn, i, mean(x), sqrt(mean(x.^2)), c);
                end
            end
            mag = median(vecnorm(a)) / median(vecnorm(g));
            resid = sqrt(mean(vecnorm(a - g).^2)) / sqrt(mean(vecnorm(g).^2));
            tc.verifyBoundAndRecord("GROOPS (TU Graz) vs Model: |log| of the median magnitude ratio", abs(log(mag)), log(1.5), "-", ...
                "measured ratio 1.119 (CI 40411254932); a frame or unit error is a factor of 1e3 or a sign");
            % Tightened from CI 40406425492 (the first run with the attitude conjugated): correlations
            % 0.983 / 1.000 / 0.968, RMS ratios 1.50 / 1.18 / 1.20 (|log| 0.404 / 0.166 / 0.180). The model's
            % along-track drag is 1.5x GROOPS's: a density and gas-surface difference, not a frame one.
            axisName = ["along-track x", "cross-track y", "radial z"];
            rhoBound = [0.05 0.01 0.1]; ratioBound = log([2 1.5 1.5]);
            for i = 1:3
                tc.verifyBoundAndRecord("GROOPS vs Model correlation shortfall, " + axisName(i), 1 - rho(i), rhoBound(i), "-", ...
                    "measured 0.017 / 0.0004 / 0.032 (CI 40406425492)");
                tc.verifyBoundAndRecord("GROOPS vs Model RMS ratio, " + axisName(i), abs(log(ratio(i))), ratioBound(i), "-", ...
                    "measured |log| 0.404 / 0.166 / 0.180 (CI 40406425492)");
            end
            tc.verifyBoundAndRecord("GROOPS vs Model RMS of the difference / RMS of GROOPS", resid, 0.6, "-", ...
                "measured 0.391 (CI 40411254932); drag dominates the along-track, radiation the cross-track");
            tc.verifyBoundAndRecord("GROOPS comparison: modelled drag share of |a| (recorded)", median(vecnorm(info.drag) ./ vecnorm(a)), 0.75, "-", ...
                "measured 0.497 (CI 40411254932): context for the along-track numbers above");
            fprintf('  magnitude ratio %.4f, RMS difference / RMS GROOPS %.4f, drag share %.4f\n', mag, resid, median(vecnorm(info.drag) ./ vecnorm(a)));
        end
        function modelAct1bReference(tc)
            % THE MODEL AGAINST THE MEASUREMENT, the way Huckfeldt et al. (2024, ASR 73, Sec. 3.5, eqs.
            % 7-10) validate theirs: ACT1B with the thruster firings cut, one bias per day and per axis,
            % scale 1, then the daily RMS of model minus calibrated ACT1B. A bias-only fit keeps model
            % errors out of the calibration parameters; a bias-plus-scale fit is recorded beside it,
            % because GRACE-C's ACT1B scale is not 1 (tools/baseline_reference.json: S_xx 0.762).
            %
            % Data: the committed 30 min of each of 2020-01-01..03 (GRACE-C, data/samples/multiday),
            % this toolbox's own state -- GNV1B, SCA1B conjugated once by io.readSCA1B -- and the
            % standard preset. THR1B is committed for 2020-01-01 only; that day's three firings in the
            % excerpt (50 ms attitude pulses at 250, 565 and 1259 s) are masked over [t - 0.5 s, t +
            % onTime + 7 s] (Huckfeldt Sec. 2.1 and Fig. 3: 500 ms before, the operational margin after)
            % widened by the CRN kernel's half span of 70.3 s on each side, since ACT1B carries the spike
            % smeared through the 1407-tap 10 Hz filter. Days 2 and 3 are unmasked; the day-1 numbers
            % with and without the mask say what that costs. ACT1B is CRN-filtered and the model is not;
            % over a 30 min arc that matters only at eclipse transitions.
            %
            % Independent expectation (stdlib Python, 2026-10-26): GROOPS against the same ACT1B at 10 s,
            % bias only, gives 2.7 / 3.1 / 1.8 nm/s^2 RMS (x, y, z), scale 0.92 / 0.85 / 0.87, and the
            % mask moves the x RMS by 4 %. Huckfeldt report 3-5 nm/s^2 along-track and never more than
            % 2 nm/s^2 cross-track and radial for GRACE-C in 2019 at solar minimum, 20 nm/s^2 at worst
            % along-track.
            %
            % MEASURED (CI 40411254932), bias only: 5.5 / 3.4 / 2.1 nm/s^2 (x, y, z); per day x 1.1, 2.8,
            % 8.7 nm/s^2. Fitted scale 1.034 / 1.026 / 0.910, per-day demeaned correlation 0.813 / 0.940 /
            % 0.955. The day-1 mask moves x by 2 %, y by 27 % (2.0 -> 1.5 nm/s^2). Along-track sits at the
            % top of Huckfeldt's 3-5 nm/s^2, cross-track and radial just above their 2: the same order as
            % GROOPS against the same ACT1B (2.7 / 3.1 / 1.8). The bounds are those values with margin.
            f = @(varargin) fullfile(tc.root, 'data', 'samples', varargin{:});
            [~, thr] = NonGravForces.io.readTHR1B(f('THR1B_2020-01-01_C_04.txt'));
            k = any(thr.fired, 1);
            tFire = reshape(thr.epoch.sec(k), 1, []); onT = max(thr.onTime_s(:, k), [], 1);
            halfCrn = (1407 - 1) / 2 * 0.1;                      % s, half span of crnKernel("ACC")
            masked = @(t) any(t(:) >= tFire - 0.5 - halfCrn & t(:) <= tFire + onT + 7 + halfCrn, 2)';
            tc.assertEqual(nnz(tFire < 1800), 3, 'three firings inside the 30 min excerpt of 2020-01-01');
            aM = cell(1, 3); aA = cell(1, 3); keep = cell(1, 3);
            for d = 1:3
                day = sprintf('2020-01-%02d', d);
                [ep, r, v] = NonGravForces.io.readGNV1B(f('multiday', "GNV1B_" + day + "_C_04.txt"), 'Eop', tc.eop);
                [epq, q] = NonGravForces.io.readSCA1B(f('multiday', "SCA1B_" + day + "_C_04.txt"));
                [epa, am] = NonGravForces.io.readACC1B(f('multiday', "ACT1B_" + day + "_C_04.txt"));
                [tf, loc] = ismember(ep.sec, epa.sec);
                tc.assertGreaterThan(nnz(tf), 1700, 'GNV1B and ACT1B share their 1 s grid over the 30 min');
                e = referenceSystem.daySecondToEpoch([ep.day(tf) ep.sec(tf)], ep.scale);
                attD = NonGravForces.geom.AttitudeProvider(epq, q);
                qk = attD.at(e);
                aM{d} = tc.M.accel(e, r(:, tf), v(:, tf), qk);
                aA{d} = am(:, loc(tf));
                keep{d} = true(1, nnz(tf));
                if d == 1, keep{d} = ~masked(ep.sec(tf)); end
            end
            tc.assertGreaterThan(nnz(~keep{1}), 300, 'the mask removes the three firings with their CRN spread');
            [rb, rs, s, rho, rDay] = localFit(aM, aA, keep);
            rb1 = localFit(aM(1), aA(1), keep(1));
            rb1u = localFit(aM(1), aA(1), {true(size(keep{1}))});
            fprintf('modelAct1bReference (2020-01-01..03 C, 3 x 30 min at 1 s, standard preset, m/s^2):\n');
            fprintf('  axis  bias-only rms  scale  bias+scale rms    corr   bias-only RMS of day 1, 2, 3   day 1 unmasked\n');
            for i = 1:3
                fprintf('  %-4d  %13.3e  %5.3f  %14.3e  %+6.3f   %9.3e %9.3e %9.3e   %9.3e\n', i, rb(i), s(i), rs(i), rho(i), rDay(i, :), rb1u(i));
            end
            axisName = ["along-track x", "cross-track y", "radial z"];
            rbBound = [10e-9 6e-9 4e-9]; rhoBound = [0.3 0.1 0.1];
            for i = 1:3
                tc.verifyBoundAndRecord("ACT1B vs Model, bias-only RMS (Huckfeldt eq. 10), " + axisName(i), rb(i), rbBound(i), "m/s^2", ...
                    "measured 5.55 / 3.37 / 2.10 nm/s^2 (CI 40411254932); Huckfeldt 2024: 3-5 x, <= 2 y and z; GROOPS 2.7 / 3.1 / 1.8");
                tc.verifyBoundAndRecord("ACT1B vs Model, bias+scale RMS over bias-only RMS, " + axisName(i), rs(i) / rb(i), 1 + 1e-12, "-", ...
                    "least squares with one more column cannot do worse: an identity");
                tc.verifyBoundAndRecord("ACT1B vs Model, |log| of the fitted scale, " + axisName(i), abs(log(s(i))), log(1.25), "-", ...
                    "measured 1.034 / 1.026 / 0.910 (CI 40411254932); GROOPS against the same ACT1B 0.92 / 0.85 / 0.87");
                tc.verifyBoundAndRecord("ACT1B vs Model correlation shortfall (per day demeaned), " + axisName(i), 1 - rho(i), rhoBound(i), "-", ...
                    "measured 0.187 / 0.060 / 0.045 (CI 40411254932); GROOPS against the same ACT1B 0.07 / 0.06 / 0.03");
                tc.verifyBoundAndRecord("ACT1B vs Model, day 1: |log| of masked over unmasked bias-only RMS, " + axisName(i), abs(log(rb1(i) / rb1u(i))), 0.5, "-", ...
                    "what the THR1B mask is worth on the one day it is committed for; measured 0.019 / 0.314 / 0.067 (CI 40411254932)");
            end
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
            % Model.rhs is GONE (LP-C WP-C9, O-4): the integrator's entry point is app.DynamicsSource, whose
            % accel() is asserted against this method in tests/tNgfIntegration.m. A handle that matched no
            % signature the integrator takes was worse than no handle at all.
            tc.verifyFalse(ismember("rhs", string(methods(tc.M))), 'Model.rhs must not come back without the interface it was supposed to serve');
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
            tc.verifyErrorId(@() mk(tc.M.options.set("geom.selfShadow", "raycast")), 'NonGravForces:Model:inputOptions');   % a table of the lut beside a ray cast
            tc.verifyClass(mk(tc.M.options.set("geom.selfShadow", "raycast").set("geom.flowShadow", "exact")), 'NonGravForces.Model');
            tc.verifyErrorId(@() mk(tc.M.options.set("thermal.internalPower_W", 100)), 'NonGravForces:Model:inputOptions');
            % the mass is the macro model's, or the caller's, never GRACE-C's by default (LP-C F10)
            o6 = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            bw = NonGravForces.geom.boxWing([3 1.9 0.7], o6);
            tc.verifyErrorId(@() NonGravForces.Model(bw, NonGravForces.Options("fast")), 'NonGravForces:Model:inputMass');
            tc.verifyErrorId(@() NonGravForces.Model(bw, NonGravForces.Options("fast"), 'Mass', -3), 'NonGravForces:Model:inputMass');
            tc.verifyEqual(NonGravForces.Model(bw, NonGravForces.Options("fast"), 'Mass', 150).mass.value_kg, 150);
            tc.verifyEqual(NonGravForces.Model(NonGravForces.geom.graceFO("D"), NonGravForces.Options("fast")).mass.value_kg, NonGravForces.const().graceFO.massD);
            % erp.backend ceres (LP-C WP-C8) evaluates through the Model with the committed EBAF excerpt, and without it refuses
            cerFiles = string(fullfile(tc.root, 'data', 'ceres', ["rsut" "rlut" "rsdt"] + "_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc"));
            oC = tc.M.options.set("erp.backend", "ceres");
            tc.verifyErrorId(@() mk(oC).accel(e, rr, vv), 'NonGravForces:Model:dependCeres');
            [aC, iC] = mk(oC, 'Ceres', NonGravForces.io.importCERES(cerFiles)).accel(e, rr, vv);
            tc.verifyTrue(all(isfinite(aC))); tc.verifyEqual(iC.backends.erp, "ceres");
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
            % WHAT THERMAL INERTIA IS WORTH, against what [H24] measured it to be worth. Their
            % Table 3 goes from instantaneous reradiation (b, 3.04 nm/s^2 residual) to thermal
            % inertia (c, 1.70), so the step buys 1.34 nm/s^2 on GRACE-FO with everything else held.
            % Ours is a MODEL difference where theirs is a RESIDUAL reduction, so the two are the
            % same order rather than the same number, and the assertion is that the option does
            % something and that the something is not enormous. A transient model that agreed with
            % steadyState to rounding would be a stub, and one that moved the force by 10 nm/s^2
            % would be wrong; both are excluded here.
            [eL, rL, vL] = tc.arc(200);
            Mtr = NonGravForces.Model(tc.macro, o.set("thermal.model", "transient"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            [~, iSt] = tc.M.accel(eL, rL, vL);
            [~, iTr] = Mtr.accel(eL, rL, vL);
            step = sqrt(mean(sum((iTr.thermal - iSt.thermal).^2, 1)));
            tc.verifyBoundAndRecord("thermal inertia vs instantaneous: the step it is worth", step, 4e-9, "m/s^2", "[H24] Table 3 (b)->(c) buys 1.34 nm/s^2 of residual on GRACE-FO; measured 0.950e-9 as a model difference, the same order");
            tc.verifyGreaterThan(step, 1e-10, 'a transient model that matched steadyState would be a stub, not a model');
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
            % Changing the density backend (msis00 -> msis21, LP-C WP-C6) changes the drag and nothing else.
            M21 = NonGravForces.Model(tc.macro, o.set("drag.density", "msis21"), 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            [a21, i21] = M21.accel(e, rr, vv);
            tc.verifyEqual(i21.backends.density, "msis21"); tc.verifyTrue(all(isfinite(a21), 'all'));
            tc.verifyEqual(i21.radiation, info.radiation, 'the radiation must not notice the density backend');
            tc.verifyEqual(i21.thermal, info.thermal, 'nor the thermal recoil');
            tc.verifyBoundAndRecord("drag msis21 / msis00 through the Model, largest departure from 1", max(abs(log(vecnorm(i21.drag) ./ vecnorm(info.drag)))), log(2), "-", "two climatologies of one thermosphere (the density ratio measured 0.88-0.99 at 490 km that day): a plausibility bound on the wiring");
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
            % 45x the cheaper (2 against 89 ms per epoch, both superseded); with the O(1) cube map, flow shadowing and
            % provenance on demand the standard preset cost 1.2 ms per epoch on the host (2026-09-15, 100 epochs; superseded
            % by tools/benchmark_reference.csv, 0.798 ms on 2026-09-18) and the fast one 0.3 ms, a ratio
            % of 0.25. The fast preset must still be the cheaper one by a clear margin: budget 3x the measured ratio.
            [e, rr, vv] = tc.arc(40);
            ratio = tc.speedRatio(@() tc.Mfast.accel(e, rr, vv), @() tc.M.accel(e, rr, vv));
            tc.verifyBoundAndRecord("fast preset vs standard preset, same arc (40 epochs)", ratio, 0.75, "-", "measured 0.25 on 2026-09-15 (it was 0.02 before WP-C3; that 88 ms of HANDOVER O-29 was per-call overhead, not physics)");
        end
    end
end

function [rb, rs, s, rho, rDay] = localFit(aM, aA, keep)
%LOCALFIT Huckfeldt eqs. 7-10 over days: one bias per day and axis (scale 1), and beside it a common
% scale with per-day biases. Inputs are 1xD cells of 3xK model, 3xK ACT1B and 1xK logical masks.
% Per-day demeaning IS the per-day bias estimate (X a column of ones in eq. 9), so both fits are closed form.
D = numel(aM); rb = zeros(3, 1); rs = zeros(3, 1); s = zeros(3, 1); rho = zeros(3, 1); rDay = zeros(3, D);
for i = 1:3
    dm = []; ym = []; xa = [];
    for d = 1:D
        m = aM{d}(i, keep{d}); a = aA{d}(i, keep{d});
        r = (m - a) - mean(m - a);
        rDay(i, d) = sqrt(mean(r.^2));
        dm = [dm, r]; ym = [ym, m - mean(m)]; xa = [xa, a - mean(a)]; %#ok<AGROW>
    end
    rb(i) = sqrt(mean(dm.^2));
    s(i) = (xa * ym') / (xa * xa');
    rs(i) = sqrt(mean((ym - s(i) * xa).^2));
    rho(i) = (xa * ym') / (norm(xa) * norm(ym));
end
end
