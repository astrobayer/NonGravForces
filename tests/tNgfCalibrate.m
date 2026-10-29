classdef tNgfCalibrate < NgfTestCase
%TNGFCALIBRATE Tests of the calibration application (N-02c): app.calibrate over timeSeriesPreparation.run.
%
%   app.calibrate composes three things that are tested on their own --
%   timeSeriesPreparation.run, app.fitCalibration and the Model -- so its
%   reference oracle is the composition: the calibration it returns must be
%   BITWISE the one app.fitCalibration gives on run's arc and the same raw
%   reference. What is new here is the Model path (per-effect evaluation along
%   the orbit, pchip onto the raw grid), the jump parameters and the
%   attribution, and the real-data integration on the 2020-01-01 GRACE-C
%   excerpt: ACT1B against the Model along GNV1B with SCA1B attitude.
%
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
    properties
        root, eop, sw, acc, orb, att, Mfast
    end
    properties (Constant)
        Quiet = 1e-30 * ones(3, 1)
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            f = @(varargin) fullfile(tc.root, varargin{:});
            tc.eop = referenceSystem.readEop(f('data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            tc.sw = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [ep, r, v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', tc.eop);
            tc.orb = struct('epoch', ep, 'r', r, 'v', v);
            [epA, q] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
            tc.att = NonGravForces.geom.AttitudeProvider(epA, q);
            [epC, a] = NonGravForces.io.readACC1B(f('data', 'samples', 'ACT1B_2020-01-01_C_04.txt'));
            names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
            units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
            aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
                "variableUnits", units, "fileName", "ACT1B_2020-01-01_C_04.txt", "product", "ACT1B");
            tc.acc = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), a.', aux);
            tc.Mfast = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("fast"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, ...
                'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function calibrateContract(tc)
            % a day after the real excerpt, so that the orbit of 2020-01-01 cannot cover it
            start = int64(7306 * 86400 - 43200 + 18) * int64(1e9);
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 2, 'Start', start);
            res = NonGravForces.app.calibrate(obs, ref, 'Fit', {'NBoot', 0});
            tc.verifyEqual(string(fieldnames(res))', ["arc" "reference" "diagnostics" "correlation" "flags" ...
                "config" "provenance" "calibration" "attribution"]);
            tc.verifyClass(res.calibration, 'NonGravForces.app.Calibration');
            tc.verifyTrue(all(ismember(["crn" "align" "qc"], res.diagnostics.stages)));
            tc.verifyClass(res.diagnostics.rawReference, 'timeSeriesPreparation.Arc');
            tc.verifyEqual(res.attribution.verdict, repmat("unavailable", 1, 3));
            C = @NonGravForces.app.calibrate;
            tc.verifyErrorId(@() C(5, ref), "NonGravForces:calibrate:inputArc");
            tc.verifyErrorId(@() C(obs, 5), "NonGravForces:calibrate:inputReference");
            tc.verifyErrorId(@() C(obs, ref, 'Config', timeSeriesPreparation.Config()), "NonGravForces:calibrate:inputConfig");
            tc.verifyErrorId(@() C(tc.acc, tc.Mfast), "NonGravForces:calibrate:inputOrbit");
            tc.verifyErrorId(@() C(obs, tc.Mfast, 'Orbit', tc.orb), "NonGravForces:calibrate:inputOrbit");   % 2020-01-01 orbit, a 2020-01-02 arc
        end
        function calibrateCampaignContract(tc)
            % N-09: an Arc array is a campaign -- run per day, one joint fit, results per day.
            obs = timeSeriesPreparation.Arc.empty; ref = obs;
            for k = 1:2
                st = int64(7305 * 86400 - 43200 + 18 + (k - 1) * 86400) * int64(1e9);
                [obs(k), ref(k)] = NonGravForces.app.syntheticPair('Days', 0.1, 'Start', st, 'Seed', k, 'Spikes', 0);
            end
            res = NonGravForces.app.calibrate(obs, ref, 'Fit', {'NBoot', 0});
            tc.verifyTrue(all(isfield(res, {'results', 'summary', 'calibration', 'used', 'fit', 'attribution'})));
            tc.verifyEqual(res.used, [true true]);
            tc.verifySize(res.calibration, [1 2]);
            tc.verifyEqual(res.calibration{2}.scale, res.calibration{1}.scale, 'one common scale');
            tc.verifyEqual(res.fit.mode, "common");
            tc.verifyErrorId(@() NonGravForces.app.calibrate(obs, ref, 'Jumps', true), "NonGravForces:calibrate:unsupportedMultiArc");
            tc.verifyErrorId(@() NonGravForces.app.calibrate(obs, [ref, ref]), "NonGravForces:calibrate:inputReference");
            tc.verifyErrorId(@() NonGravForces.app.calibrate(obs, tc.Mfast), "NonGravForces:calibrate:inputOrbit");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function calibrateReference(tc)
            % The composition is the oracle: run's arc and the raw reference through fitCalibration
            % with run's kernel give the same calibration, bit for bit.
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.1, 'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            fit = {'NBoot', 0, 'FixCrossTrack', false};
            res = NonGravForces.app.calibrate(obs, ref, 'Jumps', false, 'Fit', fit);
            K = timeSeriesPreparation.crnKernel("ACC");
            direct = NonGravForces.app.fitCalibration(res.arc, ref, 'Kernel', K, fit{:});
            tc.verifyEqual(res.calibration.scale, direct.scale);
            tc.verifyEqual(res.calibration.bias.coef, direct.bias.coef);
            % ... and on a noise-free truth it is the truth, to the rounding bound fitCalibration's
            % own reference test computes (here as a ratio of the formal cofactor)
            info = res.diagnostics.fit;
            for i = [1 3]
                q = sqrt(info.Qxx{i}(end, end)); n = nnz(info.rows(:, i));
                eRow = 4 * K.Nf * eps * sum(abs(K.F)) * (max(abs(ref.obs), [], 'all') * norm(truth.scale, inf) + max(abs(truth.bias), [], 'all'));
                bound = 2 * q * sqrt(n) * eRow + 100 * eps * res.calibration.gates.cond(i) * q * sqrt(n) * sqrt(mean(obs.obs(info.rows(:, i), i).^2));
                tc.verifyBoundAndRecord("calibrate noise-free scale axis " + i, abs(res.calibration.scale(i, i) - truth.scale(i, i)) / bound, 1, "-", "error / computed bound");
            end
        end
        function calibrateRealDay(tc)
            % Integration: ACT1B of 2020-01-01 against the Model along GNV1B with SCA1B attitude.
            res = NonGravForces.app.calibrate(tc.acc, tc.Mfast, 'Orbit', tc.orb, 'Attitude', tc.att, ...
                'Fit', {'NBoot', 50, 'FixCrossTrack', false});
            cal = res.calibration;
            tc.verifyTrue(all(isfinite(cal.scale), 'all') && all(isfinite(cal.bias.coef), 'all'));
            tc.verifyEqual(res.attribution.effects, ["drag" "radiation" "thermal"]);
            tc.verifyTrue(all(ismember(res.attribution.verdict, ["resolved" "unresolved"])));
            tc.verifyTrue(all(isnan(res.attribution.relativeCorrection(3, :))), 'the fast preset has no thermal effect, so no column');
            rows = res.diagnostics.fit.rows;
            tc.verifyGreaterThan(nnz(rows(:, 1)), 300, 'the 600 s excerpt keeps its interior after the CRN window');
            for i = 1:3
                y = tc.acc.obs(rows(:, i), i);
                tc.verifyBoundAndRecord("calibrate real day axis " + i + ": residual RMS / demeaned ACT1B RMS", ...
                    cal.uncertainty.residualRms(i) / sqrt(mean((y - mean(y)).^2)), 1 + 1e-12, "-", ...
                    "least squares with a constant in the basis cannot do worse than the mean");
            end
        end
        function calibrateCampaignRealDays(tc)
            % N-09 on real data: GRACE-C 2020-01-01..03 (30 min of each, GFZ ISDC RL04, fetched
            % 2026-10-15), ACT1B against the fast Model along each day's GNV1B with SCA1B attitude.
            % Completion and consistency, not truth: one common scale, finite per-day biases that
            % agree to a few per cent between consecutive days.
            f = @(varargin) fullfile(tc.root, varargin{:});
            names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
            units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
            accs = timeSeriesPreparation.Arc.empty; orbs = struct('epoch', {}, 'r', {}, 'v', {}); epA = []; qA = [];
            for d = 1:3
                day = sprintf('2020-01-%02d', d);
                [ep, r, v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'multiday', "GNV1B_" + day + "_C_04.txt"), 'Eop', tc.eop);
                orbs(d) = struct('epoch', ep, 'r', r, 'v', v);
                [e, q] = NonGravForces.io.readSCA1B(f('data', 'samples', 'multiday', "SCA1B_" + day + "_C_04.txt"));
                if isempty(epA), epA = e; qA = q;
                else, epA = referenceSystem.daySecondToEpoch([[epA.day; e.day], [epA.sec; e.sec]], epA.scale); qA = [qA, q]; %#ok<AGROW>
                end
                [epC, a] = NonGravForces.io.readACC1B(f('data', 'samples', 'multiday', "ACT1B_" + day + "_C_04.txt"));
                aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
                    "variableUnits", units, "fileName", "ACT1B_" + day + "_C_04.txt", "product", "ACT1B");
                accs(d) = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), a.', aux);
            end
            attP = NonGravForces.geom.AttitudeProvider(epA, qA);
            M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("fast"), 'Attitude', attP, ...
                'Eop', tc.eop, 'SpaceWeather', tc.sw, 'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));
            res = NonGravForces.app.calibrate(accs, M, 'Orbit', orbs, 'Attitude', attP, 'Fit', {'NBoot', 0});
            tc.verifyEqual(res.used, true(1, 3));
            tc.verifyEqual(res.summary.status, repmat("ok", 3, 1));
            S = res.calibration{1}.scale;
            tc.verifyTrue(all(isfinite(S), 'all'));
            tc.verifyEqual(res.calibration{3}.scale, S);
            bx = cellfun(@(c) c.bias.coef(1, 1), res.calibration);
            tc.verifyBoundAndRecord("calibrate campaign, real x bias day-to-day spread / mean", (max(bx) - min(bx)) / abs(mean(bx)), 0.2, "-", ...
                "three consecutive days of one instrument; measured 0.08 when written (-3.38, -3.32, -3.12e-7)");
        end

        function calibrateRegressorsRealDay(tc)
            % N-11a on real data: 30 min of 2020-01-01 GRACE-C, ACT1B against the fast Model, with the
            % sensor-unit temperature TSU_Y+ (AHK1B) and the three net magnetorquer currents (MAG1B) as
            % measured regressors, both on the MAG1B grid (the temperature linearly). Completion and finite
            % coefficients. TICUR is recorded, not refused: it reads absolute zero with a 3 mK quantisation
            % wobble, so the constancy check cannot tell it dead (found writing this test).
            f = @(varargin) fullfile(tc.root, 'data', 'samples', 'multiday', varargin{:});
            [epM, fld] = NonGravForces.io.readMAG1B(f('MAG1B_2020-01-01_C_04.txt'));
            tM = referenceSystem.epochTo(epM, "gpsNanoseconds");
            [tA, hk, hkAux] = gracefoReading.readL1b(f('AHK1B_2020-01-01_C_04.txt'), 'TimeFormat', "gpsns");
            nms = string(hkAux.obsNames);
            tsu = hk(:, nms == "TSU_Y+"); ok = isfinite(tsu);
            temp = interp1(double(tA(ok) - tM(1)) / 1e9, tsu(ok), double(tM - tM(1)) / 1e9, 'linear', NaN);
            cn = {'TSU_Y', 'MTQ_x', 'MTQ_y', 'MTQ_z'}; cu = struct('TSU_Y', 'C', 'MTQ_x', 'A', 'MTQ_y', 'A', 'MTQ_z', 'A');
            R = timeSeriesPreparation.Arc(tM, [temp(:), fld.torquerNet_A.'], struct("obsNames", {cn}, "variableNames", {cn}, ...
                "variableFieldNames", {cn}, "variableUnits", cu, "fileName", "MAG1B+AHK1B", "product", "REGRESSORS"));
            ticur = hk(isfinite(hk(:, nms == "TICUR")), nms == "TICUR");
            tc.verifyBoundAndRecord("calibrate real day, TICUR distance from absolute zero", max(abs(ticur + 273.15)), 0.1, "C", ...
                "a dead channel: -273.2 C on every record, varying by its 3 mK quantisation (choose live channels)");
            [ep, r, v] = NonGravForces.io.readGNV1B(f('GNV1B_2020-01-01_C_04.txt'), 'Eop', tc.eop);
            [e, q] = NonGravForces.io.readSCA1B(f('SCA1B_2020-01-01_C_04.txt'));
            attP = NonGravForces.geom.AttitudeProvider(e, q);
            [epC, a] = NonGravForces.io.readACC1B(f('ACT1B_2020-01-01_C_04.txt'));
            names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
            units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
            accR = timeSeriesPreparation.Arc(referenceSystem.epochTo(epC, "gpsNanoseconds"), a.', struct("obsNames", {names}, ...
                "variableNames", {names}, "variableFieldNames", {names}, "variableUnits", units, "fileName", "ACT1B_2020-01-01_C_04.txt", "product", "ACT1B"));
            M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("fast"), 'Attitude', attP, ...
                'Eop', tc.eop, 'SpaceWeather', tc.sw, 'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));
            orbs = struct('epoch', ep, 'r', r, 'v', v);
            res = NonGravForces.app.calibrate(accR, M, 'Orbit', orbs, 'Attitude', attP, 'Fit', {'NBoot', 0, 'Regressors', R});
            c = res.calibration;
            tc.verifyEqual(c.regressors.names, ["TSU_Y" "MTQ_x" "MTQ_y" "MTQ_z"]);
            tc.verifyTrue(all(isfinite(c.regressors.coef), 'all'));
            tc.verifyBoundAndRecord("calibrate real day, TSU_Y+ span over the 30 min", max(temp) - min(temp), 0.05, "C", ...
                "recorded: the thermally controlled sensor unit barely moves (0.022 C over the whole day)");
            out = c.apply(accR, 'Regressors', R);
            tc.verifySize(out.obs, size(accR.obs));
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function calibrateReferenceInterpolationPrecision(tc)
            % HANDOVER Part I, N-02c, the owner's question of 2026-10-15: calibrate evaluates the Model
            % at the 1 s orbit epochs and interpolates by pchip onto the 10 Hz raw grid. Measured here
            % against the Model evaluated AT 10 Hz, across a shadow entry -- the hardest case, the
            % radiation dropping in seconds -- through calibrate's own route
            % (diagnostics.rawReference), both then CRN-filtered as the fit sees them.
            % The orbit is circular at 490 km in the plane holding the Sun direction, so it crosses
            % the umbra; states at both rates are the closed form, so no orbit is interpolated.
            mu = 3.986004418e14; Re = 6378137; R = Re + 490e3; n = sqrt(mu / R^3);
            day = 58849; sec0 = 3600;                               % 2020-01-01 01:00 GPS
            s = referenceSystem.geocentricBodyState("sun", referenceSystem.daySecondToEpoch([day sec0], "GPS"));
            s = s / norm(s);
            e1 = -s; e2 = cross(s, [0; 0; 1]); e2 = e2 / norm(e2);
            phi0 = -asin(Re / R);                                   % cylindrical shadow entry
            t1 = -680:1:680; t10 = -680:0.1:680;                    % 1200 s plus the CRN margins
            r = @(t) R * (cos(phi0 + n * t) .* e1 + sin(phi0 + n * t) .* e2);
            v = @(t) R * n * (-sin(phi0 + n * t) .* e1 + cos(phi0 + n * t) .* e2);
            ep1 = referenceSystem.daySecondToEpoch([repmat(day, numel(t1), 1), sec0 + t1(:)], "GPS");
            ep10 = referenceSystem.daySecondToEpoch([repmat(day, numel(t10), 1), sec0 + t10(:)], "GPS");
            attP = NonGravForces.geom.AttitudeProvider(ep10, NonGravForces.geom.AttitudeProvider.nominal(r(t10), v(t10), true));
            M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NgfTestCase.standardUnscaled(), ...
                'Attitude', attP, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
            [a10, i10] = M.accel(ep10, r(t10), v(t10), attP.at(ep10));
            rad = vecnorm(i10.radiation);
            tc.verifyLessThan(min(rad) / max(rad), 0.5, 'the arc must reach the umbra, or it measures a sunlit arc');
            % the observation: the 10 Hz Model, CRN-filtered at the whole seconds, as ACT1B is formed
            K = timeSeriesPreparation.crnKernel("ACC");
            tns = referenceSystem.epochTo(ep10, "gpsNanoseconds");
            names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
            units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
            mk = @(t, X, prod) timeSeriesPreparation.Arc(t, X, struct("obsNames", {names}, "variableNames", {names}, ...
                "variableFieldNames", {names}, "variableUnits", units, "fileName", prod + "_2020-01-01_C_04.txt", "product", prod));
            Fd = timeSeriesPreparation.crnApply(mk(tns, a10.', "MODEL"), K, OutputEpochs = tns(1:10:end));
            keep = all(isfinite(Fd.obs), 2);
            accA = mk(Fd.t(keep), Fd.obs(keep, :), "ACT1B");
            res = NonGravForces.app.calibrate(accA, M, 'Orbit', struct('epoch', ep1, 'r', r(t1), 'v', v(t1)), ...
                'Attitude', attP, 'Fit', {'NBoot', 0});
            raw = res.diagnostics.rawReference;
            [onGrid, at] = ismember(raw.t, tns);
            tc.verifyTrue(all(onGrid), 'the raw reference lies on the 10 Hz grid');
            Fr = timeSeriesPreparation.crnApply(raw, K, OutputEpochs = accA.t).obs;
            d = Fr - accA.obs; ok = all(isfinite(d), 2);
            signal = max(abs(accA.obs), [], 'all');
            rawErr = max(abs(raw.obs - a10(:, at).'), [], 'all');
            tc.verifyBoundAndRecord("calibrate reference, pchip vs 10 Hz Model, raw grid, / signal", rawErr / signal, 0.05, "-", ...
                "recorded: unfiltered, the 1 s interpolation is visible at the shadow boundary (measured 3.9e-3)");
            tc.verifyBoundAndRecord("calibrate reference, pchip vs 10 Hz Model, CRN-filtered, / signal", max(abs(d(ok, :)), [], 'all') / signal, 1e-3, "-", ...
                "100x below the Model's own 10 % (density) error; measured 3.2e-5, i.e. 9.3e-13 m/s^2, against a filtered noise floor ESTIMATED at ~2e-11 (1e-10/rtHz over the 0.035 Hz band)");
        end

        function calibratePrecision(tc)
            % Coloured noise, spikes and a gap through run's QC: the diagonal scale within 5 bootstrap sigma.
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 7);
            res = NonGravForces.app.calibrate(obs, ref, 'Fit', {'FixCrossTrack', false, 'NBoot', 50, 'Seed', 2});
            cal = res.calibration;
            z = abs(diag(cal.scale) - diag(truth.scale)) ./ diag(cal.uncertainty.sigmaScaleBootstrap);
            tc.verifyBoundAndRecord("calibrate scale error in bootstrap sigma", max(z), 5, "sigma", "seeded noise, 50 replicates, run's QC");
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function calibrateRobustness(tc)
            % A jump in the bias, declared: the step column recovers it and leaves the scale alone.
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.1, 'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            K = timeSeriesPreparation.crnKernel("ACC");
            tj = obs.t(round(numel(obs.t) / 2));
            stepRaw = [3e-8 * double(ref.t >= tj), zeros(numel(ref.t), 2)];
            Hs = timeSeriesPreparation.crnApply(ref.derive(ref.t, stepRaw, "test.step"), K, OutputEpochs = obs.t);
            obsJ = obs.derive(obs.t, obs.obs + Hs.obs, "test.jump");
            res = NonGravForces.app.calibrate(obsJ, ref, 'Jumps', false, 'Fit', {'NBoot', 0, 'FixCrossTrack', false, 'Jumps', tj});
            cal = res.calibration;
            tc.verifyEqual(cal.bias.jumps, tj);
            tc.verifyEqual(size(cal.bias.coef), [3 3]);
            tc.verifyBoundAndRecord("calibrate recovers a declared jump", abs(cal.bias.coef(1, 3) - 3e-8) / 3e-8, 1e-6, "rel", ...
                "noise-free; the step is one more exact column");
            tc.verifyBoundAndRecord("calibrate scale beside a jump", abs(cal.scale(1, 1) - truth.scale(1, 1)), 1e-6, "-", "noise-free");
            % The default pair, with spikes and a gap, runs end to end and fits only unmasked rows.
            [o2, r2] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 3);
            res2 = NonGravForces.app.calibrate(o2, r2, 'Fit', {'NBoot', 0});
            tc.verifyFalse(any(res2.diagnostics.fit.rows(res2.arc.weights == 0, :), 'all'));
            % qc's jump events are reported in GPS ns and NOT fitted by default. The oracle for the
            % conversion is plain arithmetic: the origin is 2000-01-01 12:00 GPS (= 11:59:47 UTC), so
            % GPS ns = elapsed UTC seconds since 2000-01-01 12:00 UTC (datetime ignores leap seconds)
            % plus GPS-UTC of the date (18 s). (+5 measured exactly 13 s off on CI, run of 19038cc.)
            % The oracle passes through a double of ~6.3e8 s, so it is exact only to eps of that
            % (~120 ns); referenceSystem converts in int64. The bound is that eps, computed.
            ev = res2.flags.events;
            got = res2.diagnostics.jumps;
            if istable(ev) && height(ev) > 0 && any(ev.class == "jump")
                tj = ev.tStart(ev.class == "jump");
                if isempty(tj.TimeZone), tj.TimeZone = 'UTC'; end
                sec = seconds(tj - datetime(2000, 1, 1, 12, 0, 0, 'TimeZone', 'UTC')) + 18;
                expected = reshape(unique(round(sec * 1e9)), 1, []);
                tc.verifySize(got, size(expected));
                tc.verifyBoundAndRecord("calibrate jump epochs vs UTC->GPS arithmetic", ...
                    max(abs(double(got) - expected)), 2 * eps(max(sec)) * 1e9 + 1, "ns", ...
                    "oracle in double seconds: eps of ~6.3e8 s, plus the 1 ns rounding");
            else
                tc.verifyEmpty(got);
            end
            tc.verifyEmpty(res2.calibration.bias.jumps, 'jumps are fitted only when asked for');
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function calibrateMetamorphic(tc)
            % Two entry points, one code path: the Model reference handed back in as an Arc gives the
            % same calibration bitwise; only the attribution, which needs the effects, differs.
            fit = {'NBoot', 0};
            rm = NonGravForces.app.calibrate(tc.acc, tc.Mfast, 'Orbit', tc.orb, 'Attitude', tc.att, 'Fit', fit);
            ra = NonGravForces.app.calibrate(tc.acc, rm.diagnostics.rawReference, 'Fit', fit);
            tc.verifyEqual(ra.calibration.scale, rm.calibration.scale);
            tc.verifyEqual(ra.calibration.bias.coef, rm.calibration.bias.coef);
            tc.verifyEqual(ra.attribution.verdict, repmat("unavailable", 1, 3));
            tc.verifyEqual(rm.attribution.effects, ["drag" "radiation" "thermal"]);
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function calibrateSpeed(tc)
            % calibrate's OWN cost: its parts are run (whose change-point QC is quadratic in the worst
            % case) and fitCalibration, each with its own speed claim; what calibrate adds -- the
            % jump conversion, the merge of options, the attribution -- must stay small beside them.
            % A ratio of three timings on the same host (rule T5).
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 8);
            K = timeSeriesPreparation.crnKernel("ACC");
            cfg = set(timeSeriesPreparation.Config(obs), "run", "filterReference", true);
            calib = @() NonGravForces.app.calibrate(obs, ref, 'Fit', {'NBoot', 0});
            parts = @() NonGravForces.app.fitCalibration(timeSeriesPreparation.run(obs, Reference = ref, Config = cfg).arc, ...
                ref, 'Kernel', K, 'NBoot', 0);
            for k = 1:2, calib(); parts(); end
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; calib(); t1(k) = toc(t0);
                t0 = tic; parts(); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("calibrate / (run + fitCalibration)", median(t1) / median(t2), 1.5, "-", ...
                "calibrate adds conversions and a merge; an Arc reference has no attribution to fit");
        end
    end
end
