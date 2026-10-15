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
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
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
            [o1, r1] = NonGravForces.app.syntheticPair('Days', 0.05);
            [o2, r2] = NonGravForces.app.syntheticPair('Days', 0.1);
            f = {'Fit', {'NBoot', 0}};
            for k = 1:2, NonGravForces.app.calibrate(o1, r1, f{:}); NonGravForces.app.calibrate(o2, r2, f{:}); end
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.calibrate(o1, r1, f{:}); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.calibrate(o2, r2, f{:}); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("calibrate, 2x arc / 1x arc", median(t2) / median(t1), 4, "-", ...
                "every stage is linear or n log n in the arc; quadratic would be 4");
        end
    end
end
