classdef tNgfRevolutions < NgfTestCase
%TNGFREVOLUTIONS Tests of the per-revolution fit and its attribution (N-03c): app.perRevolution, app.attributeRevolutions.
%
%   perRevolution's truth is app.syntheticPair; attributeRevolutions' is a
%   per-revolution table built here with a planted driver law (the v3 plan's
%   WP-17 cases), and timeSeriesPreparation.regressDrivers called directly is
%   the composition oracle. The orbit drivers are checked against closed-form
%   geometry: a circular orbit whose plane holds the Sun has beta' = 0 and
%   spends asin(Re / R) / pi of a revolution in the cylindrical shadow.
%
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
    properties (Constant)
        Quiet = 1e-30 * ones(3, 1)
        Fam = struct('f', "density", 'b', "geometry")
    end
    properties
        sw
    end
    methods (TestClassSetup)
        function setup(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
        end
    end
    methods
        function T = revTable(~, s)
            %REVTABLE A perRevolution-shaped table: the given x scale, y and z not estimated.
            R = numel(s);
            T = table(zeros(R, 1, 'int64'), zeros(R, 1, 'int64'), int64(1:R)' * int64(5640e9), [s(:), nan(R, 2)], ...
                'VariableNames', {'t0', 't1', 'tMid', 'scale'});
        end
        function [T, D] = densityCase(tc, seed, R)
            %DENSITYCASE v3 WP-17: a per-revolution scale 0.9 + 0.05 z(F10.7) + noise, beside a null driver.
            rs = RandStream("twister", "Seed", seed);
            f = 70 + 10 * randn(rs, R, 1); b = 0.3 * randn(rs, R, 1);
            T = tc.revTable(0.9 + 0.05 * (f - mean(f)) / std(f) + 0.005 * randn(rs, R, 1));
            D = table(f, b);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function perRevolutionContract(tc)
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.2, 'Seed', 2);
            [T, info] = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
            tc.verifyEqual(T.Properties.VariableNames, {'window', 't0', 't1', 'tMid', 'usable', 'scale', 'sigmaScale', 'bias', 'sigmaBias', 'gate'});
            R = height(T);
            tc.verifyEqual(R, ceil(0.2 * 86400 / 5640), 'contiguous windows of one period from the first usable epoch');
            tc.verifyClass(T.tMid, 'int64');
            tc.verifyEqual(T.t1 - T.t0, repmat(int64(5640e9), R, 1));
            tc.verifySize(T.scale, [R 3]); tc.verifySize(T.gate, [R 3]);
            tc.verifyTrue(all(isnan(T.scale(:, 2))), 'S_yy is held, so it is not a per-revolution parameter');
            tc.verifyEqual(T.gate(end, :), repmat("window", 1, 3), 'the short last window is gated');
            tc.verifyTrue(all(isnan([T.scale(end, :), T.bias(end, :)])), 'a gated window is NaN, never a number');
            tc.verifyNotEmpty(char(info.messages(end)));
            tc.verifyEqual(T.Properties.UserData.sigma, "formal");
            tc.verifySize(info.calibrations, [1 R]);
            % effects add amp_<name>
            eff = struct('name', "drag", 'arc', ref);
            T2 = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0, 'Effects', eff);
            tc.verifySize(T2.amp_drag, [R 3]);
            tc.verifyGreaterThan(T2.amp_drag(1, 1), 0);
            % refusals
            tc.verifyErrorId(@() NonGravForces.app.perRevolution(5, ref), "NonGravForces:perRevolution:inputArc");
            tc.verifyErrorId(@() NonGravForces.app.perRevolution(obs, ref, 'Fit', {'BiasBasis', "linear"}), "NonGravForces:perRevolution:inputFit");
            tc.verifyErrorId(@() NonGravForces.app.perRevolution(obs, ref, 'Fit', {'Robust'}), "NonGravForces:perRevolution:inputFit");
            tc.verifyErrorId(@() NonGravForces.app.perRevolution(obs, ref, 'Effects', struct('name', "x", 'arc', 1)), "NonGravForces:perRevolution:inputEffects");
            tc.verifyErrorId(@() NonGravForces.app.perRevolution(obs.derive(obs.t, nan(size(obs.obs)), "test.nan"), ref), ...
                "NonGravForces:perRevolution:dataEmpty");
        end

        function attributeRevolutionsContract(tc)
            [T, D] = tc.densityCase(1, 40);
            [A, info] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100);
            tc.verifyEqual(A.Properties.VariableNames, {'parameter', 'verdict', 'reason', 'significant', 'r2', 'nRows'});
            tc.verifyEqual(A.parameter, "S_xx", 'only estimated scales are attributed');
            tc.verifyEqual(sort(string(fieldnames(info)))', sort(["drivers" "families" "dropped" "regression"]));
            tc.verifyEqual(info.drivers.Properties.VariableNames, {'time', 'f', 'b'});
            tc.verifyEqual(info.families.f107, "density"); tc.verifyEqual(info.families.time, "sensor");
            tc.verifyEqual(info.families.f, "density", 'a caller family extends the defaults');
            % space-weather drivers are the index file's values at tMid. 40 revolutions are 2.6 days, so
            % the daily indices are steps and some bootstrap replicates are rank deficient inside
            % regressDrivers -- documented in the unit's ACCURACY; only that warning is silenced, here.
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("MATLAB:rankDeficientMatrix"));
            T.tMid = int64(7305 * 86400 - 43200 + 18) * int64(1e9) + int64(0:39)' * int64(5640e9);
            [~, i2] = NonGravForces.app.attributeRevolutions(T, 'SpaceWeather', tc.sw, 'NBoot', 100);
            d = tc.sw.query(referenceSystem.gpsNanosecondsToEpoch(T.tMid));
            tc.verifyEqual(i2.drivers.f107, d.f107(:)); tc.verifyEqual(i2.drivers.ap, d.ap(:));
            % refusals
            tc.verifyErrorId(@() NonGravForces.app.attributeRevolutions(T(:, {'t0', 't1'})), "NonGravForces:attributeRevolutions:inputTable");
            tc.verifyErrorId(@() NonGravForces.app.attributeRevolutions(T, 'Drivers', D(1:5, :)), "NonGravForces:attributeRevolutions:inputDrivers");
            tc.verifyErrorId(@() NonGravForces.app.attributeRevolutions(T, 'Orbit', struct('epoch', 1)), "NonGravForces:attributeRevolutions:inputOrbit");
            one = T; one.scale(2:end, 1) = NaN;
            tc.verifyErrorId(@() NonGravForces.app.attributeRevolutions(one), "NonGravForces:attributeRevolutions:dataNoDriver");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function perRevolutionReference(tc)
            % The composition: each window is app.fitCalibration on that window's rows, bitwise.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.15, 'Seed', 3);
            T = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
            for w = find(T.gate(:, 1) == "pass")'
                rows = obs.t >= T.t0(w) & obs.t < T.t1(w);
                cal = NonGravForces.app.fitCalibration(obs.derive(obs.t(rows), obs.obs(rows, :), "test.window"), ref, ...
                    'BiasBasis', "constant", 'NBoot', 0, 'Seed', w);
                tc.verifyEqual(T.scale(w, [1 3]), diag(cal.scale([1 3], [1 3]))', "window " + w);
                tc.verifyEqual(T.bias(w, :), cal.bias.coef(:, 1)', "window " + w);
            end
        end

        function attributeRevolutionsReference(tc)
            % The composition: the coefficients are regressDrivers' on the same drivers, bitwise ...
            [T, D] = tc.densityCase(2, 40);
            [~, info] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100, 'Seed', 3);
            Z = [info.drivers.time, D.f, D.b];
            res = timeSeriesPreparation.regressDrivers(T.scale(:, 1), Z, Names = "S_xx", DriverNames = ["time" "f" "b"], ...
                NBoot = 100, Seed = 3, Stage = "attributeRevolutions.S_xx");
            tc.verifyEqual(info.regression{1}.coef, res.coef);
            tc.verifyEqual(info.regression{1}.sigma, res.sigma);
            % ... and the orbit drivers are closed-form geometry: a circular orbit whose plane holds the
            % Sun has beta' = 0 and spends asin(Re/R)/pi of each revolution in the cylindrical umbra
            % (the conical model differs by the penumbra, a few seconds per transit).
            mu = 3.986004418e14; Re = 6378137; R = Re + 490e3; n = sqrt(mu / R^3); P = 2 * pi / n;
            day = 58849; sec0 = 3600;
            s = referenceSystem.geocentricBodyState("sun", referenceSystem.daySecondToEpoch([day sec0], "GPS")); s = s / norm(s);
            e1 = -s; e2 = cross(s, [0; 0; 1]); e2 = e2 / norm(e2);
            t = 0:10:3 * P;
            r = R * (cos(n * t) .* e1 + sin(n * t) .* e2); v = R * n * (-sin(n * t) .* e1 + cos(n * t) .* e2);
            ep = referenceSystem.daySecondToEpoch([repmat(day, numel(t), 1), sec0 + t(:)], "GPS");
            tns = referenceSystem.epochTo(ep, "gpsNanoseconds");
            T3 = tc.revTable(0.9 + 1e-3 * (1:3)');
            T3.t0 = tns(1) + int64(0:2)' * int64(round(P * 1e9)); T3.t1 = T3.t0 + int64(round(P * 1e9)); T3.tMid = T3.t0 + int64(round(P * 5e8));
            [~, i3] = NonGravForces.app.attributeRevolutions(T3, 'Orbit', struct('epoch', ep, 'r', r, 'v', v), 'NBoot', 100);
            tc.verifyBoundAndRecord("attributeRevolutions beta' of a Sun-plane orbit", max(abs(i3.drivers.betaPrime)), 5e-3, "rad", ...
                "the Sun moves ~0.04 deg/h out of the plane fixed at the first epoch; 3 revolutions ~ 4.7 h");
            tc.verifyBoundAndRecord("attributeRevolutions eclipse fraction vs asin(Re/R)/pi", ...
                max(abs(i3.drivers.eclipseFraction - asin(Re / R) / pi)), 0.01, "-", ...
                "conical vs cylindrical shadow and a 10 s sampling: the penumbra is ~10 s of a 5650 s revolution");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function perRevolutionPrecision(tc)
            % Noise-free with a constant planted scale and bias: every passing window returns the truth;
            % the noise ASD is 1e-30, so what remains is the rounding of the filter and the solve.
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.15, 'NoiseAsd', tc.Quiet, 'BiasDrift', zeros(3, 1), ...
                'Spikes', 0, 'Gap_s', 0);
            T = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
            ok = T.gate(:, 1) == "pass";
            tc.verifyGreaterThanOrEqual(nnz(ok), 2);
            errS = max(abs(T.scale(ok, [1 3]) - diag(truth.scale([1 3], [1 3]))'), [], 'all');
            errB = max(abs(T.bias(ok, :) - truth.bias0'), [], 'all') / max(abs(truth.bias0));
            tc.verifyBoundAndRecord("perRevolution noise-free scale error", errS, 1e-9, "-", ...
                "rounding of one revolution's CRN-filtered solve; fitCalibration's computed bound is ~1e-12 on longer arcs");
            tc.verifyBoundAndRecord("perRevolution noise-free bias error, relative", errB, 1e-6, "-", ...
                "the bias is the small parameter beside a 2e-7 signal; same rounding, relative to 2e-7");
        end

        function attributeRevolutionsPrecision(tc)
            % v3 WP-17's cases over seeded realisations: a scale driven by F10.7 is "density" with the
            % F10.7 coefficient within 20 %; a pure secular drift is "sensor". A null driver sits beside
            % both, which is what the 3-sigma default is for.
            right = 0; worstCoef = 0; seeds = 1:4; R = 60;
            for s = seeds
                [T, D] = tc.densityCase(s, R);
                [A, info] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 300, 'Seed', s);
                right = right + (A.verdict == "density");
                truthCoef = 0.05 / std(D.f);
                worstCoef = max(worstCoef, abs(info.regression{1}.coef(3) / truthCoef - 1));
                rs = RandStream("twister", "Seed", 100 + s);
                Ts = tc.revTable(0.9 + 1e-3 * (1:R)' + 0.001 * randn(rs, R, 1));
                As = NonGravForces.app.attributeRevolutions(Ts, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 300, 'Seed', s);
                right = right + (As.verdict == "sensor");
            end
            tc.verifyBoundAndRecord("attributeRevolutions wrong verdicts, fraction", 1 - right / (2 * numel(seeds)), 0, "-", ...
                "4 seeds x 2 cases, 60 revolutions, 3-sigma significance");
            tc.verifyBoundAndRecord("attributeRevolutions F10.7 coefficient, relative error", worstCoef, 0.2, "-", "v3 WP-17: within 20 %");
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function perRevolutionVectorization(tc)
            % A window does not see the others: the first windows of a longer arc equal those of the
            % arc cut after them, bitwise.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.2, 'Seed', 4, 'Gap_s', 0);
            T = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
            cut = obs.t < T.t1(2);
            Tc = NonGravForces.app.perRevolution(obs.derive(obs.t(cut), obs.obs(cut, :), "test.cut"), ref, 'NBoot', 0);
            tc.verifyEqual(Tc.scale(1:2, :), T.scale(1:2, :));
            tc.verifyEqual(Tc.bias(1:2, :), T.bias(1:2, :));
        end

        function attributeRevolutionsVectorization(tc)
            % Two axes attributed together equal each attributed alone, bitwise.
            [T, D] = tc.densityCase(5, 40);
            T.scale(:, 3) = flipud(T.scale(:, 1));
            [~, both] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100);
            for i = [1 3]
                Ti = T; Ti.scale(:, setdiff(1:3, i)) = NaN;
                [~, one] = NonGravForces.app.attributeRevolutions(Ti, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100);
                tc.verifyEqual(one.regression{1}.coef, both.regression{1 + (i == 3)}.coef);
            end
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function perRevolutionRobustness(tc)
            % A window mostly lost to a gap is "window"; with the window gate lowered, the estimator's own
            % minimum-data gate refuses it and the reason is recorded -- in neither case a number.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.2, 'Seed', 5);
            y = obs.obs; T0 = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
            hole = obs.t >= T0.t0(2) & obs.t < T0.t0(2) + int64(4000e9);
            y(hole, :) = NaN;
            holed = obs.derive(obs.t, y, "test.hole");
            T = NonGravForces.app.perRevolution(holed, ref, 'NBoot', 0);
            tc.verifyEqual(T.gate(2, 1), "window"); tc.verifyTrue(all(isnan(T.scale(2, :))));
            [T2, i2] = NonGravForces.app.perRevolution(holed, ref, 'NBoot', 0, 'MinWindowFraction', 0.05);
            tc.verifyEqual(T2.gate(2, 1), "insufficientData"); tc.verifyTrue(all(isnan(T2.scale(2, :))));
            tc.verifyTrue(contains(i2.messages(2), "weight"));
            tc.verifyEqual(T2.gate(1, 1), "pass", 'the other windows are untouched');
        end

        function attributeRevolutionsRobustness(tc)
            % Collinear drivers are a verdict, not an error; a constant driver is dropped and named; gated
            % revolutions (NaN) are left out.
            [T, D] = tc.densityCase(6, 40);
            D.b2 = 2 * D.b + 1e-6 * (1:40)';
            A = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100);
            tc.verifyEqual(A.verdict, "unresolved");
            tc.verifyTrue(startsWith(A.reason, "collinear"));
            D = removevars(D, 'b2'); D.c = ones(40, 1);
            T.scale([3 7], 1) = NaN;
            [A2, info] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100);
            tc.verifyEqual(info.dropped, "c");
            tc.verifyEqual(A2.nRows, 38);
            tc.verifyEqual(A2.verdict, "density");
            % Fewer than 20 revolutions: the automatic block length refuses (found by N-05, which first
            % tried one day = 16); an explicit 'BlockLength' is accepted.
            short = T(1:16, :);
            tc.verifyErrorId(@() NonGravForces.app.attributeRevolutions(short, 'Drivers', D(1:16, :), 'NBoot', 50), ...
                "timeSeriesPreparation:blockLength:tooShort");
            A3 = NonGravForces.app.attributeRevolutions(short, 'Drivers', D(1:16, :), 'NBoot', 50, 'BlockLength', 3);
            tc.verifyEqual(A3.nRows, 14, "16 revolutions less the two gated above");
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function perRevolutionMetamorphic(tc)
            % Doubling the observation doubles every estimated parameter of the free axes, bitwise (the
            % equilibration takes a power of two out exactly).
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.15, 'Seed', 6);
            T1 = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
            T2 = NonGravForces.app.perRevolution(obs.derive(obs.t, 2 * obs.obs, "test.double"), ref, 'NBoot', 0);
            ok = T1.gate(:, 1) == "pass";
            tc.verifyEqual(T2.scale(ok, [1 3]), 2 * T1.scale(ok, [1 3]));
            tc.verifyEqual(T2.bias(ok, [1 3]), 2 * T1.bias(ok, [1 3]));
        end

        function attributeRevolutionsMetamorphic(tc)
            % A driver in other units (x 10) has its coefficient divided by 10 and the verdict unchanged.
            [T, D] = tc.densityCase(7, 40);
            [A1, i1] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', tc.Fam, 'NBoot', 100);
            D10 = D; D10.f = 10 * D.f;
            [A2, i2] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D10, 'Families', tc.Fam, 'NBoot', 100);
            tc.verifyEqual(A2.verdict, A1.verdict);
            rel = abs(10 * i2.regression{1}.coef(3) / i1.regression{1}.coef(3) - 1);
            tc.verifyBoundAndRecord("attributeRevolutions driver x10 -> coefficient / 10", rel, 1e-12, "-", "a rescaled column in a QR solve; rounding only");
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function perRevolutionSpeed(tc)
            % One fit per window: twice the revolutions, about twice the cost.
            [o1, r1] = NonGravForces.app.syntheticPair('Days', 0.15);
            [o2, r2] = NonGravForces.app.syntheticPair('Days', 0.3);
            NonGravForces.app.perRevolution(o1, r1, 'NBoot', 0); NonGravForces.app.perRevolution(o2, r2, 'NBoot', 0);
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.perRevolution(o1, r1, 'NBoot', 0); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.perRevolution(o2, r2, 'NBoot', 0); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("perRevolution, 2x arc / 1x arc", median(t2) / median(t1), 4, "-", "linear in windows would be 2");
        end

        function attributeRevolutionsSpeed(tc)
            % The bootstrap dominates and is linear in the revolutions: twice the rows, about twice the cost.
            [T1, D1] = tc.densityCase(8, 60); [T2, D2] = tc.densityCase(8, 120);
            NonGravForces.app.attributeRevolutions(T1, 'Drivers', D1, 'NBoot', 200); NonGravForces.app.attributeRevolutions(T2, 'Drivers', D2, 'NBoot', 200);
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.attributeRevolutions(T1, 'Drivers', D1, 'NBoot', 200); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.attributeRevolutions(T2, 'Drivers', D2, 'NBoot', 200); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("attributeRevolutions, 2x revolutions / 1x", median(t2) / median(t1), 4, "-", "linear would be 2");
        end
    end
end
