classdef tNgfCalibration < NgfTestCase
%TNGFCALIBRATION Tests of the calibration truth and its result object (N-02a): app.syntheticPair, app.Calibration.
%
%   The oracle for syntheticPair is the Level-1B chain it claims to imitate,
%   applied by hand: timeSeriesPreparation.crnApply on the model arc, the
%   planted scale and bias added in the order the help states, and the CRN
%   kernel's own DC gain. The kernel is symmetric (F(-n) = F(n); the taps
%   agree with their mirror images to Nf eps, asserted), so a bias linear in
%   time comes through it as dcGain * b exactly in real arithmetic, and every
%   remaining difference is the rounding of one Nf-term dot product. Every
%   bound below is computed from that, Nf * eps * sum|F| * max|raw|, and
%   never written as a literal.
%
%   The oracle for Calibration.basis("pspline") is the closed form of the
%   uniform cubic B-spline (de Boor 1978, ch. IX), written out per segment,
%   which the class does not use: it evaluates Cox-de Boor.
%
% Created 2026-10-14 by Claude (Anthropic) for Matthias Weigelt.
    properties (Constant)
        Days = 0.05          % 4320 s: room for the 140.7 s window, three spikes and the gap
        Quiet = 1e-30 * ones(3, 1)   % a noise ASD far below the rounding bound, so the chain is deterministic algebra
    end
    methods
        function bound = chainBound(~, truth, ref, factor)
            %CHAINBOUND factor * Nf * eps * sum|F| * max|raw|, the rounding of one CRN output.
            S = truth.scale;
            maxRaw = max(abs(ref.obs), [], 'all') * norm(S, inf) + max(abs(truth.bias), [], 'all') ...
                + max([abs(truth.spikes.amplitude); 0]);
            bound = factor * truth.kernel.Nf * eps * sum(abs(truth.kernel.F)) * maxRaw;
        end
        function arc = arc3(~, t, y, units)
            %ARC3 A three-column acceleration arc built the way the readers build one.
            if nargin < 4, units = 'm/s2'; end
            names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
            u = struct();
            for k = 1:3, u.(names{k}) = units; end
            aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
                "variableUnits", u, "fileName", "ACT1B_2020-01-01_C_04.txt", "product", "ACT1B");
            arc = timeSeriesPreparation.Arc(t(:), y, aux);
        end
        function cal = truthCal(~, obs, truth)
            %TRUTHCAL The calibration the truth defines: its scale and its linear bias from Start.
            bias = struct('basis', "linear", 't0', obs.t(1), 'knots', int64([]), ...
                'coef', [truth.bias0, truth.biasDrift_per_day]);
            cal = NonGravForces.app.Calibration(struct('scale', truth.scale, 'bias', bias));
        end
        function cal = curveCal(tc, slope)
            %CURVECAL splineCal(8) with a linear scale curve S(t) = S0 + slope * days on t0, tRef = t0.
            base = tc.splineCal(8);
            S0 = base.scale;
            sc = struct('basis', "linear", 't0', base.bias.t0, 'knots', int64([]), 'tRef', base.bias.t0, ...
                'coef', [S0(:), slope(:)], 'cov', {{[4e-8 1e-9; 1e-9 9e-8], [], [], [], [], [], [], [], [1e-8 0; 0 2e-8]}});
            cal = NonGravForces.app.Calibration(struct('scale', S0, 'scaleCurve', sc, 'bias', base.bias));
        end
        function [cal, R, S1] = tempCal(tc)
            %TEMPCAL regCal with a temperature-dependent scale on its 'temp' series (N-11b), reference 22 C.
            [base, R] = tc.regCal();
            S1 = [1e-3 2e-4 0; 0 0 0; -1e-4 0 -2e-3];
            cv = nan(9); est = find(S1(:) ~= 0); cv(est, est) = 1e-12 * eye(numel(est));
            st = struct('name', "temp", 'units', "C", 'reference', 22, 'coef', S1, 'cov', cv);
            cal = NonGravForces.app.Calibration(struct('scale', base.scale, 'bias', base.bias, 'regressors', base.regressors, 'scaleTemperature', st));
        end
        function [cal, R, C, m] = regCal(tc)
            %REGCAL splineCal(8) with two regressor terms and the hourly series that drives them.
            base = tc.splineCal(8);
            t = base.bias.t0 + int64(0:599)' * int64(6e9);
            z = [22 + 0.1 * sin((1:600)' / 40), 0.02 * cos((1:600)' / 13)];
            nm = {'temp', 'curr'}; u = struct('temp', 'C', 'curr', 'A');
            R = timeSeriesPreparation.Arc(t, z, struct("obsNames", {nm}, "variableNames", {nm}, "variableFieldNames", {nm}, ...
                "variableUnits", u, "fileName", "REG", "product", "REG"));
            C = [3e-9 2e-8; -1e-9 0; 5e-10 -4e-8]; m = [22 0];
            g = struct('names', ["temp" "curr"], 'units', ["C" "A"], 'mean', m, 'coef', C, 'cov', {{eye(2) * 1e-20, eye(2) * 1e-20, eye(2) * 1e-20}});
            cal = NonGravForces.app.Calibration(struct('scale', base.scale, 'bias', base.bias, 'regressors', g));
        end
        function cal = splineCal(~, P)
            %SPLINECAL A calibration with a cubic B-spline bias on uniform hourly knots (P coefficients).
            t0 = int64(631152018) * int64(1e9);
            knots = t0 + int64(-3:P) * int64(3600e9);
            coef = [linspace(1e-7, 2e-7, P); linspace(-1e-7, 0, P); 5e-8 * cos(1:P)];
            S = [0.96 1e-3 -2e-4; 5e-4 0.98 1e-3; -1e-3 2e-4 0.95];
            cal = NonGravForces.app.Calibration(struct('scale', S, ...
                'bias', struct('basis', "pspline", 't0', t0, 'knots', knots, 'coef', coef)));
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function syntheticPairContract(tc)
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', tc.Days);
            tc.verifyClass(obs, 'timeSeriesPreparation.Arc');
            tc.verifyClass(ref, 'timeSeriesPreparation.Arc');
            n = round(tc.Days * 86400);
            tc.verifySize(obs.obs, [n 3]); tc.verifySize(ref.obs, [10 * n 3]);
            tc.verifyEqual(unique(diff(obs.t)), int64(1e9));
            tc.verifyEqual(unique(diff(ref.t)), int64(1e8));
            tc.verifyEqual(obs.t, ref.t(1:10:end));
            tc.verifyEqual(obs.names, ["lin_accl_x" "lin_accl_y" "lin_accl_z"]);
            tc.verifyEqual(obs.units, repmat("m/s2", 1, 3));
            tc.verifyEqual(sort(string(fieldnames(truth)))', sort(["scale" "scaleDrift_per_day" "regressors" "regressorCoupling" "bias0" "biasDrift_per_day" ...
                "quadratic" "bias" "noise" "spikes" "gap" "kernel" "orbitPeriod_s" "seed" "provenance" "scaleTemperature" "scaleTemperatureReference"]));
            tc.verifyEqual(truth.scaleTemperature, zeros(3)); tc.verifyTrue(isnan(truth.scaleTemperatureReference));   % nothing planted by default
            tc.verifyErrorId(@() NonGravForces.app.syntheticPair('Days', tc.Days, 'ScaleTemperature', eye(3) * 1e-3), "NonGravForces:syntheticPair:inputShape");
            tc.verifySize(truth.bias, [3 n]);
            tc.verifyEqual(height(truth.spikes), 3);
            tc.verifyClass(truth.spikes.t, 'int64'); tc.verifyClass(truth.gap, 'int64');
            tc.verifySize(truth.gap, [1 2]);
            tc.verifyEqual(truth.kernel.preset, "ACC");
            tc.verifyErrorId(@() NonGravForces.app.syntheticPair('Days', tc.Days, 'Scale', eye(2)), "NonGravForces:syntheticPair:inputShape");
            tc.verifyErrorId(@() NonGravForces.app.syntheticPair('Days', tc.Days, 'Bias', [1 2]), "NonGravForces:syntheticPair:inputShape");
            tc.verifyErrorId(@() NonGravForces.app.syntheticPair('Days', tc.Days, 'ScaleDrift', [1 2 3]), "NonGravForces:syntheticPair:inputShape");
            tc.verifyEqual(truth.scaleDrift_per_day, zeros(3));
            tc.verifyEmpty(truth.regressors, 'no regressor is planted by default (N-11a)');
            tc.verifyErrorId(@() NonGravForces.app.syntheticPair('Days', 0.001), "NonGravForces:syntheticPair:inputRange");
            tc.verifyErrorId(@() NonGravForces.app.syntheticPair('Days', 0.01, 'Spikes', 40), "NonGravForces:syntheticPair:inputRange");
        end
        function calibrationContract(tc)
            cal = tc.splineCal(8);
            tc.verifyEqual(cal.convention, "observation");
            tc.verifyEqual(cal.SchemaVersion, "4");
            tc.verifyFalse(cal.isTimeVarying());
            tc.verifyTrue(all(isnan(cal.quadratic)));
            tc.verifyTrue(all(isnan(cal.covScale), 'all'));
            tc.verifySize(cal.biasAt(cal.bias.t0 + int64(0:9)' * int64(60e9)), [3 10]);
            T = cal.summaryTable();
            tc.verifySize(T, [9 + 3 * 8, 3]);
            tc.verifyEqual(T.Properties.VariableNames, {'parameter', 'value', 'sigma'});
            tc.verifyEqual(cal.invert().convention, "correction");
            good = struct('scale', eye(3), 'bias', struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'coef', zeros(3, 1)));
            bad = good; bad = rmfield(bad, 'scale');
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.scale = eye(2);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.bias.coef = zeros(3, 2);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.convention = "zarm";
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.covScale = eye(3);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.bias.basis = "fourier";
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputBasis");
            bad = good; bad.bias.basis = "hat"; bad.bias.knots = int64([5 3]); bad.bias.coef = zeros(3, 2);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputBasis");
            bad = good; bad.bias.basis = "pspline"; bad.bias.knots = int64(1:7); bad.bias.coef = zeros(3, 3);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputBasis");
            tc.verifyErrorId(@() NonGravForces.app.Calibration.basis(good.bias, 0:3), "NonGravForces:Calibration:inputBasis");
            tc.verifyErrorId(@() cal.apply(zeros(5, 3)), "NonGravForces:Calibration:inputArc");
            t = cal.bias.t0 + int64(0:4)' * int64(1e9);
            tc.verifyErrorId(@() cal.apply(tc.arc3(t, zeros(5, 3), 'm/s')), "NonGravForces:Calibration:inputArc");
            k2 = good; k2.quadratic = [0; 1e-3; 0];
            tc.verifyErrorId(@() apply(NonGravForces.app.Calibration(k2), tc.arc3(t, zeros(5, 3))), "NonGravForces:Calibration:unsupportedQuadratic");
            sing = good; sing.scale = [1 0 0; 0 1 0; 1 1 0];
            singCal = NonGravForces.app.Calibration(sing);
            tc.verifyErrorId(@() singCal.apply(tc.arc3(t, zeros(5, 3))), "NonGravForces:Calibration:numericSingular");
            tc.verifyErrorId(@() singCal.invert(), "NonGravForces:Calibration:numericSingular");
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            notCal = fullfile(folder, 'other.mat'); x = 1; save(notCal, 'x');
            tc.verifyErrorId(@() NonGravForces.app.Calibration.load(notCal), "NonGravForces:Calibration:dataSchema");
        end

        function calibrationTimeVaryingContract(tc)
            % N-03a: a scale curve -- shapes, the one-authority check of scale against the curve,
            % the refused inversion, and the drift report's fields.
            slope = zeros(3); slope(1, 1) = 0.01; slope(3, 3) = -0.02;
            cal = tc.curveCal(slope);
            tc.verifyTrue(cal.isTimeVarying());
            t = cal.bias.t0 + int64(0:4)' * int64(3600e9);
            tc.verifySize(cal.scaleAt(t), [3 3 5]);
            tc.verifySize(cal.scaleSigmaAt(t), [3 3 5]);
            sg = cal.scaleSigmaAt(t);
            tc.verifyTrue(all(isnan(sg(2, 2, :))), 'an entry without a covariance has no band');
            d = cal.drift();
            tc.verifyEqual(sort(string(fieldnames(d)))', sort(["statistic" "dof" "pValue" "significant" "alpha"]));
            tc.verifyClass(cal.isDriftSignificant(0.01), 'logical');
            T = cal.summaryTable();
            tc.verifySize(T, [9 + 3 * 8 + 2 * 2, 3]);
            tc.verifyErrorId(@() cal.invert(), "NonGravForces:Calibration:unsupportedTimeVarying");
            good = struct('scale', cal.scale, 'scaleCurve', cal.scaleCurve, 'bias', cal.bias);
            bad = good; bad.scale(1, 1) = bad.scale(1, 1) + 1e-6;
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.scaleCurve.coef = zeros(9, 3);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.scaleCurve.cov = {eye(3)};
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
            bad = good; bad.scaleCurve.basis = "fourier";
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputBasis");
            sing = good; sing.scaleCurve.coef(:, 2) = 0; sing.scaleCurve.coef([3 6 9], 2) = -cal.scale(3, :).' / 0.5;   % row 3 of S(t) vanishes at day 0.5
            singCal = NonGravForces.app.Calibration(sing);
            tc.verifyErrorId(@() singCal.apply(tc.arc3(cal.bias.t0 + int64(43200e9), zeros(1, 3))), "NonGravForces:Calibration:numericSingular");
        end

        function calibrationScaleTemperatureContract(tc)
            % N-11b: the temperature-dependent scale travels with the calibration; apply needs the temperature,
            % invert refuses, and a scale that depends on time AND temperature is refused.
            [cal, R, S1] = tc.tempCal();
            tc.verifyTrue(cal.isTemperatureDependent()); tc.verifyFalse(cal.isTimeVarying());
            S3 = cal.scaleAtTemperature([22 23]);
            tc.verifySize(S3, [3 3 2]);
            tc.verifyEqual(S3(:, :, 1), cal.scale); tc.verifyEqual(S3(:, :, 2), cal.scale + S1, 'AbsTol', 1e-15);
            t = R.t(10:20); arc = tc.arc3(t, zeros(numel(t), 3));
            tc.verifyErrorId(@() cal.apply(arc), "NonGravForces:Calibration:inputRegressors");
            other = timeSeriesPreparation.Arc(R.t, R.obs(:, 2), struct("obsNames", {{'curr'}}, "variableNames", {{'curr'}}, ...
                "variableFieldNames", {{'curr'}}, "variableUnits", struct('curr', 'A'), "fileName", "REG", "product", "REG"));
            tc.verifyErrorId(@() cal.apply(arc, 'Regressors', other), "NonGravForces:Calibration:inputRegressors");
            tc.verifyClass(cal.apply(arc, 'Regressors', R), 'timeSeriesPreparation.Arc');
            tc.verifyErrorId(@() cal.invert(), "NonGravForces:Calibration:unsupportedTemperatureScale");
            T = cal.summaryTable();
            tc.verifyEqual(sort(T.parameter(startsWith(T.parameter, "dS_")))', sort(["dS_11_dT" "dS_31_dT" "dS_12_dT" "dS_33_dT"]));
            f = [tempname '.mat']; cal.save(f); c2 = NonGravForces.app.Calibration.load(f); delete(f);
            tc.verifyEqual(c2.scaleTemperature, cal.scaleTemperature);
            tc.verifyEqual(NonGravForces.app.Calibration.SchemaVersion, "4");
            sc = struct('basis', "linear", 't0', cal.bias.t0, 'knots', int64([]), 'tRef', cal.bias.t0, 'coef', [cal.scale(:), zeros(9, 1)]);
            tc.verifyErrorId(@() NonGravForces.app.Calibration(struct('scale', cal.scale, 'bias', cal.bias, 'scaleCurve', sc, ...
                'scaleTemperature', cal.scaleTemperature)), "NonGravForces:Calibration:unsupportedTemperatureScale");
            tc.verifyErrorId(@() NonGravForces.app.Calibration(struct('scale', cal.scale, 'bias', cal.bias, ...
                'scaleTemperature', struct('name', "temp", 'reference', 22, 'coef', zeros(2)))), "NonGravForces:Calibration:inputFit");
            plain = NonGravForces.app.Calibration(struct('scale', cal.scale, 'bias', cal.bias));
            tc.verifyFalse(plain.isTemperatureDependent());
        end
        function calibrationRegressorContract(tc)
            % N-11a: measured regressor terms travel with the calibration and apply refuses without them.
            [cal, R, C] = tc.regCal();
            tc.verifyEqual(cal.regressors.names, ["temp" "curr"]);
            tc.verifyEqual(cal.regressors.units, ["C" "A"]);
            t = R.t(10:20);
            arc = tc.arc3(t, zeros(numel(t), 3));
            tc.verifyErrorId(@() cal.apply(arc), "NonGravForces:Calibration:inputRegressors");
            nm = {'temp'}; one = timeSeriesPreparation.Arc(R.t, R.obs(:, 1), struct("obsNames", {nm}, "variableNames", {nm}, ...
                "variableFieldNames", {nm}, "variableUnits", struct('temp', 'C'), "fileName", "REG", "product", "REG"));
            tc.verifyErrorId(@() cal.apply(arc, 'Regressors', one), "NonGravForces:Calibration:inputRegressors");
            tc.verifySize(cal.apply(arc, 'Regressors', R).obs, [numel(t) 3]);
            T = cal.summaryTable();
            tc.verifySize(T, [9 + 3 * 8 + 3 * 2, 3]);
            inv = cal.invert();
            tc.verifyEqual(inv.regressors.coef, -(cal.scale \ eye(3)) * C, 'RelTol', 1e-12, 'the correction direction maps the coefficients as it maps the bias');
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'reg.mat'); cal.save(f);
            tc.verifyTrue(isequaln(NonGravForces.app.Calibration.load(f), cal));
            bad = struct('scale', eye(3), 'bias', cal.bias, 'regressors', struct('names', "a", 'mean', [1 2], 'coef', zeros(3, 1)));
            tc.verifyErrorId(@() NonGravForces.app.Calibration(bad), "NonGravForces:Calibration:inputFit");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function syntheticPairReference(tc)
            % The observation IS the CRN chain on S*model + b(t): filter the model with crnApply,
            % scale it, and add dcGain * b, which is exact for a linear b and a symmetric kernel.
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', tc.Days, 'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            K = truth.kernel;
            tc.verifyEqual(K.F, flipud(K.F), 'AbsTol', K.Nf * eps * max(abs(K.F)));   % the symmetry the exactness rests on
            Hm = timeSeriesPreparation.crnApply(ref, K, OutputEpochs = obs.t);
            expected = Hm.obs * truth.scale.' + K.dcGain * truth.bias.';
            ok = all(isfinite(obs.obs), 2);
            tc.verifyEqual(ok, all(isfinite(expected), 2));
            tc.verifyGreaterThan(nnz(ok), 0.9 * numel(ok));
            err = max(abs(obs.obs(ok, :) - expected(ok, :)), [], 'all');
            tc.verifyBoundAndRecord("syntheticPair vs crnApply(S m) + dcGain b", err, tc.chainBound(truth, ref, 4), ...
                "m/s^2", "bound 4 Nf eps sum|F| max|raw|, computed");
        end
        function calibrationReference(tc)
            % Applying the planted calibration to the noise-free observation returns the filtered
            % model, up to the one approximation the help names: the bias is evaluated where the
            % chain filtered it, which leaves exactly (dcGain - 1) S \ b(t).
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', tc.Days, 'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            a = tc.truthCal(obs, truth).apply(obs);
            tc.verifyEqual(a.t, obs.t);
            Hm = timeSeriesPreparation.crnApply(ref, truth.kernel, OutputEpochs = obs.t);
            expected = Hm.obs + (truth.kernel.dcGain - 1) * (truth.scale \ truth.bias).';
            ok = all(isfinite(a.obs), 2);
            err = max(abs(a.obs(ok, :) - expected(ok, :)), [], 'all');
            tc.verifyBoundAndRecord("Calibration.apply on the truth vs the filtered model", err, ...
                cond(truth.scale) * tc.chainBound(truth, ref, 8), "m/s^2", "bound cond(S) 8 Nf eps sum|F| max|raw|");
            % ... and the pspline basis against the closed form of the uniform cubic B-spline.
            cal = tc.splineCal(8);
            t = cal.bias.t0 + int64(0:499)' * int64(36e9);       % [k4, k9): every interior segment
            B = NonGravForces.app.Calibration.basis(cal.bias, t);
            xs = double(t - cal.bias.t0) / 3600e9; fl = floor(xs); u = xs - fl;
            E = zeros(numel(t), 8);
            for r = 1:numel(t)
                E(r, fl(r) + (1:4)) = [(1 - u(r))^3, 3 * u(r)^3 - 6 * u(r)^2 + 4, -3 * u(r)^3 + 3 * u(r)^2 + 3 * u(r) + 1, u(r)^3] / 6;
            end
            tc.verifyBoundAndRecord("Calibration.basis pspline vs closed form", max(abs(B - E), [], 'all'), ...
                1e3 * eps, "-", "Cox-de Boor against the per-segment cubic, values <= 2/3");
            % hat: the identity at the knots.
            h = struct('basis', "hat", 't0', cal.bias.t0, 'knots', cal.bias.t0 + int64(0:4) * int64(3600e9), 'coef', zeros(3, 5));
            tc.verifyEqual(NonGravForces.app.Calibration.basis(h, h.knots), eye(5), 'AbsTol', 4 * eps);
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function syntheticPairPrecision(tc)
            % The spike table is the truth: the run with spikes minus the run without (same Seed,
            % so the same noise and gap) is the CRN filter of the 1 s boxcars the table lists.
            args = {'Days', tc.Days, 'NoiseAsd', tc.Quiet};
            [obsS, ref, truth] = NonGravForces.app.syntheticPair(args{:}, 'Spikes', 3);
            obs0 = NonGravForces.app.syntheticPair(args{:}, 'Spikes', 0);
            boxes = zeros(size(ref.obs));
            for s = 1:height(truth.spikes)
                i0 = find(ref.t == truth.spikes.t(s));
                boxes(i0:i0 + 9, truth.spikes.axis(s)) = truth.spikes.amplitude(s);
            end
            tc.verifyEqual(nnz(boxes), 30);
            Hs = timeSeriesPreparation.crnApply(ref.derive(ref.t, boxes, "test.spikes"), truth.kernel, OutputEpochs = obsS.t);
            d = obsS.obs - obs0.obs;
            ok = all(isfinite(d), 2);
            tc.verifyEqual(isnan(obsS.obs), isnan(obs0.obs));
            err = max(abs(d(ok, :) - Hs.obs(ok, :)), [], 'all');
            tc.verifyBoundAndRecord("syntheticPair spikes vs crnApply(boxcars)", err, tc.chainBound(truth, ref, 4), ...
                "m/s^2", "two runs, each one Nf-term dot product");
            tc.verifyGreaterThan(max(abs(Hs.obs(ok, :)), [], 'all'), 1e3 * err);   % the spikes are there to be seen
        end
        function calibrationPrecision(tc)
            % Both directions are the same algebra: the correction-direction apply of invert()
            % equals the observation-direction apply to cond(S) eps; invert twice is the identity.
            cal = tc.splineCal(8);
            t = cal.bias.t0 + int64(0:599)' * int64(30e9);
            y = [-2e-7 + 1e-8 * sin((1:600)' / 50), 1e-8 * cos((1:600)' / 40), 5e-8 * sin((1:600)' / 70)];
            arc = tc.arc3(t, y);
            a1 = cal.apply(arc).obs; a2 = cal.invert().apply(arc).obs;
            c = cond(cal.scale);
            tc.verifyBoundAndRecord("Calibration observation vs correction apply", max(abs(a1 - a2), [], 'all') / max(abs(a1), [], 'all'), ...
                16 * c * eps * (1 + max(abs(cal.biasAt(t)), [], 'all') / max(abs(y), [], 'all')), "rel", "16 cond(S) eps (1 + |b|/|y|)");
            back = cal.invert().invert();
            tc.verifyBoundAndRecord("Calibration invert twice, scale", max(abs(back.scale - cal.scale), [], 'all'), 16 * c^2 * eps, "-", "16 cond(S)^2 eps");
            tc.verifyBoundAndRecord("Calibration invert twice, bias", max(abs(back.bias.coef - cal.bias.coef), [], 'all') / max(abs(cal.bias.coef), [], 'all'), ...
                16 * c^2 * eps, "rel", "16 cond(S)^2 eps");
            B = NonGravForces.app.Calibration.basis(cal.bias, t(t < cal.bias.knots(9)));
            tc.verifyBoundAndRecord("Calibration pspline partition of unity", max(abs(sum(B, 2) - 1)), 16 * eps, "-", "four non-zero terms");
        end
    end

    methods (Test, TestTags = {'precision'})
        function calibrationTimeVaryingPrecision(tc)
            % apply against an independent per-epoch loop of 3x3 solves of y = S(t) a + b(t), with
            % S(t) written out by hand; and the Wald statistic of the linear basis against its
            % closed form, the squared slope over its variance.
            slope = zeros(3); slope(1, 1) = 0.01; slope(3, 3) = -0.02; slope(1, 2) = 1e-3;
            cal = tc.curveCal(slope);
            t = cal.bias.t0 + int64(0:239)' * int64(600e9);
            days = double(t - cal.bias.t0) / 86400e9;
            a = 1e-7 * [cos((1:240)' / 9), sin((1:240)' / 7), cos((1:240)' / 5)];
            b = cal.biasAt(t);
            y = zeros(240, 3); aRef = zeros(240, 3);
            for k = 1:240
                Sk = cal.scale + slope * days(k);
                y(k, :) = (Sk * a(k, :).' + b(:, k)).';
            end
            got = cal.apply(tc.arc3(t, y)).obs;
            for k = 1:240
                aRef(k, :) = ((cal.scale + slope * days(k)) \ (y(k, :).' - b(:, k))).';
            end
            err = max(abs(got - aRef), [], 'all') / max(abs(a), [], 'all');
            tc.verifyBoundAndRecord("Calibration apply with S(t) vs a per-epoch loop of solves", err, 100 * eps, "-", ...
                "pages vs a loop; both one LU of the same 3x3, differing only in the order S(t) is formed");
            d = cal.drift();
            W = [0.01^2 / 9e-8, 0.02^2 / 2e-8];
            tc.verifyBoundAndRecord("Calibration drift statistic vs slope^2 / var(slope)", ...
                max(abs([d.statistic(1, 1), d.statistic(3, 3)] - W) ./ W), 1e-12, "-", "pinv of a 1x1 is its reciprocal");
            tc.verifyEqual([d.dof(1, 1), d.dof(3, 3)], [1 1]);
            tc.verifyEqual(d.pValue(1, 1), erfc(sqrt(W(1) / 2)), 'RelTol', 1e-12, 'chi-square with 1 dof: P(X > W) = erfc(sqrt(W/2))');
        end

        function calibrationScaleTemperaturePrecision(tc)
            % apply against a per-epoch loop of y = S(T(t)) a + b(t) + C (z(t) - m), every series interpolated linearly.
            [cal, R, S1] = tc.tempCal();
            g = cal.regressors;
            t = R.t(1) + int64(0:239)' * int64(7.3e9);
            a = 1e-7 * [cos((1:240)' / 9), sin((1:240)' / 7), cos((1:240)' / 5)];
            b = cal.biasAt(t);
            tq = double(t - R.t(1)) / 1e9; tr = double(R.t - R.t(1)) / 1e9;
            z = [interp1(tr, R.obs(:, 1), tq), interp1(tr, R.obs(:, 2), tq)];
            y = zeros(240, 3);
            for k = 1:240
                Sk = cal.scale + S1 * (z(k, 1) - 22);
                y(k, :) = (Sk * a(k, :).' + b(:, k) + g.coef * (z(k, :) - g.mean).').';
            end
            got = cal.apply(tc.arc3(t, y), 'Regressors', R).obs;
            err = max(abs(got - a), [], 'all') / max(abs(a), [], 'all');
            tc.verifyBoundAndRecord("Calibration apply with a temperature-dependent scale vs the observation equation", err, 1e-9, "-", ...
                "one 3x3 solve per epoch on S(T(t)); ~1e-15 expected, bound as for the regressor terms");
        end
        function calibrationRegressorPrecision(tc)
            % apply against a per-epoch loop of y = S a + b(t) + C (z(t) - m), z interpolated linearly.
            [cal, R, C, m] = tc.regCal();
            t = R.t(1) + int64(0:239)' * int64(7.3e9);
            a = 1e-7 * [cos((1:240)' / 9), sin((1:240)' / 7), cos((1:240)' / 5)];
            b = cal.biasAt(t);
            tq = double(t - R.t(1)) / 1e9; tr = double(R.t - R.t(1)) / 1e9;
            z = [interp1(tr, R.obs(:, 1), tq), interp1(tr, R.obs(:, 2), tq)] - m;
            y = (cal.scale * a.' + b + C * z.').';
            got = cal.apply(tc.arc3(t, y), 'Regressors', R).obs;
            err = max(abs(got - a), [], 'all') / max(abs(a), [], 'all');
            tc.verifyBoundAndRecord("Calibration apply with regressor terms vs the observation equation", err, 1e-9, "-", ...
                "a 3x3 solve after subtracting the regressor term; the 1e-7 signal beside a 2e-7 bias leaves ~1e-15, bound generous");
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function syntheticPairVectorization(tc)
            % The arc length is the batch: a shorter arc is the start of a longer one wherever
            % both have a complete window, to the rounding of the chain (the noise is 1e-30).
            args = {'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0};
            [o1, ref1, truth] = NonGravForces.app.syntheticPair('Days', tc.Days, args{:});
            o2 = NonGravForces.app.syntheticPair('Days', 2 * tc.Days, args{:});
            n = numel(o1.t);
            tc.verifyEqual(o2.t(1:n), o1.t);
            ok = all(isfinite(o1.obs), 2);
            tc.verifyEqual(ok, all(isfinite(o2.obs(1:n, :)), 2) & ok);
            err = max(abs(o2.obs(ok, :) - o1.obs(ok, :)), [], 'all');
            tc.verifyBoundAndRecord("syntheticPair half arc vs whole arc", err, tc.chainBound(truth, ref1, 4), ...
                "m/s^2", "the same sums, the BLAS product associated by batch");
        end
        function calibrationVectorization(tc)
            % An arc applied whole equals the arc applied in two halves (not claimed bitwise:
            % the solve and the bias product are batched, CLAUDE.md pitfall).
            cal = tc.splineCal(8);
            t = cal.bias.t0 + int64(0:599)' * int64(30e9);
            y = 1e-7 * [cos((1:600)' / 30), sin((1:600)' / 45), cos((1:600)' / 60)];
            whole = cal.apply(tc.arc3(t, y)).obs;
            h1 = cal.apply(tc.arc3(t(1:300), y(1:300, :))).obs;
            h2 = cal.apply(tc.arc3(t(301:end), y(301:end, :))).obs;
            tc.verifyBoundAndRecord("Calibration.apply whole vs halves", max(abs(whole - [h1; h2]), [], 'all') / max(abs(whole), [], 'all'), ...
                16 * cond(cal.scale) * eps, "rel", "16 cond(S) eps");
            b = cal.biasAt(t);
            b1 = zeros(3, numel(t));
            for k = 1:numel(t), b1(:, k) = cal.biasAt(t(k)); end
            tc.verifyBoundAndRecord("Calibration.biasAt batch vs per epoch", max(abs(b - b1), [], 'all') / max(abs(b), [], 'all'), 8 * eps, "rel", "8 eps, four basis terms");
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function syntheticPairRobustness(tc)
            % NaN exactly where the 140.7 s window leaves the arc or meets the gap, nowhere else.
            [obs, ~, truth] = NonGravForces.app.syntheticPair('Days', tc.Days, 'Gap_s', 30);
            halfNs = int64(truth.kernel.Nh) * int64(1e8);
            t1 = obs.t(1); tn = obs.t(1) + int64(round(tc.Days * 86400 * 10) - 1) * int64(1e8);
            lo = obs.t - halfNs; hi = obs.t + halfNs;
            expected = lo < t1 | hi > tn | (hi >= truth.gap(1) & lo <= truth.gap(2));
            tc.verifyEqual(isnan(obs.obs), repmat(expected, 1, 3));
            % the gap's 29.9 s widened by 2 x 70.3 s holds 170 or 171 whole seconds, by alignment
            tc.verifyTrue(ismember(nnz(expected & lo >= t1 & hi <= tn), [170 171]));
            % Alternative input forms: a 3x1 scale is its diagonal; no spikes and no gap are empty.
            S = [0.96; 0.98; 0.95];
            a = NonGravForces.app.syntheticPair('Days', tc.Days, 'Scale', S, 'Seed', 4);
            b = NonGravForces.app.syntheticPair('Days', tc.Days, 'Scale', diag(S), 'Seed', 4);
            tc.verifyEqual(a.obs, b.obs);
            [~, ~, tr0] = NonGravForces.app.syntheticPair('Days', tc.Days, 'Spikes', 0, 'Gap_s', 0);
            tc.verifyEqual(height(tr0.spikes), 0); tc.verifyEmpty(tr0.gap);
        end
        function calibrationRobustness(tc)
            cal = tc.splineCal(8);
            % A NaN reading stays a NaN row and touches no other row.
            t = cal.bias.t0 + int64(0:9)' * int64(60e9);
            y = 1e-7 * ones(10, 3); y(4, 2) = NaN;
            a = cal.apply(tc.arc3(t, y)).obs;
            tc.verifyTrue(all(isnan(a(4, :))));
            tc.verifyTrue(all(isfinite(a([1:3 5:10], :)), 'all'));
            % The right end of the spline support is inside it (the half-open last interval).
            last = cal.bias.knots(end - 3);
            tc.verifyEqual(sum(NonGravForces.app.Calibration.basis(cal.bias, last)), 1, 'AbsTol', 16 * eps);
            % Jumps are unit steps after the smooth basis: 0 before, 1 from the jump epoch on.
            bj = struct('basis', "constant", 't0', int64(0), 'knots', int64([]), 'jumps', int64([10 20]) * int64(1e9), 'coef', zeros(3, 3));
            tc.verifyEqual(NonGravForces.app.Calibration.basis(bj, int64([0 10 15 20 25]') * int64(1e9)), ...
                [1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1]);
            bj.jumps = int64([20 10]) * int64(1e9);
            tc.verifyErrorId(@() NonGravForces.app.Calibration.basis(bj, int64(0)), "NonGravForces:Calibration:inputBasis");
            % A hat basis extrapolates linearly beyond its knots.
            h = struct('basis', "hat", 't0', int64(0), 'knots', int64([0 86400e9]), 'coef', zeros(3, 2));
            tc.verifyEqual(NonGravForces.app.Calibration.basis(h, int64(2 * 86400e9)), [-1 2], 'AbsTol', 4 * eps);
            % A file from a newer schema is read, and says so.
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'cal.mat'); cal.save(f);
            s = load(f); s.calibration.schemaVersion = string(str2double(NonGravForces.app.Calibration.SchemaVersion) + 1);   % newer than this toolbox, whatever it is
            f2 = fullfile(folder, 'cal2.mat'); save(f2, '-struct', 's');
            back = tc.verifyWarning(@() NonGravForces.app.Calibration.load(f2), "NonGravForces:Calibration:newerSchema");
            tc.verifyEqual(back.scale, cal.scale);
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function syntheticPairMetamorphic(tc)
            % The Seed is the only randomness: the same Seed is bitwise the same arcs, another is not.
            [o1, r1, tr1] = NonGravForces.app.syntheticPair('Days', tc.Days, 'Seed', 5);
            [o2, r2, tr2] = NonGravForces.app.syntheticPair('Days', tc.Days, 'Seed', 5);
            tc.verifyTrue(isequaln(o1.obs, o2.obs)); tc.verifyEqual(r1.obs, r2.obs);
            tc.verifyEqual(tr1.spikes, tr2.spikes);
            o3 = NonGravForces.app.syntheticPair('Days', tc.Days, 'Seed', 6);
            tc.verifyFalse(isequaln(o1.obs, o3.obs));
            % A bias offset d moves the observation by dcGain * d and by nothing else.
            d = [3e-8; -2e-8; 1e-8];
            o4 = NonGravForces.app.syntheticPair('Days', tc.Days, 'Seed', 5, 'Bias', tr1.bias0 + d);
            shift = o4.obs - o1.obs;
            ok = all(isfinite(shift), 2);
            err = max(abs(shift(ok, :) - tr1.kernel.dcGain * d.'), [], 'all');
            tc.verifyBoundAndRecord("syntheticPair bias offset -> dcGain * offset", err, tc.chainBound(tr1, r1, 4), "m/s^2", "two runs");
            % 'ScaleDrift' zeros IS the default pair, bitwise (N-03a adds the term only when planted).
            o5 = NonGravForces.app.syntheticPair('Days', tc.Days, 'Seed', 5, 'ScaleDrift', zeros(3));
            tc.verifyTrue(isequaln(o5.obs, o1.obs));
        end
        function calibrationMetamorphic(tc)
            cal = tc.splineCal(8);
            t = cal.bias.t0 + int64(0:99)' * int64(60e9);
            y = 1e-7 * [cos((1:100)' / 9), sin((1:100)' / 7), cos((1:100)' / 5)];
            % Linearity, bitwise: doubling the reading and the bias doubles the acceleration.
            f2 = struct('scale', cal.scale, 'bias', setfield(cal.bias, 'coef', 2 * cal.bias.coef)); %#ok<SFLD>
            a1 = cal.apply(tc.arc3(t, y)).obs;
            a2 = NonGravForces.app.Calibration(f2).apply(tc.arc3(t, 2 * y)).obs;
            tc.verifyEqual(a2, 2 * a1);
            % load(save(cal)) is cal.
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'cal.mat');
            ci = cal.invert();
            ci.save(f);
            tc.verifyTrue(isequaln(NonGravForces.app.Calibration.load(f), ci));
        end

        function calibrationTimeVaryingMetamorphic(tc)
            % A curve with zero slope applies as the constant scale does, to a few eps (pages vs
            % one matrix solve); a curve survives save/load; a schema-1 file loads as a constant.
            flat = tc.curveCal(zeros(3));
            base = tc.splineCal(8);
            t = base.bias.t0 + int64(0:99)' * int64(60e9);
            y = 1e-7 * [cos((1:100)' / 9), sin((1:100)' / 7), cos((1:100)' / 5)];
            a0 = base.apply(tc.arc3(t, y)).obs; a1 = flat.apply(tc.arc3(t, y)).obs;
            tc.verifyBoundAndRecord("Calibration flat curve vs constant scale, apply", ...
                max(abs(a1 - a0), [], 'all') / max(abs(a0), [], 'all'), 100 * eps, "-", "pages vs one 3x3 solve");
            cal = tc.curveCal(0.01 * eye(3));
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'curve.mat');
            cal.save(f);
            tc.verifyTrue(isequaln(NonGravForces.app.Calibration.load(f), cal));
            old = struct('schemaVersion', "1", 'scale', base.scale, 'covScale', base.covScale, 'bias', base.bias, ...
                'quadratic', base.quadratic, 'uncertainty', struct(), 'gates', struct(), 'convention', "observation", ...
                'provenance', base.provenance);
            calibration = old; f1 = fullfile(folder, 'schema1.mat'); save(f1, 'calibration');
            c1 = NonGravForces.app.Calibration.load(f1);
            tc.verifyFalse(c1.isTimeVarying());
            tc.verifyEqual(c1.scale, base.scale);
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function syntheticPairSpeed(tc)
            % FFT synthesis and CRN outputs both grow with the arc: twice the arc, about twice the cost.
            NonGravForces.app.syntheticPair('Days', tc.Days); NonGravForces.app.syntheticPair('Days', 2 * tc.Days);
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.syntheticPair('Days', tc.Days); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.syntheticPair('Days', 2 * tc.Days); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("syntheticPair, 2x arc / 1x arc", median(t2) / median(t1), 4, "-", ...
                "n log n would be about 2.1; quadratic 4");
        end
        function calibrationSpeed(tc)
            cal = tc.splineCal(8);
            mk = @(n) tc.arc3(cal.bias.t0 + int64(0:n - 1)' * int64(1e9), 1e-7 * ones(n, 3));
            A1 = mk(5000); A4 = mk(20000);
            for k = 1:5, cal.apply(A1); cal.apply(A4); end
            t1 = zeros(1, 3); t4 = zeros(1, 3);
            for k = 1:3
                t0 = tic; cal.apply(A1); t1(k) = toc(t0);
                t0 = tic; cal.apply(A4); t4(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("Calibration.apply, 20000 / 5000 epochs", median(t4) / median(t1), 8, "-", ...
                "linear would be 4; the basis is evaluated column by column");
        end
    end
end
