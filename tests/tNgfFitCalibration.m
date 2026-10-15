classdef tNgfFitCalibration < NgfTestCase
%TNGFFITCALIBRATION Tests of the calibration estimator (N-02b): app.fitCalibration.
%
%   The truth is app.syntheticPair: every planted number is known, and the
%   observation went through the real CRN chain. On a noise-free pair the
%   observation equation holds EXACTLY in real arithmetic (the bias is linear
%   and the kernel symmetric), so the estimate must return the truth to
%   rounding. That bound is computed, not written: each row of the system
%   carries the rounding of one Nf-term dot product, eRow = 4 Nf eps sum|F|
%   max|raw|, and a parameter moves by at most sqrt(Q_kk) * sqrt(n) * eRow
%   under it (Cauchy-Schwarz on (X'X)^-1 X' e), plus the solve's own
%   cond * eps. The noisy tests assert against the estimator's own bootstrap
%   sigma, not against a guessed absolute.
%
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
    properties (Constant)
        Quiet = 1e-30 * ones(3, 1)   % a noise ASD far below the rounding bound
        TruthS = [0.96 2e-3 0; 0 0.98 0; -1e-3 0 0.95]
    end
    methods
        function bound = rowBound(~, truth, ref)
            %ROWBOUND 4 Nf eps sum|F| max|raw|: the rounding of one filtered row.
            maxRaw = max(abs(ref.obs), [], 'all') * norm(truth.scale, inf) + max(abs(truth.bias), [], 'all') ...
                + max([abs(truth.spikes.amplitude); 0]);
            bound = 4 * truth.kernel.Nf * eps * sum(abs(truth.kernel.F)) * maxRaw;
        end
        function verifyRecovered(tc, cal, info, truth, ref, obs, i, xTrue, label)
            %VERIFYRECOVERED Every parameter of axis i within its computed rounding bound.
            [x, names] = tc.estimateOf(cal, info, i);
            Q = info.Qxx{i}; n = nnz(info.rows(:, i));
            yRms = sqrt(mean(obs.obs(info.rows(:, i), i).^2));
            q = sqrt(diag(Q));
            bound = 2 * q * sqrt(n) * tc.rowBound(truth, ref) + 100 * eps * cal.gates.cond(i) * q * sqrt(n) * yRms;
            err = abs(x - xTrue(:));
            [worst, k] = max(err ./ bound);
            tc.verifyBoundAndRecord(label + " axis " + i + " (worst " + names(k) + ")", worst, 1, "-", ...
                "error / (2 sqrt(Q n) eRow + 100 eps cond sqrt(Q n) rms y)");
        end
        function [x, names] = estimateOf(~, cal, info, i)
            %ESTIMATEOF The fitted parameters of axis i in the order of info.names{i}.
            names = info.names{i};
            x = zeros(numel(names), 1);
            ax = ["x" "y" "z"];
            for k = 1:numel(names)
                nm = names(k);
                if startsWith(nm, "b_"), x(k) = cal.bias.coef(i, str2double(extractAfter(nm, "b_")));
                elseif startsWith(nm, "S_"), x(k) = cal.scale(i, find(ax == extractAfter(nm, 3)));
                else, x(k) = cal.quadratic(i);
                end
            end
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function fitCalibrationContract(tc)
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 2);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0);
            tc.verifyClass(cal, 'NonGravForces.app.Calibration');
            tc.verifyEqual(cal.convention, "observation");
            tc.verifyEqual(cal.bias.basis, "linear");
            tc.verifySize(cal.bias.coef, [3 2]);
            tc.verifyEqual(cal.scale(2, 2), 1, 'the cross-track scale is held at its a priori by default');
            tc.verifyEqual(cal.scale(~eye(3)), zeros(6, 1), 'a diagonal fit leaves the off-diagonals at their a priori');
            tc.verifyTrue(all(isnan(cal.quadratic)));
            tc.verifyEqual(info.names{1}, ["b_1" "b_2" "S_xx"]);
            tc.verifyEqual(info.names{2}, ["b_1" "b_2"]);
            tc.verifyClass(info.results{1}, 'lsqvce.Results');
            tc.verifyEqual(info.residual.t, obs.t);
            tc.verifySize(info.rows, [numel(obs.t) 3]);
            tc.verifyEqual(sort(string(fieldnames(cal.uncertainty)))', sort(["sigmaScaleFormal" "sigmaScaleBootstrap" ...
                "inflation" "residualRms" "eiv" "noise" "nBoot" "blockLength" "stopReason"]));
            tc.verifyEqual(sort(string(fieldnames(cal.gates)))', sort(["vif" "cond" "maxAbsCorr" "weightFraction" "vifMax" "verdict"]));
            tc.verifyTrue(all(cal.uncertainty.eiv.lambda > 0 & cal.uncertainty.eiv.lambda <= 1));
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(5, ref), "NonGravForces:fitCalibration:inputArc");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'Scale', "half"), "NonGravForces:fitCalibration:inputScale");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'Scale', eye(2)), "NonGravForces:fitCalibration:inputScale");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'ScaleApriori', nan(3)), "NonGravForces:fitCalibration:inputScale");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'Kernel', 5), "NonGravForces:fitCalibration:inputKernel");
            shifted = ref.derive(ref.t + int64(5e7), ref.obs, "test.shift");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, shifted, 'Kernel', [], 'NBoot', 0), "NonGravForces:fitCalibration:inputGrid");
            tc.verifyWarning(@() NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "hat", 'BiasKnotSpacing_s', 100, 'NBoot', 0), ...
                "NonGravForces:fitCalibration:knotsBelowFilterSpan");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function fitCalibrationReference(tc)
            % Noise-free, a full S with off-diagonals and a drifting bias: the truth to rounding.
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.1, 'Scale', tc.TruthS, ...
                'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'Scale', "full", 'FixCrossTrack', false, ...
                'VifMax', Inf, 'NBoot', 0);
            for i = 1:3
                tc.verifyRecovered(cal, info, truth, ref, obs, i, [truth.bias0(i), truth.biasDrift_per_day(i), truth.scale(i, :)], ...
                    "fitCalibration full S");
            end
            % ... and with a planted K2 term, estimated per axis.
            k2 = [1e3; 0; 2e3];
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.1, 'Quadratic', k2, ...
                'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'FixCrossTrack', false, 'Quadratic', true, 'NBoot', 0);
            for i = [1 3]
                % the K2 column is demeaned over the fitted rows, and the mean went into the bias:
                % compare the scale and K2, which the demeaning leaves alone
                [x, names] = tc.estimateOf(cal, info, i);
                sel = ~startsWith(names, "b_");
                xTrue = [truth.scale(i, i); k2(i)];
                q = sqrt(diag(info.Qxx{i})); n = nnz(info.rows(:, i));
                yRms = sqrt(mean(obs.obs(info.rows(:, i), i).^2));
                bound = 2 * q(sel) * sqrt(n) * tc.rowBound(truth, ref) + 100 * eps * cal.gates.cond(i) * q(sel) * sqrt(n) * yRms;
                tc.verifyBoundAndRecord("fitCalibration K2 axis " + i, max(abs(x(sel) - xTrue) ./ bound), 1, "-", "error / computed bound");
            end
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function fitCalibrationPrecision(tc)
            % Coloured noise, spikes and a gap: every diagonal scale within 5 bootstrap sigma of the
            % truth, and the bootstrap sigma no smaller than half the formal one (the filtered
            % residual is correlated, so it should inflate, never collapse).
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 4);
            cal = NonGravForces.app.fitCalibration(obs, ref, 'FixCrossTrack', false, 'NBoot', 100, 'Seed', 1);
            sB = diag(cal.uncertainty.sigmaScaleBootstrap);
            z = abs(diag(cal.scale) - diag(truth.scale)) ./ sB;
            tc.verifyBoundAndRecord("fitCalibration scale error in bootstrap sigma", max(z), 5, "sigma", "seeded noise, 100 replicates");
            tc.verifyGreaterThan(min(diag(cal.uncertainty.inflation)), 0.5);
            tc.verifyTrue(all(isfinite(cal.covScale(logical(eye(9))))));
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function fitCalibrationRobustness(tc)
            % Spikes and a gap: the fit runs over the windows the filter left, nothing else.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 5);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0);
            tc.verifyTrue(all(isfinite(cal.scale), 'all'));
            tc.verifyFalse(any(info.rows(any(isnan(obs.obs), 2), :), 'all'), 'no NaN row of the observation is fitted');
            tc.verifyTrue(all(isnan(info.residual.obs(~info.rows))));
            % The VIF gate: a row with several model columns on collinear axes is refused, and the
            % same arc with a diagonal scale is fitted (the refusal is provoked, and paired).
            m = ref.obs;
            collinear = ref.derive(ref.t, [m(:, 1), m(:, 2), 0.5 * m(:, 1) + 1e-3 * m(:, 3)], "test.collinear");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, collinear, 'Scale', "full", 'NBoot', 0), ...
                "NonGravForces:fitCalibration:numericCollinear");
            c2 = NonGravForces.app.fitCalibration(obs, collinear, 'Scale', "diag", 'NBoot', 0);
            tc.verifyTrue(all(isfinite(c2.scale), 'all'), 'one model column per row is never gated');
            % The other noise paths: lsqvce's estimated law, and a declared ARMA law.
            ce = NonGravForces.app.fitCalibration(obs, ref, 'Noise', "estimate", 'NBoot', 0);
            tc.verifyTrue(all(isfinite(ce.scale), 'all'));
            tc.verifyEqual(ce.uncertainty.noise, "estimate");
            ca = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0, ...
                'Noise', struct('type', 'arma', 'ar', [1 -0.5], 'ma', 1, 'sigma', 1));
            tc.verifyTrue(all(isfinite(ca.scale), 'all'));
            tc.verifyEqual(ca.uncertainty.noise, "declared arma");
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function fitCalibrationMetamorphic(tc)
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 6);
            c1 = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0);
            % Three axes, three problems: changing the z observation leaves x bitwise alone.
            obsZ = obs.derive(obs.t, obs.obs + [0 0 1e-8], "test.z");
            c2 = NonGravForces.app.fitCalibration(obsZ, ref, 'NBoot', 0);
            tc.verifyEqual(c2.scale(1, :), c1.scale(1, :));
            tc.verifyEqual(c2.bias.coef(1, :), c1.bias.coef(1, :));
            tc.verifyNotEqual(c2.bias.coef(3, :), c1.bias.coef(3, :));
            % Doubling the observation doubles every estimated parameter of the free axes, bitwise:
            % the equilibration takes the factor out exactly and the solve never sees it.
            obs2 = obs.derive(obs.t, 2 * obs.obs, "test.double");
            c3 = NonGravForces.app.fitCalibration(obs2, ref, 'NBoot', 0);
            tc.verifyEqual(c3.scale([1 3], [1 3]) .* eye(2), 2 * c1.scale([1 3], [1 3]) .* eye(2));
            tc.verifyEqual(c3.bias.coef([1 3], :), 2 * c1.bias.coef([1 3], :));
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function fitCalibrationSpeed(tc)
            % Filter, gate and three solves all grow with the arc: twice the arc, about twice the cost.
            [o1, r1] = NonGravForces.app.syntheticPair('Days', 0.05);
            [o2, r2] = NonGravForces.app.syntheticPair('Days', 0.1);
            for k = 1:2
                NonGravForces.app.fitCalibration(o1, r1, 'NBoot', 0); NonGravForces.app.fitCalibration(o2, r2, 'NBoot', 0);
            end
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.fitCalibration(o1, r1, 'NBoot', 0); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.fitCalibration(o2, r2, 'NBoot', 0); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("fitCalibration, 2x arc / 1x arc", median(t2) / median(t1), 4, "-", ...
                "linear would be 2 and a fixed overhead pulls it lower; quadratic 4");
        end
    end
end
