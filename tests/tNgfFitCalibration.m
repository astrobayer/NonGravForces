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
        function [obs, ref, b0] = multiPair(tc, N, noise, varargin)
            %MULTIPAIR N consecutive days of syntheticPair: one scale, a different constant bias per day.
            obs = timeSeriesPreparation.Arc.empty; ref = obs; b0 = zeros(3, N);
            asd = [1e-10; 1e-9; 1e-10]; if ~noise, asd = tc.Quiet; end
            for k = 1:N
                st = int64(7305 * 86400 - 43200 + 18 + (k - 1) * 86400) * int64(1e9);
                b0(:, k) = [2e-7; -1e-7; 5e-8] + 1e-8 * k * [1; -2; 3];
                [obs(k), ref(k)] = NonGravForces.app.syntheticPair('Days', 0.1, 'Start', st, 'Seed', k, 'Bias', b0(:, k), ...
                    'BiasDrift', zeros(3, 1), 'NoiseAsd', asd, varargin{:});
            end
        end
        function [x, names] = estimateOf(~, cal, info, i)
            %ESTIMATEOF The fitted parameters of axis i in the order of info.names{i}.
            names = info.names{i};
            x = zeros(numel(names), 1);
            ax = ["x" "y" "z"];
            for k = 1:numel(names)
                nm = names(k);
                if startsWith(nm, "b_"), x(k) = cal.bias.coef(i, str2double(extractAfter(nm, "b_")));
                elseif startsWith(nm, "S_") && strlength(nm) > 4
                    % a scale-curve coefficient S_ij_q (N-03a): row (j-1)*3+i of scaleCurve.coef
                    j = find(ax == extractBetween(nm, 4, 4));
                    x(k) = cal.scaleCurve.coef((j - 1) * 3 + i, str2double(extractAfter(nm, 5)));
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
            % N-03b: the prediction at every row, and the span the bases cover
            tc.verifyEqual(info.prediction.t, obs.t);
            tc.verifySize(info.prediction.obs, [numel(obs.t) 3]);
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'Span', [1 2], 'NBoot', 0), "NonGravForces:fitCalibration:inputSpan");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'Span', obs.t([end 1]), 'NBoot', 0), "NonGravForces:fitCalibration:inputSpan");
        end

        function fitCalibrationTimeVaryingContract(tc)
            % N-03a: 'ScaleBasis' makes every estimated entry a curve the Calibration carries.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 2);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear", 'NBoot', 0);
            tc.verifyTrue(cal.isTimeVarying());
            tc.verifyEqual(cal.scaleCurve.basis, "linear");
            tc.verifySize(cal.scaleCurve.coef, [9 2]);
            tc.verifyEqual(info.names{1}, ["b_1" "b_2" "S_xx_1" "S_xx_2"]);
            tc.verifyEqual(find(~cellfun(@isempty, cal.scaleCurve.cov)), [1 9], 'only the estimated S_xx and S_zz carry a curve covariance');
            tc.verifyEqual(cal.scaleCurve.coef(5, :), [1 0], 'a fixed entry is its constant a priori on the linear basis');
            tc.verifyEqual(cal.scale, cal.scaleAt(cal.scaleCurve.tRef), 'AbsTol', 1e-15, 'scale is the curve at tRef');
            tc.verifyTrue(cal.scaleCurve.tRef > obs.t(1) && cal.scaleCurve.tRef < obs.t(end));
            tc.verifySize(cal.scaleSigmaAt(obs.t(1:5)), [3 3 5]);
            % the formal white covariance is refused as a basis for the drift test -- by a warning
            tc.verifyWarning(@() cal.drift(), "NonGravForces:Calibration:formalDrift");
            % a knotted scale finer than a knotted bias is refused; equal spacing is accepted
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "hat", 'BiasKnotSpacing_s', 3600, ...
                'ScaleBasis', "hat", 'ScaleKnotSpacing_s', 1800, 'NBoot', 0), "NonGravForces:fitCalibration:inputScaleBasis");
            ok = NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "hat", 'BiasKnotSpacing_s', 3600, ...
                'ScaleBasis', "hat", 'ScaleKnotSpacing_s', 3600, 'NBoot', 0);
            tc.verifyTrue(ok.isTimeVarying());
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

    methods (Test, TestTags = {'reference'})
        function fitCalibrationTimeVaryingReference(tc)
            % Noise-free with a planted linear scale drift: S(t) a is exactly in the span of the
            % columns H([1, t] .* a), so the curve returns to the same computed rounding bound.
            D = diag([0.01 0 -0.02]);
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.1, 'ScaleDrift', D, ...
                'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear", 'NBoot', 0);
            for i = [1 3]
                tc.verifyRecovered(cal, info, truth, ref, obs, i, ...
                    [truth.bias0(i), truth.biasDrift_per_day(i), truth.scale(i, i), D(i, i)], "fitCalibration S(t) linear");
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

        function fitCalibrationTimeVaryingPrecision(tc)
            % Coloured noise and a 1 %/day along-track drift (v3 WP-14's case): the recovered slope
            % within 4 bootstrap sigma, the drift flagged, and the pointwise 2-sigma band covering
            % the true S_xx(t) at more than 90 % of the fitted epochs.
            D = zeros(3); D(1, 1) = 0.01;
            [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.5, 'Seed', 4, 'ScaleDrift', D);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear", 'NBoot', 40, 'Seed', 1);
            g = cal.scaleCurve.coef(1, 2); sg = sqrt(cal.scaleCurve.cov{1}(2, 2));
            tc.verifyBoundAndRecord("fitCalibration scale drift error in bootstrap sigma", abs(g - D(1, 1)) / sg, 4, "sigma", ...
                "seeded noise, 40 replicates; the slope of a linear S_xx(t)");
            tc.verifyEqual(cal.isDriftSignificant(), logical(diag([1 0 0])), ...
                'the planted 1 %/day drift on S_xx is detected and the undrifted S_zz is not');
            t = obs.t(info.rows(:, 1));
            days = double(t - obs.t(1)) / 86400e9;
            S = squeeze(cal.scaleAt(t)); s = squeeze(cal.scaleSigmaAt(t));
            cover = mean(abs(squeeze(S(1, 1, :)) - (truth.scale(1, 1) + D(1, 1) * days)) <= 2 * squeeze(s(1, 1, :)));
            tc.verifyBoundAndRecord("fitCalibration S(t) outside its 2-sigma band, fraction of epochs", 1 - cover, 0.1, "-", ...
                "the v3 WP-14 criterion: the band covers the truth at more than 90 % of epochs");
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
            % A constant bias with S_yy held: axis y has ONE parameter, so lsqvce returns a scalar,
            % and a scalar indexed by an empty index is 1x0. N-03a broke this case and nothing caught
            % it until N-03b's cross-validation fitted a constant bias (regression test).
            cc = NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "constant", 'NBoot', 0);
            tc.verifyEqual(cc.scale(2, 2), 1);
            tc.verifySize(cc.bias.coef, [3 1]);
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

    methods (Test, TestTags = {'robustness'})
        function fitCalibrationTimeVaryingRobustness(tc)
            % Knotted scale bases over a gap and spikes: every curve finite where it is evaluated,
            % and a cubic B-spline curve on an arc shorter than its knot spacing still solves.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.1, 'Seed', 5);
            for b = ["hat" "pspline"]
                cal = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', b, 'ScaleKnotSpacing_s', 3600, 'NBoot', 0);
                tc.verifyTrue(all(isfinite(cal.scaleAt(obs.t)), 'all'), b + " curve finite over the arc");
                tc.verifyTrue(all(isfinite(cal.scaleCurve.coef), 'all'));
            end
            c = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "pspline", 'NBoot', 0);   % 43200 s > 0.1 day
            tc.verifyTrue(all(isfinite(c.scale), 'all'));
            % a calibration with a curve applies to the arc it came from without a singular page
            a = c.apply(obs);
            tc.verifyEqual(isnan(a.obs), isnan(obs.obs), 'apply adds no NaN and removes none');
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
            % The prediction and the residual are two readings of one fit: on a fitted row,
            % prediction + residual = observation, to the rounding of one subtraction each.
            [~, info] = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0);
            rows = info.rows;
            y = obs.obs(rows); back = info.prediction.obs(rows) + info.residual.obs(rows);
            tc.verifyBoundAndRecord("fitCalibration prediction + residual vs observation", max(abs(back - y)) / (4 * eps * max(abs(y))), 1, "-", ...
                "two subtractions of the same fitted value; bound 4 eps max|y|");
            % 'Span' equal to the usable span is the default, bitwise.
            live = find(any(isfinite(obs.obs), 2));
            cs = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0, 'Span', [obs.t(live(1)), obs.t(live(end))]);
            tc.verifyEqual(cs.bias.coef, c1.bias.coef);
        end

        function fitCalibrationTimeVaryingMetamorphic(tc)
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05, 'Seed', 6);
            % 'ScaleBasis' "constant" IS the default path: bitwise, scale and bias alike.
            c0 = NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 0);
            c1 = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "constant", 'NBoot', 0);
            tc.verifyEqual(c1.scale, c0.scale);
            tc.verifyEqual(c1.bias.coef, c0.bias.coef);
            tc.verifyEqual(c1.covScale, c0.covScale);
            tc.verifyFalse(c1.isTimeVarying());
            % A hat basis with one span (two knots) spans the same functions as "linear": the two
            % fits are one least-squares problem in two parameterisations, so the curves agree to
            % the solve's rounding, cond * eps of the curve.
            cl = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear", 'NBoot', 0);
            ch = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "hat", 'ScaleKnotSpacing_s', 86400, 'NBoot', 0);
            tc.verifySize(ch.scaleCurve.coef, [9 2]);
            Sl = cl.scaleAt(obs.t); Sh = ch.scaleAt(obs.t);
            bound = 100 * eps * max([cl.gates.cond([1 3]), ch.gates.cond([1 3])]) * max(abs(Sl(:)));
            tc.verifyBoundAndRecord("fitCalibration S(t), hat with one span vs linear", max(abs(Sh(:) - Sl(:))) / bound, 1, "-", ...
                "two parameterisations of one column space; bound 100 eps cond max|S|");
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

    % ===================================================================== N-09: many arcs, one scale
    methods (Test, TestTags = {'contract'})
        function fitCalibrationMultiArcContract(tc)
            [obs, ref] = tc.multiPair(3, true, 'Spikes', 0);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "constant", 'NBoot', 0);
            tc.verifyClass(cal, 'NonGravForces.app.Calibration'); tc.verifySize(cal, [1 3]);
            tc.verifyEqual(info.mode, "common");
            tc.verifyEqual(cal(2).scale, cal(1).scale, 'the scale is common: every element carries the same one');
            tc.verifyNotEqual(cal(2).bias.coef, cal(1).bias.coef, 'the biases are per day');
            tc.verifyEqual(cal(3).bias.t0, obs(3).t(1));
            tc.verifyEqual(cal(1).uncertainty.mode, "common"); tc.verifyEqual(cal(1).uncertainty.nArcs, 3);
            tc.verifySize(info.residual, [1 3]); tc.verifyEqual(info.residual{2}.t, obs(2).t);
            tc.verifyWarning(@() NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "constant", 'NBoot', 5), ...
                "NonGravForces:fitCalibration:fewArcs");
            id = "NonGravForces:fitCalibration:unsupportedMultiArc";
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "pspline", 'NBoot', 0), id);
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear", 'NBoot', 0), id);
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref, 'Jumps', obs(2).t(100), 'NBoot', 0), id);
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs, ref(1:2), 'NBoot', 0), "NonGravForces:fitCalibration:inputArc");
            tc.verifyErrorId(@() NonGravForces.app.fitCalibration(obs([2 1 3]), ref([2 1 3]), 'NBoot', 0), "NonGravForces:fitCalibration:inputArc");
        end
    end
    methods (Test, TestTags = {'reference'})
        function fitCalibrationMultiArcReference(tc)
            % Noise-free, three days, ONE planted scale and three different biases: the joint fit returns
            % them all. The scale bound is the single-arc one over all rows (the common cofactor info.Qxx
            % times the row rounding); a bias is a constant per day, recovered to the same rounding
            % relative to the signal it sits beside.
            [obs, ref, b0] = tc.multiPair(3, false, 'Spikes', 0, 'Gap_s', 0, 'Scale', tc.TruthS .* eye(3));
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "constant", 'FixCrossTrack', false, 'NBoot', 0);
            [~, ~, truth] = NonGravForces.app.syntheticPair('Days', 0.1, 'Scale', tc.TruthS .* eye(3), 'NoiseAsd', tc.Quiet, 'Spikes', 0, 'Gap_s', 0);
            n = sum(cellfun(@(r) nnz(r(:, 1)), info.rows));
            for i = 1:3
                q = sqrt(info.Qxx{i}(1, 1));
                bound = 2 * q * sqrt(n) * tc.rowBound(truth, ref(1)) + 100 * eps * cal(1).gates.cond(i) * q * sqrt(n) * 2e-7;
                tc.verifyBoundAndRecord("fitCalibration multi-arc common scale, axis " + i, abs(cal(1).scale(i, i) - tc.TruthS(i, i)) / bound, 1, "-", ...
                    "error / (2 sqrt(Q n) eRow + 100 eps cond sqrt(Q n) |signal|), n over all days");
            end
            errB = max(abs(cell2mat(arrayfun(@(c) c.bias.coef(:, 1), cal, 'UniformOutput', false)) - b0), [], 'all');
            tc.verifyBoundAndRecord("fitCalibration multi-arc per-day biases, relative to the signal", errB / 2e-7, 1e-9, "-", ...
                "rounding of the CRN-filtered solve; measured ~1e-15 when written");
        end
    end
    methods (Test, TestTags = {'precision'})
        function fitCalibrationMultiArcPrecision(tc)
            % The algebra of pre-elimination: lsqvce's joint solve with each day's bias as a local block
            % equals the DENSE joint least squares -- a block-diagonal bias design beside the common scale
            % column, solved here by QR -- at the rounding of that solve, cond * eps.
            [obs, ref] = tc.multiPair(3, true, 'Spikes', 0, 'Gap_s', 0);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "linear", 'NBoot', 0);
            K = timeSeriesPreparation.crnKernel("ACC");
            X = []; l = []; nB = 2;
            for k = 1:3
                Ha = timeSeriesPreparation.crnApply(ref(k), K, OutputEpochs = obs(k).t).obs;
                bk = struct('basis', "linear", 't0', obs(k).t(1), 'knots', int64([]), 'coef', zeros(3, 2));
                Br = NonGravForces.app.Calibration.basis(bk, ref(k).t);
                names = {'b1', 'b2'}; u = struct('b1', '-', 'b2', '-');
                bArc = timeSeriesPreparation.Arc(ref(k).t, Br, struct("obsNames", {names}, "variableNames", {names}, ...
                    "variableFieldNames", {names}, "variableUnits", u, "fileName", "DESIGN", "product", "DESIGN"));
                HB = timeSeriesPreparation.crnApply(bArc, K, OutputEpochs = obs(k).t).obs;
                r = info.rows{k}(:, 1);
                blk = zeros(nnz(r), 3 * nB); blk(:, (k - 1) * nB + (1:nB)) = HB(r, :);
                X = [X; blk, Ha(r, 1)]; l = [l; obs(k).obs(r, 1)]; %#ok<AGROW>
            end
            [Q, Rq] = qr(X, 0); xd = Rq \ (Q' * l);
            bound = 100 * eps * cond(X);
            tc.verifyBoundAndRecord("fitCalibration multi-arc scale vs dense joint LS (relative)", abs(cal(1).scale(1, 1) - xd(end)) / abs(xd(end)) / bound, 1, "-", ...
                "two solvers of one problem; bound 100 eps cond(X)");
            yb = cell2mat(arrayfun(@(c) c.bias.coef(1, :)', cal, 'UniformOutput', false));
            tc.verifyBoundAndRecord("fitCalibration multi-arc biases vs dense joint LS (relative to max bias)", ...
                max(abs(yb(:) - xd(1:end - 1))) / max(abs(xd(1:end - 1))) / bound, 1, "-", "same solve, the local block");
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function fitCalibrationVectorization(tc)
            % 'ScaleMode' "perArc" is N single-arc fits, bitwise, arc k with seed Seed + k - 1.
            [obs, ref] = tc.multiPair(3, true);
            [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'ScaleMode', "perArc", 'NBoot', 3, 'Seed', 4);
            tc.verifyEqual(info.mode, "perArc");
            for k = 1:3
                c1 = NonGravForces.app.fitCalibration(obs(k), ref(k), 'NBoot', 3, 'Seed', 4 + k - 1);
                tc.verifyEqual(cal(k).scale, c1.scale, "arc " + k);
                tc.verifyEqual(cal(k).bias.coef, c1.bias.coef, "arc " + k);
                tc.verifyEqual(cal(k).covScale, c1.covScale, "arc " + k);
            end
        end
    end
    methods (Test, TestTags = {'robustness'})
        function fitCalibrationMultiArcRobustness(tc)
            % Spikes and a gap inside every day (several sub-arcs per day), the held S_yy axis that has no
            % common column (it decouples into per-day fits), and the bootstrap over days: all finite.
            [obs, ref] = tc.multiPair(3, true);
            [cal, info] = tc.verifyWarning(@() NonGravForces.app.fitCalibration(obs, ref, 'NBoot', 4, 'Seed', 1), ...
                "NonGravForces:fitCalibration:fewArcs");
            tc.verifyTrue(all(isfinite(cal(1).scale), 'all'));
            tc.verifyTrue(all(cellfun(@(r) any(diff(find(r(:, 1))) > 1), info.rows)), 'each day really has a gap');
            tc.verifyTrue(all(isfinite(cell2mat(arrayfun(@(c) c.bias.coef(:), cal, 'UniformOutput', false))), 'all'));
            tc.verifyTrue(all(isfinite(diag(cal(1).covScale([1 9], [1 9])))), 'the bootstrap over days gives a covariance');
            tc.verifyEqual(cal(1).uncertainty.bootstrap, "days");
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function fitCalibrationMultiArcMetamorphic(tc)
            % Doubling every day's observation doubles the common scale and every bias of the free axes,
            % bitwise (the equilibration takes a power of two out exactly); and the day's own reference is
            % what matters -- passing one reference arc per day equals passing them as an array.
            [obs, ref] = tc.multiPair(3, true, 'Spikes', 0);
            c1 = NonGravForces.app.fitCalibration(obs, ref, 'BiasBasis', "constant", 'NBoot', 0);
            o2 = obs; for k = 1:3, o2(k) = obs(k).derive(obs(k).t, 2 * obs(k).obs, "test.double"); end
            c2 = NonGravForces.app.fitCalibration(o2, ref, 'BiasBasis', "constant", 'NBoot', 0);
            tc.verifyEqual(c2(1).scale([1 9]), 2 * c1(1).scale([1 9]));
            for k = 1:3, tc.verifyEqual(c2(k).bias.coef([1 3], :), 2 * c1(k).bias.coef([1 3], :)); end
        end
    end
end
