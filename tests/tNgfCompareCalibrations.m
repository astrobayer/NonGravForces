classdef tNgfCompareCalibrations < NgfTestCase
%TNGFCOMPARECALIBRATIONS Tests of the calibration model comparison (N-03b): app.compareCalibrations.
%
%   The truth is app.syntheticPair, whose default bias drifts linearly on x
%   and z and not at all on y. Blocked CV must find that order; AIC and BIC,
%   which count correlated rows as independent, must be shown NOT to (they
%   are reported as IC_unreliable, and the precision test records why).
%
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
    properties (Constant)
        Specs = {{'BiasBasis', "constant"}, {'BiasBasis', "linear"}}
        Names = ["const" "linear"]
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function compareCalibrationsContract(tc)
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.1, 'Seed', 2);
            [T, info] = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Names', tc.Names);
            tc.verifyClass(T, 'table');
            tc.verifyEqual(T.Properties.VariableNames, {'name', 'cvMse', 'cvMseSe', 'rank', 'nParams', 'nRows', 'aic', 'bic'});
            tc.verifyEqual(T.name, tc.Names');
            tc.verifySize(T.cvMse, [2 3]); tc.verifySize(T.rank, [2 3]);
            tc.verifyEqual(sort(T.rank, 1), repmat((1:2)', 1, 3), 'each axis ranks every candidate once');
            tc.verifyEqual(T.nParams(:, 1), [2; 3], 'x: bias (1 or 2) plus S_xx');
            tc.verifyTrue(all(contains(string(T.Properties.VariableDescriptions(7:8)), "IC_unreliable")));
            u = T.Properties.UserData;
            tc.verifyEqual(sort(string(fieldnames(u)))', sort(["best" "folds" "gap" "note"]));
            tc.verifySize(u.best, [1 3]);
            tc.verifyGreaterThanOrEqual(u.gap, ceil(timeSeriesPreparation.crnKernel("ACC").Tf), 'the buffer is at least the filter span');
            tc.verifyEqual(sort(string(fieldnames(info)))', sort(["foldLoss" "testRows" "gap" "calibrations"]));
            tc.verifySize(info.foldLoss{1}, [5 3]);
            tc.verifyClass(info.calibrations{2}, 'NonGravForces.app.Calibration');
            % default names
            T2 = NonGravForces.app.compareCalibrations(obs, ref, {{}}, 'Folds', 2);
            tc.verifyEqual(T2.name, "spec1");
            % refusals
            id = "NonGravForces:compareCalibrations:inputSpecs";
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {}), id);
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {'BiasBasis'}), id);
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {{'BiasBasis'}}), id);
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {{'Kernel', []}}), id);
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {{'NBoot', 10}}), id);
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Names', "one"), id);
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {{'BiasBasis', "pspline", 'BiasKnotSpacing_s', 600}}), ...
                "NonGravForces:compareCalibrations:inputFolds");
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, {{'Scale', "half"}}), ...
                "NonGravForces:fitCalibration:inputScale");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function compareCalibrationsReference(tc)
            % The composition, re-derived by hand for one interior fold: fit app.fitCalibration on
            % the rows outside the test block and its buffers, predict the test block, square and
            % average. Bitwise the fold loss compareCalibrations reports.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.1, 'Seed', 3);
            [~, info] = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs);
            k = 3;
            rows = (1:numel(obs.t))';
            tr = info.testRows(k, :);
            test = rows >= tr(1) & rows <= tr(2);
            train = rows < tr(1) - info.gap | rows > tr(2) + info.gap;
            x = obs.obs; x(~train, :) = NaN;
            live = find(any(isfinite(obs.obs), 2));
            span = [obs.t(live(1)), obs.t(live(end))];
            for c = 1:2
                [~, fi] = NonGravForces.app.fitCalibration(obs.derive(obs.t, x, "test.train"), ref, tc.Specs{c}{:}, 'NBoot', 0, 'Span', span);
                r = obs.obs(test, :) - fi.prediction.obs(test, :);
                mse = zeros(1, 3);
                for i = 1:3, e = r(isfinite(r(:, i)), i); mse(i) = mean(e.^2); end
                tc.verifyEqual(info.foldLoss{c}(k, :), mse, "fold " + k + ", candidate " + c);
            end
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function compareCalibrationsPrecision(tc)
            % The planted bias order (linear on x and z, constant on y) against a flexible spline:
            % CV must pick the truth on every axis of every seeded realisation, and AIC is recorded
            % beside it -- it picked the spline on 9 of 9 axes when this test was written, which is
            % the v3 plan's reason for reporting it as IC_unreliable.
            specs = [tc.Specs, {{'BiasBasis', "pspline", 'BiasKnotSpacing_s', 5400}}];
            names = [tc.Names, "spline"];
            truth = ["linear" "const" "linear"];
            seeds = [11 12 13];
            cvRight = 0; aicSpline = 0;
            for s = seeds
                [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', s);
                T = NonGravForces.app.compareCalibrations(obs, ref, specs, 'Names', names);
                cvRight = cvRight + nnz(T.Properties.UserData.best == truth);
                [~, ka] = min(T.aic);
                aicSpline = aicSpline + nnz(T.name(ka)' == "spline");
            end
            n = 3 * numel(seeds);
            tc.verifyBoundAndRecord("compareCalibrations CV misses of the planted bias order, fraction", 1 - cvRight / n, 0, "-", ...
                "3 seeds x 3 axes, Days 0.25, 5 folds; CV must find the truth every time");
            tc.verifyBoundAndRecord("compareCalibrations AIC picks the overfit spline, fraction", aicSpline / n, 1, "-", ...
                "recorded, not bounded: the v3 plan's claim is that it exceeds 0.5 (IC_unreliable)");
            tc.verifyGreaterThan(aicSpline / n, 0.5, 'AIC overfits on correlated residuals (v3 WP-14)');
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function compareCalibrationsVectorization(tc)
            % Candidates compared together equal candidates compared one at a time, bitwise: the folds
            % and the buffer depend on the observation only, and each candidate's fit on a fold is
            % independent of the others.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.1, 'Seed', 4);
            Tb = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs);
            for c = 1:2
                Ts = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs(c));
                tc.verifyEqual(Ts.cvMse, Tb.cvMse(c, :));
                tc.verifyEqual(Ts.cvMseSe, Tb.cvMseSe(c, :));
                tc.verifyEqual(Ts.aic, Tb.aic(c, :));
            end
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function compareCalibrationsRobustness(tc)
            % The default pair carries spikes and a mid-arc gap: every fold still scores finitely.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.1, 'Seed', 5);
            [T, info] = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Gap', 0);
            tc.verifyTrue(all(isfinite(T.cvMse), 'all'));
            tc.verifyTrue(all(cellfun(@(F) all(isfinite(F), 'all'), info.foldLoss)));
            tc.verifyEqual(T.Properties.UserData.gap, 0, 'an explicit Gap is taken as given');
            % Folds and buffer that leave a fold without training rows are refused by crossval,
            % unwrapped; the same folds with a small buffer are accepted (the refusal is provoked).
            tc.verifyErrorId(@() NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Folds', 2, 'Gap', 10000), ...
                "timeSeriesPreparation:crossval:noTraining");
            T2 = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Folds', 2, 'Gap', 200);
            tc.verifyTrue(all(isfinite(T2.cvMse), 'all'));
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function compareCalibrationsMetamorphic(tc)
            % Reordering the candidates permutes the rows, bitwise, and leaves the winners alone;
            % a 'Fit' pair a spec repeats is overridden by the spec.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.1, 'Seed', 6);
            T1 = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Names', tc.Names);
            T2 = NonGravForces.app.compareCalibrations(obs, ref, fliplr(tc.Specs), 'Names', fliplr(tc.Names));
            tc.verifyEqual(T2.cvMse, flipud(T1.cvMse));
            tc.verifyEqual(T2.aic, flipud(T1.aic));
            tc.verifyEqual(T2.Properties.UserData.best, T1.Properties.UserData.best);
            T3 = NonGravForces.app.compareCalibrations(obs, ref, tc.Specs, 'Names', tc.Names, 'Fit', {'BiasBasis', "pspline"});
            tc.verifyEqual(T3.cvMse, T1.cvMse, 'every spec sets BiasBasis, so the common one never applies');
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function compareCalibrationsSpeed(tc)
            % Each candidate costs Folds + 1 fits: two candidates, about twice one.
            [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.05);
            one = tc.Specs(1); two = tc.Specs;
            NonGravForces.app.compareCalibrations(obs, ref, one); NonGravForces.app.compareCalibrations(obs, ref, two);
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.compareCalibrations(obs, ref, one); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.compareCalibrations(obs, ref, two); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("compareCalibrations, 2 candidates / 1", median(t2) / median(t1), 4, "-", ...
                "linear in candidates would be 2; the shared buffer computation pulls it lower");
        end
    end
end
