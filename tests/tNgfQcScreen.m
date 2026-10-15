classdef tNgfQcScreen < NgfTestCase
%TNGFQCSCREEN Tests of the LP-5 quality-control units: residuals, detectImpulses, detectJumps, screen.
%
%   Oracles (none built here): the exhaustive optimal segmentation in
%   tools/pyxcheck/pelt_check.py (dynamic programming over every admissible
%   change point, no pruning) for detectJumps, planted pulses and steps at known
%   epochs with known amplitudes for the detectors, the real 2020-01-01 GRACE-C
%   THR1B day as impulse priors, and closed-form properties of the residual
%   definition. Every bound carries the number measured on 2026-09-07.
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
    properties
        root
        series      % pelt_series.csv
        oracle      % pelt_oracle.csv
        crn         % 30-sample raised-cosine pulse template
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            tc.series = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'pelt_series.csv'));
            tc.oracle = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'pelt_oracle.csv'), 'TextType', 'string');
            tc.crn = sin(pi * (0:29)' / 30).^2;
        end
    end
    methods
        function t = epochs(~, K)
            t = datetime(2020, 1, 1, 'TimeZone', 'UTC') + seconds((0:K - 1)');
        end
        function [x, t, truth] = synthetic(tc, seed)
            % 3000 white samples of 1e-9 with two CRN pulses (+8, -6), a level step (+5) and one spike (30)
            rng(seed); K = 3000; x = 1e-9 * randn(K, 1); t = tc.epochs(K);
            x(1000:1029) = x(1000:1029) + 8e-9 * tc.crn; x(2000:2029) = x(2000:2029) - 6e-9 * tc.crn;
            x(2500:end) = x(2500:end) + 5e-9; x(700) = 3e-8;
            truth = struct('pulses', [1000 2000], 'amplitudes', [8e-9 -6e-9], 'step', 2500, 'stepSize', 5e-9, 'spike', 700);
        end
        function q = speedRatio(~, fcn, base)
            fcn(); base();
            r = zeros(1, 5);
            for k = 1:5
                t1 = tic; fcn(); a = toc(t1); t2 = tic; base(); b = toc(t2); r(k) = a / b;
            end
            q = median(r);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function residualsContract(tc)
            y = randn(100, 3); m = randn(100, 3);
            [r, info] = NonGravForces.qc.residuals(y, m); tc.verifyEqual(r, y - m); tc.verifySize(info.sigmaRobust, [1 3]); tc.verifyEqual(info.whitening, "none");
            r2 = NonGravForces.qc.residuals(y, m, 'Bias', [1 2 3], 'Scale', [1.05 1 0.9]); tc.verifyEqual(r2, y - ([1 2 3] + m .* [1.05 1 0.9]), 'AbsTol', 1e-15);
            S = [1 0.1 0; 0 1 0; 0.2 0 1]; r3 = NonGravForces.qc.residuals(y, m, 'Scale', S); tc.verifyEqual(r3, y - m * S', 'AbsTol', 1e-15);
            rd = NonGravForces.qc.residuals(y, m, 'Whitening', "difference"); tc.verifyTrue(all(isnan(rd(1, :)))); tc.verifyEqual(rd(2:end, :), diff(y - m));
            rh = NonGravForces.qc.residuals(y, m, 'Whitening', "highpass", 'Window', 11); tc.verifySize(rh, [100 3]);
            f = @(varargin) NonGravForces.qc.residuals(varargin{:});
            tc.verifyErrorId(@() f(y, m(1:50, :)), 'timeSeriesPreparation:residuals:inputShape');
            tc.verifyErrorId(@() f(y, m, 'Bias', [1 2]), 'timeSeriesPreparation:residuals:inputShape');
            tc.verifyErrorId(@() f(y, m, 'Scale', eye(2)), 'timeSeriesPreparation:residuals:inputShape');
            tc.verifyErrorId(@() f(y, m, 'Whitening', "highpass"), 'timeSeriesPreparation:residuals:inputWindow');
            tc.verifyErrorId(@() f(y, m, 'Whitening', "highpass", 'Window', 10), 'timeSeriesPreparation:residuals:inputWindow');
            tc.verifyError(@() f(y, m, 'Whitening', "wavelet"), 'MATLAB:validators:mustBeMember');
        end
        function detectImpulsesContract(tc)
            [x, t] = tc.synthetic(1);
            [ev, info] = NonGravForces.qc.detectImpulses(x, 'Epoch', t, 'Response', tc.crn);
            tc.verifyEqual(ev.Properties.VariableNames, NonGravForces.qc.eventTable().Properties.VariableNames);
            tc.verifyTrue(all(ev.class == "impulse")); tc.verifyTrue(all(ev.axis == "all")); tc.verifySize(info.score, [3000 1]); tc.verifyEqual(info.templateLength, 30);
            ev3 = NonGravForces.qc.detectImpulses([x x x], 'Epoch', t, 'Response', tc.crn); tc.verifyTrue(all(ismember(ev3.axis, ["x" "y" "z"])));
            d = @(varargin) NonGravForces.qc.detectImpulses(varargin{:});
            tc.verifyErrorId(@() d(x, 'Response', tc.crn), 'timeSeriesPreparation:detectImpulses:inputEpoch');
            tNo = t; tNo.TimeZone = ''; tc.verifyErrorId(@() d(x, 'Epoch', tNo, 'Response', tc.crn), 'timeSeriesPreparation:detectImpulses:inputEpoch');
            tc.verifyErrorId(@() d(x, 'Epoch', t(1:10), 'Response', tc.crn), 'timeSeriesPreparation:detectImpulses:inputShape');
            tc.verifyErrorId(@() d([x x], 'Epoch', t, 'Response', tc.crn), 'timeSeriesPreparation:detectImpulses:inputShape');
            tc.verifyErrorId(@() d(x, 'Epoch', t), 'timeSeriesPreparation:detectImpulses:inputResponse');
            tc.verifyErrorId(@() d(x, 'Epoch', t, 'Response', ones(5, 1)), 'timeSeriesPreparation:detectImpulses:inputResponse');
            tc.verifyErrorId(@() d(x, 'Epoch', t, 'Response', tc.crn, 'Priors', table(1)), 'timeSeriesPreparation:detectImpulses:inputPriors');
            tc.verifyErrorId(@() d(x, 'Epoch', t, 'Response', tc.crn, 'Threshold', 0), 'timeSeriesPreparation:detectImpulses:inputRange');
        end
        function detectJumpsContract(tc)
            [x, t] = tc.synthetic(1);
            [ev, info] = NonGravForces.qc.detectJumps(x, 'Epoch', t);
            tc.verifyEqual(ev.Properties.VariableNames, NonGravForces.qc.eventTable().Properties.VariableNames);
            tc.verifyTrue(all(ev.class == "jump")); tc.verifyEqual(ev.tStart, ev.tEnd); tc.verifySize(info.changePoints, [1 1]); tc.verifyEqual(info.model, "level");
            evS = NonGravForces.qc.detectJumps(x, 'Epoch', t, 'Model', "slope", 'MinSegment', 20); tc.verifyTrue(istable(evS));
            j = @(varargin) NonGravForces.qc.detectJumps(varargin{:});
            tc.verifyErrorId(@() j(x), 'timeSeriesPreparation:detectJumps:inputEpoch');
            tc.verifyErrorId(@() j(x, 'Epoch', t(1:5)), 'timeSeriesPreparation:detectJumps:inputShape');
            tc.verifyErrorId(@() j([x x], 'Epoch', t), 'timeSeriesPreparation:detectJumps:inputShape');
            tc.verifyErrorId(@() j(x, 'Epoch', t, 'Penalty', -1), 'timeSeriesPreparation:detectJumps:inputRange');
            tc.verifyErrorId(@() j(x, 'Epoch', t, 'MinSegment', 1), 'timeSeriesPreparation:detectJumps:inputRange');
            tc.verifyErrorId(@() j(x, 'Epoch', t, 'Model', "slope", 'MinSegment', 2), 'timeSeriesPreparation:detectJumps:inputRange');
            tc.verifyErrorId(@() j(x, 'Epoch', t, 'Priors', table(1)), 'timeSeriesPreparation:detectJumps:inputPriors');
            tc.verifyErrorId(@() j(x(1:3), 'Epoch', t(1:3)), 'timeSeriesPreparation:detectJumps:dataEmpty');
            tc.verifyError(@() j(x, 'Epoch', t, 'Model', "quadratic"), 'MATLAB:validators:mustBeMember');
        end
        function screenContract(tc)
            [x, t] = tc.synthetic(1);
            [mask, ev, info] = NonGravForces.qc.screen(x, zeros(3000, 1), t, 'Response', tc.crn);
            tc.verifyClass(mask, 'logical'); tc.verifySize(mask, [3000 1]); tc.verifyTrue(istable(ev)); tc.verifyEqual(fieldnames(info.counts), {'outlier'; 'impulse'; 'jump'});
            tc.verifyTrue(issorted(ev.tStart));
            [~, ~, info2] = NonGravForces.qc.screen(x, zeros(3000, 1), t); tc.verifyEqual(info2.impulses, "skipped: no Response given"); tc.verifyEqual(info2.counts.impulse, 0);
            s = @(varargin) NonGravForces.qc.screen(varargin{:});
            tc.verifyErrorId(@() s(x, zeros(10, 1), t), 'timeSeriesPreparation:screen:inputShape');
            tc.verifyErrorId(@() s([x x], zeros(3000, 2), t), 'timeSeriesPreparation:screen:inputShape');
            tc.verifyErrorId(@() s(x, zeros(3000, 1), (1:3000)'), 'timeSeriesPreparation:screen:inputShape');
            tNo = t; tNo.TimeZone = ''; tc.verifyErrorId(@() s(x, zeros(3000, 1), tNo), 'timeSeriesPreparation:screen:inputEpoch');
            tc.verifyErrorId(@() s(x, zeros(3000, 1), t, 'Whitening', "highpass"), 'timeSeriesPreparation:residuals:inputWindow');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function detectJumpsReference(tc)
            % PELT equals the exhaustive optimum of pelt_check.py (level and slope, three penalties), bitwise in the indices
            x = tc.series.x; t = tc.epochs(numel(x));
            for i = 1:height(tc.oracle)
                [~, info] = NonGravForces.qc.detectJumps(x, 'Epoch', t, 'Model', tc.oracle.model(i), 'Penalty', tc.oracle.beta(i), 'MinSegment', tc.oracle.minSeg(i));
                want = str2double(split(tc.oracle.changePoints(i), ";"))';
                tc.verifyEqual(info.changePoints{1}, want, sprintf('%s, beta %g', tc.oracle.model(i), tc.oracle.beta(i)));
                tc.verifyBoundAndRecord(sprintf("PELT objective vs exhaustive optimum (%s, beta %g)", tc.oracle.model(i), tc.oracle.beta(i)), abs(info.cost - tc.oracle.cost(i)) / tc.oracle.cost(i), 1e-12, "-", "two summation orders of ~160 squares");
            end
        end
        function detectImpulsesReference(tc)
            % the real 2020-01-01 GRACE-C THR1B day as priors: every firing inside a synthetic residual span is assigned to THR1B
            thr = NonGravForces.io.readTHR1B(fullfile(tc.root, 'data', 'samples', 'THR1B_2020-01-01_C_04.txt'));
            tc.assertGreaterThan(height(thr), 0);
            t0 = thr.tStart(1) - seconds(300); K = 1800; t = t0 + seconds((0:K - 1)');
            inSpan = thr(thr.tStart >= t(1) & thr.tStart <= t(end - 30), :);
            tc.assertGreaterThan(height(inSpan), 0, 'the first thruster firing lies inside the synthetic span by construction');
            rng(3); x = 1e-9 * randn(K, 1);
            for e = 1:height(inSpan)
                q = find(t >= inSpan.tStart(e), 1); x(q:q + 29) = x(q:q + 29) + 1e-8 * tc.crn;
            end
            [ev, info] = NonGravForces.qc.detectImpulses(x, 'Epoch', t, 'Response', tc.crn, 'Priors', thr);
            tc.verifyEqual(nnz(ev.source == "THR1B"), height(inSpan), 'each planted firing is detected and assigned to its THR1B prior');
            tc.verifyEqual(height(info.unmatchedPriors), height(thr) - height(inSpan), 'priors outside the span stay unmatched');
            tc.verifyBoundAndRecord("THR1B-assigned pulse epochs vs the priors", max(abs(seconds(ev.tStart(ev.source == "THR1B") - inSpan.tStart))), 1, "s", "planted at the first sample after the prior");
        end
        function residualsReference(tc)
            % the robust scale is consistent with the Gaussian sigma (Hampel 1974 constant 1.4826): 1e5 samples
            rng(4); x = 2.5 * randn(1e5, 1);
            [~, info] = NonGravForces.qc.residuals(x, zeros(1e5, 1));
            tc.verifyBoundAndRecord("1.4826 MAD vs sigma of 1e5 Gaussian samples", abs(info.sigmaRobust - 2.5) / 2.5, 0.02, "-", "sampling error of the MAD ~ 1.2/sqrt(N) = 0.4 %; measured 0.3 %");
        end
        function screenReference(tc)
            % the real ACT1B sample (600 s, GRACE-C 2020-01-01) against a zero model: the along-track drag signal is a slow drift, no impulse, at most a few outliers
            [ep, a] = NonGravForces.io.readACC1B(fullfile(tc.root, 'data', 'samples', 'ACT1B_2020-01-01_C_04.txt'));
            K = size(a, 2); t = datetime(2020, 1, 1, 'TimeZone', 'UTC') + seconds(ep.sec(:) - ep.sec(1));
            [mask, ev, info] = NonGravForces.qc.screen(a', zeros(K, 3), t, 'Whitening', "highpass", 'Window', 61, 'Response', tc.crn);
            tc.verifySize(mask, [K 3]);
            tc.verifyBoundAndRecord("outliers flagged in the real 600-s ACT1B sample (fraction)", info.counts.outlier / (3 * K), 0.02, "-", "Hampel at 4 sigma on 1/f data: a few per 1000 expected; measured 2026-09-07");
            tc.verifyBoundAndRecord("impulses found in the real 600-s sample (recorded, no oracle)", info.counts.impulse, 20, "-", "every 30-s window could fire; a 5-sigma threshold on 1/f data with a synthetic template found 11 on 2026-09-07");
            tc.verifyTrue(all(ismember(ev.class, ["outlier" "impulse" "jump"])));
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function residualsPrecision(tc)
            % a step of height h under "difference" is one sample of height h; under "highpass" (W = 11) its centre keeps h/2 exactly for a noiseless step
            x = [zeros(50, 1); ones(50, 1)];
            rd = NonGravForces.qc.residuals(x, zeros(100, 1), 'Whitening', "difference"); tc.verifyEqual(rd(51), 1); tc.verifyEqual(rd([2:50 52:100]), zeros(98, 1));
            rh = NonGravForces.qc.residuals(x, zeros(100, 1), 'Whitening', "highpass", 'Window', 11);
            tc.verifyBoundAndRecord("highpass step response at the first new sample", abs(rh(51) - (1 - 6 / 11)), 1e-15, "-", "movmean of 11 centred: 6 of 11 samples are on the new level");
            tc.verifyEqual(rh([1:40 61:100]), zeros(80, 1), 'AbsTol', 1e-15);
        end
        function detectImpulsesPrecision(tc)
            % planted pulses: epoch to +-1 sample, amplitude within 3 sigma_a; the spike and the step do not produce plateaus of pulses
            [x, t, truth] = tc.synthetic(1);
            [ev, info] = NonGravForces.qc.detectImpulses(x, 'Epoch', t, 'Response', tc.crn);
            for k = 1:2
                q = find(abs(seconds(ev.tStart - t(truth.pulses(k)))) <= 1, 1);
                tc.assertNotEmpty(q, sprintf('pulse %d found', k));
                tc.verifyBoundAndRecord(sprintf("pulse %d amplitude error in sigma_a", k), abs(ev.amplitude(q) - truth.amplitudes(k)) / ev.sigma(q), 3, "sigma", "Gaussian noise, matched filter is the ML estimator; measured 2026-09-07");
            end
            tc.verifyBoundAndRecord("pulses reported beyond the two planted ones", height(ev) - 2, 2, "count", "the step edge and the spike may each register once with the zero-mean template; measured 2026-09-07");
            tc.verifyBoundAndRecord("sigma_a vs the analytic value sigma / (hz.h) for white noise", abs(ev.sigma(1) - 1e-9 / info.templateGain) / (1e-9 / info.templateGain), 0.1, "-", "robust scale of 3000 samples: 3 % sampling error");
        end
        function detectJumpsPrecision(tc)
            % planted 4-sigma step in white noise: epoch to +-1 sample, amplitude within 3 sigma, and no other change point
            rng(5); K = 1000; s = 1e-9; x = s * randn(K, 1); x(601:end) = x(601:end) + 4 * s; t = tc.epochs(K);
            [ev, info] = NonGravForces.qc.detectJumps(x, 'Epoch', t);
            tc.assertEqual(height(ev), 1, 'exactly one change point');
            tc.verifyBoundAndRecord("planted step epoch error", abs(info.changePoints{1} - 601), 1, "samples", "4 sigma step, BIC penalty; measured 2026-09-07");
            tc.verifyBoundAndRecord("planted step amplitude error in sigma", abs(ev.amplitude - 4 * s) / ev.sigma, 3, "sigma", "difference of two means; measured 2026-09-07");
            tc.verifyBoundAndRecord("robust sigma from first differences vs truth", abs(info.sigma - s) / s, 0.1, "-", "1000 samples, MAD sampling error ~4 %");
            % slope model: a planted kink
            y = s * randn(K, 1); y(501:end) = y(501:end) + 0.1 * s * (1:500)';
            [evS, infoS] = NonGravForces.qc.detectJumps(y, 'Epoch', t, 'Model', "slope", 'MinSegment', 10);
            tc.assertGreaterThanOrEqual(height(evS), 1);
            tc.verifyBoundAndRecord("planted kink epoch error (slope model)", min(abs(infoS.changePoints{1} - 501)), 10, "samples", "a kink of 0.1 sigma/sample: the fitted lines separate by 1 sigma after 10 samples; measured 8 (a 0.02 sigma/sample kink was off by 28)");
        end
        function screenPrecision(tc)
            % the synthetic case end to end: spike masked, both pulses masked over their 30 samples, the step reported as a jump, nothing else
            [x, t, truth] = tc.synthetic(1);
            [mask, ev, info] = NonGravForces.qc.screen(x, zeros(3000, 1), t, 'Response', tc.crn);
            tc.verifyTrue(mask(truth.spike), 'the spike is masked');
            for k = 1:2, tc.verifyTrue(all(mask(truth.pulses(k) + (1:28))), sprintf('pulse %d body is masked', k)); end
            j = ev(ev.class == "jump", :);
            tc.assertGreaterThanOrEqual(height(j), 1, 'the step is reported');
            [~, q] = min(abs(seconds(j.tStart - t(truth.step))));
            tc.verifyBoundAndRecord("step epoch error through the pipeline", abs(seconds(j.tStart(q) - t(truth.step))), 30, "s", "the step edge registers as one impulse candidate whose 30-sample window is blanked; the jump lands at its end (measured 6 s)");
            tc.verifyBoundAndRecord("step amplitude error through the pipeline", abs(j.amplitude(q) - truth.stepSize) / j.sigma(q), 3, "sigma", "measured 2026-09-07");
            tc.verifyBoundAndRecord("masked samples beyond spike and pulses", nnz(mask) - 1 - 60, 40, "samples", "Hampel false alarms at 4 sigma on 3000 Gaussian samples (~0.2) plus the step edge candidate (30); measured 2026-09-07");
            tc.verifyEqual(info.counts.jump, height(j));
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function residualsVectorization(tc)
            y = randn(200, 3); m = randn(200, 3);
            r = NonGravForces.qc.residuals(y, m, 'Bias', [1 2 3], 'Scale', [1.1 0.9 1], 'Whitening', "highpass", 'Window', 21);
            for k = 1:3, tc.verifyEqual(r(:, k), NonGravForces.qc.residuals(y(:, k), m(:, k), 'Bias', k, 'Scale', [1.1 0.9 1] * (1:3 == k)', 'Whitening', "highpass", 'Window', 21)); end
        end
        function detectImpulsesVectorization(tc)
            [x, t] = tc.synthetic(1); y = circshift(x, 137); z = -x;
            [ev3, i3] = NonGravForces.qc.detectImpulses([x y z], 'Epoch', t, 'Response', tc.crn);
            for k = 1:3
                col = [x y z]; axn = ["x" "y" "z"]; [ev1, i1] = NonGravForces.qc.detectImpulses(col(:, k), 'Epoch', t, 'Response', tc.crn);
                sel = ev3(ev3.axis == axn(k), :); sel.axis(:) = "all";
                tc.verifyEqual(sel, ev1); tc.verifyEqual(i3.score(:, k), i1.score);
            end
        end
        function detectJumpsVectorization(tc)
            rng(6); x = randn(400, 3); x(201:end, 1) = x(201:end, 1) + 3; x(101:end, 2) = x(101:end, 2) - 4; t = tc.epochs(400);
            [ev3, i3] = NonGravForces.qc.detectJumps(x, 'Epoch', t, 'Penalty', 20);
            for k = 1:3
                axn = ["x" "y" "z"]; [ev1, i1] = NonGravForces.qc.detectJumps(x(:, k), 'Epoch', t, 'Penalty', 20);
                sel = ev3(ev3.axis == axn(k), :); sel.axis(:) = "all";
                tc.verifyEqual(sel, ev1); tc.verifyEqual(i3.changePoints{k}, i1.changePoints{1}); tc.verifyEqual(i3.cost(k), i1.cost);
            end
        end
        function screenVectorization(tc)
            [x, t] = tc.synthetic(1); y = circshift(x, 500);
            [m3, ~, i3] = NonGravForces.qc.screen([x y x], zeros(3000, 3), t, 'Response', tc.crn, 'JumpPenalty', 1e-17);
            [m1, ~, i1] = NonGravForces.qc.screen(x, zeros(3000, 1), t, 'Response', tc.crn, 'JumpPenalty', 1e-17);
            [m2, ~, i2] = NonGravForces.qc.screen(y, zeros(3000, 1), t, 'Response', tc.crn, 'JumpPenalty', 1e-17);
            tc.verifyEqual(m3, [m1 m2 m1]); tc.verifyEqual(i3.residual, [i1.residual i2.residual i1.residual]);
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function residualsRobustness(tc)
            y = randn(50, 2); y(10, 1) = NaN; m = zeros(50, 2);
            [r, info] = NonGravForces.qc.residuals(y, m, 'Whitening', "highpass", 'Window', 5);
            tc.verifyTrue(isnan(r(10, 1))); tc.verifyEqual(nnz(isnan(r)), 1, 'NaN neither spreads nor vanishes'); tc.verifyTrue(all(isfinite(info.sigmaRobust)));
            [r0, info0] = NonGravForces.qc.residuals(zeros(0, 2), zeros(0, 2)); tc.verifySize(r0, [0 2]); tc.verifyTrue(all(isnan(info0.sigmaRobust)));
            rc = NonGravForces.qc.residuals(ones(50, 1), ones(50, 1)); tc.verifyEqual(rc, zeros(50, 1));
        end
        function detectImpulsesRobustness(tc)
            % NaN gaps, a constant series (no scale: no detections, no error), a template as long as the series
            [x, t] = tc.synthetic(1); x(1500:1600) = NaN;
            [ev, info] = NonGravForces.qc.detectImpulses(x, 'Epoch', t, 'Response', tc.crn);
            tc.verifyEqual(info.nanCount, 101); tc.verifyGreaterThanOrEqual(height(ev), 2);
            tc.verifyTrue(~any(ev.tStart >= t(1500) & ev.tStart <= t(1600)), 'no pulse is reported inside the gap');
            evC = NonGravForces.qc.detectImpulses(zeros(200, 1), 'Epoch', t(1:200), 'Response', tc.crn); tc.verifyEqual(height(evC), 0);
            evL = NonGravForces.qc.detectImpulses(randn(30, 1), 'Epoch', t(1:30), 'Response', tc.crn); tc.verifyTrue(istable(evL));
        end
        function detectJumpsRobustness(tc)
            % NaN samples are skipped and indices mapped back; a constant series yields no change point; a tiny series works
            rng(7); x = randn(300, 1); x(151:end) = x(151:end) + 6; x(140:160) = NaN; t = tc.epochs(300);
            [ev, info] = NonGravForces.qc.detectJumps(x, 'Epoch', t, 'Penalty', 30);
            tc.assertEqual(height(ev), 1); tc.verifyEqual(info.changePoints{1}, 161, 'the first finite sample after the gap starts the new segment');
            evC = NonGravForces.qc.detectJumps(ones(50, 1), 'Epoch', t(1:50), 'Penalty', 1); tc.verifyEqual(height(evC), 0);
            [evT, infoT] = NonGravForces.qc.detectJumps([0; 0; 5; 5], 'Epoch', t(1:4), 'Penalty', 1); tc.verifyEqual(infoT.changePoints{1}, 3); tc.verifyEqual(evT.amplitude, 5);
        end
        function screenRobustness(tc)
            % NaN in the measurement propagates into the mask as "not usable" only through the detectors' NaN policy; an empty event list is a valid result
            rng(8); x = 1e-9 * randn(500, 1); x(100:110) = NaN; t = tc.epochs(500);
            [mask, ev, info] = NonGravForces.qc.screen(x, zeros(500, 1), t, 'Response', tc.crn);
            tc.verifySize(mask, [500 1]); tc.verifyTrue(istable(ev)); tc.verifyEqual(info.counts.jump + info.counts.impulse + info.counts.outlier, height(ev));
            tc.verifyTrue(~any(mask(100:110)), 'NaN samples are not flagged (they are already unusable)');
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function residualsMetamorphic(tc)
            % swapping the roles of y and model with unit scale negates the residual; adding the same series to both leaves it unchanged
            y = randn(100, 2); m = randn(100, 2); z = randn(100, 2);
            tc.verifyEqual(NonGravForces.qc.residuals(y, m), -NonGravForces.qc.residuals(m, y));
            tc.verifyEqual(NonGravForces.qc.residuals(y + z, m + z), NonGravForces.qc.residuals(y, m), 'AbsTol', 1e-15);
            rh1 = NonGravForces.qc.residuals(y, m, 'Whitening', "highpass", 'Window', 9); rh2 = NonGravForces.qc.residuals(y + 7, m, 'Whitening', "highpass", 'Window', 9);
            tc.verifyBoundAndRecord("highpass removes a constant offset", max(abs(rh1 - rh2), [], 'all'), 1e-13, "-", "movmean of 9 samples of order 7: rounding 1e-15 x 7");
        end
        function detectImpulsesMetamorphic(tc)
            % time shift of the epochs shifts the events; scaling the series scales the amplitudes and leaves the scores; negation flips the sign
            [x, t] = tc.synthetic(1);
            [e1, i1] = NonGravForces.qc.detectImpulses(x, 'Epoch', t, 'Response', tc.crn);
            [e2, ~] = NonGravForces.qc.detectImpulses(x, 'Epoch', t + hours(5), 'Response', tc.crn);
            tc.verifyEqual(e2.tStart, e1.tStart + hours(5)); tc.verifyEqual(e2.amplitude, e1.amplitude);
            [e3, i3] = NonGravForces.qc.detectImpulses(2 * x, 'Epoch', t, 'Response', tc.crn);
            tc.verifyEqual(e3.amplitude, 2 * e1.amplitude); tc.verifyEqual(i3.score, i1.score, 'AbsTol', 1e-12);
            [e4, ~] = NonGravForces.qc.detectImpulses(-x, 'Epoch', t, 'Response', tc.crn); tc.verifyEqual(e4.amplitude, -e1.amplitude);
            % a level offset added to the series changes nothing (zero-mean template)
            [e5, i5] = NonGravForces.qc.detectImpulses(x + 3e-9, 'Epoch', t, 'Response', tc.crn);
            tc.verifyEqual(e5.tStart, e1.tStart); tc.verifyBoundAndRecord("matched-filter output under a level offset", max(abs(i5.score - i1.score)), 1e-6, "sigma", "zero-mean template: rounding of 3e-9 x sum(hz) ~ 1e-25 over sigma_c");
        end
        function detectJumpsMetamorphic(tc)
            % time shift shifts the events, affine scaling of the series scales the amplitudes, reversal mirrors the change points
            rng(9); x = randn(300, 1); x(101:200) = x(101:200) + 5; t = tc.epochs(300);
            [e1, i1] = NonGravForces.qc.detectJumps(x, 'Epoch', t, 'Penalty', 30);
            tc.assertEqual(numel(i1.changePoints{1}), 2);
            [e2, ~] = NonGravForces.qc.detectJumps(x, 'Epoch', t + days(1), 'Penalty', 30); tc.verifyEqual(e2.tStart, e1.tStart + days(1));
            [e3, i3] = NonGravForces.qc.detectJumps(3 * x + 2, 'Epoch', t, 'Penalty', 9 * 30);
            tc.verifyEqual(i3.changePoints{1}, i1.changePoints{1}); tc.verifyEqual(e3.amplitude, 3 * e1.amplitude, 'RelTol', 1e-12);
            [~, i4] = NonGravForces.qc.detectJumps(flip(x), 'Epoch', t, 'Penalty', 30);
            tc.verifyEqual(sort(302 - i4.changePoints{1}), i1.changePoints{1}, 'reversal mirrors the change points (first sample of the new segment maps to last of the old)');
        end
        function screenMetamorphic(tc)
            % time shift, affine scaling of measurement and model together, and reversal of the mask
            [x, t] = tc.synthetic(1); m = 1e-10 * randn(3000, 1);
            [k1, e1] = NonGravForces.qc.screen(x, m, t, 'Response', tc.crn);
            [k2, e2] = NonGravForces.qc.screen(x, m, t + hours(1), 'Response', tc.crn);
            tc.verifyEqual(k2, k1); tc.verifyEqual(e2.tStart, e1.tStart + hours(1)); tc.verifyEqual(e2.class, e1.class);
            [k3, e3] = NonGravForces.qc.screen(4 * x, 4 * m, t, 'Response', tc.crn);
            tc.verifyEqual(k3, k1); tc.verifyEqual(e3.amplitude, 4 * e1.amplitude, 'RelTol', 1e-12);
            [k4, e4] = NonGravForces.qc.screen(x - m, zeros(3000, 1), t, 'Response', tc.crn);
            tc.verifyEqual(k4, k1); tc.verifyEqual(e4.amplitude, e1.amplitude);
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function residualsSpeed(tc)
            y = randn(86400, 3); m = randn(86400, 3);
            q = tc.speedRatio(@() NonGravForces.qc.residuals(y, m, 'Whitening', "highpass", 'Window', 601), @() movmean(y, 601, 1, 'omitnan'));
            tc.verifyBoundAndRecord("residuals highpass on a day of 3 channels vs one movmean", q, 36, "-", "budget 3x the expected 12: one movmean plus two medians of 86400 per channel (sort-based); measured 11.3 to 14.1 across runs on the owner's host");
        end
        function detectImpulsesSpeed(tc)
            rng(10); x = randn(86400, 1); t = tc.epochs(86400);
            q = tc.speedRatio(@() NonGravForces.qc.detectImpulses(x, 'Epoch', t, 'Response', tc.crn), @() conv(x, tc.crn, 'same'));
            tc.verifyBoundAndRecord("detectImpulses on a day at 1 Hz vs one conv", q, 60, "-", "budget 3x the expected 20 (30 shifted adds, MAD, sort); 5 repeats, median");
        end
        function detectJumpsSpeed(tc)
            rng(11); x = randn(20000, 1); x(10001:end) = x(10001:end) + 3; t = tc.epochs(20000);
            q = tc.speedRatio(@() NonGravForces.qc.detectJumps(x, 'Epoch', t), @() NonGravForces.qc.detectJumps(x(1:5000), 'Epoch', t(1:5000)));
            tc.verifyBoundAndRecord("detectJumps PELT: 20000 vs 5000 samples, time ratio", q, 24, "-", "PELT is near-linear: expected ~4-8 (pruning keeps few candidates); budget 3x; 5 repeats, median");
        end
        function screenSpeed(tc)
            rng(12); x = 1e-9 * randn(20000, 1); t = tc.epochs(20000);
            q = tc.speedRatio(@() NonGravForces.qc.screen(x, zeros(20000, 1), t, 'Response', tc.crn), @() NonGravForces.qc.detectJumps(x, 'Epoch', t, 'MinSegment', 30));
            tc.verifyBoundAndRecord("screen vs detectJumps alone on 20000 samples", q, 3, "-", "budget 3x the expected ~1.05: the PELT loop dominates the pipeline (measured 2.4 s of 2.6 s at 20000 samples before the inlining, O-18), hampel 0.05 s, the rest 5 ms");
        end
    end
end
