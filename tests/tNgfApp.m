classdef tNgfApp < NgfTestCase
%TNGFAPP Tests of the application layer (LP-7): app.forceBudget.
%
%   The oracle for a reduction is the thing it reduces. app.forceBudget computes
%   no physics: with one chunk it must be BITWISE Model.budget, and that is the
%   reference test. Everything else here is about the parts Model.budget does not
%   have — chunking, several models in one table, the shares — and about the
%   arithmetic those parts add, which is a sum of squares and nothing more.
%
%   The arc is the real 2020-01-01 GRACE-C fixture, and the fast preset is used
%   wherever the test does not need the standard one: Model.accel costs 88 ms per
%   epoch in the standard preset and 2 ms in the fast (HANDOVER O-29), so a
%   200-epoch standard budget would be 18 s of gate time for nothing the fast
%   preset does not show.
%
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
    properties
        root, eop, sw, ep, r, v, att, macro, Mfast, Mstd
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
            tc.Mfast = NonGravForces.Model(tc.macro, NonGravForces.Options("fast"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, ...
                'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));
            tc.Mstd = NonGravForces.Model(tc.macro, NonGravForces.Options("standard"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw);
        end
    end
    methods
        function [e, rr, vv] = arc(tc, n)
            idx = 1:n;
            e = referenceSystem.daySecondToEpoch([tc.ep.day(idx) tc.ep.sec(idx)], tc.ep.scale);
            rr = tc.r(:, idx); vv = tc.v(:, idx);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function forceBudgetContract(tc)
            [e, rr, vv] = tc.arc(6);
            B = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, []);
            tc.verifySize(B, [4 10]);
            tc.verifyEqual(B.effect, ["drag"; "radiation"; "thermal"; "total"]);
            tc.verifyEqual(B.model, repmat("fast", 4, 1));
            tc.verifyEqual(B.Properties.VariableNames{end}, 'share');
            tc.verifyTrue(all(isfinite(B.rms_norm)));
            [~, info] = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 2);
            tc.verifyEqual(info.K, 6); tc.verifyEqual(info.chunks, 3); tc.verifyEqual(info.frame, "body");

            tc.verifyErrorId(@() NonGravForces.app.forceBudget(struct(), e, rr, vv, []), 'NonGravForces:forceBudget:inputModels');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget({}, e, rr, vv, []), 'NonGravForces:forceBudget:inputModels');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget({tc.Mfast, 7}, e, rr, vv, []), 'NonGravForces:forceBudget:inputModels');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget(tc.Mfast, e, rr(:, 1:3), vv, []), 'NonGravForces:forceBudget:inputShape');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, zeros(4, 3)), 'NonGravForces:forceBudget:inputShape');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 0), 'NonGravForces:forceBudget:inputRange');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', NaN), 'NonGravForces:forceBudget:inputRange');
            tc.verifyErrorId(@() NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Labels', ["a" "b"]), 'NonGravForces:forceBudget:inputShape');
            tc.verifyError(@() NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Frame', "rsw"), 'MATLAB:validators:mustBeMember');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function forceBudgetReference(tc)
            % THE ORACLE IS Model.budget, because this unit is a reduction over the same
            % evaluation and nothing else. One chunk covering the whole arc must therefore be
            % bitwise identical: sqrt(sum(x.^2)/K) is what mean(x.^2) computes, in the same order.
            [e, rr, vv] = tc.arc(8);
            Bm = tc.Mfast.budget(e, rr, vv, []);
            Ba = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 8);
            cols = {'rms_x', 'rms_y', 'rms_z', 'max_x', 'max_y', 'max_z', 'rms_norm'};
            worst = 0;
            for c = 1:numel(cols)
                worst = max(worst, max(abs(Ba.(cols{c}) - Bm.(cols{c}))));
            end
            tc.verifyBoundAndRecord("forceBudget vs Model.budget, one chunk", worst, 0, "m/s^2", ...
                "the same reduction over the same evaluation: bitwise or the unit computes physics of its own");
            tc.verifyEqual(Ba.effect, Bm.effect);
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function forceBudgetPrecision(tc)
            % Chunking changes the association of the sum of squares and nothing else, so the RMS
            % may move at the rounding level and the maxima may not move at all.
            [e, rr, vv] = tc.arc(40);
            B1 = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 40);
            B7 = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 7);
            rel = max(abs(B7.rms_norm - B1.rms_norm) ./ max(B1.rms_norm, realmin));
            tc.verifyBoundAndRecord("forceBudget RMS, 40 epochs in chunks of 7 vs one chunk", rel, 1e-14, "-", ...
                "K squares summed in one pass or in six: association only");
            mxDiff = max(abs([B7.max_x - B1.max_x; B7.max_y - B1.max_y; B7.max_z - B1.max_z]));
            tc.verifyBoundAndRecord("forceBudget maxima under chunking", mxDiff, 0, "m/s^2", ...
                "a maximum of maxima is exact whatever the grouping");
            % The share column is formed from rms_norm, so it must reproduce the ratio exactly.
            tc.verifyBoundAndRecord("forceBudget share vs rms_norm ratio", ...
                max(abs(B1.share - B1.rms_norm / B1.rms_norm(end))), 0, "-", "one division, quoted back");
        end
    end

    % ================================================================ vectorization
    methods (Test, TestTags = {'vectorization'})
        function forceBudgetVectorization(tc)
            % Two models in one call must give exactly what two calls give, one each: stacking is
            % bookkeeping and may not touch a number.
            [e, rr, vv] = tc.arc(6);
            Bboth = NonGravForces.app.forceBudget({tc.Mfast, tc.Mstd}, e, rr, vv, [], 'Chunk', 6);
            Bf = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 6);
            Bs = NonGravForces.app.forceBudget(tc.Mstd, e, rr, vv, [], 'Chunk', 6);
            tc.verifyEqual(height(Bboth), 8);
            tc.verifyEqual(Bboth.model, [repmat("fast", 4, 1); repmat("standard", 4, 1)]);
            worst = max([max(abs(Bboth.rms_norm(1:4) - Bf.rms_norm)); max(abs(Bboth.rms_norm(5:8) - Bs.rms_norm))]);
            tc.verifyBoundAndRecord("forceBudget, two models stacked vs one each", worst, 0, "m/s^2", ...
                "stacking is bookkeeping; a difference here would mean the models see each other");
            % Two models of the same preset get distinguishable labels rather than a silent clash.
            Bdup = NonGravForces.app.forceBudget({tc.Mfast, tc.Mfast}, e, rr, vv, [], 'Chunk', 6);
            tc.verifyEqual(unique(Bdup.model), ["fast#1"; "fast#2"]);
        end
    end

    % ==================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function forceBudgetRobustness(tc)
            [e, rr, vv] = tc.arc(5);
            % A single epoch is a legal arc: the RMS is the magnitude and the maximum is the same.
            e1 = referenceSystem.daySecondToEpoch([tc.ep.day(1) tc.ep.sec(1)], tc.ep.scale);
            B1 = NonGravForces.app.forceBudget(tc.Mfast, e1, tc.r(:, 1), tc.v(:, 1), []);
            tc.verifyBoundAndRecord("forceBudget, K = 1: rms equals max", ...
                max(abs([B1.rms_x - B1.max_x; B1.rms_y - B1.max_y; B1.rms_z - B1.max_z])), 0, "m/s^2", ...
                "one sample: the root mean square of one number is its absolute value");
            % A chunk larger than the arc is not an error; it is one chunk.
            [~, i2] = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 1000);
            tc.verifyEqual(i2.chunks, 1);
            % The fast preset sets thermal.model = "none" (Options D.5), so that row must be
            % EXACTLY zero and the shares must stay finite rather than becoming 0/0.
            Bf = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, []);
            tc.verifyEqual(Bf.rms_norm(Bf.effect == "thermal"), 0);
            tc.verifyEqual(Bf.share(Bf.effect == "thermal"), 0);
            tc.verifyTrue(all(isfinite(Bf.share)));
            % Explicit quaternions and the provider's own must agree; passing q is not a second path.
            q = tc.att.at(e);
            Bq = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, q, 'Chunk', 5);
            Bn = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 5);
            tc.verifyBoundAndRecord("forceBudget, explicit q vs the attitude provider", ...
                max(abs(Bq.rms_norm - Bn.rms_norm)), 0, "m/s^2", ...
                "the provider is asked for the same epochs, so the quaternions are the same numbers");
        end
    end

    % =================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function forceBudgetMetamorphic(tc)
            [e, rr, vv] = tc.arc(6);
            B = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 6);
            % 1. Every RMS is bounded by its own maximum, per axis. A reduction that broke this
            %    would be mixing rows or axes.
            tc.verifyTrue(all(B.rms_x <= B.max_x + eps(B.max_x)) && all(B.rms_y <= B.max_y + eps(B.max_y)) ...
                && all(B.rms_z <= B.max_z + eps(B.max_z)));
            % 2. The share of the total row is exactly one, whatever the arc.
            tc.verifyEqual(B.share(B.effect == "total"), 1);
            % 3. The frame is a rotation of the same accelerations, so the NORM of the per-axis RMS
            %    is not invariant (the rotation varies over the arc) but the total maximum magnitude
            %    cannot grow: rotating cannot lengthen a vector.
            Bg = NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 6, 'Frame', "gcrs");
            nb = norm([B.max_x(end), B.max_y(end), B.max_z(end)]);
            ng = norm([Bg.max_x(end), Bg.max_y(end), Bg.max_z(end)]);
            tc.verifyBoundAndRecord("forceBudget, gcrs vs body maximum magnitude ratio", ng / nb, 1.6, "-", ...
                "per-axis maxima can fall on different epochs after a rotation, so the norm of the maxima may grow; it cannot grow without bound");
            % 4. Scaling the gas force by cdScale scales that row and leaves the radiation row
            %    alone. Radiation, not thermal: the fast preset has no thermal recoil at all, so a
            %    thermal invariance here would assert 0 == 0 and measure nothing.
            M2 = NonGravForces.Model(tc.macro, NonGravForces.Options("fast", "drag.cdScale", 2), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, ...
                'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));
            B2 = NonGravForces.app.forceBudget(M2, e, rr, vv, [], 'Chunk', 6);
            d1 = B.rms_norm(B.effect == "drag"); d2 = B2.rms_norm(B2.effect == "drag");
            p1 = B.rms_norm(B.effect == "radiation"); p2 = B2.rms_norm(B2.effect == "radiation");
            tc.assertGreaterThan(p1, 0, 'the radiation row must carry signal or this invariance measures nothing');
            tc.verifyBoundAndRecord("forceBudget: cdScale 2 doubles the drag row", abs(d2 - 2 * d1) / (2 * d1), 1e-15, "-", ...
                "cdScale multiplies the per-panel gas force and nothing else");
            tc.verifyBoundAndRecord("forceBudget: cdScale leaves the radiation row alone", abs(p2 - p1) / p1, 1e-15, "-", ...
                "the effects are separate reductions over separate per-panel arrays");
        end
    end

    % ========================================================================= speed
    methods (Test, TestTags = {'speed'})
        function forceBudgetSpeed(tc)
            % Chunking exists to bound the memory, not to be free. What must be true is that it
            % does not cost a MULTIPLE of the unchunked run: the evaluation dominates and the
            % accumulation is four sums of three numbers per chunk.
            [e, rr, vv] = tc.arc(20);
            NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 20);   % warm
            t = zeros(1, 3); tc2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 20); t(k) = toc(t0);
                t0 = tic; NonGravForces.app.forceBudget(tc.Mfast, e, rr, vv, [], 'Chunk', 5);  tc2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("forceBudget, 4 chunks / 1 chunk over 20 epochs", median(tc2) / median(t), 4.0, "-", ...
                "four Model.accel calls instead of one on the same 20 epochs; the per-call overhead is what this measures, and a ratio near 4 would mean the evaluation does not vectorise over epochs at all");
        end
    end
end
