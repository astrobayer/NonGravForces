classdef tNgfApp < NgfTestCase
%TNGFAPP Tests of the application layer (LP-7): app.forceBudget, app.densityInversion, app.calibrate.
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

    % ============================================================ densityInversion
    methods
        function [e, rr, vv, qq, a, ai] = invScene(tc, n, scale)
            %INVSCENE A model-generated measurement whose gas part is scaled by a known factor.
            % The truth is therefore known exactly, which is the only situation in which an
            % inversion can be tested at all: against real data the answer is unknown by definition.
            idx = 1:n;
            e = referenceSystem.daySecondToEpoch([tc.ep.day(idx) tc.ep.sec(idx)], tc.ep.scale);
            rr = tc.r(:, idx); vv = tc.v(:, idx);
            qq = tc.att.at(e);
            [~, ai] = tc.Mfast.accel(e, rr, vv, qq);
            a = ai.radiation + ai.thermal + scale * ai.drag;
        end
    end
    methods (Test, TestTags = {'contract'})
        function densityInversionContract(tc)
            [e, rr, vv, qq, a] = tc.invScene(6, 1);
            [T, info] = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a);
            tc.verifyEqual(height(T), 6);
            tc.verifyEqual(string(T.Properties.VariableNames), ...
                ["rho_kg_m3" "rho_model_kg_m3" "scale" "dragProjected_m_s2" "dragModel_m_s2" "residualPerp_m_s2" "usable"]);
            tc.verifyTrue(all(T.usable));
            tc.verifyEqual(info.K, 6); tc.verifyEqual(info.usable, 6);
            tc.verifyTrue(istimetable(T));

            D = @NonGravForces.app.densityInversion;
            tc.verifyErrorId(@() D(struct(), e, rr, vv, qq, a), 'NonGravForces:densityInversion:inputModel');
            tc.verifyErrorId(@() D(tc.Mfast, e, rr(:, 1:3), vv, qq, a), 'NonGravForces:densityInversion:inputShape');
            tc.verifyErrorId(@() D(tc.Mfast, e, rr, vv, qq, a(:, 1:3)), 'NonGravForces:densityInversion:inputShape');
            tc.verifyErrorId(@() D(tc.Mfast, e, rr, vv, zeros(4, 3), a), 'NonGravForces:densityInversion:inputShape');
            bad = a; bad(1) = NaN;
            tc.verifyErrorId(@() D(tc.Mfast, e, rr, vv, qq, bad), 'NonGravForces:densityInversion:numericNaN');
            tc.verifyErrorId(@() D(tc.Mfast, e, rr, vv, qq, a, 'Chunk', 0), 'NonGravForces:densityInversion:inputRange');
            tc.verifyErrorId(@() D(tc.Mfast, e, rr, vv, qq, a, 'MinDrag', -1), 'NonGravForces:densityInversion:inputRange');
            tc.verifyError(@() D(tc.Mfast, e, rr, vv, qq, a, 'Frame', "rsw"), 'MATLAB:validators:mustBeMember');
        end
    end
    methods (Test, TestTags = {'reference'})
        function densityInversionReference(tc)
            % THE ROUND TRIP IS THE ORACLE. A model's own acceleration, fed back in, must return
            % the model's own density: the inversion undoes exactly what the forward model did, and
            % anything else would mean the two disagree about what the gas force is.
            [e, rr, vv, qq, a, ai] = tc.invScene(12, 1);
            T = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a);
            tc.verifyBoundAndRecord("density round trip: the model's own acceleration returns its own density", ...
                max(abs(T.rho_kg_m3 - ai.rho(:)) ./ ai.rho(:)), 1e-14, "-", ...
                "the gas force is exactly linear in rho, so the inversion is a division and the round trip is exact to rounding");
            tc.verifyBoundAndRecord("density round trip: the recovered scale is one", ...
                max(abs(T.scale - 1)), 1e-14, "-", "same statement, read off the scale column");
            % A KNOWN scale must come back, not just unity: a test that only ever recovers 1 would
            % pass for an inversion that returns the model's density and ignores the measurement.
            for s = [0.5, 1.7, 3.25]
                [e2, r2, v2, q2, a2, i2] = tc.invScene(8, s);
                T2 = NonGravForces.app.densityInversion(tc.Mfast, e2, r2, v2, q2, a2);
                tc.verifyBoundAndRecord("density inversion recovers a planted scale of " + s, ...
                    max(abs(T2.scale - s)) / s, 1e-13, "-", ...
                    "the measurement's gas part was scaled by a known factor and nothing else was touched");
                tc.verifyBoundAndRecord("planted scale " + s + ": the perpendicular residual stays at zero", ...
                    max(T2.residualPerp_m_s2) / max(max(abs(i2.drag), [], 'all'), realmin), 1e-13, "-", ...
                    "scaling along the drag direction cannot create a component across it");
            end
        end
    end
    methods (Test, TestTags = {'precision'})
        function densityInversionPrecision(tc)
            % Chunking must not move the answer: every epoch is inverted on its own, so a different
            % grouping is a different loop over the same independent divisions.
            [e, rr, vv, qq, a] = tc.invScene(40, 1.3);
            T1 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a, 'Chunk', 40);
            T7 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a, 'Chunk', 7);
            tc.verifyBoundAndRecord("densityInversion, chunks of 7 vs one chunk", ...
                max(abs(T7.rho_kg_m3 - T1.rho_kg_m3)), 0, "kg/m^3", ...
                "each epoch is an independent division; unlike an RMS there is nothing to associate differently");
            % The scale and the density must be consistent with the model density to the last bit.
            tc.verifyBoundAndRecord("rho equals scale times the model density", ...
                max(abs(T1.rho_kg_m3 - T1.scale .* T1.rho_model_kg_m3)), 0, "kg/m^3", ...
                "one multiplication, reported twice");
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function densityInversionVectorization(tc)
            % One epoch at a time must equal the whole arc, epoch by epoch.
            [e, rr, vv, qq, a] = tc.invScene(10, 0.8);
            T = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a);
            worst = 0;
            for k = 1:10
                e1 = referenceSystem.daySecondToEpoch([tc.ep.day(k) tc.ep.sec(k)], tc.ep.scale);
                T1 = NonGravForces.app.densityInversion(tc.Mfast, e1, rr(:, k), vv(:, k), qq(:, k), a(:, k));
                worst = max(worst, abs(T1.rho_kg_m3 - T.rho_kg_m3(k)));
            end
            % NOT bitwise, and the bound said 0 until the runner measured 1.01e-28 kg/m^3 on
            % densities of order 1e-13, i.e. 1e-15 relative. The inversion itself has no
            % cross-epoch term, but Model.accel underneath it does not take the same route for
            % K = 1 as for K = 10 -- vectorised products and norms associate by batch size, which
            % this project has now met three times (pagemtimes in Model.evaluate, the chunked RMS
            % in forceBudget, and here). THE RULE: never claim bitwise across a change in batch
            % size. Claim it only across a change in grouping that leaves the arithmetic alone.
            tc.verifyBoundAndRecord("densityInversion, one epoch at a time vs the whole arc", ...
                worst / max(T.rho_kg_m3), 1e-13, "-", ...
                "the inversion has no cross-epoch term; what is left is Model.accel associating differently at K = 1 than at K = 10");
        end
    end
    methods (Test, TestTags = {'robustness'})
        function densityInversionRobustness(tc)
            [e, rr, vv, qq, a] = tc.invScene(8, 1);
            % A modelled drag below MinDrag is refused, not divided by: usable false and rho NaN.
            T = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a, 'MinDrag', 1);
            tc.verifyFalse(any(T.usable));
            tc.verifyTrue(all(isnan(T.rho_kg_m3)) && all(isnan(T.scale)));
            tc.verifyTrue(all(isfinite(T.rho_model_kg_m3)), 'the model density is still reported where the inversion is refused');

            % A NEGATIVE scale is returned, not clipped. Clipping would turn a visible failure -
            % a measurement pointing against the modelled drag - into a plausible density floor.
            [e2, r2, v2, q2, aNeg] = tc.invScene(8, -0.5);
            T2 = NonGravForces.app.densityInversion(tc.Mfast, e2, r2, v2, q2, aNeg);
            tc.verifyTrue(all(T2.usable));
            tc.verifyBoundAndRecord("a measurement against the modelled drag returns a negative scale", ...
                max(abs(T2.scale + 0.5)) / 0.5, 1e-13, "-", ...
                "the sign is data, not an error: clipping it at zero would hide the case it reports");
            tc.verifyTrue(all(T2.rho_kg_m3 < 0), 'a negative scale must carry through to the density it implies');

            % An error that is NOT along the drag direction lands in the perpendicular residual and
            % changes the density only by its projection. That is the whole point of reporting it.
            perpDir = null(reshape(a(:, 1) / norm(a(:, 1)), 1, 3));
            bump = 0.1 * norm(a(:, 1)) * repmat(perpDir(:, 1), 1, 8);
            T3 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a + bump);
            tc.verifyTrue(all(T3.residualPerp_m_s2 > 0));
            tc.verifyTrue(all(isfinite(T3.rho_kg_m3)));
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function densityInversionMetamorphic(tc)
            [e, rr, vv, qq, a, ai] = tc.invScene(10, 1);
            T1 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a);
            % 1. Doubling the measurement's gas part doubles the inverted density, exactly.
            T2 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, ai.radiation + ai.thermal + 2 * ai.drag);
            tc.verifyBoundAndRecord("doubling the measured gas part doubles the density", ...
                max(abs(T2.rho_kg_m3 - 2 * T1.rho_kg_m3) ./ (2 * T1.rho_kg_m3)), 1e-13, "-", ...
                "linearity in rho, seen from the other side");
            % 2. An error in the modelled radiation moves the density by its projection on the drag
            %    direction and by exactly that: this is the failure mode the help warns about, and
            %    it is asserted rather than described.
            dRad = 0.05 * ai.drag;                              % a radiation error along the drag
            T3 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a + dRad);
            tc.verifyBoundAndRecord("a radiation error along the drag direction moves rho by exactly its size", ...
                max(abs((T3.rho_kg_m3 - T1.rho_kg_m3) ./ T1.rho_kg_m3 - 0.05)), 1e-12, "-", ...
                "nothing separates the two along that direction, which is why residualPerp is reported beside every value");
            % 3. An error ACROSS the drag direction does not move rho at all.
            cross1 = cross(ai.drag, repmat([0; 0; 1], 1, size(ai.drag, 2)));
            cross1 = 0.2 * cross1 ./ max(vecnorm(cross1), realmin) .* vecnorm(ai.drag);
            T4 = NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a + cross1);
            tc.verifyBoundAndRecord("an error across the drag direction leaves rho untouched", ...
                max(abs(T4.rho_kg_m3 - T1.rho_kg_m3) ./ T1.rho_kg_m3), 1e-12, "-", ...
                "the projection is what the inversion reads, so an orthogonal error can only reach residualPerp");
            tc.verifyTrue(all(T4.residualPerp_m_s2 > T1.residualPerp_m_s2));
        end
    end
    methods (Test, TestTags = {'speed'})
        function densityInversionSpeed(tc)
            % The cost is Model.accel and nothing else; the inversion itself is a division per
            % epoch. So chunking must cost what the extra Model.accel calls cost and no more.
            [e, rr, vv, qq, a] = tc.invScene(20, 1);
            NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a, 'Chunk', 20);
            t1 = zeros(1, 3); t4 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a, 'Chunk', 20); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.densityInversion(tc.Mfast, e, rr, vv, qq, a, 'Chunk', 5);  t4(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("densityInversion, 4 chunks / 1 chunk over 20 epochs", ...
                median(t4) / median(t1), 4.0, "-", ...
                "four Model.accel calls instead of one; the divisions are free by comparison");
        end
    end

    % =================================================================== calibrate
    methods
        function [e, aM, y, truth] = calScene(tc, K, biasKind, noise)
            %CALSCENE A measurement built from KNOWN parameters, so the answer is known exactly.
            % An estimator can only be tested where the truth is; on real data it is not, which is
            % why the real-arc check further down asserts an order of magnitude and not a value.
            idx = 1:K;
            e = referenceSystem.daySecondToEpoch([tc.ep.day(idx) tc.ep.sec(idx)], tc.ep.scale);
            [~, ai] = tc.Mfast.accel(e, tc.r(:, idx), tc.v(:, idx), tc.att.at(e));
            aM = ai.drag + ai.radiation + ai.thermal;
            t = (tc.ep.sec(idx) - tc.ep.sec(1))';
            truth = struct('scale', diag([1.05 0.96 1.01]), 'b0', [1e-6; -3e-6; 2e-7]);
            switch biasKind
                case "constant", b = repmat(truth.b0, 1, K);
                otherwise
                    truth.rate = [2e-10; -5e-10; 1e-10];        % per second
                    b = truth.b0 + truth.rate .* t(:)';
            end
            truth.bias = b;
            rng(20260909);                                       % fixed: a test that drifts with
            y = b + truth.scale * aM + noise * randn(3, K);      % the seed is not a measurement
        end
    end
    methods (Test, TestTags = {'contract'})
        function calibrateContract(tc)
            [e, aM, y] = tc.calScene(40, "constant", 0);
            [cal, info] = NonGravForces.app.calibrate(y, aM, e);
            tc.verifySize(cal.scale, [3 3]); tc.verifySize(cal.bias, [3 40]);
            tc.verifySize(cal.sigma0, [3 1]); tc.verifyEqual(numel(cal.covariance), 3);
            tc.verifyEqual(cal.dof, repmat(40 - 2, 3, 1));
            tc.verifyEqual(info.nBasis, 1); tc.verifyEqual(info.nScale, 1);
            tc.verifyEqual(cal.names, ["bias0", "scale_ii"]);
            tc.verifyTrue(all(isfinite(cal.condition)));

            C = @NonGravForces.app.calibrate;
            tc.verifyErrorId(@() C(y(:, 1:5), aM, e), 'NonGravForces:calibrate:inputShape');
            tc.verifyErrorId(@() C(y, aM, e, 'Weights', ones(2, 40)), 'NonGravForces:calibrate:inputShape');
            tc.verifyErrorId(@() C(y, aM, e, 'Mask', true(2, 40)), 'NonGravForces:calibrate:inputShape');
            tc.verifyErrorId(@() C(y, aM, e, 'Weights', -ones(3, 40)), 'NonGravForces:calibrate:inputRange');
            tc.verifyErrorId(@() C(y, aM, e, 'Bias', "hat", 'KnotInterval_s', 0), 'NonGravForces:calibrate:inputRange');
            tc.verifyErrorId(@() C(y, aM, e, 'Events', 1e9), 'NonGravForces:calibrate:inputRange');
            bad = y; bad(1, 3) = NaN;
            tc.verifyErrorId(@() C(bad, aM, e), 'NonGravForces:calibrate:numericNaN');
            % ...unless the offending sample is masked out, which is the point of the mask.
            m = true(3, 40); m(1, 3) = false;
            tc.verifyWarningFree(@() C(bad, aM, e, 'Mask', m));
            % Fewer usable samples than unknowns is refused, not solved.
            tc.verifyErrorId(@() C(y(:, 1:2), aM(:, 1:2), ...
                referenceSystem.daySecondToEpoch([tc.ep.day(1:2) tc.ep.sec(1:2)], tc.ep.scale), ...
                'Bias', "hat", 'KnotInterval_s', 1), 'NonGravForces:calibrate:numericSingular');
        end
    end
    methods (Test, TestTags = {'reference'})
        function calibrateReference(tc)
            % PLANTED PARAMETERS ARE THE ORACLE. The model is linear, so a noise-free arc must
            % return the truth to the conditioning of the QR and to nothing else.
            [e, aM, y, truth] = tc.calScene(60, "constant", 0);
            cal = NonGravForces.app.calibrate(y, aM, e, 'Scale', "diag");
            tc.verifyBoundAndRecord("calibrate recovers a planted diagonal scale", ...
                max(abs(diag(cal.scale) - diag(truth.scale)) ./ diag(truth.scale)), 1e-10, "-", ...
                "noise-free and linear: what is left is the conditioning of the equilibrated solve");
            tc.verifyBoundAndRecord("calibrate recovers a planted constant bias", ...
                max(abs(cal.bias(:, 1) - truth.b0) ./ abs(truth.b0)), 1e-9, "-", ...
                "the bias and the scale are separable here because the modelled acceleration varies over the arc");
            tc.verifyBoundAndRecord("calibrate residuals vanish on an exact fit", ...
                max(abs(cal.residual), [], 'all') / max(abs(y), [], 'all'), 1e-12, "-", ...
                "the truth is in the column space, so the projection removes all of it");
            % A drifting bias needs a basis that can hold it: "constant" cannot, "linear" can, and
            % the difference is what the residual reports. Both halves asserted.
            [e2, aM2, y2, truth2] = tc.calScene(60, "linear", 0);
            calC = NonGravForces.app.calibrate(y2, aM2, e2, 'Bias', "constant");
            calL = NonGravForces.app.calibrate(y2, aM2, e2, 'Bias', "linear");
            tc.verifyBoundAndRecord("a linear basis recovers a drifting bias", ...
                max(abs(calL.bias - truth2.bias), [], 'all') / max(abs(truth2.bias), [], 'all'), 1e-9, "-", ...
                "the drift is in the column space of [1, tau]");
            tc.verifyBoundAndRecord("a constant basis cannot, and the residual says so", ...
                max(abs(calC.residual), [], 'all') / max(abs(calL.residual), [], 'all'), 1e12, "-", ...
                "recorded, not bounded from above by a guess: the point is that it is orders larger, and the run report carries the number");
            tc.verifyGreaterThan(max(abs(calC.residual), [], 'all'), 1e3 * max(abs(calL.residual), [], 'all'));
        end
    end
    methods (Test, TestTags = {'precision'})
        function calibratePrecision(tc)
            % THE 1/sqrt(K) LAW HOLDS FOR THE BIAS AND NOT FOR THE SCALE, and the first version of
            % this test asserted it for the scale and measured 8.3 where it expected 2. The reason
            % is worth keeping: a constant bias has a constant design column, so its formal sigma
            % is exactly sigma0/sqrt(K); the scale's column is the modelled acceleration, whose
            % spread over the arc is what determines it, and 240 epochs of a GRACE orbit do not
            % carry four times the spread of 60. So the law is asserted where it is a law, and the
            % scale's ratio is RECORDED beside it rather than bounded by a wish.
            sig = 2e-9;
            [e1, a1, y1] = tc.calScene(60, "constant", sig);
            [e2, a2, y2] = tc.calScene(240, "constant", sig);
            b1 = NonGravForces.app.calibrate(y1, a1, e1, 'Scale', "none");
            b2 = NonGravForces.app.calibrate(y2, a2, e2, 'Scale', "none");
            ratio = b1.sigma(:, 1) ./ b2.sigma(:, 1);
            tc.verifyBoundAndRecord("formal bias sigma falls as sqrt(K): 240 vs 60 epochs", ...
                max(abs(ratio - 2)) / 2, 0.15, "-", ...
                "a constant column gives sigma0/sqrt(K) exactly; the only spread left is in sigma0 itself, about 5 per cent at 59 degrees of freedom");
            c1 = NonGravForces.app.calibrate(y1, a1, e1);
            c2 = NonGravForces.app.calibrate(y2, a2, e2);
            tc.verifyBoundAndRecord("formal scale sigma, 60 epochs / 240 epochs (recorded, not a law)", ...
                max(c1.sigma(:, end) ./ c2.sigma(:, end)), 20, "-", ...
                "determined by the spread of the modelled acceleration over the arc, not by the sample count; the bound only catches a runaway");
            tc.verifyBoundAndRecord("sigma0 recovers the planted noise level", ...
                max(abs(c2.sigma0 - sig)) / sig, 0.2, "-", ...
                "a-posteriori sigma from a 240-sample arc: the chi-square spread at 237 degrees of freedom is about 9 per cent");
            tc.verifyBoundAndRecord("the scale estimate stays within three formal sigmas of the truth", ...
                max(abs(diag(c2.scale) - [1.05; 0.96; 1.01]) ./ c2.sigma(:, end)), 3, "-", ...
                "unbiasedness, tested the only way one realisation allows");
        end
    end

    methods (Test, TestTags = {'vectorization'})
        function calibrateVectorization(tc)
            % THE PROBLEM SEPARATES BY AXIS, and that is a claim about the estimator, not an
            % implementation detail: perturbing one axis of the measurement may not move another
            % axis's parameters at all. If it did, the three solves would be sharing something.
            [e, aM, y] = tc.calScene(50, "constant", 1e-9);
            c0 = NonGravForces.app.calibrate(y, aM, e, 'Scale', "diag");
            yp = y; yp(2, :) = yp(2, :) + 1e-6;
            cp = NonGravForces.app.calibrate(yp, aM, e, 'Scale', "diag");
            tc.verifyBoundAndRecord("perturbing axis y leaves the x and z parameters untouched", ...
                max(abs([cp.scale(1,1) - c0.scale(1,1); cp.scale(3,3) - c0.scale(3,3); ...
                         cp.bias([1 3], 1) - c0.bias([1 3], 1)])), 0, "-", ...
                "a diagonal scale makes the three axes independent problems; a shared factorisation would leak here");
            % With a FULL scale they are no longer independent, by construction, and that must show.
            cf0 = NonGravForces.app.calibrate(y, aM, e, 'Scale', "full");
            cfp = NonGravForces.app.calibrate(yp, aM, e, 'Scale', "full");
            tc.verifyBoundAndRecord("a full scale still leaves the OTHER rows untouched", ...
                max(abs(cfp.scale([1 3], :) - cf0.scale([1 3], :)), [], 'all'), 0, "-", ...
                "row i of S is fitted to row i of the measurement, so only row 2 may move");
            tc.verifyGreaterThan(max(abs(cfp.scale(2, :) - cf0.scale(2, :))), 0);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function calibrateRobustness(tc)
            % An outlier must not carry the fit. Huber IRLS is asserted against the ordinary
            % solution on the same data: both are computed, and the robust one must stay near the
            % truth while the ordinary one does not.
            [e, aM, y, truth] = tc.calScene(120, "constant", 1e-9);
            yBad = y; yBad(1, 40:44) = yBad(1, 40:44) + 5e-6;    % five gross outliers on one axis
            cOls = NonGravForces.app.calibrate(yBad, aM, e);
            cRob = NonGravForces.app.calibrate(yBad, aM, e, 'Robust', "huber");
            eOls = abs(cOls.scale(1,1) - truth.scale(1,1));
            eRob = abs(cRob.scale(1,1) - truth.scale(1,1));
            tc.verifyBoundAndRecord("Huber IRLS beats the ordinary solve on five gross outliers", ...
                eRob / max(eOls, realmin), 1, "-", ...
                "the ratio is the claim; the absolute errors depend on where the outliers fall and are recorded in the run log");
            tc.verifyGreaterThan(cRob.iterations(1), 1, 'the robust path must actually iterate');
            % A masked gap is not a zero: the estimate must match a call on the trimmed arc.
            m = true(3, 120); m(:, 40:44) = false;
            cMask = NonGravForces.app.calibrate(yBad, aM, e, 'Mask', m);
            keep = [1:39, 45:120];
            eTrim = referenceSystem.daySecondToEpoch([tc.ep.day(keep) tc.ep.sec(keep)], tc.ep.scale);
            cTrim = NonGravForces.app.calibrate(yBad(:, keep), aM(:, keep), eTrim);
            tc.verifyBoundAndRecord("masking samples equals removing them", ...
                max(abs(diag(cMask.scale) - diag(cTrim.scale)) ./ diag(cTrim.scale)), 1e-12, "-", ...
                "the masked samples never enter a sum, so the two problems are the same problem");
            % A weight of zero is a mask.
            w = ones(3, 120); w(:, 40:44) = 0;
            cW = NonGravForces.app.calibrate(yBad, aM, e, 'Weights', w);
            tc.verifyBoundAndRecord("a zero weight equals a mask", ...
                max(abs(diag(cW.scale) - diag(cMask.scale))), 0, "-", "both drop the sample before the design matrix is built");
            % MORE UNKNOWNS THAN SAMPLES is refused, not fitted. A 2 s knot interval was tried
            % here first and did NOT refuse, correctly: at this sampling it still leaves more
            % samples than knots, so it is estimable and merely useless. 0.5 s puts the knot count
            % above the sample count whatever the fixture's rate, which is the guard's actual job.
            tc.verifyErrorId(@() NonGravForces.app.calibrate(y, aM, e, 'Bias', "hat", 'KnotInterval_s', 0.5), ...
                'NonGravForces:calibrate:numericSingular');
            % The condition guard is a separate branch and gets its own provocation, paired with
            % the call that must succeed (rule T6).
            tc.verifyErrorId(@() NonGravForces.app.calibrate(y, aM, e, 'MaxCondition', 1), ...
                'NonGravForces:calibrate:numericSingular');
            tc.verifyWarningFree(@() NonGravForces.app.calibrate(y, aM, e, 'MaxCondition', 1e12));
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function calibrateMetamorphic(tc)
            [e, aM, y, truth] = tc.calScene(60, "constant", 0);
            c = NonGravForces.app.calibrate(y, aM, e);
            % 1. Adding a constant to the measurement moves the bias by it and the scale not at all.
            c2 = NonGravForces.app.calibrate(y + 7e-6, aM, e);
            tc.verifyBoundAndRecord("a constant added to the measurement lands entirely in the bias", ...
                max(abs((c2.bias - c.bias) - 7e-6), [], 'all') / 7e-6, 1e-9, "-", ...
                "the constant is in the bias basis, so the projection puts it there and nowhere else");
            tc.verifyBoundAndRecord("...and leaves the scale alone", ...
                max(abs(diag(c2.scale) - diag(c.scale))) / max(diag(c.scale)), 1e-9, "-", ...
                "the scale column is orthogonalised against the constant by the solve");
            % 2. Doubling the modelled acceleration halves the scale and leaves the bias.
            c3 = NonGravForces.app.calibrate(y, 2 * aM, e);
            tc.verifyBoundAndRecord("doubling the model halves the scale", ...
                max(abs(diag(c3.scale) - diag(c.scale) / 2) ./ (diag(c.scale) / 2)), 1e-9, "-", ...
                "the product S a is what is observed, so only the product is determined");
            % 3. A declared jump is recovered where a constant basis alone would smear it.
            yJ = y; yJ(:, 31:end) = yJ(:, 31:end) + [4e-6; -2e-6; 1e-6];
            tJump = tc.ep.sec(31) - tc.ep.sec(1);
            cJ = NonGravForces.app.calibrate(yJ, aM, e, 'Events', tJump);
            tc.verifyBoundAndRecord("a declared jump is recovered", ...
                max(abs(cJ.jump - [4e-6; -2e-6; 1e-6]) ./ abs([4e-6; -2e-6; 1e-6])), 1e-9, "-", ...
                "the step is in the column space once it is declared");
            tc.verifyBoundAndRecord("...and the scale survives the jump", ...
                max(abs(diag(cJ.scale) - diag(truth.scale)) ./ diag(truth.scale)), 1e-9, "-", ...
                "a jump absorbed by its own column cannot bias the scale");
        end
    end
    methods (Test, TestTags = {'speed'})
        function calibrateSpeed(tc)
            % Three QR solves of a few unknowns: the cost must be linear in the arc length, because
            % the design matrix grows in rows only. Four times the epochs, not sixteen times the cost.
            [e1, a1, y1] = tc.calScene(60, "constant", 1e-9);
            [e2, a2, y2] = tc.calScene(240, "constant", 1e-9);
            NonGravForces.app.calibrate(y2, a2, e2);
            t1 = zeros(1, 3); t2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.app.calibrate(y1, a1, e1); t1(k) = toc(t0);
                t0 = tic; NonGravForces.app.calibrate(y2, a2, e2); t2(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("calibrate, 240 epochs / 60 epochs", median(t2) / median(t1), 8, "-", ...
                "linear in rows would be 4 and a fixed overhead pulls it below; quadratic would be 16");
        end
    end

end
