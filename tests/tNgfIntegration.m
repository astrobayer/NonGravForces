classdef tNgfIntegration < NgfTestCase
%TNGFINTEGRATION Tests of app.DynamicsSource: the integrator interface and its native C evaluator (LP-C, WP-C9).
%
%   THE ORACLE IS Model.accel, which is normative (CLAUDE.md sect. 2). The native
%   evaluator ports the whole chain to C, so every claim here compares the two at
%   the TABLE NODES, where the time tables are exact by construction, and records
%   what the interpolation costs between them. The second oracle is gjToolbox
%   itself: its frozen checkNativeEntry, and a propagated arc whose two bindings
%   (MATLAB accel, C accel with the C stepper) must give one trajectory.
%
%   This class needs gjToolbox on the path (gjToolboxSetup) and the rhs MEX built;
%   it FAILS, never skips, when either is missing.
% Created 2026-09-16 by Claude Opus 5 (Anthropic).
    properties
        root, eop, sw, ep, r, v, att, macro, tai, cer, eph
        combinedScale                  % the standard preset's density scale, built once per class (scaleOf)
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            % exist() does not resolve package-qualified names; which() does, and it is what
            % buildMex uses to find gjToolbox's csrc in the first place.
            tc.assertNotEmpty(which('gjToolbox.DynamicsSource'), ...
                'gjToolbox is not on the path: run gjToolboxSetup from the odeSuite checkout (the app layer serves its frozen interface)');
            rep = NonGravForces.buildMex('Check', true, 'Target', "rhs");
            tc.assertTrue(rep.upToDate, 'the rhs MEX must be built and current before this class runs (ciGate builds it)');
            f = @(varargin) fullfile(tc.root, varargin{:});
            tc.eop = referenceSystem.readEop(f('data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            tc.sw = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [tc.ep, tc.r, tc.v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', tc.eop);
            [eps2, qq] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
            tc.att = NonGravForces.geom.AttitudeProvider(eps2, qq);
            tc.macro = NonGravForces.geom.graceFO("C");
            tc.cer = NonGravForces.io.importCERES(string(f('data', 'ceres', ["rsut" "rlut" "rsdt"] + "_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc")));
            % The precise preset's own sun.source is rs.de440, so the kernel is a requirement of this
            % class, not an extra. It is fetched by setupNonGravForces (FetchKernels) and never committed.
            kernel = f('data', 'kernels', 'de440s.bsp');
            tc.assertTrue(isfile(kernel), 'data/kernels/de440s.bsp is missing: run setupNonGravForces so the precise preset can be evaluated as shipped');
            tc.eph = referenceSystem.readEphemeris(string(kernel));
            e1 = referenceSystem.daySecondToEpoch([tc.ep.day(1) tc.ep.sec(1)], tc.ep.scale);
            tc.tai = referenceSystem.epochTo(referenceSystem.convertTimeScale(e1, "TAI"), "datetime");
        end
    end
    methods
        function M = modelFor(tc, o)
            %MODELFOR A Model of the sample day with everything the presets can ask for -- since 2026-10-29 that
            % includes the standard preset's combined density scale, from the committed Swarm and CLS excerpts.
            extra = {};
            if o.drag.densityScale == "combined"
                if isempty(tc.combinedScale), tc.combinedScale = tc.scaleOf("combined"); end
                extra = {'DensityScale', tc.combinedScale};
            end
            M = NonGravForces.Model(tc.macro, o, 'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, ...
                'Ceres', tc.cer, 'Ephemeris', tc.eph, 'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3), extra{:});
        end
        function o = optionsFor(~, name)
            %OPTIONSFOR The seven configurations the parity is claimed over, INCLUDING the precise preset
            % exactly as shipped: since 2026-09-16 the DE440 kernel is fetched at setup, so sun.source
            % rs.de440 is evaluated here rather than substituted away (HANDOVER O-7, O-15).
            switch name
                case "precise"
                    o = NonGravForces.Options("precise");
                case "walker2014"
                    % Since O-42 (2026-09-17) the precise preset takes a FIXED alpha, so walker2014 --
                    % accommodation enum value 3 in the C evaluator -- is reached by no preset. Without
                    % this configuration it would leave the parity claim entirely, which is the same
                    % blindness the hwm14 case below was added to prevent. The mass is Mehta et al.
                    % (2022)'s, the value the preset used to ship.
                    o = NonGravForces.Options("precise").set("drag.alpha", []) ...
                        .set("drag.accommodation", "walker2014").set("drag.surfaceMass_amu", 263.3223);
                case "sesam"
                    % WP-E3: accommodation enum value 4, and the ONLY model whose pressure needs the
                    % flow. If the C read the thermodynamic pressure instead of eq. (15)'s ram
                    % pressure both paths would still run and would differ by a factor of 140.
                    o = NonGravForces.Options("precise").set("drag.alpha", []).set("drag.accommodation", "sesam");
                case "hwm14"
                    % pilinski2010 rides along here because the C accommodation enum gained a value in
                    % LP-D WP-D1 and an enum that MATLAB and C number differently is invisible in every
                    % other test: both paths would run, and quietly evaluate different models.
                    o = NonGravForces.Options("standard").set("drag.wind", "hwm14").set("drag.gsi", "maxwell") ...
                        .set("drag.specularFraction", 0.1).set("erp.backend", "pointSource").set("srp.lunarEclipse", true) ...
                        .set("drag.alpha", []).set("drag.accommodation", "pilinski2010");
                case "pilinski2011"
                    % msis00 is the backend both Pilinski fits were made with, so this is the pairing
                    % their K is actually valid for (the precise preset's msis21 is not).
                    o = NonGravForces.Options("standard").set("drag.alpha", []).set("drag.accommodation", "pilinski2011");
                case "standardSteady"
                    % the standard preset with the steady-state walls it had until 2026-10-31: for the tests that
                    % compare the native path with Model.accel THROUGH THE SOURCE ('Native', false), which the
                    % transient walls cannot serve (they need the native temperature table)
                    o = NonGravForces.Options("standard").set("thermal.model", "steadyState");
                otherwise
                    o = NonGravForces.Options(name);
            end
        end
        function M = scaledModel(tc, kind, varargin)
            %SCALEDMODEL The standard preset with drag.densityScale of the given kind, from the committed Swarm and
            % CLS excerpts; extra arguments go to env.densityScale (e.g. 'MinRecords' to force the fallback).
            S = tc.scaleOf(kind, varargin{:});
            % steady walls: these tests are about the density scale, and one of them propagates an arc ('Native' false)
            M = NonGravForces.Model(tc.macro, NonGravForces.Options("standard").set("drag.densityScale", kind).set("thermal.model", "steadyState"), 'Attitude', tc.att, ...
                'Eop', tc.eop, 'SpaceWeather', tc.sw, 'Ceres', tc.cer, 'Ephemeris', tc.eph, 'DensityScale', S);
        end
        function S = scaleOf(tc, kind, varargin)
            f = @(varargin) fullfile(tc.root, varargin{:});
            files = string(f('data', 'samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
            cls = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
            switch kind
                case "swarm"
                    S = NonGravForces.env.densityScale(files, tc.sw, varargin{:});
                case "driver"
                    S = NonGravForces.env.densityDriverScale(tc.sw, cls, varargin{:});
                case "combined"
                    S = NonGravForces.env.densityScale(files, tc.sw, varargin{:}, 'Fallback', NonGravForces.env.densityDriverScale(tc.sw, cls));
            end
        end
        function src = preparedFor(tc, M, varargin)
            %PREPAREDFOR The native source of M on the sample day; transient walls (the standard preset since
            % 2026-10-31) are integrated along the sample day's own GNV1B, as 'ThermalOrbit' requires.
            extra = {};
            if M.options.thermal.model == "transient", extra = {'ThermalOrbit', tc.orbitOf(1:numel(tc.ep.sec))}; end
            src = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60, extra{:}, varargin{:});
            src = src.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
        end
        function [src, M] = sourceFor(tc, name, varargin)
            M = tc.modelFor(tc.optionsFor(name));
            src = tc.preparedFor(M, varargin{:});
        end
        function orb = orbitOf(tc, idx)
            %ORBITOF The sample day's GNV1B (GCRS) over the given records, as 'ThermalOrbit' takes it.
            orb = struct('epoch', referenceSystem.daySecondToEpoch([reshape(tc.ep.day(idx), [], 1) reshape(tc.ep.sec(idx), [], 1)], tc.ep.scale), ...
                'r', tc.r(:, idx), 'v', tc.v(:, idx));
        end
        function a = modelAccel(tc, M, k)
            %MODELACCEL Model.accel at record k -- for transient walls along the WHOLE sample orbit, as the native
            % table integrates them, and then record k: one epoch alone would start the walls at their steady state.
            if M.options.thermal.model == "transient"
                orb = tc.orbitOf(1:numel(tc.ep.sec));
                A = M.accel(orb.epoch, orb.r, orb.v, [], 'Frame', "gcrs");
                a = A(:, k);
                return
            end
            epk = referenceSystem.daySecondToEpoch([tc.ep.day(k) tc.ep.sec(k)], tc.ep.scale);
            a = M.accel(epk, tc.r(:, k), tc.v(:, k), [], 'Frame', "gcrs");
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
        function DynamicsSourceContract(tc)
            M = tc.modelFor(tc.optionsFor("standardSteady"));   % the refusals below test other arguments than the walls
            src = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60);
            tc.verifyEqual(src.StateSize, 3); tc.verifyTrue(src.Vectorized); tc.verifyFalse(src.IsPrepared);
            tc.verifyErrorId(@() src.accel(0, tc.r(:, 1), tc.v(:, 1)), 'gjToolbox:dynamics:notPrepared');
            % WP-E2: erp.backend ceresClimatology has no native evaluator and must REFUSE, because
            % find() on a name the C enum lacks returns empty and the packed value would become 0,
            % which is "none" -- a Model that computes Earth radiation silently becoming one that
            % does not, on this path only. The parity test cannot catch what it cannot construct,
            % so the refusal is the check (HANDOVER O-44).
            Mc = tc.modelFor(NgfTestCase.standardUnscaled().set("erp.backend", "ceresClimatology"));
            tc.verifyErrorId(@() NonGravForces.app.DynamicsSource(Mc, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60), ...
                'NonGravForces:DynamicsSource:unsupportedOption');
            % The same rule for drag.density dtm2020 (2026-10-27, O-52): the C evaluator links no DTM2020.
            Md = tc.modelFor(NgfTestCase.standardUnscaled().set("drag.density", "dtm2020"));
            tc.verifyErrorId(@() NonGravForces.app.DynamicsSource(Md, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60), ...
                'NonGravForces:DynamicsSource:unsupportedOption');
            % The observed density scale was refused here from 2026-10-27 to 2026-10-29; the native path now
            % carries it as a daily table (O-19), so every kind constructs and prepares.
            for kind = ["swarm" "driver" "combined"]
                Ms = tc.scaledModel(kind);
                ss = NonGravForces.app.DynamicsSource(Ms, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60);
                ss = ss.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
                tc.verifyTrue(ss.IsPrepared, kind);
                ss.release();
            end
            ctx = struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 7, 'Epoch', tc.tai);
            src = src.prepare(ctx);
            tc.verifyTrue(src.IsPrepared);
            tc.verifyEqual(sort(string(fieldnames(src.NativeEntry))'), sort(["AbiVersion" "EntryPtr" "StateSize" "MaxColumns"]));
            tc.verifyEqual(src.NativeEntry.MaxColumns, 7); tc.verifyEqual(src.NativeEntry.StateSize, 3);
            tc.verifyWarningFree(@() gjToolbox.DynamicsSource.checkNativeEntry(src.NativeEntry, 3));   % gjToolbox's own rule
            tc.verifyEqual(src.ParameterSpec.Count, 0); tc.verifyEqual(src.ParameterSpec.Id, "NonGravForces:standard");
            tc.verifyFalse(src.hasGradient());
            tc.verifyEqual(numel(src.BreakpointKinds), numel(src.Breakpoints));
            A = src.accel(180, tc.r(:, [181 182]), tc.v(:, [181 182]));
            tc.verifySize(A, [3 2]); tc.verifyTrue(all(isfinite(A), 'all'));
            tc.verifyErrorId(@() src.accel(180, tc.r(1:2, 181), tc.v(1:2, 181)), 'NonGravForces:DynamicsSource:inputModel');
            src = src.release();
            tc.verifyFalse(src.IsPrepared); tc.verifyEmpty(src.NativeEntry);
            % refusals at construction, each paired with the construction that must succeed (T6)
            D = @(varargin) NonGravForces.app.DynamicsSource(varargin{:});
            tc.verifyErrorId(@() D(struct(), 'Epoch', tc.tai), 'NonGravForces:DynamicsSource:inputModel');
            tc.verifyErrorId(@() D(M), 'NonGravForces:DynamicsSource:inputEpoch');
            tc.verifyErrorId(@() D(M, 'Epoch', datetime(2020, 1, 1, 'TimeZone', 'UTC')), 'NonGravForces:DynamicsSource:inputEpoch');
            tc.verifyErrorId(@() D(M, 'Epoch', tc.tai, 'Parameters', "alpha"), 'NonGravForces:DynamicsSource:inputModel');
            tc.verifyErrorId(@() D(M, 'Epoch', tc.tai, 'Gravity', struct()), 'NonGravForces:DynamicsSource:inputGravity');
            % thermal.model transient (2026-10-30): served by a temperature table integrated along 'ThermalOrbit',
            % so it is refused without that orbit, refused off the native path, and refused at prepare when the
            % orbit does not cover the tables -- each beside the construction that must succeed.
            Mt = tc.modelFor(NonGravForces.Options("standard").set("thermal.model", "transient"));
            orbAll = tc.orbitOf(1:numel(tc.ep.sec));
            tc.verifyErrorId(@() D(Mt, 'Epoch', tc.tai), 'NonGravForces:DynamicsSource:inputModel');
            tc.verifyErrorId(@() D(Mt, 'Epoch', tc.tai, 'Native', false, 'ThermalOrbit', orbAll), 'NonGravForces:DynamicsSource:inputModel');
            sShort = D(Mt, 'Epoch', tc.tai, 'Margin_s', 30, 'ThermalOrbit', tc.orbitOf(1:200));
            tc.verifyErrorId(@() sShort.prepare(struct('Span', [120 420], 'StencilSpan', [120 420])), 'NonGravForces:DynamicsSource:dataRange');
            sFull = D(Mt, 'Epoch', tc.tai, 'Margin_s', 30, 'ThermalOrbit', orbAll);
            sFull = sFull.prepare(struct('Span', [120 420], 'StencilSpan', [120 420]));
            tc.verifyTrue(sFull.IsPrepared);
            sFull.release();
            % a context epoch that disagrees, and a native path the options cannot serve
            s2 = D(M, 'Epoch', tc.tai);
            tc.verifyErrorId(@() s2.prepare(struct('Span', [0 60], 'StencilSpan', [0 60], 'Epoch', tc.tai + seconds(1))), 'NonGravForces:DynamicsSource:inputEpoch');
            s3 = D(tc.modelFor(tc.optionsFor("standardSteady").set("geom.selfShadow", "raycast").set("geom.flowShadow", "exact")), 'Epoch', tc.tai, 'Margin_s', 30);
            tc.verifyErrorId(@() s3.prepare(struct('Span', [120 420], 'StencilSpan', [120 420])), 'NonGravForces:DynamicsSource:inputModel');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function DynamicsSourceReference(tc)
            % Model.accel is the reference and the native evaluator is the port: at a table node the two
            % are two implementations of the same expressions, and nothing else differs. Seven
            % configurations, chosen to reach every branch the native path serves: exponential density with
            % constantCd and no shadowing (fast); msis00, sentman, Knocke tiles, the cube-map shadow and the
            % steady thermal model (standard); msis21, CERES, SOLAARS-CF and a 32-cell table (precise);
            % HWM14 winds, maxwell, the point source and the lunar shadow (hwm14). The last two exist
            % only to reach the accommodation enum: MATLAB and C number it independently, so a value no
            % configuration selects would let both paths run and quietly evaluate different models.
            worstNode = 0; worstBetween = 0; worstBetweenTransient = 0;
            names = ["fast" "standard" "precise" "hwm14" "pilinski2011" "walker2014" "sesam"];
            for name = names
                [src, M] = tc.sourceFor(name);
                for k = [121 151 181 241 301 361 391 421]
                    t = double(tc.ep.sec(k) - tc.ep.sec(1));
                    a1 = src.accel(t, tc.r(:, k), tc.v(:, k));
                    a0 = tc.modelAccel(M, k);
                    rel = norm(a1 - a0) / norm(a0);
                    if mod(t, 60) == 0
                        worstNode = max(worstNode, rel);
                    elseif M.options.thermal.model == "transient"
                        worstBetweenTransient = max(worstBetweenTransient, rel);
                    else
                        worstBetween = max(worstBetween, rel);
                    end
                end
                src.release();
            end
            tc.verifyBoundAndRecord("native evaluator vs Model.accel at the table nodes (7 configurations)", worstNode, 1e-13, "-", ...
                "accumulated rounding of two implementations of one expression set, written in the same association order where it decides the digits; measured 2.504e-14 on 2026-09-16");
            tc.verifyBoundAndRecord("native evaluator vs Model.accel between the table nodes (60 s tables)", worstBetween, 1e-9, "-", ...
                "the cubic tables' interpolation error, not the model's: measured 8.124e-11 at 60 s spacing, and it falls as the fourth power of the step");
            % Since 2026-10-31 standard, hwm14 and pilinski2011 carry transient walls, whose temperatures the native
            % path ALSO reads from a cubic table: between nodes that table's error comes on top, and it is the
            % temperatures' (an eclipse edge is not smooth), not the time tables'. CI 41152904075 measured 1.564e-05
            % for it against the 1e-9 claimed for the time tables alone; the two are now recorded apart.
            tc.verifyBoundAndRecord("native evaluator vs Model.accel between the table nodes, transient walls (60 s tables)", worstBetweenTransient, 1e-4, "-", ...
                "the wall-temperature table's cubic between nodes; measured 1.564e-05 (CI 41152904075); shrink TableStep_s where it matters");
        end
    end

    methods (Test, TestTags = {'reference'})
        function DynamicsSourceThermalTableReference(tc)
            % thermal.model transient on the native path (2026-10-30, the owner's ruling that standard goes
            % transient and stays integrable): prepare integrates the walls ONCE along 'ThermalOrbit' with
            % Model.accel, and the C reads the temperatures by time. On that orbit and at a table node (an orbit
            % epoch) the native path and Model.accel of the SAME orbit evaluate the same temperatures, so they
            % are again two implementations of one expression set; the thermal recoil is compared on its own
            % through the diagnostics, so a parity carried by the drag cannot pass for one of the walls.
            M = tc.modelFor(NonGravForces.Options("standard").set("thermal.model", "transient"));
            orb = tc.orbitOf(1:numel(tc.ep.sec));
            src = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60, 'ThermalOrbit', orb);
            src = src.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
            aG = M.accel(orb.epoch, orb.r, orb.v, [], 'Frame', "gcrs");
            [~, infoB] = M.accel(orb.epoch, orb.r, orb.v, [], 'Frame', "body");
            worstNode = 0; worstThermal = 0; worstBetween = 0;
            for k = [121 151 181 241 301 361 391 421]
                t = double(tc.ep.sec(k) - tc.ep.sec(1));
                a1 = src.accel(t, tc.r(:, k), tc.v(:, k));
                rel = norm(a1 - aG(:, k)) / norm(aG(:, k));
                if mod(t, 60) == 0
                    worstNode = max(worstNode, rel);
                    d = src.diagnostics(t, tc.r(:, k), tc.v(:, k));
                    worstThermal = max(worstThermal, norm(d.thermal - infoB.thermal(:, k)) / norm(infoB.thermal(:, k)));
                else
                    worstBetween = max(worstBetween, rel);
                end
            end
            src.release();
            tc.verifyBoundAndRecord("native transient table vs Model.accel at the table nodes, on the reference orbit", worstNode, 1e-13, "-", ...
                "the same parity as the seven configurations: the temperatures are the same numbers, the recoil the same expression");
            tc.verifyBoundAndRecord("native transient thermal recoil vs Model's at the table nodes", worstThermal, 1e-13, "rel", ...
                "kernel.thermalRecoil and its C port on identical temperatures");
            tc.verifyBoundAndRecord("native transient table vs Model.accel between the nodes (60 s tables)", worstBetween, 1e-2, "-", ...
                "the cubic reads the wall temperatures as smooth between nodes; recorded, not claimed (shrink TableStep_s where it matters)");
        end
    end

    methods (Test, TestTags = {'reference'})
        function DynamicsSourceDensityScaleReference(tc)
            % drag.densityScale on the native path (2026-10-29, O-19): the daily table against Model.accel at the
            % table nodes, for Swarm's bracket (interpolated in height in C), the driver scale (one factor per
            % day) and a combined scale forced onto its fallback (MinRecords above any day's count). The factor
            % itself is compared through the diagnostics, so a parity that holds only because the drag is small
            % against the radiation cannot pass for one that holds on the scale.
            cases = {{"swarm"}, {"driver"}, {"combined", 'MinRecords', 2000}};
            worst = 0; worstF = 0;
            for c = 1:numel(cases)
                M = tc.scaledModel(cases{c}{:});
                src = tc.preparedFor(M);
                for k = [121 181 241 301 361 421]
                    t = double(tc.ep.sec(k) - tc.ep.sec(1));
                    a1 = src.accel(t, tc.r(:, k), tc.v(:, k));
                    a0 = tc.modelAccel(M, k);
                    worst = max(worst, norm(a1 - a0) / norm(a0));
                    d = src.diagnostics(t, tc.r(:, k), tc.v(:, k));
                    epk = referenceSystem.daySecondToEpoch([tc.ep.day(k) tc.ep.sec(k)], tc.ep.scale);
                    [~, info] = M.accel(epk, tc.r(:, k), tc.v(:, k), [], 'Frame', "body");
                    worstF = max(worstF, abs(d.densityScale / info.densityScale - 1));
                end
                src.release();
            end
            tc.verifyBoundAndRecord("native evaluator vs Model.accel at the table nodes, drag.densityScale swarm / driver / combined", worst, 1e-13, "-", ...
                "the same parity as the seven configurations: the scale is one multiplication after the density");
            tc.verifyBoundAndRecord("native density scale vs Model's", worstF, 1e-14, "rel", ...
                "Swarm's bracket: the same log-linear expression in the same order; driver days: the factor itself, computed by S.at");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function DynamicsSourcePrecision(tc)
            % The two PATHS of one source: the native entry and Model.accel, at a node, with everything else
            % identical. And the gravity source is added, not folded in: with it the acceleration is the sum.
            [srcN, M] = tc.sourceFor("standardSteady");   % 'Native' false cannot serve transient walls
            srcM = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60, 'Native', false);
            srcM = srcM.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
            k = 181; t = double(tc.ep.sec(k) - tc.ep.sec(1));
            aN = srcN.accel(t, tc.r(:, k), tc.v(:, k));
            aM = srcM.accel(t, tc.r(:, k), tc.v(:, k));
            tc.verifyBoundAndRecord("the source's two paths at a table node", norm(aN - aM) / norm(aM), 1e-12, "-", "the C evaluator against Model.accel through the same source");
            % the diagnostics of the native pass reproduce the Model's own environment
            d = srcN.diagnostics(t, tc.r(:, k), tc.v(:, k));
            epk = referenceSystem.daySecondToEpoch([tc.ep.day(k) tc.ep.sec(k)], tc.ep.scale);
            [~, info] = M.accel(epk, tc.r(:, k), tc.v(:, k), [], 'Frame', "body");
            tc.verifyBoundAndRecord("native density vs Model's", abs(d.rho / info.rho - 1), 1e-13, "rel", "the same backend called with the same drivers");
            tc.verifyBoundAndRecord("native relative velocity vs Model's", norm(d.vrel_body - info.vrel_body) / norm(info.vrel_body), 1e-12, "-", "the co-rotation from the same rotation tables");
            tc.verifyBoundAndRecord("native shadow factor vs Model's", abs(d.nu - info.nu), 1e-12, "-", "the same eclipse model");
            tc.verifyEqual(d.mass_kg, info.mass_kg, 'the mass table is exact at a node');
            srcN.release(); srcM.release();
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function DynamicsSourceVectorization(tc)
            % An n-by-m block equals m single-column calls, BITWISE: the evaluator loops over columns and
            % carries nothing between them, which is what makes the integrator's ensembles trustworthy.
            src = tc.sourceFor("standard");
            k = [151 181 241 301];
            t = double(tc.ep.sec(181) - tc.ep.sec(1));
            A = src.accel(t, tc.r(:, k), tc.v(:, k));
            for j = 1:numel(k)
                tc.verifyEqual(A(:, j), src.accel(t, tc.r(:, k(j)), tc.v(:, k(j))), sprintf('column %d batched vs alone', j));
            end
            src.release();
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function DynamicsSourceRobustness(tc)
            src = tc.sourceFor("standard");
            % outside the tables the evaluator REFUSES; it does not extrapolate the rotation or the drivers
            tc.verifyErrorId(@() src.accel(100000, tc.r(:, 181), tc.v(:, 181)), 'NonGravForces:rhsGateway:dataRange');
            % a released handle is refused by name, not dereferenced
            released = src.release();
            tc.verifyEmpty(released.NativeEntry);
            % two sources alive at once keep separate contexts
            s1 = tc.sourceFor("standard");
            s2 = tc.sourceFor("fast");
            a1 = s1.accel(180, tc.r(:, 181), tc.v(:, 181));
            a2 = s2.accel(180, tc.r(:, 181), tc.v(:, 181));
            tc.verifyTrue(all(isfinite([a1; a2])) && norm(a1 - a2) > 0, 'two prepared sources must not share a context');
            s1.release(); s2.release();
            % the attitude must cover the span plus the margin: a wider margin than the excerpt is refused
            M = tc.modelFor(tc.optionsFor("standardSteady"));   % the 'Native' false sources below cannot serve transient walls
            wide = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 3600);
            tc.verifyErrorId(@() wide.prepare(struct('Span', [120 420], 'StencilSpan', [120 420])), 'NonGravForces:DynamicsSource:dataRange');
            % a gravity source without a native entry cannot be added by the C evaluator, and says so
            g = gjToolbox.FunctionHandleSource(@(t, R, V) zeros(size(R)), 'StateSize', 3, 'Vectorized', "on");
            withG = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'Gravity', g);
            tc.verifyErrorId(@() withG.prepare(struct('Span', [120 420], 'StencilSpan', [120 420])), 'NonGravForces:DynamicsSource:inputGravity');
            % ... and the MATLAB path adds exactly that source's acceleration
            withGm = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'Gravity', g, 'Native', false);
            withGm = withGm.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
            plain = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'Native', false);
            plain = plain.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
            tc.verifyEqual(withGm.accel(180, tc.r(:, 181), tc.v(:, 181)), plain.accel(180, tc.r(:, 181), tc.v(:, 181)), ...
                'a zero gravity source adds zero, bitwise');
            withGm.release(); plain.release();
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function DynamicsSourceMetamorphicFlowTable(tc)
            % The reference test's parity covers the flow-shadow table only if BOTH paths read it: the standard
            % preset's native acceleration agrees with Model.accel at a node, and the exact flow shadow moves the
            % result by far more than the two paths differ.
            [src, Mt] = tc.sourceFor("standard");
            Me = tc.modelFor(NonGravForces.Options("standard").set("geom.flowShadow", "exact"));
            k = 181; t = double(tc.ep.sec(k) - tc.ep.sec(1));
            aN = src.accel(t, tc.r(:, k), tc.v(:, k));
            src.release();
            aT = tc.modelAccel(Mt, k); aE = tc.modelAccel(Me, k);
            tc.verifyBoundAndRecord("native vs Model.accel, both on the flow-shadow table, at a node", norm(aN - aT) / norm(aT), 1e-13, "-", ...
                "the C trilinear is MATLAB's in the same order");
            tc.verifyGreaterThan(norm(aT - aE), 100 * norm(aN - aT));
        end
        function DynamicsSourceMetamorphic(tc)
            % 1. cdScale multiplies the gas part and nothing else, through the native evaluator.
            [s1, M] = tc.sourceFor("standard");
            M2 = tc.modelFor(M.options.set("drag.cdScale", 2));
            s2 = tc.preparedFor(M2);
            d1 = s1.diagnostics(180, tc.r(:, 181), tc.v(:, 181));
            d2 = s2.diagnostics(180, tc.r(:, 181), tc.v(:, 181));
            tc.verifyBoundAndRecord("native drag doubles with cdScale", norm(d2.drag - 2 * d1.drag) / norm(d1.drag), 1e-15, "-", "an exact power of two through a linear factor");
            tc.verifyEqual(d2.radiation, d1.radiation, 'the radiation must not notice the drag scale');
            s1.release(); s2.release();
            % 2. A propagated arc: the two bindings are one trajectory. gjToolbox promises agreement to a
            % stated band of ulps, not bitwise, because the two stepping kernels sum differently; what is
            % asserted here is that the SOURCE is the same on both sides.
            warnState = warning('off', 'gjToolbox:dynamics:fdGradient');
            restore = onCleanup(@() warning(warnState));
            opt = gjToolbox.Options('Epoch', tc.tai, 'Step', 5, 'Order', 8, 'Startup', "onesided");
            % The arc has NO gravity (this toolbox has none to give), so the state flies a straight line and
            % climbs ~80 km in 300 s -- out of the standard preset's Swarm bracket (446-513 km on the day),
            % where env.densityScale rightly refuses (CI 40884016269). The driver scale, one factor per day
            % and no height dependence, keeps a density scale on both bindings without that domain limit.
            Ma = tc.scaledModel("driver");
            srcV = NonGravForces.app.DynamicsSource(Ma, 'Epoch', tc.tai, 'Margin_s', 30, 'Native', false);
            srcN = NonGravForces.app.DynamicsSource(Ma, 'Epoch', tc.tai, 'Margin_s', 30);
            dynV = gjToolbox.Dynamics(srcV, 'Binding', "vectorized", 'ValidateGradient', false, 'FiniteDifferenceFallback', false);
            dynN = gjToolbox.Dynamics(srcN, 'Binding', "native", 'ValidateGradient', false, 'FiniteDifferenceFallback', false);
            r0 = tc.r(:, 121); v0 = tc.v(:, 121);
            solV = gjToolbox.propagate(dynV, [120 420], r0, v0, opt);
            solN = gjToolbox.propagate(dynN, [120 420], r0, v0, opt);
            tc.verifyEqual(solN.Diagnostics.ForceEvals.total, solV.Diagnostics.ForceEvals.total, 'the two bindings must cost the same evaluations');
            drift = norm(solN.Y(1:3, end) - solV.Y(1:3, end));
            tc.verifyBoundAndRecord("propagated arc: native binding vs vectorized binding over 300 s", drift, 1e-6, "m", ...
                "one model through two stepping kernels; the acceleration differs by 1e-14 relative and the arc is 300 s (measured 0 m)");
        end
    end

    methods (Test, TestTags = {'metamorphic'})
        function DynamicsSourceDensityScaleMetamorphic(tc)
            % The gas force is linear in the densities, so the native drag with the driver scale is the unscaled
            % native drag times the day's factor, and radiation and thermal do not notice it.
            M0 = tc.modelFor(NgfTestCase.standardUnscaled());
            M1 = tc.scaledModel("driver");
            s0 = tc.preparedFor(M0); s1 = tc.preparedFor(M1);
            d0 = s0.diagnostics(180, tc.r(:, 181), tc.v(:, 181));
            d1 = s1.diagnostics(180, tc.r(:, 181), tc.v(:, 181));
            tc.verifyEqual(d0.densityScale, 1, 'without drag.densityScale the factor is exactly 1');
            tc.verifyBoundAndRecord("native drag with the driver scale / (factor x unscaled native drag)", norm(d1.drag - d1.densityScale * d0.drag) / norm(d1.drag), 1e-14, "-", ...
                "a linear factor through the per-species force sums; fixed accommodation does not read the density");
            tc.verifyEqual(d1.radiation, d0.radiation, 'the radiation must not notice the density scale');
            tc.verifyEqual(d1.thermal, d0.thermal, 'the thermal recoil must not notice the density scale');
            s0.release(); s1.release();
        end
    end

    methods (Test, TestTags = {'robustness'})
        function DynamicsSourceDensityScaleRobustness(tc)
            % A Swarm-only scale on days without a factor prepares (the table marks them missing) and is REFUSED
            % when a step reaches one, as Model.accel refuses it; the combined scale's fallback serves the same day.
            Ms = tc.scaledModel("swarm", 'MinRecords', 2000);
            src = tc.preparedFor(Ms);
            tc.verifyErrorId(@() src.accel(180, tc.r(:, 181), tc.v(:, 181)), 'NonGravForces:rhsGateway:dataRange');
            src.release();
            epk = referenceSystem.daySecondToEpoch([tc.ep.day(181) tc.ep.sec(181)], tc.ep.scale);
            tc.verifyErrorId(@() Ms.accel(epk, tc.r(:, 181), tc.v(:, 181), [], 'Frame', "gcrs"), 'NonGravForces:densityScale:dataDay');
            Mc = tc.scaledModel("combined", 'MinRecords', 2000);
            src = tc.preparedFor(Mc);
            tc.verifyTrue(all(isfinite(src.accel(180, tc.r(:, 181), tc.v(:, 181)))), 'the fallback serves the day Swarm cannot');
            src.release();
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function DynamicsSourceSpeed(tc)
            % The native entry exists to be cheap in an inner loop: one evaluation against the same source's
            % MATLAB path, both timed inside one repeat on this host (rule T5, a ratio and never an absolute).
            [srcN, M] = tc.sourceFor("standardSteady");   % 'Native' false cannot serve transient walls
            srcM = NonGravForces.app.DynamicsSource(M, 'Epoch', tc.tai, 'Margin_s', 30, 'TableStep_s', 60, 'Native', false);
            srcM = srcM.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
            r1 = tc.r(:, 181); v1 = tc.v(:, 181);
            q = tc.speedRatio(@() srcN.accel(180, r1, v1), @() srcM.accel(180, r1, v1));
            tc.verifyBoundAndRecord("native accel / MATLAB accel, one state", q, 0.08, "-", ...
                "measured 0.0234 (43x) on 2026-09-16, R2026a, standard preset; budget 3x the measured ratio");
            srcN.release(); srcM.release();
        end
    end
end
