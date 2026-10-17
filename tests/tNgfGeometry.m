classdef tNgfGeometry < NgfTestCase
%TNGFGEOMETRY Tests of the LP-4 geom units: boxWing, selfShadow, thermalState, AttitudeProvider, checkAttitudeSense.
%
%   Oracles (none built here): the pure-Python ray cast tools/pyxcheck/selfshadow_check.py
%   on the GRACE-FO polygons (160 x 160 samples), an analytic two-plate occlusion, the real
%   2020-01-01 GRACE-C attitude and state (data/samples), Wöske (2019) eq. 5.17/5.18 in
%   closed form (energy balance, linearised exponential relaxation), and
%   referenceSystem's own quaternion algebra for the nominal attitude round trip.
%   Every bound carries the number measured on 2026-09-07.
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
    properties
        root
        mm          % graceFO("C")
        oracle      % selfshadow_oracle.csv
        au
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            tc.mm = NonGravForces.geom.graceFO("C");
            tc.oracle = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'selfshadow_oracle.csv'), 'VariableNamingRule', 'preserve');
            tc.au = NonGravForces.const().au;
        end
    end
    methods
        function [ep, r, v, q] = sampleState(tc)
            eop = referenceSystem.readEop(fullfile(tc.root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            [ep, r, v] = NonGravForces.io.readGNV1B(fullfile(tc.root, 'data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
            [~, q] = NonGravForces.io.readSCA1B(fullfile(tc.root, 'data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
        end
        function [d, fo] = oracleDirections(tc)
            d = [tc.oracle.dx tc.oracle.dy tc.oracle.dz]'; fo = table2array(tc.oracle(:, 5:12))';
        end
        function two = twoPlates(~, shift)
            % lower unit square at z = 0 (normal +z), upper unit square at z = 0.5 shifted by `shift` in x
            pl = repmat(struct('name',"",'kind',"plate",'area_m2',1,'normal',[0;0;1],'centroid_m',[0;0;0],'absVis',0.3,'specVis',0.3,'diffVis',0.4, ...
                'emisIR',0.8,'specIR',0.1,'diffIR',0.1,'alpha',0.85,'group',"bus",'vertices_m',zeros(3,0),'heatCapacity_J_m2K',1000), 1, 2);
            pl(1).name = "lower"; pl(1).vertices_m = [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5; 0 0 0 0];
            pl(2).name = "upper"; pl(2).centroid_m = [shift; 0; 0.5]; pl(2).vertices_m = [shift-0.5 shift+0.5 shift+0.5 shift-0.5; -0.5 -0.5 0.5 0.5; 0.5 0.5 0.5 0.5];
            two = NonGravForces.geom.MacroModel(pl, 'Name', "twoPlates");
        end
        function mm = tpmcWall(~)
            % The geometry of tools/pyxcheck/tpmc_check.py: a 1e-4 m target plate at the origin facing +z, and a two-sided
            % wall in the plane x = 0.5 spanning y in [-2, 2], z in [0, 0.5].
            h = 0.5e-4;
            pl = repmat(struct('name',"",'kind',"plate",'area_m2',1,'normal',[0;0;1],'centroid_m',[0;0;0],'absVis',0.3,'specVis',0.3,'diffVis',0.4, ...
                'emisIR',0.8,'specIR',0.1,'diffIR',0.1,'group',"bus",'vertices_m',zeros(3,0),'heatCapacity_J_m2K',1000), 1, 2);
            pl(1).name = "target"; pl(1).area_m2 = (2*h)^2; pl(1).vertices_m = [-h h h -h; -h -h h h; 0 0 0 0];
            pl(2).name = "wall"; pl(2).area_m2 = 2; pl(2).normal = [-1; 0; 0]; pl(2).centroid_m = [0.5; 0; 0.25];
            pl(2).vertices_m = [0.5 0.5 0.5 0.5; -2 -2 2 2; 0 0.5 0.5 0];
            mm = NonGravForces.geom.MacroModel(pl, 'Name', "tpmcWall");
        end
        function mm = farPlates(~)
            % Two unit plates 100 m apart in one plane: neither obstructs the other in any front-hemisphere direction.
            pl = repmat(struct('name',"",'kind',"plate",'area_m2',1,'normal',[0;0;1],'centroid_m',[0;0;0],'absVis',0.3,'specVis',0.3,'diffVis',0.4, ...
                'emisIR',0.8,'specIR',0.1,'diffIR',0.1,'group',"bus",'vertices_m',zeros(3,0),'heatCapacity_J_m2K',1000), 1, 2);
            pl(1).vertices_m = [-0.5 0.5 0.5 -0.5; -0.5 -0.5 0.5 0.5; 0 0 0 0];
            pl(2).centroid_m = [100; 0; 0]; pl(2).vertices_m = pl(1).vertices_m + [100; 0; 0];
            mm = NonGravForces.geom.MacroModel(pl, 'Name', "farPlates");
        end
        function T = temperatureFor(tc, s, U)
            % the gas temperature at which atomic oxygen (species 2) has speed ratio s at relative speed U
            m = tc.c().speciesMass(2) / tc.c().NA; T = U^2 * m / (2 * tc.c().kB * s^2);
        end
        function k = c(~)
            k = NonGravForces.const();
        end
        function src = constantSource(~, dirBody, P)
            src = struct('dir', reshape(dirBody / norm(dirBody), 3, 1, 1), 'P', P, 'band', "vis");
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
        function flowShadowTableContract(tc)
            % Method "table": the table itself (both cones), the struct it returns, and what it refuses.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            tb = NonGravForces.geom.flowShadow(tc.mm, [], [], 'Method', "table", 'Lut', lut);
            tc.verifyEqual(tb.shaded, find(~all(lut.f == 1, 2))');
            tc.verifySize(tb.plus, [24 * numel(tb.shaded), tb.n^2 * tb.ns]); tc.verifySize(tb.minus, size(tb.plus));
            vb = [-1 -0.99; 0.02 0.1; 0.01 -0.05]; vb = 7600 * vb ./ vecnorm(vb); T = [950 1000];
            [S, info] = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "table", 'Lut', lut);
            tc.verifySize(S.number, [8 2 8]); tc.verifySize(S.normal, [8 2 8]); tc.verifySize(S.tangential, [8 2 8]);
            tc.verifyClass(S.tableFallback, 'logical'); tc.verifyEqual(S.tableFallback, false(1, 2));
            tc.verifyEqual(info.method, "table"); tc.verifyEqual(info.fallbackEpochs, 0);
            Sb = NonGravForces.geom.flowShadow(tc.mm, vb, [], 'Method', "table", 'Lut', lut);
            tc.verifyEmpty(Sb.number); tc.verifyEqual(Sb.beam, S.beam); tc.verifyEqual(Sb.tableFallback, false(1, 2));
            Sl = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut);
            tc.verifyEmpty(Sl.tableFallback); tc.verifyEqual(Sl.beam, S.beam);
            g = @(varargin) NonGravForces.geom.flowShadow(varargin{:});
            tc.verifyErrorId(@() g(tc.mm, vb, T, 'Method', "table"), 'NonGravForces:flowShadow:inputLut');
            tc.verifyErrorId(@() g(tc.mm, vb, T, 'Method', "table", 'Lut', lut, 'Integrals', true), 'NonGravForces:flowShadow:unsupportedIntegrals');
            tc.verifyErrorId(@() g(tc.mm, vb(1:2, :), T, 'Method', "table", 'Lut', lut), 'NonGravForces:flowShadow:inputShape');
            tc.verifyErrorId(@() g(tc.mm, vb, [950 1000 1050], 'Method', "table", 'Lut', lut), 'NonGravForces:flowShadow:inputShape');
            tc.verifyErrorId(@() g(tc.mm, vb, [950 -1], 'Method', "table", 'Lut', lut), 'NonGravForces:flowShadow:inputRange');
        end
        function flowShadowContract(tc)
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 4, 'Samples', 10);
            vb = 7600 * [-1 -0.99; 0.02 0.1; 0.01 -0.05]; vb = 7600 * vb ./ vecnorm(vb); T = [950 1000];
            [S, info] = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut);
            tc.verifyEqual(S.kind, "flow"); tc.verifySize(S.beam, [8 2]); tc.verifySize(S.number, [8 2 8]);
            tc.verifySize(S.normal, [8 2 8]); tc.verifySize(S.tangential, [8 2 8]); tc.verifySize(S.speedRatio, [8 2]);
            tc.verifyTrue(all([S.beam(:); S.number(:); S.normal(:)] >= 0 & [S.beam(:); S.number(:); S.normal(:)] <= 1));
            tc.verifyEqual(info.method, "lut"); tc.verifyEqual(info.quadrature, [12 16]); tc.verifyEqual(info.directionsPerEpoch, 8 * 6 * 12 * 16);
            tc.verifyEqual(info.provenance.unit, "geom.flowShadow");
            Sb = NonGravForces.geom.flowShadow(tc.mm, vb, [], 'Method', "lut", 'Lut', lut);
            tc.verifyEqual(Sb.beam, S.beam); tc.verifyEmpty(Sb.number); tc.verifyEmpty(Sb.speedRatio);
            [~, infoI] = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut, 'Integrals', true);
            tc.verifySize(infoI.integrals.CpInc, [8 2 8]);
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]), 'array', struct('vis', [0.65 0.05 0.30], 'ir', [0.81 0.03 0.16]));
            bwA = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Array', struct('area_m2', 2, 'normal', [0; 0; -1], 'centroid_m', [0; 0; -1]));
            g = @(varargin) NonGravForces.geom.flowShadow(varargin{:});
            tc.verifyErrorId(@() g(bwA, vb, T), 'NonGravForces:flowShadow:inputMacro');
            tc.verifyErrorId(@() g(tc.mm, vb(1:2, :), T, 'Method', "lut", 'Lut', lut), 'NonGravForces:flowShadow:inputShape');
            tc.verifyErrorId(@() g(tc.mm, vb, [950 1000 1050], 'Method', "lut", 'Lut', lut), 'NonGravForces:flowShadow:inputShape');
            tc.verifyErrorId(@() g(tc.mm, vb, T, 'Method', "lut", 'Lut', lut, 'Quadrature', [0 16]), 'NonGravForces:flowShadow:inputShape');
            tc.verifyErrorId(@() g(tc.mm, [0 1; 0 0; 0 0], T, 'Method', "lut", 'Lut', lut), 'NonGravForces:flowShadow:inputRange');
            tc.verifyErrorId(@() g(tc.mm, vb, [950 -1], 'Method', "lut", 'Lut', lut), 'NonGravForces:flowShadow:inputRange');
            tc.verifyErrorId(@() g(tc.mm, vb, T, 'Method', "lut"), 'NonGravForces:flowShadow:inputLut');
        end
        function boxWingContract(tc)
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]), 'array', struct('vis', [0.65 0.05 0.30], 'ir', [0.81 0.03 0.16]));
            bw = NonGravForces.geom.boxWing([3 1.9 0.7], o);
            tc.verifyClass(bw, 'NonGravForces.geom.MacroModel'); tc.verifyEqual(bw.n, 6); tc.verifyTrue(bw.hasVertices);
            tc.verifyEqual(sum(bw.A), 2 * (1.9*0.7 + 3*0.7 + 3*1.9), 'AbsTol', 1e-12); tc.verifyEqual(bw.silhouetteDefect(), zeros(3, 1), 'AbsTol', 1e-12);
            arr = struct('area_m2', 2, 'normal', [0; 0; -1], 'centroid_m', [0; 0; -1], 'name', "wing");
            bwA = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Array', arr, 'Name', "Sat", 'HeatCapacity_J_m2K', [3000 6000]);
            tc.verifyEqual(bwA.n, 7); tc.verifyFalse(bwA.hasVertices); tc.verifyEqual(bwA.name, "Sat"); tc.verifyEqual(bwA.panels(7).name, "wing");
            tc.verifyEqual(bwA.heatCapacity, [3000*ones(1,6) 6000]); tc.verifyEqual(bwA.optVis(:, 7), [0.65; 0.05; 0.30]);
            tc.verifyErrorId(@() NonGravForces.geom.boxWing([3 -1 0.7], o), 'NonGravForces:boxWing:inputDims');
            tc.verifyErrorId(@() NonGravForces.geom.boxWing([3 1 0.7], struct('vis', 1)), 'NonGravForces:boxWing:inputOptics');
            tc.verifyErrorId(@() NonGravForces.geom.boxWing([3 1 0.7], struct('bus', o.bus), 'Array', arr), 'NonGravForces:boxWing:inputOptics');
            tc.verifyErrorId(@() NonGravForces.geom.boxWing([3 1 0.7], o, 'Array', struct('area_m2', 2, 'normal', [0; 0; -2], 'centroid_m', [0; 0; -1])), 'NonGravForces:boxWing:inputArray');
            tc.verifyErrorId(@() NonGravForces.geom.boxWing([3 1 0.7], struct('bus', struct('vis', [0.5 0.5 0.5], 'ir', o.bus.ir))), 'NonGravForces:MacroModel:inputEnergy');
        end
        function selfShadowContract(tc)
            [d, ~] = tc.oracleDirections();
            [f, info] = NonGravForces.geom.selfShadow(tc.mm, d);
            tc.verifySize(f, [8 size(d, 2)]); tc.verifyTrue(all(f(:) >= 0 & f(:) <= 1)); tc.verifyEqual(info.method, "raycast"); tc.verifyEqual(info.samples, 40);
            % OCCLUSION ONLY (LP-C, F1): a plate facing away keeps a geometric value, it is not forced to 0
            lit = (tc.mm.N' * d) > 0; tc.verifyTrue(all(f(lit) > 0), 'lit plates are at least partly exposed');
            tc.verifyTrue(any(f(~lit) > 0), 'back-facing plates carry their mirror value, not a forced 0');
            [fl, infoL] = NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'LutResolution', 4);
            tc.verifySize(fl, [8 size(d, 2)]); tc.verifyEqual(infoL.lut.resolution, 4); tc.verifySize(infoL.lut.f, [8 6*25]); tc.verifySize(infoL.lut.nodes, [3 6*25]);
            tab = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 4); tc.verifyEqual(tab.f, infoL.lut.f);
            [~, infoNoLut] = NonGravForces.geom.selfShadow(tc.mm, d); tc.verifyEmpty(infoNoLut.lut); tc.verifyEqual(infoNoLut.provenance.unit, "geom.selfShadow");
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]), 'array', struct('vis', [0.65 0.05 0.30], 'ir', [0.81 0.03 0.16]));
            bwA = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Array', struct('area_m2', 2, 'normal', [0; 0; -1], 'centroid_m', [0; 0; -1]));
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(bwA, d), 'NonGravForces:selfShadow:inputMacro');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, [1; 1; 0]), 'NonGravForces:selfShadow:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, []), 'NonGravForces:selfShadow:inputShape');
            for bad = {2.5, 0, 300, [4 4]}
                tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'LutResolution', bad{1}), 'NonGravForces:selfShadow:inputResolution');
            end
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'LutResolution', 5, 'Lut', infoL.lut), 'NonGravForces:selfShadow:inputLut');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.twoPlates(0.5), [0; 0; 1], 'Method', "lut", 'Lut', infoL.lut), 'NonGravForces:selfShadow:inputLut');
        end
        function thermalStateContract(tc)
            [ep, r, ~, q] = tc.sampleState(); ep5 = referenceSystem.daySecondToEpoch([ep.day(1:5) ep.sec(1:5)], "GPS");
            src = NonGravForces.env.radiationSources(ep5, r(:, 1:5), q(:, 1:5), NonGravForces.Options("standard"));
            [T, info] = NonGravForces.geom.thermalState(tc.mm, src);
            tc.verifySize(T, [8 5]); tc.verifyTrue(all(T(:) >= 0 & T(:) < 400)); tc.verifyEqual(info.model, "steadyState");
            tc.verifySize(info.absorbed_W_m2, [8 5]); tc.verifySize(info.emitted_W_m2, [8 5]); tc.verifyEqual(info.provenance.unit, "geom.thermalState");
            [Tt, it] = NonGravForces.geom.thermalState(tc.mm, src, 'Model', "transient", 'Epoch', ep5, 'T0_K', T(:, 1)');
            tc.verifySize(Tt, [8 5]); tc.verifyEqual(Tt(:, 1), T(:, 1)); tc.verifySize(it.substeps, [1 4]); tc.verifySize(it.timeConstant_s, [8 5]);
            share = zeros(1, 8); share(7) = 0.9; share([1 2 4 6]) = 0.025;
            [Ti, ii] = NonGravForces.geom.thermalState(tc.mm, src, 'InternalPower_W', 250, 'InternalShare', share);
            tc.verifyTrue(all(Ti(7, :) > T(7, :))); tc.verifyEqual(sum(ii.internal_W_m2(:, 1) .* tc.mm.A'), 250, 'AbsTol', 1e-9);
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(struct(), src), 'NonGravForces:thermalState:inputMacro');
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(tc.mm, struct('dir', 1)), 'NonGravForces:thermalState:inputSources');
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(tc.mm, src, 'InternalPower_W', 100), 'NonGravForces:thermalState:inputInternal');
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(tc.mm, src, 'Exposed', ones(8, 4)), 'NonGravForces:thermalState:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(tc.mm, src, 'Model', "transient", 'Epoch', ep5), 'NonGravForces:thermalState:inputShape');
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(NonGravForces.geom.boxWing([3 1.9 0.7], o), src, 'Model', "transient", 'Epoch', ep5, 'T0_K', 300*ones(1, 6)), 'NonGravForces:thermalState:inputMacro');
        end
        function attitudeProviderContract(tc)
            [ep, r, v, q] = tc.sampleState();
            [epS] = ep;   % SCA1B and GNV1B share the 1 Hz epochs of the sample
            ap = NonGravForces.geom.AttitudeProvider(epS, q);
            tc.verifyEqual(ap.nodeCount, 600); tc.verifyEqual(ap.method, "slerp"); tc.verifyEqual(ap.scale, "GPS"); tc.verifyEqual(ap.convention.quaternion, "scalarFirst");
            tc.verifyEqual(ap.span, [58849, 58849 + 599/86400]);
            epQ = referenceSystem.daySecondToEpoch([58849*ones(3, 1) [0; 10.5; 599]], "GPS");
            [qi, info] = ap.at(epQ); tc.verifySize(qi, [4 3]); tc.verifyEqual(qi(:, [1 3]), q(:, [1 600])); tc.verifyFalse(info.converted);
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider(epS, q(1:3, :)), 'NonGravForces:AttitudeProvider:inputNodes');
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider(struct('x', 1), q), 'NonGravForces:AttitudeProvider:inputNodes');
            tc.verifyErrorId(@() ap.at(struct('x', 1)), 'NonGravForces:AttitudeProvider:inputEpoch');
            tc.verifyError(@() ap.at(referenceSystem.daySecondToEpoch([58849 700], "GPS")), 'referenceSystem:interpolateQuaternions:queryOutOfRange');
            tc.verifyError(@() NonGravForces.geom.AttitudeProvider(epS, q, 'Method', "linear"), 'MATLAB:validators:mustBeMember');
            qn = NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:4), v(:, 1:4), false); tc.verifySize(qn, [4 4]); tc.verifyLessThan(max(abs(vecnorm(qn) - 1)), 1e-14);
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:4), v(:, 1:3), false), 'NonGravForces:AttitudeProvider:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider.nominal([7e6; 0; 0], [7e3; 0; 0], true), 'NonGravForces:AttitudeProvider:inputShape');
        end
        function checkAttitudeSenseContract(tc)
            [~, r, v, q] = tc.sampleState();
            [ok, rep] = NonGravForces.geom.checkAttitudeSense(q, r, v);
            tc.verifyTrue(ok); tc.verifySize(rep.nadirAngle_deg, [1 600]); tc.verifySize(rep.nadirAngleConjugate_deg, [1 600]); tc.verifySize(rep.alongTrackAngle_deg, [1 600]);
            tc.verifyEqual(rep.xSense, "againstVelocity"); tc.verifyEqual(rep.K, 600); tc.verifyEqual(rep.tolerance_deg, 5); tc.verifyEqual(rep.provenance.unit, "geom.checkAttitudeSense");
            [okc, repc] = NonGravForces.geom.checkAttitudeSense(referenceSystem.invertQuaternion(q), r, v);
            tc.verifyFalse(okc, 'the conjugated series must fail'); tc.verifyGreaterThan(repc.maxNadirAngle_deg, 30);
            tc.verifyErrorId(@() NonGravForces.geom.checkAttitudeSense(q, r(:, 1:10), v), 'NonGravForces:checkAttitudeSense:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.checkAttitudeSense([q(:, 1:599) [NaN; 0; 0; 0]], r, v), 'NonGravForces:checkAttitudeSense:inputShape');
            tc.verifyError(@() NonGravForces.geom.checkAttitudeSense(q, r, v, 'NadirTolerance_deg', -1), 'MATLAB:validators:mustBePositive');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function flowShadowTableReference(tc)
            % The table against the computation it tabulates (Method lut), over GRACE-FO's stream: yaw +-4 deg, pitch
            % +-1 deg about both +x and -x, s_O 7..10. Measured 2026-10-17 on 81952 epochs of ten real days (C and D):
            % 1.2e-3 at most, which moved the drag by 6.7e-13 m/s^2 of 2.1e-8. Bound 4x.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            rng(11); K = 400;
            sg = [-ones(1, K / 2), ones(1, K / 2)];
            u = [sg; tand(4 * (2 * rand(1, K) - 1)); tand(2 * rand(1, K) - 1)]; vb = 7600 * u ./ vecnorm(u);
            T = arrayfun(@(s) tc.temperatureFor(s, 7600), 7 + 3 * rand(1, K));
            St = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "table", 'Lut', lut);
            Se = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut);
            tc.verifyFalse(any(St.tableFallback));
            err = max(abs([St.number(:) - Se.number(:); St.normal(:) - Se.normal(:); St.tangential(:) - Se.tangential(:)]));
            tc.verifyBoundAndRecord("flowShadow table vs Method lut, every ratio, GRACE-FO cone (400 epochs)", err, 5e-3, "-", ...
                "trilinear at 0.5 deg and 0.5 in s_O; measured 1.2e-3 on 81952 real epochs (2026-10-17)");
            tc.verifyEqual(St.beam, Se.beam);
        end
        function flowShadowReference(tc)
            % Oracle: tools/pyxcheck/tpmc_check.py, a standard-library test-particle Monte Carlo of the incident number,
            % normal-momentum and tangential-momentum fluxes on a tiny plate screened by a wall, written from a physics
            % statement by an agent that did not open this toolbox (T9). The model is the flow ratios times Sentman's
            % incident closed forms. Measured 2026-09-15 at quadrature [96 128] with an exact ray cast: worst |z| = 2.58
            % over the components the Monte Carlo resolves (standard error below 10 % of the value). The occlusion is a
            % step in direction, which the quadrature resolves to its node spacing: at [24 32] the same cases missed by up
            % to 20 standard errors and converged as the nodes were refined, so the fine quadrature is what isolates the
            % physics from the discretisation. Bound |z| <= 4. Where the Monte Carlo saw too few unoccluded molecules to
            % measure a value, the model must lie below the oracle's detection bound instead.
            tp = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'tpmc_oracle.csv'));
            mm = tc.tpmcWall(); U = 7500; worst = 0;
            for i = 1:height(tp)
                vh = [tp.vhat_x(i); tp.vhat_y(i); tp.vhat_z(i)];
                S = NonGravForces.geom.flowShadow(mm, U * vh, tc.temperatureFor(tp.s(i), U), 'Samples', 2, 'Quadrature', [96 128]);
                g = vh(3); l = sqrt(max(0, 1 - g^2)); s = S.speedRatio(2); P = exp(-(g*s)^2)/s; Z = erfc(-g*s);
                model = [S.number(1, 1, 2) * (P/(2*sqrt(pi)) + g*Z/2), S.normal(1, 1, 2) * (g*P/sqrt(pi) + (g^2 + 1/(2*s^2))*Z), ...
                         S.tangential(1, 1, 2) * l * (P/sqrt(pi) + g*Z)];
                ref = [tp.CN_occ(i) tp.CP_occ(i) tp.CT_occ(i)]; se = [tp.CN_occ_se(i) tp.CP_occ_se(i) tp.CT_occ_se(i)];
                resolvable = se > 0 & se < 0.1 * abs(ref);
                if tp.elevation_deg(i) == 90, resolvable(3) = false; end        % no tangential stream: the oracle's t is a fallback
                if any(resolvable)
                    worst = max(worst, max(abs(model(resolvable) - ref(resolvable)) ./ se(resolvable)));
                end
                unresolved = ~resolvable & [true true tp.elevation_deg(i) ~= 90];
                floorBound = ref + 3 * max(se, tp.CN_free(i) / tp.N(i));
                tc.verifyTrue(all(model(unresolved) <= floorBound(unresolved)), sprintf('case %d: below the oracle''s detection bound', i));
            end
            tc.verifyBoundAndRecord("flowShadow vs test-particle Monte Carlo, worst |z| (resolved components)", worst, 4, "-", "measured 2.58 on 2026-09-15; 18 cases, s = 2, 5, 8, elevation -10 to 90 deg");
        end
        function boxWingReference(tc)
            % The bus of a 1 m cube is the unit cube: six unit areas, face centres at +-1/2, closed silhouette.
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            bw = NonGravForces.geom.boxWing([1 1 1], o);
            tc.verifyEqual(bw.A, ones(1, 6)); tc.verifyEqual(sort(abs(bw.C(:))), [zeros(12, 1); 0.5*ones(6, 1)]);
            tc.verifyEqual(bw.N, [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1]);
        end
        function selfShadowReference(tc)
            % Oracle: the pure-Python ray cast (160 x 160 per plate). Measured 2026-09-07: 40 samples 1.25e-2, 100 samples 3.8e-3;
            % the sampling error bound is ~2/m at a shadow edge: 5e-2 and 2e-2. Bounds 3e-2 and 1e-2.
            % The oracle writes 0 for a plate facing away (the pre-LP-C convention); occlusion is compared where the
            % plate faces the direction, which is the only place any consumer weights it.
            [d, fo] = tc.oracleDirections();
            lit = (tc.mm.N' * d) > 0;
            f40 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 40); f100 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 100);
            tc.verifyBoundAndRecord("selfShadow raycast 40 vs Python 160 (lit plates)", max(abs(f40(lit) - fo(lit))), 3e-2, "-", "measured 1.25e-2");
            tc.verifyBoundAndRecord("selfShadow raycast 100 vs Python 160 (lit plates)", max(abs(f100(lit) - fo(lit))), 1e-2, "-", "measured 3.75e-3");
            % which plate is shadowed and by how much is the physics: only the nadir plate, at low elevation, by the overhangs
            shadowed = any(f100 < 0.999 & lit, 2);
            tc.verifyEqual(find(shadowed)', 7, 'only the Nadir plate is partly shadowed on GRACE-FO');
            tc.verifyLessThan(f100(7, 3), 0.75, 'at 85 deg from nadir the overhang shadows about 30 % of the nadir plate');
        end
        function thermalStateReference(tc)
            % Oracle: Stefan-Boltzmann in closed form. One plate normal to the Sun at 1 au with alpha 0.34, eps 0.62:
            % T = (0.34 * 1361 / (0.62 sigma))^(1/4).
            sigma = NonGravForces.const().sigmaSB; c = referenceSystem.iersConstants().speedOfLight;
            src = tc.constantSource([1; 0; 0], 1361 / c);
            T = NonGravForces.geom.thermalState(tc.mm, src);
            tc.verifyBoundAndRecord("steady-state temperature vs closed form", abs(T(1) - (0.34 * 1361 / (0.62 * sigma))^(1/4)), 1e-9, "K", "one fourth root");
            tc.verifyEqual(T(2:end), zeros(7, 1), 'unlit plates without internal power sit at 0 K by construction');
        end
        function attitudeProviderReference(tc)
            % Oracle: the real 2020-01-01 attitude. The nominal Earth-pointing frame with +X against the velocity is
            % within 1 deg of the measured SCA1B attitude (measured 0.84 deg mean, 0.86 max); with +X along it, 180 deg.
            [~, r, v, q] = tc.sampleState();
            Rm = referenceSystem.rotationFromQuaternion(q);
            for sense = [false true]
                Rn = referenceSystem.rotationFromQuaternion(NonGravForces.geom.AttitudeProvider.nominal(r, v, sense));
                dR = pagemtimes(Rn, 'transpose', Rm, 'none');
                ang = acosd(max(-1, min(1, (dR(1,1,:) + dR(2,2,:) + dR(3,3,:) - 1) / 2)));
                if ~sense
                    tc.verifyBoundAndRecord("nominal attitude (+X against v) vs SCA1B", max(ang), 3, "deg", "pointing offsets ~1 deg; measured 0.86");
                else
                    tc.verifyGreaterThan(min(ang), 170, 'the other sense is a half turn off');
                end
            end
        end
        function checkAttitudeSenseReference(tc)
            % The measured GRACE-C attitude: nadir within 1 deg, flying with +X against the velocity (JPL D-56935: GRACE-C
            % trails GRACE-D with its K-band horn towards it on 2020-01-01, measured in tests/tNgfIo.m as 179.2 deg).
            [~, r, v, q] = tc.sampleState();
            [ok, rep] = NonGravForces.geom.checkAttitudeSense(q, r, v);
            tc.verifyTrue(ok);
            tc.verifyBoundAndRecord("checkAttitudeSense max nadir angle, real day", rep.maxNadirAngle_deg, 3, "deg", "measured 0.86");
            tc.verifyGreaterThan(min(rep.nadirAngleConjugate_deg), 30);
            tc.verifyGreaterThan(mean(rep.alongTrackAngle_deg), 170);
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function thermalStateOcclusionPrecision(tc)
            % F2 (LP-C): the absorbed flux with per-source occlusion is the flux radiationPanelForce forms from the same
            % sources and the same occlusion, source by source; and a panel facing away from the Sun keeps its Earth heating.
            [ep, r, ~, q] = tc.sampleState(); ep5 = referenceSystem.daySecondToEpoch([ep.day(1:5) ep.sec(1:5)], "GPS");
            src = NonGravForces.env.radiationSources(ep5, r(:, 1:5), q(:, 1:5), NonGravForces.Options("standard"));
            [~, K] = size(src.P); M = size(src.dir, 2);
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 8, 'Samples', 20);
            X = reshape(NonGravForces.geom.selfShadow(tc.mm, reshape(src.dir, 3, []), 'Method', "lut", 'Lut', lut), tc.mm.n, M, K);
            X = permute(X, [1 3 2]);                                                     % N x K x M
            [~, it] = NonGravForces.geom.thermalState(tc.mm, src, 'Exposed', X);
            [~, ir] = NonGravForces.kernel.radiationPanelForce(tc.mm, src, 'Shadow', X);
            tc.verifyBoundAndRecord("thermalState absorbed flux vs radiationPanelForce's, per-source occlusion", ...
                max(abs(it.absorbed_W_m2 - ir.absorbed_W_m2), [], 'all') / max(ir.absorbed_W_m2, [], 'all'), 1e-13, "-", ...
                "the same product a P c gamma o summed over 577 sources in two association orders");
            % the nadir radiator (plate 7) sees the Earth on every epoch: with the Sun's occlusion alone applied to every
            % source (the pre-LP-C composition) its Earth flux would follow the Sun's side; per source it never vanishes
            tc.verifyTrue(all(it.absorbed_W_m2(7, :) > 50), 'the nadir plate absorbs Earth radiation at every epoch');
            sunOnly = repmat(X(:, :, 1), 1, 1, M);
            [~, itOld] = NonGravForces.geom.thermalState(tc.mm, src, 'Exposed', sunOnly);
            tc.verifyBoundAndRecord("nadir absorbed flux: per-source minus Sun-only occlusion (F2 regression)", ...
                -min(it.absorbed_W_m2(7, :) - itOld.absorbed_W_m2(7, :)), 0, "W/m^2", "per-source occlusion never absorbs less than the Sun column forced onto the Earth tiles");
        end
        function flowShadowTablePrecision(tc)
            % At its own nodes the table returns what it stored, and what it stored is Method lut at that node: a node
            % read back differs from the exact evaluation there only by the rounding of the node's coordinates.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            tb = NonGravForces.geom.flowShadow(tc.mm, [], [], 'Method', "table", 'Lut', lut);
            % interior nodes: a node ON the table's edge comes back from v/|v| an ulp outside it and is, correctly,
            % computed exactly instead (flowShadowTableRobustness)
            ix = [1 12 23 3]; iy = [1 12 23 20]; is = [1 7 13 2]; sg = [1 -1 1 -1];
            u = [sg; tb.t0 + ix * tb.h; tb.t0 + iy * tb.h]; u = u ./ vecnorm(u);
            T = arrayfun(@(s) tc.temperatureFor(s, 7600), tb.s0 + is * tb.hs);
            St = NonGravForces.geom.flowShadow(tc.mm, 7600 * u, T, 'Method', "table", 'Lut', lut);
            Se = NonGravForces.geom.flowShadow(tc.mm, 7600 * u, T, 'Method', "lut", 'Lut', lut);
            tc.verifyFalse(any(St.tableFallback));
            err = max(abs([St.number(:) - Se.number(:); St.normal(:) - Se.normal(:); St.tangential(:) - Se.tangential(:)]));
            tc.verifyBoundAndRecord("flowShadow table at its nodes vs Method lut there", err, 1e-10, "-", ...
                "node coordinates rebuilt from v/|v| carry a few ulp, times a slope below 1e2 per unit of tan");
        end
        function flowShadowPrecision(tc)
            % Without occlusion the directional integrals ARE Sentman's incident coefficients (number flux, normal and
            % tangential momentum): an identity of kinetic theory, so the comparison measures the quadrature alone.
            % Measured 2026-09-15 at the default [12 16]: worst 1.6e-6 in coefficient units over s = 1..15 and incidence
            % 0..180 deg (a ram plate's coefficient is ~2). Bound 5e-6. On two plates that cannot see each other every
            % ratio is exactly 1.
            mm = tc.farPlates(); U = 7500; worst = 0;
            lut = NonGravForces.geom.selfShadow(mm, [], 'Method', "lut", 'LutResolution', 2, 'Samples', 4);
            for s = [1 2 5 8 15]
                for thdeg = [0 30 60 85 90 95 120 150 180]
                    th = deg2rad(thdeg); vb = U * [sin(th); 0; cos(th)];
                    [S, info] = NonGravForces.geom.flowShadow(mm, vb, tc.temperatureFor(s, U), 'Method', "lut", 'Lut', lut, 'Integrals', true);
                    g = cos(th); l = sin(th); so = S.speedRatio(2); P = exp(-(g*so)^2)/so; Z = erfc(-g*so);
                    ref = [P/(2*sqrt(pi)) + g*Z/2, g*P/sqrt(pi) + (g^2 + 1/(2*so^2))*Z, l*(P/sqrt(pi) + g*Z)];
                    got = [info.integrals.CN(1, 1, 2), info.integrals.CpInc(1, 1, 2), info.integrals.CtauInc(1, 1, 2)];
                    worst = max(worst, max(abs(got - ref)));
                    tc.verifyTrue(all(S.number(:) == 1 & S.normal(:) == 1 & S.tangential(:) == 1), 'nothing obstructs: every ratio exactly 1');
                end
            end
            tc.verifyBoundAndRecord("flowShadow integrals vs Sentman incident closed forms", worst, 5e-6, "-", "measured 1.6e-6 on 2026-09-15 at [12 16]; [8 16] gives 3e-4");
        end
        function boxWingPrecision(tc)
            % Polygon corners reproduce the face areas exactly (products of two edge lengths, bitwise).
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            bw = NonGravForces.geom.boxWing([3.1225 1.944 0.721], o);
            for k = 1:6
                V = bw.V{k}; a = norm(cross(V(:, 2) - V(:, 1), V(:, 4) - V(:, 1)));
                tc.verifyEqual(a, bw.A(k), 'AbsTol', 4 * eps(bw.A(k)));
            end
        end
        function selfShadowPrecision(tc)
            % Analytic two-plate case: a unit square half covered by a parallel unit square shifted by 1/2 has exposed
            % fraction exactly 1/2 for a beam along +z; with an even sample count the grid straddles the edge exactly.
            two = tc.twoPlates(0.5);
            f = NonGravForces.geom.selfShadow(two, [0; 0; 1], 'Samples', 40);
            tc.verifyBoundAndRecord("selfShadow two plates half covered", abs(f(1) - 0.5), 1e-12, "-", "edge on a grid line: exact");
            tc.verifyEqual(f(2), 1);
            f3 = NonGravForces.geom.selfShadow(two, [0; 0; 1], 'Samples', 39);
            tc.verifyBoundAndRecord("selfShadow two plates, odd grid", abs(f3(1) - 0.5), 1/39 + 1e-12, "-", "one straddling row");
        end
        function thermalStatePrecision(tc)
            % Transient against the linearised relaxation of a 1 K excursion on the Front plate (constant source):
            % T - Ts = exp(-t/tau), tau = C/(4 eps sigma Ts^3). The nonlinearity of T^4 for 1 K at 300 K is ~1 %;
            % RK4 with dt <= tau/10 adds < 1e-5. Bound 5e-2 K (measured see last_run). Untouched plates stay exactly.
            c = referenceSystem.iersConstants().speedOfLight; sigma = NonGravForces.const().sigmaSB;
            src1 = tc.constantSource([1; 0; 0], 1361 / c);
            K = 400; epK = referenceSystem.daySecondToEpoch([58849*ones(K, 1), (0:K-1)' * 5], "GPS");
            srcK = struct('dir', repmat(src1.dir, 1, 1, K), 'P', repmat(src1.P, 1, K), 'band', src1.band);
            Ts = NonGravForces.geom.thermalState(tc.mm, src1); T0 = Ts; T0(1) = Ts(1) + 1;
            [Tt, it] = NonGravForces.geom.thermalState(tc.mm, srcK, 'Model', "transient", 'Epoch', epK, 'T0_K', T0');
            tau = tc.mm.heatCapacity(1) / (4 * tc.mm.optIR(1, 1) * sigma * Ts(1)^3);
            Tlin = Ts(1) + exp(-(0:K-1) * 5 / tau);
            tc.verifyBoundAndRecord("transient vs linearised relaxation (1 K)", max(abs(Tt(1, :) - Tlin)), 5e-3, "K", "T^4 nonlinearity of a 1 K step; measured 1.1e-3");
            tc.verifyBoundAndRecord("transient: steady plates stay put", max(abs(Tt(2:end, :) - Ts(2:end)), [], 'all'), 1e-12, "K", "fixed point of RK4");
            tc.verifyTrue(all(it.substeps >= 1));
        end
        function attitudeProviderPrecision(tc)
            % A slerp query exactly at a node returns the node bitwise (referenceSystem reproduces nodes to 1e-16 and the
            % node value is the input), and the nominal attitude round-trips through its rotation matrix to 1e-15.
            [ep, r, v, q] = tc.sampleState();
            ap = NonGravForces.geom.AttitudeProvider(ep, q);
            idx = [1 7 300 600];
            qi = ap.at(referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], "GPS"));
            tc.verifyBoundAndRecord("AttitudeProvider at the nodes", max(abs(qi - q(:, idx)), [], 'all'), 2e-16, "-", "renormalisation of a unit quaternion");
            qn = NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:50), v(:, 1:50), false);
            Rn = referenceSystem.rotationFromQuaternion(qn);
            zb = reshape(pagemtimes(Rn, [0; 0; 1]), 3, 50);
            tc.verifyBoundAndRecord("nominal attitude +Z is nadir", max(abs(zb + r(:, 1:50) ./ vecnorm(r(:, 1:50))), [], 'all'), 1e-14, "-", "rotation round trip");
        end
        function checkAttitudeSensePrecision(tc)
            % The nominal frame points +Z exactly at nadir: the angle is rounding only. With atan2 (not acos) the
            % rounding is 1e-14 deg (acos would give sqrt(eps) ~ 1e-6 deg, measured 1.2e-6 before the change).
            [~, r, v] = tc.sampleState();
            qn = NonGravForces.geom.AttitudeProvider.nominal(r, v, true);
            [ok, rep] = NonGravForces.geom.checkAttitudeSense(qn, r, v);
            tc.verifyTrue(ok);
            tc.verifyBoundAndRecord("checkAttitudeSense nadir angle of the nominal frame", rep.maxNadirAngle_deg, 1e-9, "deg", "atan2 of cross and dot: rounding only");
            tc.verifyBoundAndRecord("checkAttitudeSense along-track angle of the nominal frame", max(rep.alongTrackAngle_deg), 1, "deg", "flight-path angle of the near-circular orbit (< 0.2 deg)");
            tc.verifyEqual(rep.xSense, "alongVelocity");
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function flowShadowTableVectorization(tc)
            % Element by element, inside and outside the table alike: a batch equals its epochs one at a time, bitwise.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            u = [-1 -1 1 -1 0.2 -1; 0.03 -0.08 0.05 0.4 1 0; -0.01 0.02 0 0.1 0 0.004]; vb = 7600 * u ./ vecnorm(u);
            T = [900 1000 800 950 1100 4000];
            S = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "table", 'Lut', lut);
            for k = 1:size(vb, 2)
                Sk = NonGravForces.geom.flowShadow(tc.mm, vb(:, k), T(k), 'Method', "table", 'Lut', lut);
                tc.verifyEqual(Sk.number, S.number(:, k, :)); tc.verifyEqual(Sk.normal, S.normal(:, k, :));
                tc.verifyEqual(Sk.tangential, S.tangential(:, k, :)); tc.verifyEqual(Sk.tableFallback, S.tableFallback(k));
            end
        end
        function flowShadowVectorization(tc)
            % five epochs at once equal five single epochs, bitwise: every sum runs along a dimension whose length does
            % not depend on the number of epochs
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 4, 'Samples', 10);
            rng(11); vb = [-7600 + 50 * randn(1, 5); 300 * randn(1, 5); 200 * randn(1, 5)]; T = 900 + 100 * rand(1, 5);
            S = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut);
            for k = 1:5
                S1 = NonGravForces.geom.flowShadow(tc.mm, vb(:, k), T(k), 'Method', "lut", 'Lut', lut);
                tc.verifyEqual(S.number(:, k, :), S1.number); tc.verifyEqual(S.normal(:, k, :), S1.normal);
                tc.verifyEqual(S.tangential(:, k, :), S1.tangential); tc.verifyEqual(S.beam(:, k), S1.beam);
            end
        end
        function selfShadowVectorization(tc)
            [d, ~] = tc.oracleDirections();
            fAll = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 30);
            for k = 1:size(d, 2)
                tc.verifyEqual(fAll(:, k), NonGravForces.geom.selfShadow(tc.mm, d(:, k), 'Samples', 30));
            end
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 4, 'Samples', 20);
            flAll = NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'Lut', lut);
            for k = 1:size(d, 2)
                tc.verifyEqual(flAll(:, k), NonGravForces.geom.selfShadow(tc.mm, d(:, k), 'Method', "lut", 'Lut', lut));
            end
        end
        function thermalStateVectorization(tc)
            [ep, r, ~, q] = tc.sampleState(); idx = [1 50 200 600];
            epI = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], "GPS");
            src = NonGravForces.env.radiationSources(epI, r(:, idx), q(:, idx), NonGravForces.Options("standard"));
            T = NonGravForces.geom.thermalState(tc.mm, src);
            for k = 1:4
                s1 = struct('dir', src.dir(:, :, k), 'P', src.P(:, k), 'band', src.band);
                tc.verifyEqual(T(:, k), NonGravForces.geom.thermalState(tc.mm, s1));
            end
        end
        function attitudeProviderVectorization(tc)
            [ep, r, v, q] = tc.sampleState();
            ap = NonGravForces.geom.AttitudeProvider(ep, q);
            secs = [0.25 12.5 300.75 598.5];
            qAll = ap.at(referenceSystem.daySecondToEpoch([58849*ones(4, 1) secs'], "GPS"));
            for k = 1:4
                tc.verifyEqual(qAll(:, k), ap.at(referenceSystem.daySecondToEpoch([58849 secs(k)], "GPS")));
            end
            qn = NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:5), v(:, 1:5), false);
            for k = 1:5
                tc.verifyEqual(qn(:, k), NonGravForces.geom.AttitudeProvider.nominal(r(:, k), v(:, k), false));
            end
        end
        function checkAttitudeSenseVectorization(tc)
            [~, r, v, q] = tc.sampleState(); idx = [1 100 600];
            [~, rep] = NonGravForces.geom.checkAttitudeSense(q(:, idx), r(:, idx), v(:, idx));
            for k = 1:3
                [~, r1] = NonGravForces.geom.checkAttitudeSense(q(:, idx(k)), r(:, idx(k)), v(:, idx(k)));
                tc.verifyEqual(rep.nadirAngle_deg(k), r1.nadirAngle_deg); tc.verifyEqual(rep.alongTrackAngle_deg(k), r1.alongTrackAngle_deg);
            end
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function flowShadowTableRobustness(tc)
            % Outside the table -- a stream 20 deg off the axis, side-on, or a speed ratio beyond [5, 12] -- the epoch is
            % computed exactly, bitwise as Method lut computes it, and flagged; the epochs inside are not disturbed.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            u = [-1 -1 0.2 -1 1; 0.02 tand(20) 1 0 -0.01; 0 0 0 0 0.01]; vb = 7600 * u ./ vecnorm(u);
            T = [950 950 950 tc.temperatureFor(4, 7600) tc.temperatureFor(13, 7600)];
            S = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "table", 'Lut', lut);
            tc.verifyEqual(S.tableFallback, [false true true true true]);
            out = S.tableFallback;
            E = NonGravForces.geom.flowShadow(tc.mm, vb(:, out), T(out), 'Method', "lut", 'Lut', lut);
            tc.verifyEqual(S.number(:, out, :), E.number); tc.verifyEqual(S.normal(:, out, :), E.normal);
            tc.verifyEqual(S.tangential(:, out, :), E.tangential);
            S1 = NonGravForces.geom.flowShadow(tc.mm, vb(:, 1), T(1), 'Method', "table", 'Lut', lut);
            tc.verifyEqual(S1.number, S.number(:, 1, :));
            tc.verifyTrue(all(S.number(:) >= 0 & S.number(:) <= 1 & S.normal(:) >= 0 & S.normal(:) <= 1));
        end
        function flowShadowRobustness(tc)
            % the stream exactly along a plate normal (no tangential direction), exactly edge-on, exactly leeward, a very cold
            % and a very hot gas, and no epochs at all: finite ratios in [0, 1], and exactly 1 on GRACE-FO plates nothing screens
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 4, 'Samples', 10);
            vb = 7600 * [1 0 0 -1 0; 0 0 1 0 0; 0 1 0 0 -1]; T = [900 900 50 1e5 900];
            S = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut);
            v = [S.number(:); S.normal(:)];
            tc.verifyTrue(all(isfinite(v) & v >= 0 & v <= 1), 'number and normal ratios have non-negative weights: [0, 1]');
            tc.verifyTrue(all(isfinite(S.tangential(:))), 'the tangential ratio is signed and finite (above 1 when counter-flow is screened)');
            clear1 = find(all(lut.f == 1, 2))';
            tc.verifyNotEmpty(clear1);
            tc.verifyTrue(all(S.number(clear1, :, :) == 1, 'all') && all(S.tangential(clear1, :, :) == 1, 'all'), 'a plate whose table row is 1 everywhere is exactly 1');
            S0 = NonGravForces.geom.flowShadow(tc.mm, zeros(3, 0), zeros(1, 0), 'Method', "lut", 'Lut', lut);
            tc.verifySize(S0.number, [8 0 8]); tc.verifySize(S0.beam, [8 0]);
        end
        function boxWingRobustness(tc)
            % A flat box (thin plate), a tiny box, an array without a name, and a polygon check on a hand-built model.
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]), 'array', struct('vis', [0.65 0.05 0.30], 'ir', [0.81 0.03 0.16]));
            thin = NonGravForces.geom.boxWing([2 1 0.001], o); tc.verifyEqual(thin.n, 6); tc.verifyTrue(thin.hasVertices);
            tiny = NonGravForces.geom.boxWing([1e-3 1e-3 1e-3], o); tc.verifyEqual(sum(tiny.A), 6e-6, 'AbsTol', 1e-20);
            bwA = NonGravForces.geom.boxWing([1 1 1], o, 'Array', struct('area_m2', 1, 'normal', [1; 0; 0], 'centroid_m', [2; 0; 0])); tc.verifyEqual(bwA.panels(7).name, "array1");
            % MacroModel refuses a clockwise polygon, a non-planar one, an area mismatch and a stray centroid
            p = thin.panels; p(1).vertices_m = p(1).vertices_m(:, [1 4 3 2]);
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputVertices');
            p = thin.panels; p(1).vertices_m(1, 1) = p(1).vertices_m(1, 1) + 0.01;
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputVertices');
            p = thin.panels; p(1).area_m2 = 2 * p(1).area_m2;
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputVertices');
            p = thin.panels; p(1).centroid_m = p(1).centroid_m + [0; 0.2; 0];
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputVertices');
            p = thin.panels; p(1).heatCapacity_J_m2K = -1;
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputPanels');
        end
        function selfShadowRobustness(tc)
            % Exactly edge-on beams, a beam along the body axis, and the two-plate model with no overlap and with full overlap.
            % Along +x nothing obstructs any plate: the Front plate faces it, the Rear plate's mirror direction is -x, and
            % the six edge-on plates take their front-side limit (1e-6 rad out of plane), which clears the Front plate's edge.
            fx = NonGravForces.geom.selfShadow(tc.mm, [1; 0; 0]); tc.verifyEqual(fx, ones(8, 1), 'along the body axis nothing is occluded');
            % Along +y the nadir plate's front-side limit runs just below its own plane into the starboard overhang, and
            % the port inner flank looks straight across under the body at the starboard one: both fully occluded.
            f = NonGravForces.geom.selfShadow(tc.mm, [0; 1; 0]); tc.verifyTrue(all(f(:) >= 0 & f(:) <= 1));
            tc.verifyEqual(f([7 6]), [0; 0], 'grazing +y: the overhang occludes the nadir plate and the far inner flank');
            noOverlap = tc.twoPlates(3); f0 = NonGravForces.geom.selfShadow(noOverlap, [0; 0; 1]); tc.verifyEqual(f0, [1; 1]);
            full = tc.twoPlates(0); f1 = NonGravForces.geom.selfShadow(full, [0; 0; 1]); tc.verifyEqual(f1, [0; 1], 'fully covered lower plate');
            fb = NonGravForces.geom.selfShadow(full, [0; 0; -1]); tc.verifyEqual(fb, f1, 'from below both plates take their mirror direction +z');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, [0; 0; NaN]), 'NonGravForces:selfShadow:inputShape');
        end
        function thermalStateRobustness(tc)
            % No sources: 0 K steady, and a transient cooling from 300 K decays monotonically; an exposed-fraction array scales
            % the absorbed flux; NaN pressures refuse.
            none = struct('dir', zeros(3, 0, 1), 'P', zeros(0, 1), 'band', string.empty(1, 0));
            tc.verifyEqual(NonGravForces.geom.thermalState(tc.mm, none), zeros(8, 1));
            K = 20; epK = referenceSystem.daySecondToEpoch([58849*ones(K, 1), (0:K-1)' * 60], "GPS");
            noneK = struct('dir', zeros(3, 0, K), 'P', zeros(0, K), 'band', string.empty(1, 0));
            Tc = NonGravForces.geom.thermalState(tc.mm, noneK, 'Model', "transient", 'Epoch', epK, 'T0_K', 300 * ones(1, 8));
            tc.verifyTrue(all(diff(Tc, 1, 2) < 0, 'all')); tc.verifyTrue(all(Tc(:) > 0));
            c = referenceSystem.iersConstants().speedOfLight; src = tc.constantSource([1; 0; 0], 1361 / c);
            [T1, i1] = NonGravForces.geom.thermalState(tc.mm, src); [Th, ih] = NonGravForces.geom.thermalState(tc.mm, src, 'Exposed', 0.5 * ones(8, 1));
            tc.verifyEqual(ih.absorbed_W_m2, 0.5 * i1.absorbed_W_m2); tc.verifyEqual(Th(1), T1(1) * 0.5^0.25, 'AbsTol', 1e-9);
            bad = src; bad.P = NaN;
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(tc.mm, bad), 'NonGravForces:thermalState:numericNaN');
            share = ones(1, 8) / 8;
            tc.verifyErrorId(@() NonGravForces.geom.thermalState(tc.mm, src, 'InternalPower_W', 10, 'InternalShare', share(1:7)), 'NonGravForces:thermalState:inputInternal');
        end
        function attitudeProviderRobustness(tc)
            % Two nodes suffice; a query in another time scale is converted; squad and slerp agree at the nodes and
            % differ between them only at the interpolation level for this slowly turning attitude.
            [ep, ~, ~, q] = tc.sampleState();
            ep2 = referenceSystem.daySecondToEpoch([ep.day([1 600]) ep.sec([1 600])], "GPS");
            ap2 = NonGravForces.geom.AttitudeProvider(ep2, q(:, [1 600])); tc.verifyEqual(ap2.nodeCount, 2);
            qm = ap2.at(referenceSystem.daySecondToEpoch([58849 300], "GPS")); tc.verifyLessThan(abs(norm(qm) - 1), 1e-14);
            ap = NonGravForces.geom.AttitudeProvider(ep, q); aps = NonGravForces.geom.AttitudeProvider(ep, q, 'Method', "squad");
            epQ = referenceSystem.daySecondToEpoch([58849*ones(7, 1) (0.5:99:595)'], "GPS");
            [qU, infoU] = ap.at(referenceSystem.convertTimeScale(epQ, "UTC"));
            tc.verifyTrue(infoU.converted);
            % the 18 s UTC round trip moves the interval fraction by rounding only (measured: not bitwise, below 1e-12)
            tc.verifyBoundAndRecord("AttitudeProvider UTC vs GPS query", max(abs(qU - ap.at(epQ)), [], 'all'), 1e-12, "-", "epoch arithmetic rounding");
            tc.verifyBoundAndRecord("squad vs slerp between 1 Hz nodes", max(abs(aps.at(epQ) - ap.at(epQ)), [], 'all'), 1e-6, "-", "orbital rate 1.1e-3 rad/s: cubic vs linear over 0.5 s");
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider(ep2, q(:, [1 2 3])), 'NonGravForces:AttitudeProvider:inputNodes');
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider(referenceSystem.daySecondToEpoch([58849 0], "GPS"), q(:, 1)), 'NonGravForces:AttitudeProvider:inputNodes');
        end
        function checkAttitudeSenseRobustness(tc)
            % A tolerance below the real pointing offset turns ok false without an error; a single epoch works; mixed sense reported.
            [~, r, v, q] = tc.sampleState();
            [ok, rep] = NonGravForces.geom.checkAttitudeSense(q, r, v, 'NadirTolerance_deg', 0.5);
            tc.verifyFalse(ok); tc.verifyGreaterThan(rep.maxNadirAngle_deg, 0.5);
            [ok1, rep1] = NonGravForces.geom.checkAttitudeSense(q(:, 1), r(:, 1), v(:, 1)); tc.verifyTrue(ok1); tc.verifyEqual(rep1.K, 1);
            qmix = [q(:, 1:2) NonGravForces.geom.AttitudeProvider.nominal(r(:, 3), v(:, 3), true)];
            [~, repm] = NonGravForces.geom.checkAttitudeSense(qmix, r(:, 1:3), v(:, 1:3)); tc.verifyEqual(repm.xSense, "mixed");
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function flowShadowTableMetamorphic(tc)
            % Switching geom.flowShadow from exact to table may move the DRAG and nothing else: radiation and thermal
            % recoil are bitwise unchanged, and the drag moves, but by no more than the table's error carries.
            [ep, r, v, q] = tc.sampleState();
            eop = referenceSystem.readEop(fullfile(tc.root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            k = 1:20:numel(ep.sec);
            e = referenceSystem.daySecondToEpoch([ep.day(k) ep.sec(k)], ep.scale);
            Mt = NonGravForces.Model(tc.mm, NonGravForces.Options("standard"), 'Eop', eop, 'SpaceWeather', sw);
            Me = NonGravForces.Model(tc.mm, NonGravForces.Options("standard").set("geom.flowShadow", "exact"), 'Eop', eop, 'SpaceWeather', sw);
            [~, it] = Mt.accel(e, r(:, k), v(:, k), q(:, k)); [~, ie] = Me.accel(e, r(:, k), v(:, k), q(:, k));
            tc.verifyEqual(it.radiation, ie.radiation); tc.verifyEqual(it.thermal, ie.thermal);
            tc.verifyFalse(any(it.flowShadow.tableFallback)); tc.verifyEmpty(ie.flowShadow.tableFallback);
            d = max(abs(it.drag - ie.drag), [], 'all');
            tc.verifyGreaterThan(d, 0);
            tc.verifyBoundAndRecord("Model drag, flowShadow table vs exact (sample arc)", d / max(vecnorm(ie.drag)), 1e-4, "rel", ...
                "measured 6.7e-13 m/s^2 of 4.5e-8 on real days (1.5e-5); bound ~7x");
        end
        function flowShadowMetamorphic(tc)
            % (1) the beam field IS selfShadow along v_rel, bitwise; (2) mirroring the stream in y swaps starboard and port;
            % (3) as the gas cools (s -> inf) every ratio converges on the beam value: the error falls as s rises.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 8, 'Samples', 20);
            vb = [-7000; 2500; 1200]; U = norm(vb);
            S = NonGravForces.geom.flowShadow(tc.mm, vb, 900, 'Method', "lut", 'Lut', lut);
            tc.verifyEqual(S.beam, NonGravForces.geom.selfShadow(tc.mm, vb / U, 'Method', "lut", 'Lut', lut));
            Sm = NonGravForces.geom.flowShadow(tc.mm, vb .* [1; -1; 1], 900, 'Method', "lut", 'Lut', lut);
            swap = [1 2 5 6 3 4 7 8];
            tc.verifyBoundAndRecord("flowShadow y-mirror swaps starboard and port", max(abs(Sm.normal(swap, :, :) - S.normal), [], 'all'), 1e-12, "-", ...
                "the cube map is symmetric in y to the ray cast's sampling, which is symmetric; the flow frame is not, so rounding enters");
            err = zeros(1, 2); sList = [20 400];
            for i = 1:2
                Si = NonGravForces.geom.flowShadow(tc.mm, vb, tc.temperatureFor(sList(i), U), 'Method', "lut", 'Lut', lut);
                err(i) = max(abs(Si.normal(:, 1, 2) - Si.beam(:, 1)));
            end
            tc.verifyLessThan(err(2), err(1), 'the ratios approach the parallel beam as the thermal spread closes');
            tc.verifyBoundAndRecord("flowShadow normal ratio vs beam at s = 400", err(2), 2e-2, "-", "spread 1/s = 2.5e-3 rad across a cube cell of 0.2 rad with a value jump <= 1");
        end
        function boxWingMetamorphic(tc)
            % Scaling every dimension by s scales every area by s^2 and every centroid by s; permuting the dimensions
            % permutes the faces.
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            a = NonGravForces.geom.boxWing([2 1 0.5], o); b = NonGravForces.geom.boxWing(2 * [2 1 0.5], o);
            tc.verifyBoundAndRecord("boxWing area scaling", max(abs(b.A - 4 * a.A)), 1e-14, "m^2", "products of two lengths");
            tc.verifyEqual(b.C, 2 * a.C); tc.verifyEqual(b.N, a.N);
            p = NonGravForces.geom.boxWing([1 0.5 2], o);
            tc.verifyEqual(sort(p.A), sort(a.A));
        end
        function selfShadowMetamorphic(tc)
            % (1) sampling convergence: the 40- and 100-sample casts differ by less than the sampling bound 2/40;
            % (2) the table read at its own nodes returns the stored values (bilinear weights (1,0,0,0)); (3) the table
            % error decreases with resolution; (4) mirror symmetry of GRACE-FO in y swaps starboard and port;
            % (5) a back-facing direction returns exactly the value of its mirror direction in the plate's plane.
            [d, ~] = tc.oracleDirections();
            f40 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 40); f100 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 100);
            tc.verifyBoundAndRecord("selfShadow 40 vs 100 samples", max(abs(f40(:) - f100(:))), 2 / 40, "-", "one edge row of the coarser grid");
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 8, 'Samples', 30);
            fv = NonGravForces.geom.selfShadow(tc.mm, lut.nodes, 'Method', "lut", 'Lut', lut);
            tc.verifyTrue(all(isfinite(fv(:))), 'the interpolation must return a number for every node direction');
            tc.verifyBoundAndRecord("cube map at its nodes vs the table", max(abs(fv - lut.f), [], 'all'), 1e-12, "-", ...
                "a node is read back through two arctangents; one ulp in the face angle moves the weight by 1e-16, times a jump of at most 1 between neighbours");
            rng(7); dr = randn(3, 40); dr = dr ./ vecnorm(dr); fr = NonGravForces.geom.selfShadow(tc.mm, dr, 'Samples', 30);
            err = zeros(1, 3); res = [4 8 16];
            for i = 1:3, fl = NonGravForces.geom.selfShadow(tc.mm, dr, 'Method', "lut", 'LutResolution', res(i), 'Samples', 30); err(i) = mean(abs(fl(:) - fr(:))); end
            tc.verifyTrue(all(diff(err) < 0), 'the mean table error decreases with resolution');
            tc.verifyBoundAndRecord("cube map 16 mean error vs ray cast (40 random directions)", err(3), 1e-2, "-", ...
                "measured 3.1e-3 on 2026-09-15 (the 642-vertex icosphere it replaces measured 2.6e-2: the zeroed back faces smeared the horizon)");
            tc.verifyBoundAndRecord("cube map 16 / cube map 4 mean error ratio", err(3) / err(1), 0.5, "-", "first-order convergence in the cell size gives 0.25; measured 0.11");
            dm = d .* [1; -1; 1]; fm = NonGravForces.geom.selfShadow(tc.mm, dm, 'Samples', 40);
            swap = [1 2 5 6 3 4 7 8];
            tc.verifyEqual(fm(swap, :), f40, 'AbsTol', 1e-12, 'mirror in y swaps starboard and port');
            k = 7; nk = tc.mm.N(:, k); back = dr(:, (nk' * dr) < -1e-3); mir = back - 2 * nk .* (nk' * back);
            fb = NonGravForces.geom.selfShadow(tc.mm, back, 'Samples', 30); fmir = NonGravForces.geom.selfShadow(tc.mm, mir, 'Samples', 30);
            tc.verifyBoundAndRecord("back-facing direction vs its mirror (nadir plate)", max(abs(fb(k, :) - fmir(k, :))), 0, "-", ...
                "the same rays are cast: the mirror of the mirror is the direction itself to rounding, and the occlusion test is discrete");
        end
        function thermalStateMetamorphic(tc)
            % Energy balance in steady state (emitted = absorbed + internal) exactly; doubling every pressure raises every
            % steady temperature by 2^(1/4); a transient run from the steady state stays there.
            [ep, r, ~, q] = tc.sampleState(); idx = 1:20:600;
            epI = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], "GPS");
            src = NonGravForces.env.radiationSources(epI, r(:, idx), q(:, idx), NonGravForces.Options("standard"));
            share = zeros(1, 8); share(7) = 0.9; share([1 2 4 6]) = 0.025;
            [~, info] = NonGravForces.geom.thermalState(tc.mm, src, 'InternalPower_W', 250, 'InternalShare', share);
            balance = abs(info.emitted_W_m2 - info.absorbed_W_m2 - info.internal_W_m2) ./ (info.absorbed_W_m2 + info.internal_W_m2);
            tc.verifyBoundAndRecord("steady-state energy balance", max(balance(:)), 1e-13, "rel", "fourth root and fourth power");
            src2 = src; src2.P = 2 * src.P; T2 = NonGravForces.geom.thermalState(tc.mm, src2);
            T1 = NonGravForces.geom.thermalState(tc.mm, src);
            tc.verifyBoundAndRecord("steady temperature scales with P^(1/4)", max(abs(T2(:) - 2^0.25 * T1(:))), 1e-9, "K", "");
            Tt = NonGravForces.geom.thermalState(tc.mm, struct('dir', repmat(src.dir(:, :, 1), 1, 1, 5), 'P', repmat(src.P(:, 1), 1, 5), 'band', src.band), ...
                'Model', "transient", 'Epoch', referenceSystem.daySecondToEpoch([58849*ones(5, 1) (0:4)' * 10], "GPS"), 'T0_K', T1(:, 1)');
            tc.verifyBoundAndRecord("transient from the steady state stays", max(abs(Tt - T1(:, 1)), [], 'all'), 1e-9, "K", "fixed point");
        end
        function attitudeProviderMetamorphic(tc)
            % (1) the nominal attitude of a state rotated by R is the rotated nominal attitude (R q); (2) reversing the
            % velocity sense flips the +X choice: nominal(r, v, true) == nominal(r, -v, false) up to sign; (3) a provider
            % built from q reproduces the nominal frame it was built from at every node.
            [ep, r, v, ~] = tc.sampleState(); n = 30;
            qn = NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:n), v(:, 1:n), false);
            Rz = referenceSystem.rotationFromQuaternion([cos(0.4); 0; 0; sin(0.4)]);
            qr = NonGravForces.geom.AttitudeProvider.nominal(Rz * r(:, 1:n), Rz * v(:, 1:n), false);
            Rn = referenceSystem.rotationFromQuaternion(qn); Rr = referenceSystem.rotationFromQuaternion(qr);
            tc.verifyBoundAndRecord("nominal attitude rotates with the state", max(abs(pagemtimes(Rz, Rn) - Rr), [], 'all'), 1e-14, "-", "one rotation");
            qa = NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:n), v(:, 1:n), true); qb = NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:n), -v(:, 1:n), false);
            tc.verifyBoundAndRecord("sense flag vs reversed velocity", max(abs(abs(sum(qa .* qb, 1)) - 1)), 1e-14, "-", "same rotation up to the quaternion sign");
            epN = referenceSystem.daySecondToEpoch([ep.day(1:n) ep.sec(1:n)], "GPS");
            ap = NonGravForces.geom.AttitudeProvider(epN, qn);
            tc.verifyBoundAndRecord("provider reproduces its own nominal nodes", max(abs(ap.at(epN) - qn), [], 'all'), 1e-15, "-", "renormalisation of unit quaternions");
        end
        function checkAttitudeSenseMetamorphic(tc)
            % Rotating the whole scene (state and attitude) by R leaves every angle; swapping the sense flag of a nominal
            % attitude turns alongVelocity into againstVelocity with along-track angles summing to 180 deg.
            [~, r, v, q] = tc.sampleState(); n = 50;
            Rz = referenceSystem.rotationFromQuaternion([cos(0.7); sin(0.7) * [0.6; 0.8; 0]]);
            % rotate the scene: BODY->GCRS' = Rz * BODY->GCRS, built through the matrices (convention-free)
            qz = referenceSystem.quaternionFromRotation(pagemtimes(Rz, referenceSystem.rotationFromQuaternion(q(:, 1:n))));
            [~, a] = NonGravForces.geom.checkAttitudeSense(q(:, 1:n), r(:, 1:n), v(:, 1:n));
            [~, b] = NonGravForces.geom.checkAttitudeSense(qz, Rz * r(:, 1:n), Rz * v(:, 1:n));
            tc.verifyBoundAndRecord("sense angles invariant under a scene rotation", max(abs(a.nadirAngle_deg - b.nadirAngle_deg)), 1e-9, "deg", "two rotations and a quaternion round trip");
            [~, c1] = NonGravForces.geom.checkAttitudeSense(NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:n), v(:, 1:n), true), r(:, 1:n), v(:, 1:n));
            [~, c2] = NonGravForces.geom.checkAttitudeSense(NonGravForces.geom.AttitudeProvider.nominal(r(:, 1:n), v(:, 1:n), false), r(:, 1:n), v(:, 1:n));
            tc.verifyEqual(c1.xSense, "alongVelocity"); tc.verifyEqual(c2.xSense, "againstVelocity");
            tc.verifyBoundAndRecord("along-track angles of the two senses sum to 180", max(abs(c1.alongTrackAngle_deg + c2.alongTrackAngle_deg - 180)), 1e-9, "deg", "");
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function flowShadowTableSpeed(tc)
            % The table read against the computation it replaces, 200 epochs, after the session's table is built.
            % Measured 2026-10-17: 2-4 us against 450 us per epoch. Budget 5x the measured 0.01.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            % GRACE-FO's stream (yaw +-4 deg, pitch +-1 deg): every epoch inside the table, so this times the read, not
            % the exact fallback an epoch outside it gets
            rng(6); u = [-ones(1, 200); tand(4 * (2 * rand(1, 200) - 1)); tand(2 * rand(1, 200) - 1)];
            vb = 7600 * u ./ vecnorm(u); T = 900 + 100 * rand(1, 200);
            S = NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "table", 'Lut', lut);
            tc.assertFalse(any(S.tableFallback));
            q = tc.speedRatio(@() NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "table", 'Lut', lut), ...
                              @() NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut));
            tc.verifyBoundAndRecord("flowShadow table / Method lut time ratio (200 epochs)", q, 0.05, "-", "measured ~0.01 on 2026-10-17; budget 5x");
        end
        function flowShadowSpeed(tc)
            % flowShadow for 32 epochs against the radiation look-up the same Model epoch needs (577 source directions x 32):
            % the gas quadrature must stay within an order of magnitude of the radiation occlusion it sits beside.
            % Measured 2026-09-15: see last_run; budget 3x the measured ratio.
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut");
            rng(5); vb = [-7600 + 50 * randn(1, 32); 300 * randn(1, 32); 200 * randn(1, 32)]; T = 900 + 100 * rand(1, 32);
            d = randn(3, 577 * 32); d = d ./ vecnorm(d);
            q = tc.speedRatio(@() NonGravForces.geom.flowShadow(tc.mm, vb, T, 'Method', "lut", 'Lut', lut), ...
                              @() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'Lut', lut));
            tc.verifyBoundAndRecord("flowShadow (32 epochs) / radiation occlusion look-up (577 x 32)", q, 16, "-", "measured 5.3 on 2026-09-15; budget 3x");
        end
        function checkAttitudeSenseSpeed(tc)
            % Its whole cost is two rotations of the same quaternion series plus two atan2 angles,
            % so it must stay within a small multiple of one rotationFromQuaternion over the same K.
            [~, r, v, q] = tc.sampleState();
            base = @() referenceSystem.rotationFromQuaternion(q);
            fcn = @() NonGravForces.geom.checkAttitudeSense(q, r, v);
            ratio = tc.speedRatio(fcn, base);
            tc.verifyBoundAndRecord("checkAttitudeSense vs one rotationFromQuaternion over 600 epochs", ratio, 12, "-", "budget 3x the expected 4: two rotations, one quaternion conjugation and two atan2 angles; 5 repeats, median");
        end
        function selfShadowSpeed(tc)
            % A table query against the ray cast of the same directions. The cube map is O(1) per direction (two
            % arctangents and four reads); measured 2026-09-15: 3.6 us per direction against 1.5 ms per direction for the
            % 40-sample cast, ratio 2.4e-3. Budget 3x: 1e-2. (The icosphere it replaced scanned every face per direction.)
            rng(3); d = randn(3, 600); d = d ./ vecnorm(d);
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 16, 'Samples', 40);
            q = tc.speedRatio(@() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'Lut', lut), @() NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 40));
            tc.verifyBoundAndRecord("selfShadow cube-map query / ray cast time ratio (600 directions)", q, 1e-2, "-", "O(1) look-up; measured 2.4e-3");
        end
        function thermalStateSpeed(tc)
            % Transient over the steady state for the 600-epoch sample: the RK4 loop with dt <= tau/10 substeps.
            % Measured 2026-09-07: 0.02 s for 400 epochs; bound 20x the steady state.
            [ep, r, ~, q] = tc.sampleState();
            src = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard"));
            T0 = NonGravForces.geom.thermalState(tc.mm, struct('dir', src.dir(:, :, 1), 'P', src.P(:, 1), 'band', src.band));
            ratio = tc.speedRatio(@() NonGravForces.geom.thermalState(tc.mm, src, 'Model', "transient", 'Epoch', ep, 'T0_K', T0'), @() NonGravForces.geom.thermalState(tc.mm, src));
            tc.verifyBoundAndRecord("thermalState transient / steady time ratio (600 epochs)", ratio, 20, "-", "RK4 substeps");
        end
        function attitudeProviderSpeed(tc)
            % at() over referenceSystem.interpolateQuaternions called directly for 5000 queries: the wrapper adds a
            % scale check and a struct field. Measured 2026-09-07: 5.6 in a cold probe (the three-output call of the
            % dependency computes its diagnostics); bound 8, recorded per run.
            [ep, ~, ~, q] = tc.sampleState(); ap = NonGravForces.geom.AttitudeProvider(ep, q);
            epQ = referenceSystem.daySecondToEpoch([58849*ones(5000, 1), linspace(0, 599, 5000)'], "GPS");
            ratio = tc.speedRatio(@() ap.at(epQ), @() referenceSystem.interpolateQuaternions(q, ep, epQ));
            tc.verifyBoundAndRecord("AttitudeProvider.at / interpolateQuaternions time ratio", ratio, 8, "-", "three-output form of the dependency");
        end
    end
end
