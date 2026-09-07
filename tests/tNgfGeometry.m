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
            lit = (tc.mm.N' * d) > 0; tc.verifyTrue(all(f(~lit) == 0), 'unlit plates carry 0'); tc.verifyTrue(all(f(lit) > 0), 'lit plates are at least partly exposed');
            [fl, infoL] = NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'LutResolution', 42);
            tc.verifySize(fl, [8 size(d, 2)]); tc.verifyEqual(infoL.lut.resolution, 42); tc.verifySize(infoL.lut.f, [8 42]); tc.verifySize(infoL.lut.verts, [3 42]); tc.verifySize(infoL.lut.faces, [3 80]);
            tab = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 42); tc.verifyEqual(tab.f, infoL.lut.f);
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]), 'array', struct('vis', [0.65 0.05 0.30], 'ir', [0.81 0.03 0.16]));
            bwA = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Array', struct('area_m2', 2, 'normal', [0; 0; -1], 'centroid_m', [0; 0; -1]));
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(bwA, d), 'NonGravForces:selfShadow:inputMacro');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, [1; 1; 0]), 'NonGravForces:selfShadow:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, []), 'NonGravForces:selfShadow:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'LutResolution', 100), 'NonGravForces:selfShadow:inputResolution');
            tc.verifyErrorId(@() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'LutResolution', 162, 'Lut', infoL.lut), 'NonGravForces:selfShadow:inputLut');
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
            [d, fo] = tc.oracleDirections();
            f40 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 40); f100 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 100);
            tc.verifyBoundAndRecord("selfShadow raycast 40 vs Python 160", max(abs(f40(:) - fo(:))), 3e-2, "-", "measured 1.25e-2");
            tc.verifyBoundAndRecord("selfShadow raycast 100 vs Python 160", max(abs(f100(:) - fo(:))), 1e-2, "-", "measured 3.8e-3");
            % which plate is shadowed and by how much is the physics: only the nadir plate, at low elevation, by the overhangs
            shadowed = any(f100 < 0.999 & f100 > 0, 2);
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
        function selfShadowVectorization(tc)
            [d, ~] = tc.oracleDirections();
            fAll = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 30);
            for k = 1:size(d, 2)
                tc.verifyEqual(fAll(:, k), NonGravForces.geom.selfShadow(tc.mm, d(:, k), 'Samples', 30));
            end
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 42, 'Samples', 20);
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
            % Grazing incidence, exactly edge-on beams, a direction with the whole body lit only on one plate, and the
            % two-plate model with no overlap and with full overlap.
            f = NonGravForces.geom.selfShadow(tc.mm, [0; 1; 0]); tc.verifyTrue(all(f(:) >= 0 & f(:) <= 1)); tc.verifyEqual(f(1), 0); tc.verifyEqual(f(7), 0, 'edge-on: the nadir plate is unlit');
            fx = NonGravForces.geom.selfShadow(tc.mm, [1; 0; 0]); tc.verifyEqual(fx, [1; 0; 0; 0; 0; 0; 0; 0], 'from +x only the Front plate is lit and fully exposed');
            noOverlap = tc.twoPlates(3); f0 = NonGravForces.geom.selfShadow(noOverlap, [0; 0; 1]); tc.verifyEqual(f0, [1; 1]);
            full = tc.twoPlates(0); f1 = NonGravForces.geom.selfShadow(full, [0; 0; 1]); tc.verifyEqual(f1, [0; 1], 'fully covered lower plate');
            fb = NonGravForces.geom.selfShadow(full, [0; 0; -1]); tc.verifyEqual(fb, [0; 0], 'from below both normals face away');
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
            % (2) LUT at its own vertices reproduces the ray cast bitwise; (3) LUT error decreases with resolution;
            % (4) mirror symmetry of GRACE-FO in y: f(dy -> -dy) swaps starboard and port plates bitwise.
            [d, ~] = tc.oracleDirections();
            f40 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 40); f100 = NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 100);
            tc.verifyBoundAndRecord("selfShadow 40 vs 100 samples", max(abs(f40(:) - f100(:))), 2 / 40, "-", "one edge row of the coarser grid");
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 162, 'Samples', 30);
            fv = NonGravForces.geom.selfShadow(tc.mm, lut.verts(:, 1:40), 'Method', "lut", 'Lut', lut);
            tc.verifyBoundAndRecord("LUT at its vertices vs the table", max(abs(fv - lut.f(:, 1:40)), [], 'all'), 1e-14, "-", "barycentric weights (1,0,0)");
            rng(7); dr = randn(3, 40); dr = dr ./ vecnorm(dr); fr = NonGravForces.geom.selfShadow(tc.mm, dr, 'Samples', 30);
            err = zeros(1, 3); res = [42 162 642];
            for i = 1:3, fl = NonGravForces.geom.selfShadow(tc.mm, dr, 'Method', "lut", 'LutResolution', res(i), 'Samples', 30); err(i) = mean(abs(fl(:) - fr(:))); end
            tc.verifyTrue(all(diff(err) < 0), 'the mean LUT error decreases with resolution');
            tc.verifyBoundAndRecord("LUT 642 mean error vs ray cast (40 random directions)", err(3), 5e-2, "-", "measured see last_run; a shadow edge sharper than 6 deg is not resolved");
            tc.verifyBoundAndRecord("LUT 642 / LUT 42 mean error ratio", err(3) / err(1), 1, "-", "must improve");
            dm = d .* [1; -1; 1]; fm = NonGravForces.geom.selfShadow(tc.mm, dm, 'Samples', 40);
            swap = [1 2 5 6 3 4 7 8];
            tc.verifyEqual(fm(swap, :), f40, 'AbsTol', 1e-12, 'mirror in y swaps starboard and port');
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
            % A LUT query over the 600-epoch sample against the ray cast of the same directions: the table must pay off.
            % Measured 2026-09-07 (50 directions): ray cast 0.06 s, LUT 642 query 0.02 s; bound: LUT query <= ray cast.
            rng(3); d = randn(3, 600); d = d ./ vecnorm(d);
            lut = NonGravForces.geom.selfShadow(tc.mm, [], 'Method', "lut", 'LutResolution', 642, 'Samples', 40);
            q = tc.speedRatio(@() NonGravForces.geom.selfShadow(tc.mm, d, 'Method', "lut", 'Lut', lut), @() NonGravForces.geom.selfShadow(tc.mm, d, 'Samples', 40));
            tc.verifyBoundAndRecord("selfShadow LUT query / ray cast time ratio (600 directions)", q, 1, "-", "the table exists to be cheaper");
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
