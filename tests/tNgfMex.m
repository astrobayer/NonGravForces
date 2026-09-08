classdef tNgfMex < NgfTestCase
%TNGFMEX Parity of the compiled single-epoch kernel against the MATLAB units it ports (LP-6b, WP-17b).
%
%   MATLAB IS NORMATIVE. mex/ngfKernel.c is a transcription of
%   +NonGravForces/+kernel/*.m; it is verified against them here and against the
%   frozen Python mirror in mex/test/ngfKernelParity.c, which the CI audit job
%   runs without MATLAB. Where the two disagree the MATLAB is right.
%
%   TWO ORACLES, DELIBERATELY DIFFERENT IN KIND. The reference test uses the
%   frozen mirror fixtures (independent code, rule T9); the precision test uses
%   the MATLAB kernel chain directly, which is not independent but IS the thing
%   the port must reproduce, and it is the only place the 1e-14 claim of the
%   level prompt can be measured. Neither replaces the other.
%
%   The mirror fixtures were exported from geom once and agree with the live
%   geom.graceFO to 1e-15 (asserted by tNgfMirror/mirrorFixtureFidelity), so the
%   reference test below uses the live model and its bounds carry that round
%   trip. Only the graceFO_C cases are replayed; the fixture's synthetic cube is
%   already covered against the MATLAB by tNgfMirror.
%
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
    properties
        root
        mm            % geom.graceFO("C")
        species       % 8x1 number densities, the oxygen-dominated 490 km case
        rho           % 1x1 mass density re-formed from the species
        eop, sw, ep, r, v, att, Mfast, expo   % the real-arc fixture for the end-to-end check
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            rep = NonGravForces.buildMex('Check', true, 'Target', "kernel");
            tc.assertTrue(rep.present, ...
                'the single-epoch kernel MEX must be built before this class runs (ciGate builds it)');
            tc.assertTrue(rep.upToDate, ...
                'the kernel MEX is stale against mex/ngfKernel.c -> run NonGravForces.buildMex');
            tc.mm = NonGravForces.geom.graceFO("C");
            tc.species = [2.0e13; 8.0e13; 4.0e11; 1.0e10; 2.0e9; 5.0e10; 1.0e12; 1.0e12];
            c = NonGravForces.const();
            tc.rho = sum(tc.species .* (c.speciesMass(:) / c.NA));
            % The real-arc fixture, for the one test that closes the loop against Model.accel.
            f = @(varargin) fullfile(tc.root, varargin{:});
            tc.eop = referenceSystem.readEop(f('data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            tc.sw = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [tc.ep, tc.r, tc.v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', tc.eop);
            [eps2, qq] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
            tc.att = NonGravForces.geom.AttitudeProvider(eps2, qq);
            tc.expo = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
            tc.Mfast = NonGravForces.Model(tc.mm, NonGravForces.Options("fast"), ...
                'Attitude', tc.att, 'Eop', tc.eop, 'SpaceWeather', tc.sw, 'Exponential', tc.expo);
        end
    end
    methods
        function [gas, src, Tp] = scene(tc, seed)
            %SCENE One epoch: gas stream, three radiation sources, panel temperatures.
            N = tc.mm.n;
            v = [7590; 120; -45] + seed * [3; -7; 2];
            gas = struct('rho', tc.rho, 'species', tc.species, 'T_K', 950 + 5 * seed, ...
                         'vrel_body', v, 'shadow', 0.5 + 0.5 * cos((1:N)' + seed).^2);
            d = [1 -0.19518001458970663 0.0984135662610246; ...
                 0  0.09759000729485331 -0.14762034939153687; ...
                 0 -0.9759000729485331 -0.9841356626102459];
            d = d ./ vecnorm(d);
            src = struct('dir', d, 'P', [4.5606e-6; 1.83e-7; 7.4e-7], ...
                         'band', ["vis" "vis" "ir"], ...
                         'shadow', 0.25 + 0.25 * mod((1:N)' + (0:2), 3));
            Tp = 250 + 7 * (1:N)' + seed;
        end
        function [F, tau] = matlabChain(tc, gas, src, Tp, cdScale, crScale)
            %MATLABCHAIN The four MATLAB kernel units, composed exactly as Model.evaluate composes them.
            N = tc.mm.n; M = size(src.dir, 2);
            Fd = NonGravForces.kernel.dragPanelForce(tc.mm, gas.rho, gas.species, gas.T_K, ...
                gas.vrel_body, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', gas.shadow);
            if M > 0
                s1 = struct('dir', reshape(src.dir, 3, M, 1), 'P', reshape(src.P, M, 1), 'band', src.band);
                Fr = NonGravForces.kernel.radiationPanelForce(tc.mm, s1, ...
                    'Shadow', reshape(src.shadow, N, 1, M));
            else
                Fr = zeros(3, N, 1);
            end
            if isempty(Tp)
                Ft = zeros(3, N, 1);
            else
                Ft = NonGravForces.kernel.thermalRecoil(tc.mm, Tp);
            end
            F = zeros(3, 3); tau = zeros(3, 3);
            pages = {cdScale * Fd, crScale * Fr, Ft};
            for j = 1:3
                [F(:, j), tau(:, j)] = NonGravForces.kernel.assemble(pages{j}, tc.mm);
            end
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function singleEpochContract(tc)
            [gas, src, Tp] = tc.scene(0);
            N = tc.mm.n;
            args = {'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            % The call that must succeed, paired with every refusal below (rule T6).
            [F, tau, info] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:});
            tc.verifySize(F, [3 3]); tc.verifySize(tau, [3 3]);
            tc.verifySize(info.Fk, [3 N 3]); tc.verifySize(info.gamma, [N 1]);
            tc.verifySize(info.speedRatio, [8 1]); tc.verifySize(info.absorbed_W_m2, [N 1]);
            tc.verifyEqual(info.backend, "mex");
            tc.verifyTrue(all(isfinite(F(:))) && all(isfinite(tau(:))));

            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(struct(), gas, src, Tp, args{:}), 'NonGravForces:singleEpoch:inputMacro');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp), 'NonGravForces:singleEpoch:inputModel');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, 'Model', "cll"), 'NonGravForces:singleEpoch:inputModel');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, 'Model', "sentman"), 'NonGravForces:singleEpoch:inputModel');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, 'Model', "constantCd"), 'NonGravForces:singleEpoch:inputModel');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, 'Model', "maxwell", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:singleEpoch:inputModel');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, rmfield(gas, 'shadow'), src, Tp, args{:}), 'NonGravForces:singleEpoch:inputShape');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, rmfield(src, 'band'), Tp, args{:}), 'NonGravForces:singleEpoch:inputSources');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp(1:end-1), args{:}), 'NonGravForces:singleEpoch:inputShape');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:}, 'CentreOfMass', [0 0 0]), 'NonGravForces:singleEpoch:inputShape');

            bad = gas; bad.vrel_body = [0; 0; 0];
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, bad, src, Tp, args{:}), 'NonGravForces:singleEpoch:inputRange');
            bad = gas; bad.rho = -1;
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, bad, src, Tp, args{:}), 'NonGravForces:singleEpoch:inputRange');
            bad = gas; bad.shadow(1) = 1.5;
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, bad, src, Tp, args{:}), 'NonGravForces:singleEpoch:inputRange');
            bad = gas; bad.T_K = NaN;
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, bad, src, Tp, args{:}), 'NonGravForces:singleEpoch:numericNaN');
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:}, 'Alpha', 1.5), 'NonGravForces:singleEpoch:inputRange');
            bads = src; bads.dir(:, 1) = 2 * bads.dir(:, 1);
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, bads, Tp, args{:}), 'NonGravForces:singleEpoch:inputSources');
            bads = src; bads.band(1) = "uv";
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, bads, Tp, args{:}), 'NonGravForces:singleEpoch:inputSources');
            bads = src; bads.P(1) = -1;
            tc.verifyErrorId(@() NonGravForces.kernel.singleEpoch(tc.mm, gas, bads, Tp, args{:}), 'NonGravForces:singleEpoch:inputRange');
        end
        function singleEpochGatewayContract(tc)
            % The gas-surface model reaches the C as a code, and only three codes exist. The gateway
            % must be unreachable except through the unit above: the private folder is what enforces
            % that, and this asserts the folder is where the unit looks for it.
            gw = fullfile(tc.root, '+NonGravForces', '+kernel', 'private', "kernelGateway." + mexext);
            tc.verifyTrue(isfile(gw), 'the kernel gateway must live in +kernel/private, reachable from kernel.singleEpoch alone');
            tc.verifyEqual(exist('kernelGateway', 'file'), 0, 'the gateway must NOT be visible outside the kernel package');
            [gas, src, Tp] = tc.scene(0);
            % constantCd needs neither species nor a temperature, and must not read them.
            g2 = gas; g2.species = []; g2.T_K = 950;
            F1 = NonGravForces.kernel.singleEpoch(tc.mm, g2, src, Tp, 'Model', "constantCd", 'Cd', 2.2);
            g3 = g2; g3.T_K = 1400;
            F2 = NonGravForces.kernel.singleEpoch(tc.mm, g3, src, Tp, 'Model', "constantCd", 'Cd', 2.2);
            tc.verifyBoundAndRecord("constantCd ignores the gas temperature", max(abs(F1 - F2), [], 'all'), 0, "N", "the model has no speed ratio, so T cannot enter");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function singleEpochMirrorReference(tc)
            % The frozen Python mirror (WP-33) is the independent authority. Each unit's cases are
            % replayed through the composite by silencing the other two effects: a zero source list
            % and no panel temperature leave only the gas force, and a zero species vector leaves
            % only the radiation and the recoil.
            p = @(f) fullfile(tc.root, 'tools', 'pyxcheck', f);
            noSrc = struct('dir', zeros(3, 0), 'P', zeros(0, 1), 'band', strings(1, 0), 'shadow', []);
            zeroGas = struct('rho', 0, 'species', zeros(8, 1), 'T_K', 950, ...
                             'vrel_body', [7600; 0; 0], 'shadow', []);

            % --- drag
            T = readtable(p("mirror_drag.csv"), 'TextType', 'string');
            T = T(T.body == "graceFO_C" & T.model == "sentman", :);
            cases = unique(T.caseId, 'stable'); worst = 0;
            for c = cases'
                r = T(T.caseId == c, :);
                sp = [r.d1(1) r.d2(1) r.d3(1) r.d4(1) r.d5(1) r.d6(1) r.d7(1) r.d8(1)]';
                cst = NonGravForces.const();
                gas = struct('rho', sum(sp .* (cst.speciesMass(:) / cst.NA)), 'species', sp, ...
                             'T_K', r.t_k(1), 'vrel_body', [r.vx(1); r.vy(1); r.vz(1)], 'shadow', r.shadow);
                F = NonGravForces.kernel.singleEpoch(tc.mm, gas, noSrc, [], ...
                    'Model', "sentman", 'Alpha', r.alpha(1), 'Tw_K', r.tw_k(1));
                want = sum([r.fx r.fy r.fz]', 2);
                worst = max(worst, max(abs(F(:, 1) - want)) / max(max(abs(want)), realmin));
                tc.verifyEqual(F(:, 2:3), zeros(3, 2), 'a silenced effect must contribute exactly nothing');
            end
            tc.verifyBoundAndRecord("singleEpoch gas force vs the mirror (" + numel(cases) + " cases)", worst, 1e-12, "-", ...
                "the C, the MATLAB and the mirror are three derivations of the same closed form in double precision");

            % --- radiation
            T = readtable(p("mirror_radiation.csv"), 'TextType', 'string');
            T = T(T.body == "graceFO_C", :);
            cases = unique(T.caseId, 'stable'); worstR = 0; worstA = 0;
            for c = cases'
                r = T(T.caseId == c, :); N = height(r);
                d = zeros(3, 3); P = zeros(3, 1); band = strings(1, 3); sh = zeros(N, 3);
                for m = 1:3
                    fld = @(suffix) char("s" + m + suffix);      % dynamic field, then index (a call result cannot be indexed)
                    d(:, m) = [r.(fld('dx'))(1); r.(fld('dy'))(1); r.(fld('dz'))(1)];
                    P(m) = r.(fld('p'))(1); band(m) = r.(fld('band'))(1); sh(:, m) = r.(fld('shadow'));
                end
                src = struct('dir', d, 'P', P, 'band', band, 'shadow', sh);
                [F, ~, info] = NonGravForces.kernel.singleEpoch(tc.mm, zeroGas, src, [], ...
                    'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
                want = sum([r.fx r.fy r.fz]', 2);
                worstR = max(worstR, max(abs(F(:, 2) - want)) / max(max(abs(want)), realmin));
                worstA = max(worstA, max(abs(info.absorbed_W_m2 - r.absorbed)) / max(max(abs(r.absorbed)), realmin));
            end
            tc.verifyBoundAndRecord("singleEpoch radiation force vs the mirror (" + numel(cases) + " cases)", worstR, 1e-12, "-", "independently derived optical model");
            tc.verifyBoundAndRecord("singleEpoch absorbed flux vs the mirror", worstA, 1e-12, "-", "a E gamma summed over sources");

            % --- thermal
            T = readtable(p("mirror_thermal.csv"), 'TextType', 'string');
            T = T(T.body == "graceFO_C", :);
            cases = unique(T.caseId, 'stable'); worstT = 0;
            for c = cases'
                r = T(T.caseId == c, :);
                F = NonGravForces.kernel.singleEpoch(tc.mm, zeroGas, noSrc, r.temp_k, ...
                    'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
                want = sum([r.fx r.fy r.fz]', 2);
                worstT = max(worstT, max(abs(F(:, 3) - want)) / max(max(abs(want)), realmin));
            end
            tc.verifyBoundAndRecord("singleEpoch thermal recoil vs the mirror (" + numel(cases) + " cases)", worstT, 1e-13, "-", "Stefan-Boltzmann and one hemispherical factor: only rounding may differ");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function singleEpochMatlabParityPrecision(tc)
            % THE CLAIM OF THE LEVEL PROMPT: the compiled path reproduces the MATLAB units to 1e-14.
            % Relative to the largest component of each column, because the three effects differ by
            % two orders of magnitude and an absolute bound would only measure the largest.
            worstF = 0; worstT = 0; worstK = 0;
            for seed = 0:4
                [gas, src, Tp] = tc.scene(seed);
                cd = 1 + 0.01 * seed; cr = 1 - 0.005 * seed;
                [F, tau, info] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, ...
                    'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'CdScale', cd, 'CrScale', cr);
                [Fm, taum] = tc.matlabChain(gas, src, Tp, cd, cr);
                for j = 1:3
                    sF = max(max(abs(Fm(:, j))), realmin);
                    sT = max(max(abs(taum(:, j))), realmin);
                    worstF = max(worstF, max(abs(F(:, j) - Fm(:, j))) / sF);
                    worstT = max(worstT, max(abs(tau(:, j) - taum(:, j))) / sT);
                end
                % The per-panel arrays too, or a compensating pair of errors in the reduction
                % would pass the totals and hide in the panels.
                Fd = NonGravForces.kernel.dragPanelForce(tc.mm, gas.rho, gas.species, gas.T_K, ...
                    gas.vrel_body, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', gas.shadow);
                k = abs(info.Fk(:, :, 1) - cd * reshape(Fd, 3, tc.mm.n));
                worstK = max(worstK, max(k, [], 'all') / max(max(abs(cd * Fd), [], 'all'), realmin));
            end
            tc.verifyBoundAndRecord("singleEpoch force vs the MATLAB kernel chain (5 epochs)", worstF, 1e-14, "-", "the same closed forms in the same association order; only rounding may differ");
            tc.verifyBoundAndRecord("singleEpoch torque vs the MATLAB kernel chain", worstT, 1e-14, "-", "cross products and sums in the same plate order");
            tc.verifyBoundAndRecord("singleEpoch per-panel gas force vs the MATLAB kernel", worstK, 1e-14, "-", "the reduction cannot hide a per-panel error if the panels are compared too");

            % The diagnostics must be the same numbers, not merely consistent ones.
            [gas, src, Tp] = tc.scene(1);
            [~, ~, info] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            [~, di] = NonGravForces.kernel.dragPanelForce(tc.mm, gas.rho, gas.species, gas.T_K, ...
                gas.vrel_body, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', gas.shadow);
            tc.verifyBoundAndRecord("singleEpoch gamma vs the MATLAB kernel", max(abs(info.gamma - di.gamma)), 1e-15, "-", "one dot product per plate");
            tc.verifyBoundAndRecord("singleEpoch speed ratio vs the MATLAB kernel", max(abs(info.speedRatio - di.speedRatio) ./ di.speedRatio), 1e-15, "-", "U / sqrt(2 k T / m) per species");
            tc.verifyBoundAndRecord("singleEpoch mixture Cp vs the MATLAB kernel", max(abs(info.Cp - di.Cp) ./ max(abs(di.Cp), realmin)), 1e-13, "-", "mass-weighted sum over eight species in the same order");
        end
        function singleEpochEndToEndPrecision(tc)
            % THE TWO PATHS MUST AGREE, and this is where that is checked at the level a user sees.
            % Owner ruling of 2026-09-08: the vectorised MATLAB path and the pointwise C path are
            % BOTH first class; neither replaces the other. Two paths that are both kept are two
            % paths that must be kept in step, so this test composes the environment exactly as
            % Model.evaluate composes it, hands one epoch to the compiled kernel, and requires the
            % acceleration back to be the one Model.accel returns for that epoch.
            %
            % The fast preset is used deliberately: it sets geom.selfShadow "none" and
            % thermal.model "none", so the shadow arrays are all ones and the recoil is zero, and
            % what is left under test is the composition itself rather than two more interpolators.
            o = tc.Mfast.options;
            idx = 7;                                            % an arbitrary epoch of the real arc
            ep = referenceSystem.daySecondToEpoch([tc.ep.day(idx) tc.ep.sec(idx)], tc.ep.scale);
            rr = tc.r(:, idx); vv = tc.v(:, idx);
            q = tc.att.at(ep);
            [aRef, ai] = tc.Mfast.accel(ep, rr, vv, q);

            % --- the same environment, through the same public units, in the same order
            vrel = NonGravForces.env.wind(o.drag.wind, ep, rr, vv, 'Eop', tc.eop);
            R = referenceSystem.rotationFromQuaternion(q);
            vrelBody = R' * vrel;
            Ritrs = referenceSystem.rotationItrsToGcrs(ep, 'Eop', tc.eop);
            lla = referenceSystem.cartesianToGeodetic((Ritrs' * rr)', 'Ellipsoid', "WGS84")';
            [rho, sp, T] = NonGravForces.env.density("exponential", ep, lla, [], 'Exponential', tc.expo);
            src = NonGravForces.env.radiationSources(ep, rr, q, o);
            gas = struct('rho', rho, 'species', sp, 'T_K', T(2), 'vrel_body', vrelBody, 'shadow', []);

            [F, ~, info] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, [], ...
                'Model', o.drag.gsi, 'Cd', o.drag.cd, 'CdScale', o.drag.cdScale, 'CrScale', o.srp.crScale);
            m = ai.mass_kg;
            got = [F(:, 1) / m, F(:, 2) / m, F(:, 3) / m];
            want = [ai.drag, ai.radiation, ai.thermal];
            scale = max(max(abs(want), [], 'all'), realmin);
            tc.verifyBoundAndRecord("compiled single epoch vs the vectorised path, per effect", ...
                max(abs(got - want), [], 'all') / scale, 1e-13, "-", ...
                "the same environment, the same panel loop, the same reduction: only the association of the mass division and the C rounding may differ");
            tc.verifyBoundAndRecord("compiled single epoch vs the vectorised path, total acceleration", ...
                max(abs(sum(got, 2) - aRef)) / max(max(abs(aRef)), realmin), 1e-13, "-", ...
                "Model.accel adds three quotients, this adds three forces and divides once");
            % The relative velocity is the input most easily got wrong (co-rotation, not omega x r),
            % so it is asserted against the Model's own report rather than assumed to match.
            tc.verifyBoundAndRecord("the relative velocity handed to the C equals the Model's", ...
                max(abs(vrelBody - ai.vrel_body)), 0, "m/s", "both come from env.wind through the same rotation");
            tc.verifyEqual(numel(info.gamma), tc.mm.n);
        end
    end

    % ================================================================ vectorization
    methods (Test, TestTags = {'vectorization'})
        function singleEpochVectorization(tc)
            % The compiled path is single-epoch by construction, so the claim is the other way
            % round: K calls to it must reproduce the K-vectorised MATLAB kernel column by column.
            % That is the property an integrator relies on when it swaps one for the other.
            K = 4; N = tc.mm.n;
            vrel = zeros(3, K); shadow = zeros(N, K); Tp = zeros(N, K);
            for k = 1:K
                [g, ~, t] = tc.scene(k - 1);
                vrel(:, k) = g.vrel_body; shadow(:, k) = g.shadow; Tp(:, k) = t;
            end
            gasT = 950 + 5 * (0:K-1);
            Fbatch = NonGravForces.kernel.dragPanelForce(tc.mm, repmat(tc.rho, 1, K), ...
                repmat(tc.species, 1, K), gasT, vrel, 'Model', "sentman", 'Alpha', 0.85, ...
                'Tw_K', 300, 'Shadow', shadow);
            Tbatch = NonGravForces.kernel.thermalRecoil(tc.mm, Tp);
            worst = 0;
            noSrc = struct('dir', zeros(3, 0), 'P', zeros(0, 1), 'band', strings(1, 0), 'shadow', []);
            for k = 1:K
                gas = struct('rho', tc.rho, 'species', tc.species, 'T_K', gasT(k), ...
                             'vrel_body', vrel(:, k), 'shadow', shadow(:, k));
                [~, ~, info] = NonGravForces.kernel.singleEpoch(tc.mm, gas, noSrc, Tp(:, k), ...
                    'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
                s = max(max(abs(Fbatch(:, :, k)), [], 'all'), realmin);
                worst = max(worst, max(abs(info.Fk(:, :, 1) - Fbatch(:, :, k)), [], 'all') / s);
                sT = max(max(abs(Tbatch(:, :, k)), [], 'all'), realmin);
                worst = max(worst, max(abs(info.Fk(:, :, 3) - Tbatch(:, :, k)), [], 'all') / sT);
            end
            tc.verifyBoundAndRecord("K single-epoch calls vs the K-vectorised MATLAB kernel", worst, 1e-14, "-", "epoch k of a batch and the same epoch alone are the same physics");
        end
    end

    % ==================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function singleEpochRobustness(tc)
            N = tc.mm.n;
            [gas, src, Tp] = tc.scene(0);
            args = {'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};

            % No sources at all: legal, and the radiation column is exactly zero, not NaN.
            noSrc = struct('dir', zeros(3, 0), 'P', zeros(0, 1), 'band', strings(1, 0), 'shadow', []);
            [F, tau, info] = NonGravForces.kernel.singleEpoch(tc.mm, gas, noSrc, Tp, args{:});
            tc.verifyEqual(F(:, 2), zeros(3, 1)); tc.verifyEqual(tau(:, 2), zeros(3, 1));
            tc.verifyEqual(info.absorbed_W_m2, zeros(N, 1));
            tc.verifyTrue(all(isfinite(F(:))));

            % No panel temperatures: the thermal column is exactly zero and nothing else moves.
            [F2, ~, i2] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, [], args{:});
            tc.verifyEqual(F2(:, 3), zeros(3, 1));
            tc.verifyEqual(i2.emitted_W_m2, zeros(N, 1));

            % A plate at 0 K emits nothing and does not poison its neighbours.
            Tz = Tp; Tz(1) = 0;
            [~, ~, i3] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tz, args{:});
            tc.verifyEqual(i3.emitted_W_m2(1), 0);
            tc.verifyTrue(all(i3.emitted_W_m2(2:end) > 0));

            % Fully shadowed: every panel force vanishes, and the reduction stays finite.
            gz = gas; gz.shadow = zeros(N, 1);
            sz = src; sz.shadow = zeros(N, 3);
            [F4, ~, i4] = NonGravForces.kernel.singleEpoch(tc.mm, gz, sz, [], args{:});
            tc.verifyBoundAndRecord("singleEpoch under total shadow", max(abs(F4), [], 'all'), 0, "N", "a zero fraction multiplies every panel force, both effects");
            tc.verifyTrue(all(isfinite(i4.Fk(:))));

            % Zero density: the gas force vanishes but gamma and theta are still reported, because
            % a caller screening a pass needs them whether or not there was any gas.
            gz = gas; gz.species = zeros(8, 1); gz.rho = 0;
            [F5, ~, i5] = NonGravForces.kernel.singleEpoch(tc.mm, gz, src, Tp, args{:});
            tc.verifyEqual(F5(:, 1), zeros(3, 1));
            tc.verifyTrue(all(isfinite(i5.gamma)) && all(isfinite(i5.theta_rad)));
            tc.verifyEqual(i5.Cp, zeros(N, 1));

            % The stream reversed: a different set of plates faces the flow, so the magnitude is not
            % predictable, but the gas force must still oppose the motion. Sentman's shear is lift and
            % tilts it off the stream by a fraction of a degree; the projection on the stream cannot
            % change sign, because a free-molecular gas cannot push a body forward.
            gl = gas; gl.vrel_body = -gas.vrel_body; gl.shadow = ones(N, 1);
            [F6, ~, i6] = NonGravForces.kernel.singleEpoch(tc.mm, gl, noSrc, [], args{:});
            tc.verifyTrue(all(isfinite(F6(:))) && all(isfinite(i6.Fk(:))));
            cosAngle = dot(F6(:, 1), gl.vrel_body) / (norm(F6(:, 1)) * norm(gl.vrel_body));
            tc.verifyBoundAndRecord("reversed stream: cosine between the gas force and the velocity", ...
                cosAngle, 0, "-", "a dissipative force has a non-positive projection on the velocity, whichever plates face the flow");
        end
    end

    % =================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function singleEpochMetamorphic(tc)
            [gas, src, Tp] = tc.scene(2);
            args = {'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            F1 = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:});

            % 1. The gas force is linear in the number densities: double them, double it.
            g2 = gas; g2.species = 2 * gas.species; g2.rho = 2 * gas.rho;
            F2 = NonGravForces.kernel.singleEpoch(tc.mm, g2, src, Tp, args{:});
            tc.verifyBoundAndRecord("gas force doubles with the density", ...
                max(abs(F2(:, 1) - 2 * F1(:, 1))) / max(abs(2 * F1(:, 1))), 1e-15, "-", ...
                "rho enters every species term linearly; the coefficients do not depend on it");

            % 2. The radiation force is linear in the pressures.
            s3 = src; s3.P = 3 * src.P;
            F3 = NonGravForces.kernel.singleEpoch(tc.mm, gas, s3, Tp, args{:});
            tc.verifyBoundAndRecord("radiation force triples with the pressure", ...
                max(abs(F3(:, 2) - 3 * F1(:, 2))) / max(abs(3 * F1(:, 2))), 1e-15, "-", ...
                "P multiplies the whole optical expression");

            % 3. The recoil follows T^4 exactly.
            F4 = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, 2 * Tp, args{:});
            tc.verifyBoundAndRecord("thermal recoil scales as T^4", ...
                max(abs(F4(:, 3) - 16 * F1(:, 3))) / max(abs(16 * F1(:, 3))), 1e-15, "-", ...
                "Stefan-Boltzmann with an exact power of two");

            % 4. CdScale and CrScale reach their own column and no other.
            F5 = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:}, 'CdScale', 1.5, 'CrScale', 0.5);
            tc.verifyBoundAndRecord("CdScale multiplies the gas column alone", ...
                max([abs(F5(:, 1) - 1.5 * F1(:, 1)) ./ max(abs(1.5 * F1(:, 1))); ...
                     abs(F5(:, 3) - F1(:, 3)) ./ max(abs(F1(:, 3)))]), 1e-15, "-", ...
                "the scales are applied per effect, where Model.evaluate applies them");
            tc.verifyBoundAndRecord("CrScale multiplies the radiation column alone", ...
                max(abs(F5(:, 2) - 0.5 * F1(:, 2))) / max(abs(0.5 * F1(:, 2))), 1e-15, "-", ...
                "same, for the radiation column");

            % 5. The torque reference point shifts the torque by r x F and nothing else.
            shift = [0.1; -0.2; 0.05];
            [Fa, taua] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:});
            [Fb, taub] = NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:}, ...
                'CentreOfMass', tc.mm.centreOfMass + shift);
            tc.verifyEqual(Fb, Fa, 'moving the reference point cannot change a force');
            predicted = taua + cross(repmat(-shift, 1, 3), Fa);
            tc.verifyBoundAndRecord("torque shifts by -d x F when the reference point moves", ...
                max(abs(taub - predicted), [], 'all') / max(abs(taua), [], 'all'), 1e-13, "-", ...
                "the lever arms all change by the same vector, so the sum changes by that vector crossed into the total force");
        end
    end

    % ========================================================================= speed
    methods (Test, TestTags = {'speed'})
        function singleEpochSpeed(tc)
            % WHAT THIS MEASURES AND WHY. HANDOVER D.6 justifies the compiled path with a MATLAB
            % call overhead of 5-15 microseconds, but LP-6's own measurement (O-29) puts Model.accel
            % at 88 ms per epoch in the standard preset and 2 ms in the fast one, where 15 us is
            % 0.02 % and 0.75 %. The premise therefore needs a number, and this is it: the ratio of
            % the compiled single-epoch kernel to the MATLAB kernel chain that computes exactly the
            % same three effects. Nothing here claims anything about Model.accel, whose cost is
            % dominated by geom.selfShadow and not by this kernel.
            [gas, src, Tp] = tc.scene(0);
            args = {'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            reps = 200;
            NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:});   % warm the path
            tc.matlabChain(gas, src, Tp, 1, 1);
            tMex = zeros(1, 5); tMat = zeros(1, 5);
            for trial = 1:5
                t0 = tic;
                for i = 1:reps, NonGravForces.kernel.singleEpoch(tc.mm, gas, src, Tp, args{:}); end
                tMex(trial) = toc(t0) / reps;
                t0 = tic;
                for i = 1:reps, tc.matlabChain(gas, src, Tp, 1, 1); end
                tMat(trial) = toc(t0) / reps;
            end
            % Medians, not means: one scheduling hiccup in 200 calls should not decide a bound (T5).
            mex_s = median(tMex); mat_s = median(tMat);
            tc.verifyBoundAndRecord("single-epoch kernel: compiled / MATLAB chain (median of 5 x 200)", ...
                mex_s / mat_s, 1.0, "-", ...
                sprintf("compiled %.1f us, MATLAB %.1f us per epoch on this host; a value above 1 is a FINDING, not a flake: it says the port buys nothing and D.6 must be rewritten", 1e6 * mex_s, 1e6 * mat_s));
            % NOT a claim about this host, and the 500 us written here on 2026-09-08 was a guess with
            % nothing behind it: the CI runner measured 901.8 us and the test went red for no reason
            % but the guess (rule T5 forbids exactly this, and rule A.4 records the same mistake being
            % made once before). Decomposed rather than widened: the assertion above is the ratio, and
            % this ceiling sits an order of magnitude above the measured value so that a pathological
            % rebuild -- a MEX linked without -O, a marshalling path that starts copying per call --
            % still turns the gate red while a slower machine does not.
            tc.verifyBoundAndRecord("single-epoch compiled kernel, absolute cost per epoch", 1e6 * mex_s, 10000, "us", ...
                "smoke ceiling, not a claim; 901.8 us measured on the CI runner 2026-09-08, MATLAB chain 4149.6 us");
        end
    end
end
