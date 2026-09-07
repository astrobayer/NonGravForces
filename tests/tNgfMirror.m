classdef tNgfMirror < NgfTestCase
%TNGFMIRROR Parity of the MATLAB kernel against the independent Python mirror (WP-33).
%
%   The oracle is tools/pyxcheck/ngfkernel/: a standard-library Python mirror of the
%   frozen kernel, derived from the physics by agents that were denied the MATLAB
%   sources, the existing Sentman oracle, the kernel tests and the records, so it
%   cannot have inherited this implementation's mistakes (rule T9). The hardest
%   quantity, the gas-surface coefficients, was derived three times by three routes
%   and the routes reconciled numerically before the module was accepted.
%
%   tools/pyxcheck/kernel_mirror_check.py evaluates the mirror over a fixed case list
%   and writes one CSV per unit holding the inputs AND the mirror's answer. This class
%   feeds the same inputs to the MATLAB kernel and compares. Plate geometry and optics
%   are input data exported once from geom (the ruling mirrors the kernel only), and
%   mirrorFixtureFidelity checks that export still matches the live geom.
%
%   AGREEMENT CLAIMED. Both sides are double precision but evaluate different closed
%   forms, so the claim is relative agreement on the scale of each case, not bitwise
%   equality. Quantities far below the scale of their case (a leeward plate at speed
%   ratio 8 carries e-30 of the windward force) are required to be negligible on both
%   sides rather than relatively equal: relative error is meaningless there and
%   asserting it would only measure which form cancels first.
%
% Created 2026-09-07 by Claude Opus 5 (Anthropic).
    properties
        root
        plates      % table, exported plate geometry and optics
        bodies      % table, centre of mass and inertia
        macro       % containers.Map body name -> MacroModel rebuilt from the fixture
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            p = @(f) fullfile(tc.root, 'tools', 'pyxcheck', f);
            for f = ["kernel_mirror_plates.csv", "kernel_mirror_bodies.csv", "mirror_gsi.csv", ...
                     "mirror_drag.csv", "mirror_radiation.csv", "mirror_thermal.csv", ...
                     "mirror_assemble.csv", "mirror_torques.csv"]
                tc.assertTrue(isfile(p(f)), sprintf('mirror fixture %s is missing -> python tools/pyxcheck/kernel_mirror_check.py', f));
            end
            tc.plates = readtable(p("kernel_mirror_plates.csv"), 'TextType', 'string');
            tc.bodies = readtable(p("kernel_mirror_bodies.csv"), 'TextType', 'string');
            tc.macro = containers.Map();
            for b = 1:height(tc.bodies)
                name = tc.bodies.body(b);
                tc.macro(char(name)) = tc.buildMacro(name);
            end
        end
    end
    methods
        function mm = buildMacro(tc, name)
            %BUILDMACRO Rebuild a MacroModel from the exported fixture so both sides see one geometry.
            rows = tc.plates(tc.plates.body == name, :);
            br = tc.bodies(tc.bodies.body == name, :);
            n = height(rows);
            pl = repmat(struct('name', "", 'kind', "plate", 'area_m2', 0, 'normal', [0; 0; 1], 'centroid_m', [0; 0; 0], ...
                'absVis', 0, 'specVis', 0, 'diffVis', 0, 'emisIR', 0, 'specIR', 0, 'diffIR', 0, 'alpha', 0.85, 'group', "bus"), 1, n);
            for k = 1:n
                pl(k).name = rows.name(k);
                pl(k).area_m2 = rows.area_m2(k);
                pl(k).normal = [rows.nx(k); rows.ny(k); rows.nz(k)];
                pl(k).centroid_m = [rows.cx(k); rows.cy(k); rows.cz(k)];
                pl(k).absVis = rows.absVis(k); pl(k).specVis = rows.specVis(k); pl(k).diffVis = rows.diffVis(k);
                pl(k).emisIR = rows.emisIR(k); pl(k).specIR = rows.specIR(k); pl(k).diffIR = rows.diffIR(k);
            end
            I = [br.I11 br.I12 br.I13; br.I21 br.I22 br.I23; br.I31 br.I32 br.I33];
            mm = NonGravForces.geom.MacroModel(pl, 'Name', name, 'Inertia', I, 'CentreOfMass', [br.cmx; br.cmy; br.cmz]);
        end
        function e = relErr(~, got, want, scale)
            %RELERR Error of got against want, on the given scale; NaN-free, empty-safe.
            if isempty(want), e = 0; return; end
            s = max(scale, realmin);
            e = max(abs(got(:) - want(:))) / s;
        end
        function [big, small] = splitByScale(~, want, floorAbs)
            big = abs(want) >= floorAbs; small = ~big;
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function mirrorFixtureFidelity(tc)
            % The exported plate table must still be what geom produces, or every case below is stale.
            live = NonGravForces.geom.graceFO("C");
            fx = tc.macro('graceFO_C');
            tc.verifyEqual(fx.n, live.n);
            tc.verifyBoundAndRecord("fixture areas vs live geom.graceFO", max(abs(fx.A(:) - live.A(:))), 0, "m^2", "exported from the same object; must be bitwise");
            tc.verifyBoundAndRecord("fixture normals vs live geom.graceFO", max(abs(fx.N(:) - live.N(:))), 1e-15, "-", "one CSV round trip through 17 significant digits");
            tc.verifyBoundAndRecord("fixture centroids vs live geom.graceFO", max(abs(fx.C(:) - live.C(:))), 1e-15, "m", "one CSV round trip");
            tc.verifyBoundAndRecord("fixture optics vs live geom.graceFO", max(abs([fx.optVis(:) - live.optVis(:); fx.optIR(:) - live.optIR(:)])), 1e-15, "-", "one CSV round trip");
            tc.verifyBoundAndRecord("fixture inertia vs live geom.graceFO", max(abs(fx.inertia(:) - live.inertia(:))), 1e-9, "kg m^2", "printed to 17 digits from Table 4");
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function gsiCoefficientsMirror(tc)
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'mirror_gsi.csv'), 'TextType', 'string');
            n = height(T); cp = zeros(n, 1); ct = zeros(n, 1);
            for i = 1:n
                switch T.model(i)
                    case "constantCd"
                        [cp(i), ct(i)] = NonGravForces.kernel.gsiCoefficients("constantCd", [], T.theta(i), [], [], [], 'Cd', T.cd(i));
                    case "maxwell"
                        [cp(i), ct(i)] = NonGravForces.kernel.gsiCoefficients("maxwell", T.s(i), T.theta(i), T.alpha(i), T.tw_k(i), T.tinf_k(i), 'SpecularFraction', T.specular(i));
                    otherwise
                        [cp(i), ct(i)] = NonGravForces.kernel.gsiCoefficients("sentman", T.s(i), T.theta(i), T.alpha(i), T.tw_k(i), T.tinf_k(i));
                end
            end
            % Coefficients above 1e-12 are compared relatively; below it both sides must simply be negligible.
            for pair = {{cp, T.cp, "cp"}, {ct, T.ctau, "ctau"}}
                [got, want, nm] = pair{1}{:};
                big = abs(want) >= 1e-12;
                tc.verifyBoundAndRecord("gsiCoefficients " + nm + " vs the mirror (relative, " + nnz(big) + " of " + n + " cases)", ...
                    max(abs(got(big) - want(big)) ./ abs(want(big))), 1e-12, "-", "two independently derived closed forms in double precision");
                tc.verifyBoundAndRecord("gsiCoefficients " + nm + " vs the mirror (" + nnz(~big) + " negligible cases)", ...
                    max([abs(got(~big)); abs(want(~big)); 0]), 1e-11, "-", "leeward plates: both sides must be negligible, relative agreement is meaningless there");
            end
        end
        function dragPanelForceMirror(tc)
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'mirror_drag.csv'), 'TextType', 'string');
            dens = [T.d1 T.d2 T.d3 T.d4 T.d5 T.d6 T.d7 T.d8];
            cases = unique(T.caseId, 'stable'); worst = 0; worstNeg = 0;
            for c = cases'
                r = T(T.caseId == c, :); mm = tc.macro(char(r.body(1)));
                sp = dens(find(T.caseId == c, 1), :)';
                rho = sum(sp .* (NonGravForces.const().speciesMass(:) / NonGravForces.const().NA));
                v = [r.vx(1); r.vy(1); r.vz(1)]; f = r.shadow;
                switch r.model(1)
                    case "constantCd"
                        F = NonGravForces.kernel.dragPanelForce(mm, rho, [], [], v, 'Model', "constantCd", 'Cd', r.cd(1), 'Shadow', f);
                    case "maxwell"
                        F = NonGravForces.kernel.dragPanelForce(mm, rho, sp, r.t_k(1), v, 'Model', "maxwell", 'Alpha', r.alpha(1), 'Tw_K', r.tw_k(1), 'SpecularFraction', r.specular(1), 'Shadow', f);
                    otherwise
                        F = NonGravForces.kernel.dragPanelForce(mm, rho, sp, r.t_k(1), v, 'Model', "sentman", 'Alpha', r.alpha(1), 'Tw_K', r.tw_k(1), 'Shadow', f);
                end
                want = [r.fx r.fy r.fz]'; got = reshape(F, 3, []);
                scale = max(abs(want(:)));
                per = abs(got - want) / max(scale, realmin);
                big = max(abs(want), [], 1) >= 1e-12 * scale;      % plates that carry a real share of this case
                if any(big), worst = max(worst, max(per(:, big), [], 'all')); end
                if any(~big), worstNeg = max(worstNeg, max(abs(got(:, ~big)), [], 'all') / max(scale, realmin)); end
            end
            tc.verifyBoundAndRecord("dragPanelForce vs the mirror (" + numel(cases) + " cases, relative to the case scale)", worst, 1e-12, "-", "independent mixture rule and plate geometry, double precision");
            tc.verifyBoundAndRecord("dragPanelForce, plates negligible in their case", worstNeg, 1e-11, "-", "leeward plates: negligible on both sides");
        end
        function radiationPanelForceMirror(tc)
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'mirror_radiation.csv'), 'TextType', 'string');
            cases = unique(T.caseId, 'stable'); worst = 0; worstAbs = 0;
            for c = cases'
                r = T(T.caseId == c, :); mm = tc.macro(char(r.body(1))); N = height(r);
                d = zeros(3, 3, 1); P = zeros(3, 1); band = strings(1, 3); sh = zeros(N, 1, 3);
                for m = 1:3
                    fld = @(suffix) char("s" + m + suffix);
                    d(:, m, 1) = [r.(fld('dx'))(1); r.(fld('dy'))(1); r.(fld('dz'))(1)];
                    P(m, 1) = r.(fld('p'))(1); band(m) = r.(fld('band'))(1);
                    sh(:, 1, m) = r.(fld('shadow'));
                end
                src = struct('dir', d, 'P', P, 'band', band);
                [F, info] = NonGravForces.kernel.radiationPanelForce(mm, src, 'Shadow', sh);
                want = [r.fx r.fy r.fz]'; got = reshape(F, 3, []);
                scale = max(abs(want(:)));
                worst = max(worst, max(abs(got - want), [], 'all') / max(scale, realmin));
                wa = r.absorbed; ga = info.absorbed_W_m2(:);
                worstAbs = max(worstAbs, max(abs(ga - wa)) / max(max(abs(wa)), realmin));
            end
            tc.verifyBoundAndRecord("radiationPanelForce vs the mirror (" + numel(cases) + " cases, relative to the case scale)", worst, 1e-12, "-", "independently derived optical model, double precision");
            tc.verifyBoundAndRecord("radiationPanelForce absorbed flux vs the mirror", worstAbs, 1e-12, "-", "a E gamma summed over sources");
        end
        function thermalRecoilMirror(tc)
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'mirror_thermal.csv'), 'TextType', 'string');
            cases = unique(T.caseId, 'stable'); worst = 0;
            for c = cases'
                r = T(T.caseId == c, :); mm = tc.macro(char(r.body(1)));
                F = NonGravForces.kernel.thermalRecoil(mm, r.temp_k);
                want = [r.fx r.fy r.fz]'; got = reshape(F, 3, []);
                worst = max(worst, max(abs(got - want), [], 'all') / max(max(abs(want(:))), realmin));
            end
            tc.verifyBoundAndRecord("thermalRecoil vs the mirror (" + numel(cases) + " cases)", worst, 1e-13, "-", "Stefan-Boltzmann and one hemispherical factor: only rounding may differ");
        end
        function assembleMirror(tc)
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'mirror_assemble.csv'), 'TextType', 'string');
            cases = unique(T.caseId, 'stable'); wF = 0; wT = 0;
            for c = cases'
                r = T(T.caseId == c, :); mm = tc.macro(char(r.body(1)));
                fk = reshape([r.inx r.iny r.inz]', 3, [], 1);
                [F, tau] = NonGravForces.kernel.assemble(fk, mm, 'CentreOfMass', [r.rcmx(1); r.rcmy(1); r.rcmz(1)]);
                wantF = [r.Fx(1); r.Fy(1); r.Fz(1)]; wantT = [r.Tx(1); r.Ty(1); r.Tz(1)];
                wF = max(wF, max(abs(F - wantF)) / max(max(abs(wantF)), realmin));
                wT = max(wT, max(abs(tau - wantT)) / max(max(abs(wantT)), realmin));
            end
            tc.verifyBoundAndRecord("assemble force vs the mirror (" + numel(cases) + " cases)", wF, 1e-13, "-", "sums of 6 to 8 terms in two orders");
            tc.verifyBoundAndRecord("assemble torque vs the mirror", wT, 1e-12, "-", "cross products then sums, two orders");
        end
        function bodyTorquesMirror(tc)
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'mirror_torques.csv'), 'TextType', 'string');
            worst = 0; nz = 0;
            for i = 1:height(T)
                mm = tc.macro(char(T.body(i)));
                tau = NonGravForces.kernel.bodyTorques(mm, [T.rx(i); T.ry(i); T.rz(i)], ...
                    'GravityGradient', logical(T.gravityGradient(i)), 'Magnetic', logical(T.magnetic(i)), ...
                    'B_body', [T.bx(i); T.by(i); T.bz(i)], 'Dipole_Am2', [T.mx(i); T.my(i); T.mz(i)], 'GM', T.gm(i));
                want = [T.tx(i); T.ty(i); T.tz(i)];
                s = max(abs(want));
                if s > 0
                    worst = max(worst, max(abs(tau - want)) / s); nz = nz + 1;
                else
                    tc.verifyEqual(tau, zeros(3, 1), sprintf('case %d selects no torque term', T.caseId(i)));
                end
            end
            tc.verifyBoundAndRecord("bodyTorques vs the mirror (" + nz + " non-zero of " + height(T) + " cases)", worst, 1e-12, "-", "independently derived gravity-gradient and dipole forms");
        end
    end
end
