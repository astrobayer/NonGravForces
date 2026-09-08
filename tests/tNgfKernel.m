classdef tNgfKernel < NgfTestCase
%TNGFKERNEL Tests of the LP-5 kernel units: gsiCoefficients, accommodation, dragPanelForce, radiationPanelForce, thermalRecoil, assemble, bodyTorques.
%
%   Oracles (none built here): the pure-Python kinetic-theory integration
%   tools/pyxcheck/sentman_check.py (Simpson quadrature of the momentum flux
%   of a drifting Maxwellian, no erf/exp closed form), the hyperthermal and
%   static-gas limits in closed form, Doornbos' (2011) published drag and lift
%   forms as the projections of the normal/tangential pair, the Lambert sphere
%   Cr = 1 + 4/9 delta and the absorbing and specular spheres, the Stefan-
%   Boltzmann law with CODATA constants, geom.thermalState's steady state, and
%   the gravity-gradient closed form of Wertz (1978) for a diagonal inertia.
%   Every bound carries the number measured on 2026-09-07.
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
    properties
        root
        mm          % graceFO("C")
        cube        % 1 m boxWing, diagonal inertia
        oracle      % sentman_oracle.csv
        c
    end
    methods (TestClassSetup)
        function setup(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            tc.mm = NonGravForces.geom.graceFO("C");
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            tc.cube = NonGravForces.geom.boxWing([1 1 1], o, 'Name', "cube");
            tc.oracle = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'sentman_oracle.csv'));
            tc.c = NonGravForces.const();
        end
    end
    methods
        function sph = sphere(~, n, vis)
            % unit sphere of n tangent plates on a Fibonacci lattice, equal areas 4 pi / n
            k = (0:n-1)'; z = 1 - 2 * (k + 0.5) / n; phi = k * pi * (3 - sqrt(5));
            nrm = [sqrt(1 - z.^2) .* cos(phi), sqrt(1 - z.^2) .* sin(phi), z]';
            pl = repmat(struct('name', "", 'kind', "plate", 'area_m2', 4 * pi / n, 'normal', [0; 0; 1], 'centroid_m', [0; 0; 1], ...
                'absVis', vis(1), 'specVis', vis(2), 'diffVis', vis(3), 'emisIR', 0.9, 'specIR', 0.05, 'diffIR', 0.05, 'alpha', 0.85, 'group', "bus"), 1, n);
            for i = 1:n, pl(i).name = "p" + i; pl(i).normal = nrm(:, i); pl(i).centroid_m = nrm(:, i); end
            sph = NonGravForces.geom.MacroModel(pl, 'Name', "sphere");
        end
        function [rho, sp, T] = gas(tc, K)
            % an atomic-oxygen dominated thermosphere at ~490 km: number densities per species, 1xK scaling
            f = linspace(1, 1.5, K);
            sp = zeros(8, K); sp(1, :) = 2e13 * f; sp(2, :) = 8e13 * f; sp(3, :) = 4e11 * f; sp(6, :) = 5e10 * f; sp(7, :) = 1e12 * f; sp(8, :) = 1e12 * f;
            T = 950 + 50 * f;
            rho = sum(sp .* (tc.c.speciesMass(:) / tc.c.NA), 1);
        end
        function src = pointSource(~, d, P, band)
            d = d / norm(d);
            src = struct('dir', reshape(d, 3, 1, 1), 'P', P, 'band', band);
        end
        function s = scaledSource(~, src, factor)
            %SCALEDSOURCE the same source list with every pressure scaled, for a numerical derivative
            s = src; s.P = src.P * factor;
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
        function partialsContract(tc)
            [rho, sp, T] = tc.gas(3);
            vb = 7600 * [1 0.99 0.98; 0.02 0.05 0.1; -0.01 0 0.02];
            vb = vb ./ vecnorm(vb) * 7600;
            p = @(varargin) NonGravForces.kernel.partials(varargin{:});
            args = {'Rho', rho, 'Species', sp, 'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            [dF, info] = p("drag", "alpha", tc.mm, args{:});
            tc.verifySize(dF, [3 tc.mm.n 3]); tc.verifyEqual(info.what, "drag"); tc.verifyEqual(info.wrt, "alpha");
            tc.verifyEqual(info.N, tc.mm.n); tc.verifyEqual(info.K, 3); tc.verifyEqual(info.unit, "N");
            tc.verifySize(p("drag", "speed", tc.mm, args{:}), [3 tc.mm.n 3]);
            tc.verifySize(p("drag", "densityScale", tc.mm, args{:}), [3 tc.mm.n 3]);
            src = tc.pointSource([0.3; -0.4; sqrt(0.75)], 4.56e-6, "vis");
            tc.verifySize(p("radiation", "pressureScale", tc.mm, 'Sources', src), [3 tc.mm.n 1]);
            Tp = 280 + zeros(tc.mm.n, 3);
            tc.verifySize(p("thermal", "temperature", tc.mm, 'Temperature_K', Tp), [3 tc.mm.n 3]);
            tc.verifyErrorId(@() p("drag", "alpha", struct(), args{:}), 'NonGravForces:partials:inputMacro');
            tc.verifyErrorId(@() p("drag", "r", tc.mm, args{:}), 'NonGravForces:partials:unsupportedWrt');
            tc.verifyErrorId(@() p("radiation", "alpha", tc.mm, 'Sources', src), 'NonGravForces:partials:unsupportedWrt');
            tc.verifyErrorId(@() p("drag", "alpha", tc.mm), 'NonGravForces:partials:inputMissing');
            tc.verifyErrorId(@() p("radiation", "pressureScale", tc.mm), 'NonGravForces:partials:inputMissing');
            tc.verifyErrorId(@() p("thermal", "temperature", tc.mm), 'NonGravForces:partials:inputMissing');
            tc.verifyErrorId(@() p("drag", "alpha", tc.mm, 'Rho', rho, 'Vrel_body', vb(:, 1:2), 'Model', "sentman"), 'NonGravForces:partials:inputShape');
            tc.verifyError(@() p("lift", "alpha", tc.mm, args{:}), 'MATLAB:validators:mustBeMember');
        end
        function gsiCoefficientsContract(tc)
            [Cp, Ct, info] = NonGravForces.kernel.gsiCoefficients("sentman", [5 8], [0; pi/4; pi/2], 0.85, 300, 1000);
            tc.verifySize(Cp, [3 2]); tc.verifySize(Ct, [3 2]); tc.verifySize(info.CD, [3 2]); tc.verifyEqual(info.model, "sentman");
            tc.verifyGreaterThanOrEqual(Cp(:), 0); tc.verifyGreaterThanOrEqual(Ct(:), 0);
            [Cp2, Ct2] = NonGravForces.kernel.gsiCoefficients("constantCd", [], [0 pi/3 2], [], [], [], 'Cd', 2.2);
            tc.verifyEqual(Cp2, 2.2 * max(cos([0 pi/3 2]), 0) .* cos([0 pi/3 2])); tc.verifyEqual(Ct2(3), 0);
            g = @(varargin) NonGravForces.kernel.gsiCoefficients(varargin{:});
            tc.verifyErrorId(@() g("cll", 8, 0, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:unsupportedModel');
            tc.verifyErrorId(@() g("schamberg", 8, 0, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:unsupportedModel');
            tc.verifyError(@() g("dria", 8, 0, 0.85, 300, 1000), 'MATLAB:validators:mustBeMember');
            tc.verifyErrorId(@() g("constantCd", [], 0, [], [], []), 'NonGravForces:gsiCoefficients:inputCd');
            tc.verifyErrorId(@() g("constantCd", [], 0, [], [], [], 'Cd', -1), 'NonGravForces:gsiCoefficients:inputCd');
            tc.verifyErrorId(@() g("maxwell", 8, 0, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:inputSpecular');
            tc.verifyErrorId(@() g("maxwell", 8, 0, 0.85, 300, 1000, 'SpecularFraction', 1.5), 'NonGravForces:gsiCoefficients:inputSpecular');
            tc.verifyErrorId(@() g("sentman", [5 8 9], [0 1], 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:inputShape');
            tc.verifyErrorId(@() g("sentman", 0, 0, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:inputRange');
            tc.verifyErrorId(@() g("sentman", 8, -0.1, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:inputRange');
            tc.verifyErrorId(@() g("sentman", 8, 0, 1.2, 300, 1000), 'NonGravForces:gsiCoefficients:inputRange');
            tc.verifyErrorId(@() g("sentman", 8, 0, 0.85, 0, 1000), 'NonGravForces:gsiCoefficients:inputRange');
            tc.verifyErrorId(@() g("sentman", NaN, 0, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:numericNaN');
            tc.verifyErrorId(@() g("sentman", 8, NaN, 0.85, 300, 1000), 'NonGravForces:gsiCoefficients:numericNaN');
        end
        function accommodationContract(tc)
            sp = zeros(8, 4); T = 900 * ones(1, 4);
            [a, info] = NonGravForces.kernel.accommodation("fixed", sp, T, 'Alpha', 0.85);
            tc.verifyEqual(a, 0.85 * ones(1, 4)); tc.verifyEqual(info.K, 4); tc.verifyEqual(info.model, "fixed");
            tc.verifyEqual(NonGravForces.kernel.accommodation("fixed", sp, T, 'Alpha', [0.8 0.9 1 0]), [0.8 0.9 1 0]);
            f = @(varargin) NonGravForces.kernel.accommodation(varargin{:});
            tc.verifyErrorId(@() f("sesam", sp, T), 'NonGravForces:accommodation:notImplemented');
            tc.verifyErrorId(@() f("langmuir", sp, T), 'NonGravForces:accommodation:notImplemented');
            tc.verifyErrorId(@() f("fixed", sp, T), 'NonGravForces:accommodation:inputAlpha');
            tc.verifyErrorId(@() f("fixed", sp, T, 'Alpha', 1.1), 'NonGravForces:accommodation:inputAlpha');
            tc.verifyErrorId(@() f("fixed", sp, T, 'Alpha', [0.8 0.9]), 'NonGravForces:accommodation:inputAlpha');
            tc.verifyErrorId(@() f("fixed", sp(1:7, :), T, 'Alpha', 0.85), 'NonGravForces:accommodation:inputShape');
            tc.verifyErrorId(@() f("fixed", sp, T', 'Alpha', 0.85), 'NonGravForces:accommodation:inputShape');
            tc.verifyErrorId(@() f("fixed", sp, T, 'Alpha', NaN), 'NonGravForces:accommodation:numericNaN');
            tc.verifyErrorId(@() f("fixed", sp, [900 NaN 900 900], 'Alpha', 0.85), 'NonGravForces:accommodation:numericNaN');
        end
        function dragPanelForceContract(tc)
            [rho, sp, T] = tc.gas(3); v = [7600 7600 7500; 10 -30 0; -5 0 20];
            [F, info] = NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            tc.verifySize(F, [3 tc.mm.n 3]); tc.verifySize(info.gamma, [tc.mm.n 3]); tc.verifySize(info.speedRatio, [8 3]);
            tc.verifyEqual(info.model, "sentman"); tc.verifyEqual(info.rhoSpecies_kg_m3, rho, 'RelTol', 1e-15);
            Fc = NonGravForces.kernel.dragPanelForce(tc.mm, rho, [], [], v, 'Model', "constantCd", 'Cd', 2.2);
            tc.verifySize(Fc, [3 tc.mm.n 3]);
            Fm = NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, v, 'Model', "maxwell", 'Alpha', 0.85, 'Tw_K', 300, 'SpecularFraction', 0.2, 'Shadow', 0.5 * ones(tc.mm.n, 3));
            tc.verifySize(Fm, [3 tc.mm.n 3]);
            d = @(varargin) NonGravForces.kernel.dragPanelForce(varargin{:});
            tc.verifyErrorId(@() d(struct(), rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputMacro');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v, 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputModel');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v, 'Model', "cll", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputModel');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputModel');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v(:, 1:2), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputShape');
            tc.verifyErrorId(@() d(tc.mm, rho, sp(1:7, :), T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputShape');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', ones(2, 3)), 'NonGravForces:dragPanelForce:inputShape');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', 2 * ones(tc.mm.n, 3)), 'NonGravForces:dragPanelForce:inputRange');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, zeros(3, 3), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:inputRange');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', -1), 'NonGravForces:dragPanelForce:inputRange');
            tc.verifyErrorId(@() d(tc.mm, rho, sp, [T(1:2) NaN], v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300), 'NonGravForces:dragPanelForce:numericNaN');
            tc.verifyErrorId(@() d(tc.mm, rho, [], [], v, 'Model', "constantCd"), 'NonGravForces:gsiCoefficients:inputCd');
        end
        function radiationPanelForceContract(tc)
            d = [0.3; -0.4; sqrt(0.75)]; src = struct('dir', cat(2, reshape(repmat(d, 1, 2), 3, 1, 2), reshape(repmat(-d, 1, 2), 3, 1, 2)), 'P', [4.56e-6 4.5e-6; 1e-7 1.1e-7], 'band', ["vis" "ir"]);
            [F, info] = NonGravForces.kernel.radiationPanelForce(tc.mm, src);
            tc.verifySize(F, [3 tc.mm.n 2]); tc.verifyEqual(info.M, 2); tc.verifyEqual(info.K, 2); tc.verifySize(info.absorbed_W_m2, [tc.mm.n 2]);
            F2 = NonGravForces.kernel.radiationPanelForce(tc.mm, src, 'Shadow', ones(tc.mm.n, 2, 2)); tc.verifyEqual(F2, F);
            r = @(varargin) NonGravForces.kernel.radiationPanelForce(varargin{:});
            tc.verifyErrorId(@() r(struct(), src), 'NonGravForces:radiationPanelForce:inputMacro');
            tc.verifyErrorId(@() r(tc.mm, struct('dir', src.dir, 'P', src.P)), 'NonGravForces:radiationPanelForce:inputSources');
            tc.verifyErrorId(@() r(tc.mm, struct('dir', src.dir, 'P', src.P(1, :), 'band', src.band)), 'NonGravForces:radiationPanelForce:inputSources');
            tc.verifyErrorId(@() r(tc.mm, struct('dir', src.dir, 'P', src.P, 'band', ["vis" "uv"])), 'NonGravForces:radiationPanelForce:inputSources');
            tc.verifyErrorId(@() r(tc.mm, struct('dir', 2 * src.dir, 'P', src.P, 'band', src.band)), 'NonGravForces:radiationPanelForce:inputSources');
            tc.verifyErrorId(@() r(tc.mm, struct('dir', src.dir, 'P', -src.P, 'band', src.band)), 'NonGravForces:radiationPanelForce:inputRange');
            tc.verifyErrorId(@() r(tc.mm, src, 'Shadow', ones(tc.mm.n, 2)), 'NonGravForces:radiationPanelForce:inputShape');
            tc.verifyErrorId(@() r(tc.mm, src, 'Shadow', 2 * ones(tc.mm.n, 2, 2)), 'NonGravForces:radiationPanelForce:inputRange');
            bad = src; bad.P(1) = NaN;
            tc.verifyErrorId(@() r(tc.mm, bad), 'NonGravForces:radiationPanelForce:numericNaN');
        end
        function thermalRecoilContract(tc)
            T = 280 + 20 * rand(tc.mm.n, 4);
            [F, info] = NonGravForces.kernel.thermalRecoil(tc.mm, T);
            tc.verifySize(F, [3 tc.mm.n 4]); tc.verifySize(info.emitted_W_m2, [tc.mm.n 4]); tc.verifySize(info.emittedPower_W, [1 4]);
            t = @(varargin) NonGravForces.kernel.thermalRecoil(varargin{:});
            tc.verifyErrorId(@() t(struct(), T), 'NonGravForces:thermalRecoil:inputMacro');
            tc.verifyErrorId(@() t(tc.mm, T(1:3, :)), 'NonGravForces:thermalRecoil:inputShape');
            tc.verifyErrorId(@() t(tc.mm, -T), 'NonGravForces:thermalRecoil:inputRange');
            T(2, 3) = NaN; tc.verifyErrorId(@() t(tc.mm, T), 'NonGravForces:thermalRecoil:numericNaN');
        end
        function assembleContract(tc)
            Fk = randn(3, tc.mm.n, 5);
            [F, tau, info] = NonGravForces.kernel.assemble(Fk, tc.mm);
            tc.verifySize(F, [3 5]); tc.verifySize(tau, [3 5]); tc.verifySize(info.lever_m, [3 tc.mm.n]); tc.verifyEqual(info.K, 5);
            [~, tau2] = NonGravForces.kernel.assemble(Fk, tc.mm, 'CentreOfMass', [0.1; 0; 0]); tc.verifySize(tau2, [3 5]);
            a = @(varargin) NonGravForces.kernel.assemble(varargin{:});
            tc.verifyErrorId(@() a(Fk, struct()), 'NonGravForces:assemble:inputMacro');
            tc.verifyErrorId(@() a(Fk(:, 1:3, :), tc.mm), 'NonGravForces:assemble:inputShape');
            tc.verifyErrorId(@() a(Fk, tc.mm, 'CentreOfMass', [1 2 3]), 'NonGravForces:assemble:inputShape');
            Fk(1) = Inf; tc.verifyErrorId(@() a(Fk, tc.mm), 'NonGravForces:assemble:numericNaN');
        end
        function bodyTorquesContract(tc)
            r = 6.87e6 * [1 0 0.6; 0 1 0.8; 0 0 0]; B = 3e-5 * randn(3, 3); m = [0.1; -0.2; 0.3];
            [tau, info] = NonGravForces.kernel.bodyTorques(tc.mm, r, 'GravityGradient', true, 'Magnetic', true, 'B_body', B, 'Dipole_Am2', m);
            tc.verifySize(tau, [3 3]); tc.verifyEqual(tau, info.gravityGradient + info.magnetic); tc.verifyTrue(info.enabled.gravityGradient);
            [tau0, info0] = NonGravForces.kernel.bodyTorques(tc.mm, r); tc.verifyEqual(tau0, zeros(3, 3)); tc.verifyFalse(info0.enabled.magnetic);
            b = @(varargin) NonGravForces.kernel.bodyTorques(varargin{:});
            tc.verifyErrorId(@() b(struct(), r, 'GravityGradient', true), 'NonGravForces:bodyTorques:inputMacro');
            noI = NonGravForces.geom.boxWing([1 1 1], struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15])));
            tc.verifyErrorId(@() b(noI, r, 'GravityGradient', true), 'NonGravForces:bodyTorques:inputMacro');
            tc.verifyErrorId(@() b(tc.mm, r(1:2, :)), 'NonGravForces:bodyTorques:inputShape');
            tc.verifyErrorId(@() b(tc.mm, r, 'Magnetic', true), 'NonGravForces:bodyTorques:inputMagnetic');
            tc.verifyErrorId(@() b(tc.mm, r, 'Magnetic', true, 'B_body', B(:, 1:2), 'Dipole_Am2', m), 'NonGravForces:bodyTorques:inputShape');
            tc.verifyErrorId(@() b(tc.mm, zeros(3, 1), 'GravityGradient', true), 'NonGravForces:bodyTorques:inputRange');
            tc.verifyErrorId(@() b(tc.mm, r, 'GM', -1), 'NonGravForces:bodyTorques:inputRange');
            tc.verifyErrorId(@() b(tc.mm, [r(:, 1:2) [NaN; 0; 0]]), 'NonGravForces:bodyTorques:numericNaN');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function partialsReference(tc)
            % Richardson-extrapolated central differences of the force itself: an independent
            % numerical route to the same derivative. A complex step is impossible here and the
            % toolbox says so rather than pretending: MATLAB's erf and erfc raise notFullReal on
            % a complex argument and abs/vecnorm discard the imaginary part (measured 2026-09-08).
            [rho, sp, T] = tc.gas(3);
            vb = 7600 * [1 0.99 0.98; 0.02 0.05 0.1; -0.01 0 0.02]; vb = vb ./ vecnorm(vb) * 7600;
            rich = @(f, h) (4 * (f(h) - f(-h)) / (2 * h) - (f(2*h) - f(-2*h)) / (4 * h)) / 3;
            rel = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));
            args = {'Rho', rho, 'Species', sp, 'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            dA = NonGravForces.kernel.partials("drag", "alpha", tc.mm, args{:});
            fa = @(d) NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, vb, 'Model', "sentman", 'Alpha', 0.85 + d, 'Tw_K', 300);
            tc.verifyBoundAndRecord("d(drag)/d(alpha) vs Richardson central differences", rel(dA, rich(fa, 1e-4)), 1e-9, "-", "the extrapolation's own truncation dominates; measured 1e-12");
            U0 = vecnorm(vb);
            dU = NonGravForces.kernel.partials("drag", "speed", tc.mm, args{:});
            fu = @(d) NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, vb .* ((U0 + d) ./ U0), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            tc.verifyBoundAndRecord("d(drag)/d(speed) vs Richardson central differences", rel(dU, rich(fu, 1e-2)), 1e-9, "-", "measured 7e-11");
            dR = NonGravForces.kernel.partials("drag", "densityScale", tc.mm, args{:});
            fr = @(d) NonGravForces.kernel.dragPanelForce(tc.mm, rho * (1 + d), sp * (1 + d), T, vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            tc.verifyBoundAndRecord("d(drag)/d(densityScale) vs Richardson", rel(dR, rich(fr, 1e-5)), 1e-9, "-", "measured 2e-11");
            src = tc.pointSource([0.3; -0.4; sqrt(0.75)], 4.56e-6, "vis");
            dP = NonGravForces.kernel.partials("radiation", "pressureScale", tc.mm, 'Sources', src);
            fp = @(d) NonGravForces.kernel.radiationPanelForce(tc.mm, tc.scaledSource(src, 1 + d));
            tc.verifyBoundAndRecord("d(radiation)/d(pressureScale) vs Richardson", rel(dP, rich(fp, 1e-5)), 1e-9, "-", "measured 8e-11");
            Tp = 280 + 20 * cos((1:tc.mm.n)') + zeros(tc.mm.n, 3);
            dT = NonGravForces.kernel.partials("thermal", "temperature", tc.mm, 'Temperature_K', Tp);
            ft = @(d) NonGravForces.kernel.thermalRecoil(tc.mm, Tp + d);
            tc.verifyBoundAndRecord("d(thermal)/d(temperature) vs Richardson", rel(dT, rich(ft, 1e-3)), 1e-9, "-", "measured 2e-11");
        end
        function gsiCoefficientsReference(tc)
            % kinetic-theory integration (Simpson, no erf) vs the closed forms; oracle error < 1e-10
            o = tc.oracle; Tinf = 1000;
            [Cp, Ct, info] = NonGravForces.kernel.gsiCoefficients("sentman", o.s, deg2rad(o.theta_deg), o.alpha, o.twRatio * Tinf, Tinf);
            Cps = NonGravForces.kernel.gsiCoefficients("maxwell", o.s, deg2rad(o.theta_deg), o.alpha, o.twRatio * Tinf, Tinf, 'SpecularFraction', 1);
            tc.verifyBoundAndRecord("sentman Cp vs kinetic-theory oracle (105 cases)", max(abs(Cp - o.Cp)), 1e-12, "-", "oracle Simpson error ~1e-13, coefficients up to 7.8; measured 4e-14");
            tc.verifyBoundAndRecord("sentman Ctau vs oracle", max(abs(Ct - o.Ctau)), 1e-12, "-", "measured 1.2e-14");
            tc.verifyBoundAndRecord("drag projection Cp g + Ctau l vs oracle CD", max(abs(info.CD - o.CD)), 1e-12, "-", "measured 4e-14");
            tc.verifyBoundAndRecord("lift projection Cp l - Ctau g vs oracle CL", max(abs(info.CL - o.CL)), 1e-12, "-", "measured 1e-14");
            tc.verifyBoundAndRecord("specular (maxwell sigma = 1) Cp vs oracle", max(abs(Cps - o.Cp_spec)), 1e-12, "-", "measured 5e-14");
            % the drag projection equals Doornbos' (2011) published C_D form: P/sqrt(pi) + gamma Q Z + gamma/2 vre (gamma sqrt(pi) Z + P)
            g = cos(deg2rad(o.theta_deg)); P = exp(-(g .* o.s).^2) ./ o.s; Z = 1 + erf(g .* o.s); Q = 1 + 1 ./ (2 * o.s.^2);
            vr = sqrt(0.5 * (1 + o.alpha .* (2 * o.twRatio ./ o.s.^2 - 1)));
            CDdoornbos = P / sqrt(pi) + g .* Q .* Z + 0.5 * g .* vr .* (g * sqrt(pi) .* Z + P);
            tc.verifyBoundAndRecord("Cp g + Ctau l vs Doornbos C_D closed form", max(abs(info.CD - CDdoornbos)), 1e-14, "-", "algebraic identity, rounding only");
            % the predecessor's defect D1 (shear divided by s, not s^2) is not reproduced: at s = 8 the wrong form is off by ~s
            i = find(o.s == 8 & o.theta_deg == 60 & o.alpha == 0.85);
            wrong = sin(deg2rad(60)) * (1 ./ o.s(i) + g(i) .* Z(i));    % P/sqrt(pi) replaced by 1/s: the v1 grouping
            tc.verifyGreaterThan(abs(wrong - o.Ctau(i)) / o.Ctau(i), 0.05, 'the D1 regression must be distinguishable from the correct value');
        end
        function dragPanelForceReference(tc)
            % Woeske (2019) eq. 5.20 with the projected area as reference: the constantCd plate sum reproduces it
            v = [7600; 30; -20]; rho = 1e-12;
            F = NonGravForces.kernel.dragPanelForce(tc.mm, rho, [], [], v, 'Model', "constantCd", 'Cd', 2.2);
            Ft = NonGravForces.kernel.assemble(F, tc.mm);
            Aproj = sum(tc.mm.A(:) .* max(tc.mm.N' * v / norm(v), 0));
            want = -0.5 * rho * 2.2 * Aproj * norm(v) * v;
            tc.verifyBoundAndRecord("constantCd plate sum vs eq. 5.20", norm(Ft - want) / norm(want), 1e-14, "-", "8-term sum, rounding only; measured 1e-16");
            % Sentman on the GRACE-FO panels at s ~ 8: the effective drag coefficient on the projected area lies in the published 3.2-3.8 band (Doornbos 2011 Fig. 4.8 for GRACE with Sentman)
            [rho, sp, T] = tc.gas(1);
            Fs = NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            Fst = NonGravForces.kernel.assemble(Fs, tc.mm);
            cdEff = norm(Fst) / (0.5 * rho * norm(v)^2 * Aproj);
            tc.verifyGreaterThan(cdEff, 3.2); tc.verifyLessThan(cdEff, 3.8);
            tc.verifyBoundAndRecord("GRACE-FO Sentman Cd on the projected area (deviation from 3.5)", abs(cdEff - 3.5), 0.3, "-", "published band 3.2-3.8 at s = 7.9; measured 3.59");
        end
        function radiationPanelForceReference(tc)
            % sphere of tangent plates: absorbing Cr = 1, Lambert Cr = 1 + 4/9 delta, specular Cr = 1 (standard results)
            n = 4000; d = [0.3; -0.4; sqrt(0.75)]; P = 4.56e-6; F0 = P * pi;
            cr = @(vis) -dot(NonGravForces.kernel.assemble(NonGravForces.kernel.radiationPanelForce(tc.sphere(n, vis), tc.pointSource(d, P, "vis")), tc.sphere(n, vis)), d) / F0;
            tc.verifyBoundAndRecord("absorbing sphere Cr - 1", abs(cr([1 0 0]) - 1), 1e-4, "-", "Fibonacci lattice of 4000 plates, discretisation ~1/n; measured 3.8e-5");
            tc.verifyBoundAndRecord("Lambert sphere Cr - (1 + 4/9)", abs(cr([0 0 1]) - (1 + 4 / 9)), 1e-4, "-", "measured 3.8e-5");
            tc.verifyBoundAndRecord("specular sphere Cr - 1", abs(cr([0 1 0]) - 1), 1e-4, "-", "measured 1e-6");
        end
        function thermalRecoilReference(tc)
            % Stefan-Boltzmann with CODATA 2018: a black 1 m^2 plate at 300 K emits 459.30 W and recoils with (2/3) P / c
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [1 0 0]));
            bw = NonGravForces.geom.boxWing([1 1 1], o);
            [F, info] = NonGravForces.kernel.thermalRecoil(bw, 300 * ones(6, 1));
            tc.verifyBoundAndRecord("emitted flux of a black plate at 300 K vs sigma T^4", abs(info.emitted_W_m2(1) - 5.670374419e-8 * 300^4), 1e-9, "W/m^2", "CODATA 2018 sigma, exact arithmetic");
            tc.verifyBoundAndRecord("recoil of the +x plate vs (2/3) sigma T^4 A / c", abs(norm(F(:, 1, 1)) - (2 / 3) * 459.2999 / 299792458), 1e-11, "N", "published 459.3 W/m^2 to 4 digits");
            tc.verifyBoundAndRecord("cube at uniform temperature: net recoil", norm(sum(F(:, :, 1), 2)), 1e-20, "N", "opposite faces cancel exactly");
        end
        function bodyTorquesReference(tc)
            % Wertz (1978) eq. 17-30: tau_x = 3 GM / r^3 (I3 - I2) sin(phi) cos(phi) for r in the y-z plane and a diagonal inertia
            I = diag([100 500 700]);
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            bw = NonGravForces.geom.boxWing([1 1 1], o);
            pl = bw.panels; mmI = NonGravForces.geom.MacroModel(pl, 'Name', "diag", 'Inertia', I);
            GM = referenceSystem.iersConstants().gmEarth; rr = 6.87e6; phi = 0.3;
            r = rr * [0; sin(phi); cos(phi)];
            tau = NonGravForces.kernel.bodyTorques(mmI, r, 'GravityGradient', true);
            want = [3 * GM / rr^3 * (700 - 500) * sin(phi) * cos(phi); 0; 0];
            tc.verifyBoundAndRecord("gravity-gradient torque vs Wertz eq. 17-30", norm(tau - want) / norm(want), 1e-14, "-", "closed form, rounding only");
            tc.verifyBoundAndRecord("gravity-gradient torque magnitude at 490 km", abs(norm(want) - 2.0824e-4) / 2.0824e-4, 1e-3, "-", "hand calculation: 3 GM/r^3 = 3.688e-6 s^-2, x 200 kg m^2 x sin(0.3) cos(0.3) = 2.0824e-4 N m");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function partialsPrecision(tc)
            % The linear parameters have exact derivatives: the force itself, bitwise.
            [rho, sp, T] = tc.gas(2);
            vb = 7600 * [1 0.99; 0.02 0.05; -0.01 0]; vb = vb ./ vecnorm(vb) * 7600;
            args = {'Rho', rho, 'Species', sp, 'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            F = NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            dR = NonGravForces.kernel.partials("drag", "densityScale", tc.mm, args{:});
            tc.verifyBoundAndRecord("d(drag)/d(densityScale) minus the force", max(abs(dR - F), [], 'all'), 0, "N", "the force is linear in the density: the derivative IS the force, bitwise");
            src = tc.pointSource([0; 0; 1], 4.56e-6, "vis");
            Fr = NonGravForces.kernel.radiationPanelForce(tc.mm, src);
            dP = NonGravForces.kernel.partials("radiation", "pressureScale", tc.mm, 'Sources', src);
            tc.verifyBoundAndRecord("d(radiation)/d(pressureScale) minus the force", max(abs(dP - Fr), [], 'all'), 0, "N", "as above, bitwise");
            Tp = 300 + zeros(tc.mm.n, 2);
            Ft = NonGravForces.kernel.thermalRecoil(tc.mm, Tp);
            dT = NonGravForces.kernel.partials("thermal", "temperature", tc.mm, 'Temperature_K', Tp);
            tc.verifyBoundAndRecord("d(thermal)/d(T) minus 4F/T", max(abs(dT - 4 * Ft / 300), [], 'all'), 1e-30, "N/K", "a T^4 law differentiated once");
        end
        function gsiCoefficientsPrecision(tc)
            % hyperthermal limit at s = 1e4: Cp -> 2 g^2, Ctau -> 2 g l up to the re-emission sqrt(pi) sqrt(Tw/(Tinf s^2)) g and 1/(2 s^2)
            th = linspace(0, pi / 2, 9); s = 1e4; Tw = 300; Tinf = 1000;
            [Cp, Ct] = NonGravForces.kernel.gsiCoefficients("sentman", s, th, 1, Tw, Tinf);
            bound = (sqrt(pi) * sqrt(Tw / (Tinf * s^2)) + 1 / s^2) * (1 + 1e-9);
            tc.verifyBoundAndRecord("hyperthermal Cp - 2 g^2", max(abs(Cp - 2 * cos(th).^2)), bound, "-", "re-emission term sqrt(pi) sqrt(Tw/Tinf)/s = 9.71e-5 plus 1/s^2, exact to 1e-9 after the alpha bracket rewrite (was 4e-9 off); measured 9.709e-5");
            tc.verifyBoundAndRecord("hyperthermal Ctau - 2 g l", max(abs(Ct - 2 * cos(th) .* sin(th))), 1e-4, "-", "P/sqrt(pi) l = 5.6e-5 at s = 1e4; measured 5.64e-5");
            % static gas (s -> 0), head-on, Tw = Tinf: the plate feels the full gas pressure n k T = rho cbar^2 / 2, i.e. Cp (1/2) s^2 -> 1/2
            s0 = 1e-3; Cp0 = NonGravForces.kernel.gsiCoefficients("sentman", s0, 0, 1, 1000, 1000);
            tc.verifyBoundAndRecord("static-gas pressure limit Cp s^2/2 - 1/2", abs(Cp0 * 0.5 * s0^2 - 0.5), 2 * s0, "-", "first-order term in s (bulk drift), oracle 0.501 at s = 1e-3");
            % head-on at s = 8, alpha = 1, Tw = 300, Tinf = 1000: 2 + 1/64 + sqrt(pi) sqrt(0.6/64)/... = 2.136977 (help example)
            Cp8 = NonGravForces.kernel.gsiCoefficients("sentman", 8, 0, 1, 300, 1000);
            tc.verifyBoundAndRecord("head-on s = 8 example value", abs(Cp8 - (2 + 1 / 64 + sqrt(pi) * sqrt(0.5 * 0.6 / 64))), 1e-15, "-", "closed form in the help text");
        end
        function dragPanelForcePrecision(tc)
            % one plate, one species: the force equals (1/2) rho U^2 A [Cp (-n) + Ctau t] with the plate's own theta
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
            bw = NonGravForces.geom.boxWing([1 1 1], o);                % +x plate is panel 1
            sp = zeros(8, 1); sp(2) = 1e14; T = 1000; U = 7600; th = 0.4;
            v = U * [cos(th); sin(th); 0];
            m = tc.c.speciesMass(2) / tc.c.NA; rho = sp(2) * m; s = U / sqrt(2 * tc.c.kB * T / m);
            [F, info] = NonGravForces.kernel.dragPanelForce(bw, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            [Cp, Ct] = NonGravForces.kernel.gsiCoefficients("sentman", s, th, 0.85, 300, T);
            tdir = (-v / U + cos(th) * [1; 0; 0]) / sin(th);
            want = 0.5 * rho * U^2 * 1 * (-Cp * [1; 0; 0] + Ct * tdir);
            tc.verifyBoundAndRecord("+x plate force vs closed form", norm(F(:, 1, 1) - want) / norm(want), 1e-14, "-", "rounding of a 3-term product; measured 2e-16");
            sWant = 7600 / sqrt(2 * 1.380649e-23 * 1000 / (15.9994e-3 / 6.02214076e23));   % CODATA 2018, restated here
            tc.verifyBoundAndRecord("speed ratio of atomic oxygen at 1000 K, 7.6 km/s", abs(info.speedRatio(2) - sWant) / sWant, 1e-12, "-", "U / sqrt(2 k T / m_O) = 7.4548 from CODATA constants; rounding");
            % an edge-on plate (gamma = 0) still feels the thermal pressure 1/(2 s^2) and the re-emission, plus the shear P/sqrt(pi) along the stream -v
            F3 = F(:, 5, 1);                                            % +z plate: n = [0 0 1], v in the x-y plane
            [Cp0, Ct0] = NonGravForces.kernel.gsiCoefficients("sentman", s, pi / 2, 0.85, 300, T);
            want3 = 0.5 * rho * U^2 * (-Cp0 * [0; 0; 1] + Ct0 * (-v / U));
            tc.verifyBoundAndRecord("edge-on plate force vs closed form", norm(F3 - want3) / norm(want3), 1e-14, "-", "gamma = 0 exactly: Cp = 1/(2 s^2) + vre P / 2, Ctau = P / sqrt(pi); rounding");
            tc.verifyBoundAndRecord("edge-on plate: normal over tangential force", abs(Cp0 / Ct0 - (1 / (2 * s^2) + 0.5 * sqrt(0.5 * (0.15 + 0.85 * 2 * 300 / (T * s^2))) / s) / (1 / (s * sqrt(pi)))), 1e-12, "-", "the two closed forms at theta = pi/2");
        end
        function radiationPanelForcePrecision(tc)
            % single plate closed forms: absorber -P A g d, mirror -2 P A g^2 n, Lambert -P A g (d + 2/3 n)
            d = [0.6; 0; 0.8]; P = 4.56e-6; g = 0.8;
            mk = @(vis) NonGravForces.geom.boxWing([1 1 1], struct('bus', struct('vis', vis, 'ir', [0.62 0.23 0.15])));
            src = tc.pointSource(d, P, "vis");
            Fa = NonGravForces.kernel.radiationPanelForce(mk([1 0 0]), src); Fs = NonGravForces.kernel.radiationPanelForce(mk([0 1 0]), src); Fd = NonGravForces.kernel.radiationPanelForce(mk([0 0 1]), src);
            n = [0; 0; 1];                                              % +z plate is panel 5 of boxWing
            tc.verifyBoundAndRecord("absorbing plate", norm(Fa(:, 5) + P * g * d), 1e-21, "N", "exact expression, rounding");
            tc.verifyBoundAndRecord("specular plate", norm(Fs(:, 5) + 2 * P * g^2 * n), 1e-21, "N", "exact expression, rounding");
            tc.verifyBoundAndRecord("Lambert plate", norm(Fd(:, 5) + P * g * (d + 2 / 3 * n)), 1e-21, "N", "exact expression, rounding");
            tc.verifyBoundAndRecord("unlit -z plate", norm(Fa(:, 6)), 0, "N", "gamma < 0 contributes exactly nothing");
            % absorbed flux reported for thermalState: a E gamma
            [~, info] = NonGravForces.kernel.radiationPanelForce(mk([1 0 0]), src);
            tc.verifyBoundAndRecord("absorbed flux a E gamma", abs(info.absorbed_W_m2(5) - P * 299792458 * g), 1e-9, "W/m^2", "E = P c; rounding");
        end
        function thermalRecoilPrecision(tc)
            % steady-state thermalState + thermalRecoil: recoil = (2/3) absorbed power / c on every lit plate
            src = tc.pointSource([1; 0.2; 0.1], 4.56e-6, "vis");
            T = NonGravForces.geom.thermalState(tc.mm, src, 'Model', "steadyState");
            [F, info] = NonGravForces.kernel.thermalRecoil(tc.mm, T);
            [~, rinfo] = NonGravForces.kernel.radiationPanelForce(tc.mm, src);
            lit = rinfo.absorbed_W_m2(:, 1) > 0;
            tc.verifyBoundAndRecord("steady state: emitted equals absorbed on lit plates", max(abs(info.emitted_W_m2(lit) ./ rinfo.absorbed_W_m2(lit) - 1)), 1e-12, "-", "T = (q/(eps sigma))^(1/4) then eps sigma T^4: two roundings");
            mag = vecnorm(reshape(F(:, :, 1), 3, []))';
            tc.verifyBoundAndRecord("recoil = (2/3) absorbed A / c", max(abs(mag(lit) ./ ((2 / 3) * rinfo.absorbed_W_m2(lit) .* tc.mm.A(lit)' / 299792458) - 1)), 1e-12, "-", "unit normals (O-2 closed 2026-09-07: N renormalised); was 7e-7 with printed normals");
        end
        function assemblePrecision(tc)
            % dumbbell: equal and opposite forces on the +x and -x plates of the cube -> F = 0, tau = 2 L f (lever 0.5 m each)
            Fk = zeros(3, 6); Fk(:, 1) = [0; 1; 0]; Fk(:, 2) = [0; -1; 0];
            [F, tau] = NonGravForces.kernel.assemble(Fk, tc.cube);
            tc.verifyEqual(F, zeros(3, 1)); tc.verifyEqual(tau, [0; 0; 1]);
            % centroid of the +x plate is [0.5 0 0]; force [0 0 1] there gives tau = r x F = [0 -0.5 0]
            Fk = zeros(3, 6); Fk(:, 1) = [0; 0; 1];
            [~, tau1] = NonGravForces.kernel.assemble(Fk, tc.cube); tc.verifyEqual(tau1, [0; -0.5; 0]);
            [~, tau2] = NonGravForces.kernel.assemble(Fk, tc.cube, 'CentreOfMass', [0.5; 0; 0]); tc.verifyEqual(tau2, zeros(3, 1));
        end
        function bodyTorquesPrecision(tc)
            % magnetic torque m x B in closed form; a principal axis gives zero gravity-gradient torque
            B = [1e-5; -2e-5; 3e-5]; m = [0.3; 0.2; -0.1];
            [tau, info] = NonGravForces.kernel.bodyTorques(tc.cube, [7e6; 0; 0], 'Magnetic', true, 'B_body', B, 'Dipole_Am2', m);
            tc.verifyEqual(tau, cross(m, B)); tc.verifyEqual(info.gravityGradient, zeros(3, 1));
            I = diag([100 500 700]); mmI = NonGravForces.geom.MacroModel(tc.cube.panels, 'Name', "diag", 'Inertia', I);
            tg = NonGravForces.kernel.bodyTorques(mmI, 6.87e6 * eye(3), 'GravityGradient', true);
            tc.verifyEqual(tg, zeros(3, 3), 'principal axes give no gravity-gradient torque');
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function partialsVectorization(tc)
            [rho, sp, T] = tc.gas(4);
            vb = 7600 * [1 0.99 0.98 0.97; 0.02 0.05 0.1 0.15; -0.01 0 0.02 0.03];
            vb = vb ./ vecnorm(vb) * 7600;
            for wrt = ["alpha", "speed", "densityScale"]
                dAll = NonGravForces.kernel.partials("drag", wrt, tc.mm, 'Rho', rho, 'Species', sp, ...
                    'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
                for k = 1:4
                    d1 = NonGravForces.kernel.partials("drag", wrt, tc.mm, 'Rho', rho(k), 'Species', sp(:, k), ...
                        'T_K', T(k), 'Vrel_body', vb(:, k), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
                    tc.verifyEqual(dAll(:, :, k), d1, sprintf('%s, epoch %d batched vs alone', wrt, k));
                end
            end
        end
        function gsiCoefficientsVectorization(tc)
            s = [3 8 12]; th = [0.1; 1.0; 2.0; 3.0]; a = [0.7 0.85 1];
            [Cp, Ct, info] = NonGravForces.kernel.gsiCoefficients("sentman", s, th, a, 300, [900 1000 1100]);
            for i = 1:4
                for j = 1:3
                    [cp1, ct1, i1] = NonGravForces.kernel.gsiCoefficients("sentman", s(j), th(i), a(j), 300, 900 + 100 * (j - 1));
                    tc.verifyEqual(Cp(i, j), cp1); tc.verifyEqual(Ct(i, j), ct1); tc.verifyEqual(info.CD(i, j), i1.CD);
                end
            end
        end
        function accommodationVectorization(tc)
            sp = rand(8, 5) * 1e13; T = 900 + 10 * (1:5); a = [0.8 0.85 0.9 0.95 1];
            al = NonGravForces.kernel.accommodation("fixed", sp, T, 'Alpha', a);
            for k = 1:5, tc.verifyEqual(al(k), NonGravForces.kernel.accommodation("fixed", sp(:, k), T(k), 'Alpha', a(k))); end
        end
        function dragPanelForceVectorization(tc)
            [rho, sp, T] = tc.gas(6); v = 7600 * [ones(1, 6); 0.01 * (1:6); -0.005 * (1:6)]; f = 0.5 + 0.5 * rand(tc.mm.n, 6);
            F = NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', f);
            Fc = NonGravForces.kernel.dragPanelForce(tc.mm, rho, [], [], v, 'Model', "constantCd", 'Cd', 2.3);
            for k = 1:6
                tc.verifyEqual(F(:, :, k), NonGravForces.kernel.dragPanelForce(tc.mm, rho(k), sp(:, k), T(k), v(:, k), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', f(:, k)));
                tc.verifyEqual(Fc(:, :, k), NonGravForces.kernel.dragPanelForce(tc.mm, rho(k), [], [], v(:, k), 'Model', "constantCd", 'Cd', 2.3));
            end
        end
        function radiationPanelForceVectorization(tc)
            K = 5; d1 = randn(3, K); d1 = d1 ./ vecnorm(d1); d2 = randn(3, K); d2 = d2 ./ vecnorm(d2);
            src = struct('dir', cat(2, reshape(d1, 3, 1, K), reshape(d2, 3, 1, K)), 'P', [4.5e-6 * ones(1, K); 1e-7 * (1:K)], 'band', ["vis" "ir"]);
            f = rand(tc.mm.n, K, 2);
            F = NonGravForces.kernel.radiationPanelForce(tc.mm, src, 'Shadow', f);
            for k = 1:K
                s1 = struct('dir', src.dir(:, :, k), 'P', src.P(:, k), 'band', src.band);
                tc.verifyEqual(F(:, :, k), NonGravForces.kernel.radiationPanelForce(tc.mm, s1, 'Shadow', f(:, k, :)));
            end
        end
        function thermalRecoilVectorization(tc)
            T = 250 + 60 * rand(tc.mm.n, 7);
            F = NonGravForces.kernel.thermalRecoil(tc.mm, T);
            for k = 1:7, tc.verifyEqual(F(:, :, k), NonGravForces.kernel.thermalRecoil(tc.mm, T(:, k))); end
        end
        function assembleVectorization(tc)
            Fk = 1e-6 * randn(3, tc.mm.n, 9);
            [F, tau] = NonGravForces.kernel.assemble(Fk, tc.mm);
            for k = 1:9, [f1, t1] = NonGravForces.kernel.assemble(Fk(:, :, k), tc.mm); tc.verifyEqual(F(:, k), f1); tc.verifyEqual(tau(:, k), t1); end
        end
        function bodyTorquesVectorization(tc)
            r = 6.9e6 * (randn(3, 6) ./ vecnorm(randn(3, 6))); B = 3e-5 * randn(3, 6); m = [0.1; -0.2; 0.3];
            tau = NonGravForces.kernel.bodyTorques(tc.mm, r, 'GravityGradient', true, 'Magnetic', true, 'B_body', B, 'Dipole_Am2', m);
            for k = 1:6, tc.verifyEqual(tau(:, k), NonGravForces.kernel.bodyTorques(tc.mm, r(:, k), 'GravityGradient', true, 'Magnetic', true, 'B_body', B(:, k), 'Dipole_Am2', m)); end
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function partialsRobustness(tc)
            [rho, sp, T] = tc.gas(2);
            vb = 7600 * [1 0.99; 0.02 0.05; -0.01 0]; vb = vb ./ vecnorm(vb) * 7600;
            % constantCd has no accommodation coefficient, so that derivative is exactly zero
            dA = NonGravForces.kernel.partials("drag", "alpha", tc.mm, 'Rho', rho, 'Species', sp, 'T_K', T, ...
                'Vrel_body', vb, 'Model', "constantCd", 'Cd', 2.2);
            tc.verifyEqual(dA, zeros(3, tc.mm.n, 2), 'constantCd carries no alpha');
            % and its speed derivative is exactly 2F/U, the dynamic pressure alone
            F = NonGravForces.kernel.dragPanelForce(tc.mm, rho, [], [], vb, 'Model', "constantCd", 'Cd', 2.2);
            dU = NonGravForces.kernel.partials("drag", "speed", tc.mm, 'Rho', rho, 'Species', sp, 'T_K', T, ...
                'Vrel_body', vb, 'Model', "constantCd", 'Cd', 2.2);
            tc.verifyBoundAndRecord("constantCd d/dU minus 2F/U", max(abs(dU - 2 * F ./ reshape(vecnorm(vb), 1, 1, 2)), [], 'all'), 1e-30, "N/(m/s)", "no coefficient depends on the speed there");
            % an empty species column must not produce NaN
            spz = sp; spz(:, 2) = 0;
            d0 = NonGravForces.kernel.partials("drag", "alpha", tc.mm, 'Rho', [rho(1) 0], 'Species', spz, 'T_K', T, ...
                'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            tc.verifyTrue(all(isfinite(d0(:))));
            % a panel at 0 K emits nothing and its derivative is finite
            dz = NonGravForces.kernel.partials("thermal", "temperature", tc.mm, 'Temperature_K', zeros(tc.mm.n, 1));
            tc.verifyTrue(all(isfinite(dz(:))));
        end
        function gsiCoefficientsRobustness(tc)
            % leeward plates: exponentially small, finite, non-negative; theta = pi exactly; extreme s finite
            [Cp, Ct] = NonGravForces.kernel.gsiCoefficients("sentman", [0.01 8 1e6], [pi/2 + 0.1; 2.5; pi], 0.85, 300, 1000);
            tc.verifyTrue(all(isfinite([Cp(:); Ct(:)]))); tc.verifyGreaterThanOrEqual([Cp(:); Ct(:)], 0, 'erfc keeps leeward coefficients non-negative (1 + erf cancelled below zero, measured)');
            tc.verifyLessThan(Cp(3, 2), 1e-20, 'leeward at s = 8 is exponentially small');
            tc.verifyBoundAndRecord("theta = pi: Ctau at the sin(pi) rounding level", max(Ct(3, :)), 1e-14, "-", "sin(pi) = 1.2e-16 times P/sqrt(pi) = 56 at s = 0.01");
            Cd0 = NonGravForces.kernel.gsiCoefficients("constantCd", [], [2 3 pi], [], [], [], 'Cd', 2.2); tc.verifyEqual(Cd0, zeros(1, 3), 'constantCd on leeward plates is exactly zero');
            [CpE, CtE] = NonGravForces.kernel.gsiCoefficients("sentman", zeros(1, 0), zeros(1, 0), 0.85, 300, 1000); tc.verifySize(CpE, [1 0]); tc.verifySize(CtE, [1 0]);
        end
        function accommodationRobustness(tc)
            a = NonGravForces.kernel.accommodation("fixed", zeros(8, 0), zeros(1, 0), 'Alpha', 0.85); tc.verifySize(a, [1 0]);
            a1 = NonGravForces.kernel.accommodation("fixed", zeros(8, 3), [1 1 1], 'Alpha', [0 1 0.5]); tc.verifyEqual(a1, [0 1 0.5], 'boundary values 0 and 1 are admissible');
        end
        function dragPanelForceRobustness(tc)
            % a species column of zeros is skipped (no NaN from s = U/0), an exactly edge-on and an exactly head-on plate give finite forces
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15])); bw = NonGravForces.geom.boxWing([1 1 1], o);
            sp = zeros(8, 2); sp(2, :) = [1e14 0]; sp(3, :) = [0 1e12]; T = [1000 1000]; v = [7600 7600; 0 0; 0 0];
            rho = sum(sp .* (tc.c.speciesMass(:) / tc.c.NA), 1);
            [F, info] = NonGravForces.kernel.dragPanelForce(bw, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            tc.verifyTrue(all(isfinite(F(:)))); tc.verifyEqual(info.gamma(1, :), [1 1]); tc.verifyEqual(info.gamma(3, :), [0 0]);
            tc.verifyEqual(F(2:3, 1, :), zeros(2, 1, 2), 'head-on plate: normal force only');
            tc.verifyBoundAndRecord("leeward -x plate force at s = 7.45", abs(F(1, 2, 1)), 1e-28, "N", "exp(-s^2) = 1e-24 times (1/2) rho U^2 A = 8e-5; measured 3.2e-30");
            FK0 = NonGravForces.kernel.dragPanelForce(bw, zeros(1, 0), zeros(8, 0), zeros(1, 0), zeros(3, 0), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300); tc.verifySize(FK0, [3 6 0]);
        end
        function radiationPanelForceRobustness(tc)
            % grazing incidence (gamma = 0) and zero pressure give exact zeros; a source exactly behind contributes nothing
            src = struct('dir', reshape([1 0 -1; 0 0 0; 0 1 0], 3, 3, 1), 'P', [4.56e-6; 0; 4.56e-6], 'band', ["vis" "vis" "ir"]);
            [F, info] = NonGravForces.kernel.radiationPanelForce(tc.cube, src);
            tc.verifyEqual(F(:, 5, 1), zeros(3, 1), '+z plate is edge-on to source 1 and 3 and source 2 has P = 0');
            tc.verifyEqual(info.litCount(1), 1); tc.verifyEqual(info.litCount(2), 1);
            F0 = NonGravForces.kernel.radiationPanelForce(tc.cube, struct('dir', zeros(3, 1, 0), 'P', zeros(1, 0), 'band', "vis")); tc.verifySize(F0, [3 6 0]);
        end
        function thermalRecoilRobustness(tc)
            F = NonGravForces.kernel.thermalRecoil(tc.cube, [1e-3 * ones(6, 1) 1e4 * ones(6, 1)]);
            tc.verifyTrue(all(isfinite(F(:)))); tc.verifyEqual(sum(F(:, :, 2), 2), zeros(3, 1), 'AbsTol', 1e-20);
            F0 = NonGravForces.kernel.thermalRecoil(tc.cube, zeros(6, 0)); tc.verifySize(F0, [3 6 0]);
        end
        function assembleRobustness(tc)
            [F, tau] = NonGravForces.kernel.assemble(zeros(3, tc.mm.n, 0), tc.mm); tc.verifySize(F, [3 0]); tc.verifySize(tau, [3 0]);
            [F1, tau1] = NonGravForces.kernel.assemble(zeros(3, tc.mm.n), tc.mm); tc.verifyEqual(F1, zeros(3, 1)); tc.verifyEqual(tau1, zeros(3, 1));
        end
        function bodyTorquesRobustness(tc)
            tau = NonGravForces.kernel.bodyTorques(tc.mm, zeros(3, 0), 'GravityGradient', true); tc.verifySize(tau, [3 0]);
            tau1 = NonGravForces.kernel.bodyTorques(tc.mm, [7e6; 0; 0], 'Magnetic', true, 'B_body', zeros(3, 1), 'Dipole_Am2', [1; 1; 1]); tc.verifyEqual(tau1, zeros(3, 1));
            tau2 = NonGravForces.kernel.bodyTorques(tc.mm, [7e6; 0; 0], 'GravityGradient', true, 'GM', 1); tc.verifyTrue(all(isfinite(tau2)));
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function partialsMetamorphic(tc)
            [rho, sp, T] = tc.gas(2);
            vb = 7600 * [1 0.99; 0.02 0.05; -0.01 0]; vb = vb ./ vecnorm(vb) * 7600;
            a1 = {'Rho', rho, 'Species', sp, 'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            a2 = {'Rho', 2 * rho, 'Species', 2 * sp, 'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            for wrt = ["alpha", "speed", "densityScale"]
                d1 = NonGravForces.kernel.partials("drag", wrt, tc.mm, a1{:});
                d2 = NonGravForces.kernel.partials("drag", wrt, tc.mm, a2{:});
                tc.verifyBoundAndRecord("d(drag)/d(" + wrt + ") doubles with the density", ...
                    max(abs(d2 - 2 * d1), [], 'all') / max(abs(d1), [], 'all'), 1e-15, "-", "an exact power of two through a linear factor");
            end
            f = 0.5 * ones(tc.mm.n, 2);
            ds = NonGravForces.kernel.partials("drag", "alpha", tc.mm, a1{:}, 'Shadow', f);
            d0 = NonGravForces.kernel.partials("drag", "alpha", tc.mm, a1{:});
            tc.verifyBoundAndRecord("half the exposure halves d/d(alpha)", max(abs(ds - 0.5 * d0), [], 'all') / max(abs(d0), [], 'all'), 1e-15, "-", "the shadow multiplies the same scale");
        end
        function gsiCoefficientsMetamorphic(tc)
            % maxwell with sigma = 0 is sentman bitwise; sigma = 1 has no re-emission (alpha has no effect); the drag/lift projections invert
            s = [2 8]; th = [0.2; 1.1; 2.3]; [Cp, Ct, info] = NonGravForces.kernel.gsiCoefficients("sentman", s, th, 0.85, 300, 1000);
            [Cm, Tm] = NonGravForces.kernel.gsiCoefficients("maxwell", s, th, 0.85, 300, 1000, 'SpecularFraction', 0); tc.verifyEqual(Cm, Cp); tc.verifyEqual(Tm, Ct);
            C1 = NonGravForces.kernel.gsiCoefficients("maxwell", s, th, 0.85, 300, 1000, 'SpecularFraction', 1); C2 = NonGravForces.kernel.gsiCoefficients("maxwell", s, th, 0.1, 900, 1000, 'SpecularFraction', 1);
            tc.verifyEqual(C1, C2, 'fully specular: the wall temperature and alpha drop out');
            g = cos(th); l = sin(th);
            tc.verifyBoundAndRecord("normal/tangential from drag/lift round trip", max(abs(info.CD .* g + info.CL .* l - Cp), [], 'all'), 1e-15, "-", "rotation of a 2-vector, rounding");
            % the coefficients do not depend on the mass scaling of Tw and Tinf separately, only on Tw/Tinf and s
            Ca = NonGravForces.kernel.gsiCoefficients("sentman", 8, 0.5, 0.85, 300, 1000); Cb = NonGravForces.kernel.gsiCoefficients("sentman", 8, 0.5, 0.85, 600, 2000);
            tc.verifyEqual(Ca, Cb, 'only the ratio Tw/Tinf enters');
        end
        function accommodationMetamorphic(tc)
            sp = zeros(8, 4); T = 900 * ones(1, 4);
            tc.verifyEqual(NonGravForces.kernel.accommodation("fixed", sp, T, 'Alpha', 0.85), NonGravForces.kernel.accommodation("fixed", sp, T, 'Alpha', 0.85 * ones(1, 4)), 'scalar and expanded Alpha agree bitwise');
        end
        function dragPanelForceMetamorphic(tc)
            % density scaling by 2 doubles the force exactly; reversing v on the symmetric cube reverses the force; shadow 0 kills a plate
            o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15])); bw = NonGravForces.geom.boxWing([1 1 1], o);
            [rho, sp, T] = tc.gas(2); v = [7600 7600; 300 -300; -100 100];
            F1 = NonGravForces.kernel.dragPanelForce(bw, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            F2 = NonGravForces.kernel.dragPanelForce(bw, 2 * rho, 2 * sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            tc.verifyEqual(F2, 2 * F1, 'linear in the density, exact for a power of two');
            Fr = NonGravForces.kernel.dragPanelForce(bw, rho, sp, T, -v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            [Ft, ~] = NonGravForces.kernel.assemble(F1, bw); [Frt, ~] = NonGravForces.kernel.assemble(Fr, bw);
            tc.verifyBoundAndRecord("cube: F(-v) = -F(v)", norm(Frt + Ft) / norm(Ft), 1e-14, "-", "opposite plates swap roles; rounding");
            f = ones(6, 2); f(1, :) = 0;
            Fs = NonGravForces.kernel.dragPanelForce(bw, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300, 'Shadow', f);
            tc.verifyEqual(Fs(:, 1, :), zeros(3, 1, 2)); tc.verifyEqual(Fs(:, 2:6, :), F1(:, 2:6, :));
        end
        function radiationPanelForceMetamorphic(tc)
            % superposition: two sources together equal the sum of each alone; a duplicated source doubles the force exactly
            d1 = [0.6; 0; 0.8]; d2 = [-0.3; 0.5; sqrt(1 - 0.34)];
            s1 = tc.pointSource(d1, 4.56e-6, "vis"); s2 = tc.pointSource(d2, 1.2e-7, "ir");
            both = struct('dir', cat(2, s1.dir, s2.dir), 'P', [s1.P; s2.P], 'band', ["vis" "ir"]);
            F1 = NonGravForces.kernel.radiationPanelForce(tc.mm, s1); F2 = NonGravForces.kernel.radiationPanelForce(tc.mm, s2); Fb = NonGravForces.kernel.radiationPanelForce(tc.mm, both);
            tc.verifyBoundAndRecord("superposition of two sources", max(abs(Fb - (F1 + F2)), [], 'all') / max(abs(Fb), [], 'all'), 1e-15, "-", "one addition per source, rounding");
            twice = struct('dir', cat(2, s1.dir, s1.dir), 'P', [s1.P; s1.P], 'band', ["vis" "vis"]);
            tc.verifyEqual(NonGravForces.kernel.radiationPanelForce(tc.mm, twice), 2 * F1, 'a duplicated source doubles exactly');
        end
        function thermalRecoilMetamorphic(tc)
            T = 250 + 60 * rand(tc.mm.n, 3);
            F1 = NonGravForces.kernel.thermalRecoil(tc.mm, T); F2 = NonGravForces.kernel.thermalRecoil(tc.mm, 2 * T);
            tc.verifyEqual(F2, 16 * F1, 'T^4 scaling, exact for a power of two');
            % the net recoil of a body at uniform temperature is the emissivity-weighted silhouette defect: zero for the closed cube
            Fc = NonGravForces.kernel.thermalRecoil(tc.cube, 300 * ones(6, 1)); tc.verifyEqual(sum(Fc, 2), zeros(3, 1), 'AbsTol', 1e-22);
        end
        function assembleMetamorphic(tc)
            % concatenating panel arrays equals adding them; moving the reference point changes the torque by -delta x F
            Fa = 1e-6 * randn(3, tc.mm.n, 4); Fb = 1e-6 * randn(3, tc.mm.n, 4);
            [Fs, taus] = NonGravForces.kernel.assemble(Fa + Fb, tc.mm); [F1, tau1] = NonGravForces.kernel.assemble(Fa, tc.mm); [F2, tau2] = NonGravForces.kernel.assemble(Fb, tc.mm);
            tc.verifyBoundAndRecord("assemble is additive", max(abs([Fs - (F1 + F2); taus - (tau1 + tau2)]), [], 'all') / max(abs(taus), [], 'all'), 1e-14, "-", "8-term sums in two orders, rounding");
            delta = [0.2; -0.1; 0.05];
            [~, tauD] = NonGravForces.kernel.assemble(Fa, tc.mm, 'CentreOfMass', tc.mm.centreOfMass + delta);
            tc.verifyBoundAndRecord("torque transfer tau' = tau - delta x F", max(abs(tauD - (tau1 - cross(repmat(delta, 1, 4), F1))), [], 'all') / max(abs(tau1), [], 'all'), 1e-14, "-", "rounding");
        end
        function bodyTorquesMetamorphic(tc)
            r = 6.9e6 * [0.3; 0.5; sqrt(1 - 0.34)]; B = [1e-5; -2e-5; 3e-5]; m = [0.3; 0.2; -0.1];
            t1 = NonGravForces.kernel.bodyTorques(tc.mm, r, 'GravityGradient', true); t2 = NonGravForces.kernel.bodyTorques(tc.mm, r, 'GravityGradient', true, 'GM', 2 * referenceSystem.iersConstants().gmEarth);
            tc.verifyEqual(t2, 2 * t1, 'linear in GM, exact for a power of two');
            t3 = NonGravForces.kernel.bodyTorques(tc.mm, 2 * r, 'GravityGradient', true);
            tc.verifyBoundAndRecord("gravity gradient scales as r^-3", norm(t3 - t1 / 8) / norm(t1), 1e-14, "-", "rounding");
            ta = NonGravForces.kernel.bodyTorques(tc.mm, r, 'Magnetic', true, 'B_body', B, 'Dipole_Am2', m); tb = NonGravForces.kernel.bodyTorques(tc.mm, r, 'Magnetic', true, 'B_body', m, 'Dipole_Am2', B);
            tc.verifyEqual(ta, -tb, 'm x B = -(B x m)');
            tc.verifyEqual(dot(ta, m), 0, 'AbsTol', 1e-20); tc.verifyEqual(dot(ta, B), 0, 'AbsTol', 1e-20);
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function partialsSpeed(tc)
            % One derivative costs about one force evaluation: the coefficients are computed once more.
            [rho, sp, T] = tc.gas(200);
            vb = 7600 * [ones(1, 200); 0.02 * cos(1:200); 0.01 * sin(1:200)];
            vb = vb ./ vecnorm(vb) * 7600;
            args = {'Rho', rho, 'Species', sp, 'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300};
            base = @() NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            fcn = @() NonGravForces.kernel.partials("drag", "alpha", tc.mm, args{:});
            tc.verifyBoundAndRecord("d(drag)/d(alpha) vs one force evaluation", tc.speedRatio(fcn, base), 9, "-", "budget 3x the expected 3: the derivative evaluates the force and the coefficients again; 5 repeats, median");
        end
        function gsiCoefficientsSpeed(tc)
            % 1e5 elements at once vs 100 scalar calls: the array form must be at least 20x cheaper per element (measured ~200x)
            th = rand(1, 1e5) * pi; th1 = th(1:100);
            fcn = @() NonGravForces.kernel.gsiCoefficients("sentman", 8, th, 0.85, 300, 1000);
            base = @() arrayfun(@(t) NonGravForces.kernel.gsiCoefficients("sentman", 8, t, 0.85, 300, 1000), th1);
            q = tc.speedRatio(fcn, base);                 % time(1e5 array) / time(100 scalars)
            tc.verifyBoundAndRecord("gsiCoefficients: 1e5-element array vs 100 scalar calls, time ratio", q, 50, "-", "budget 3x the expected 1000/60; N = 1e5, 5 repeats, median");
        end
        function dragPanelForceSpeed(tc)
            [rho, sp, T] = tc.gas(1000); v = 7600 * [ones(1, 1000); 0.01 * randn(1, 1000); 0.01 * randn(1, 1000)];
            fcn = @() NonGravForces.kernel.dragPanelForce(tc.mm, rho, sp, T, v, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            base = @() NonGravForces.kernel.dragPanelForce(tc.mm, rho(1), sp(:, 1), T(1), v(:, 1), 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
            q = tc.speedRatio(fcn, base);
            tc.verifyBoundAndRecord("dragPanelForce: K = 1000 batched vs one epoch, time ratio", q, 100, "-", "budget 3x the expected ~30; 8 species x 8 plates; 5 repeats, median");
        end
        function radiationPanelForceSpeed(tc)
            K = 1000; M = 51; d = randn(3, M, K); d = d ./ vecnorm(d, 2, 1);
            src = struct('dir', d, 'P', 1e-7 * rand(M, K), 'band', ["vis" repmat("ir", 1, M - 1)]);
            s1 = struct('dir', d(:, :, 1), 'P', src.P(:, 1), 'band', src.band);
            q = tc.speedRatio(@() NonGravForces.kernel.radiationPanelForce(tc.mm, src), @() NonGravForces.kernel.radiationPanelForce(tc.mm, s1));
            tc.verifyBoundAndRecord("radiationPanelForce: K = 1000 x 51 sources vs one epoch, time ratio", q, 300, "-", "budget 3x the expected ~100 (the source loop is not batched over K); 5 repeats, median");
        end
        function thermalRecoilSpeed(tc)
            T = 280 + 20 * rand(tc.mm.n, 10000);
            q = tc.speedRatio(@() NonGravForces.kernel.thermalRecoil(tc.mm, T), @() NonGravForces.kernel.thermalRecoil(tc.mm, T(:, 1)));
            tc.verifyBoundAndRecord("thermalRecoil: K = 1e4 vs one epoch, time ratio", q, 60, "-", "budget 3x the expected ~20; 5 repeats, median");
        end
        function assembleSpeed(tc)
            Fk = randn(3, tc.mm.n, 10000);
            q = tc.speedRatio(@() NonGravForces.kernel.assemble(Fk, tc.mm), @() NonGravForces.kernel.assemble(Fk(:, :, 1), tc.mm));
            tc.verifyBoundAndRecord("assemble: K = 1e4 vs one epoch, time ratio", q, 60, "-", "budget 3x the expected ~20; 5 repeats, median");
        end
        function bodyTorquesSpeed(tc)
            r = 6.9e6 * (randn(3, 10000) ./ vecnorm(randn(3, 10000))); B = 3e-5 * randn(3, 10000);
            q = tc.speedRatio(@() NonGravForces.kernel.bodyTorques(tc.mm, r, 'GravityGradient', true, 'Magnetic', true, 'B_body', B, 'Dipole_Am2', [0.1; 0.2; 0.3]), ...
                @() NonGravForces.kernel.bodyTorques(tc.mm, r(:, 1), 'GravityGradient', true, 'Magnetic', true, 'B_body', B(:, 1), 'Dipole_Am2', [0.1; 0.2; 0.3]));
            tc.verifyBoundAndRecord("bodyTorques: K = 1e4 vs one epoch, time ratio", q, 60, "-", "budget 3x the expected ~20; 5 repeats, median");
        end
    end
end
