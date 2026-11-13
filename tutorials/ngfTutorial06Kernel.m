%% Tutorial 6 - Inside the kernel: one panel, one gas, one beam
% Tutorial 1 called |Model.accel| and got an acceleration. This tutorial opens the box
% underneath it: the KERNEL layer, where every force is computed panel by panel from three
% ingredients - a gas stream (drag), a list of radiation beams (Sun, Earth albedo, Earth
% infrared) and the panel temperatures (thermal recoil). The kernel knows no orbit, no time and
% no atmosphere model: everything it receives is already in the body frame, in SI units, at
% one or many epochs. That is what makes it testable on its own, and what lets the toolbox keep
% a compiled copy of it (|kernel.singleEpoch|) that must agree with the MATLAB units to
% rounding.
%
% *What you will learn*
%
% * how the normal and tangential plate coefficients of |kernel.gsiCoefficients| are defined,
%   and what they do edge-on and leeward (sentman, maxwell, constantCd);
% * why the speed ratio is computed per SPECIES and never from a mean molar mass;
% * what |kernel.accommodation| returns (fixed, pilinski2010, walker2014, sesam) and when it
%   is extrapolating;
% * how |dragPanelForce|, |radiationPanelForce| and |thermalRecoil| produce per-panel forces,
%   and how |assemble| reduces ONE per-panel array to a force and a torque;
% * the body torques of |kernel.bodyTorques| (gravity gradient);
% * that the compiled |kernel.singleEpoch| reproduces |Model.accel| at one epoch;
% * how |kernel.partials| compares with finite differences of the force.
%
% *Requirements:* none - one epoch of the committed 2020-01-01 GRACE-C day. Step 9 uses the
% compiled kernel (|NonGravForces.buildMex('Target', "kernel")|; |ciGate| builds it) and says
% so if it is missing. Run time: about ten seconds.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - One epoch of the real orbit: the gas the satellite meets
% The kernel needs the gas in the body frame. We build it exactly as |Model| does (its help
% lists the order): the satellite velocity relative to the co-rotating atmosphere from
% |env.wind| (GCRS, m/s), rotated into the body frame with the SCA1B attitude; the geodetic
% position on WGS-84 from the ITRS position; and NRLMSISE-00 (|env.density("msis00", ...)|)
% for the total mass density |rho| (kg/m^3), the eight number densities |sp| (1/m^3, in the
% order of |const.speciesName|) and the temperatures |T| (K; |T(2)| is the gas temperature at
% the satellite). No observed density scale is applied here - it multiplies every species alike
% and changes nothing in what this tutorial shows.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
[ep, r, vel] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
macro = NonGravForces.geom.graceFO("C");
c = NonGravForces.const();

k = 7;                                                 % an arbitrary epoch of the arc
e1 = referenceSystem.daySecondToEpoch([ep.day(k) ep.sec(k)], ep.scale);
r1 = r(:, k); v1 = vel(:, k); q1 = att.at(e1);
R = referenceSystem.rotationFromQuaternion(q1);         % body -> GCRS
vrel = NonGravForces.env.wind("corotation", e1, r1, v1, 'Eop', eop);
vb = R' * vrel;                                        % relative velocity, body frame
Ritrs = referenceSystem.rotationItrsToGcrs(e1, 'Eop', eop);
lla = referenceSystem.cartesianToGeodetic((Ritrs' * r1)', 'Ellipsoid', "WGS84")';
[rho, sp, T] = NonGravForces.env.density("msis00", e1, lla, sw);
Tg = T(2);
U = norm(vb);
fprintf('epoch %d: height %.1f km, |v_rel| %.1f m/s, v_rel/|v_rel| (body) = [%.4f %.4f %.4f]\n', ...
    k, lla(3) / 1e3, U, vb / U);
fprintf('rho %.3e kg/m^3, gas temperature %.1f K\n', rho, Tg);
mj = c.speciesMass(:) / c.NA;                          % kg per particle
fprintf('%-12s %12s %10s\n', 'species', 'n (1/m^3)', 'mass frac');
for j = 1:8
    fprintf('%-12s %12.3e %10.4f\n', c.speciesName(j), sp(j), sp(j) * mj(j) / sum(sp .* mj));
end

%% Step 2 - The speed ratio, species by species
% Free-molecular flow is governed by the speed ratio s = U / sqrt(2 k T / m): the satellite
% speed over the most probable thermal speed of the molecules. Each species has its OWN m and
% therefore its own s. A light gas (He, H) is slow compared with its own thermal motion and hits
% even a leeward plate; atomic oxygen at s ~ 7-8 behaves almost like a beam. Treating the gas
% as "air" with one mean molar mass gets s wrong for every species at once (the predecessor's
% defect D1b; CLAUDE.md section 7 quotes -24 % on a plate at 100 deg incidence).
s = U ./ sqrt(2 * c.kB * Tg ./ mj);                    % 8x1
mMean = sum(sp .* mj) / sum(sp);                       % number-weighted mean mass
sMean = U / sqrt(2 * c.kB * Tg / mMean);
fprintf('%-12s %8s\n', 'species', 's');
for j = 1:8, fprintf('%-12s %8.3f\n', c.speciesName(j), s(j)); end
fprintf('mean molar mass %.2f g/mol -> one "mean" speed ratio %.3f\n', mMean * c.NA * 1e3, sMean);

%% Step 3 - The plate coefficients: normal and tangential, three models
% |kernel.gsiCoefficients(model, s, theta, alpha, Tw_K, Tinf_K)| returns, per unit plate area and
% in the (1/2) rho U^2 convention, the NORMAL coefficient |Cp| (pushes the plate inward, along
% -n) and the TANGENTIAL coefficient |Ctau| (drags it along the tangential part of the stream).
% theta is the angle between the outward normal and the relative velocity: 0 is head-on,
% pi/2 edge-on, beyond pi/2 leeward. The drag and lift of the plate are projections of that
% pair, |info.CD = Cp cos(theta) + Ctau sin(theta)| and |info.CL|. Sentman's published C_D and
% C_L are drag and lift, NOT normal and tangential - CLAUDE.md section 7 records how long that
% took to get right.
%
% * |"sentman"|: diffuse re-emission with energy accommodation |alpha| at wall temperature
%   |Tw_K|; the standard preset, alpha 0.85 (March et al. 2021), Tw 300 K;
% * |"maxwell"|: a |'SpecularFraction'| reflects like a mirror, the rest as sentman;
% * |"constantCd"|: the textbook Cp = Cd max(gamma, 0) gamma, Ctau = Cd max(gamma, 0) sin(theta) -
%   no lift, nothing on a leeward plate; the fast preset (Cd 2.2).
%
% The EDGE-ON limit is the check that the convention is the right one: at gamma = 0 Sentman's
% plate still feels the thermal pressure 1/(2 s^2) plus the re-emission vr/(2 s), and a shear
% 1/(s sqrt(pi)) along the stream; constantCd gives exactly 0.
sO = s(2);                                             % atomic oxygen carries most of the drag
alpha = 0.85; Tw = 300;
th = linspace(0, pi, 181)';
[CpS, CtS, iS] = NonGravForces.kernel.gsiCoefficients("sentman", sO, th, alpha, Tw, Tg);
[CpM, CtM] = NonGravForces.kernel.gsiCoefficients("maxwell", sO, th, alpha, Tw, Tg, 'SpecularFraction', 0.1);
[CpC, CtC] = NonGravForces.kernel.gsiCoefficients("constantCd", [], th, [], [], [], 'Cd', 2.2);
vr = sqrt(0.5 * ((1 - alpha) + alpha * 2 * Tw / (Tg * sO^2)));
[Cp0, Ct0] = NonGravForces.kernel.gsiCoefficients("sentman", sO, pi / 2, alpha, Tw, Tg);
fprintf('edge-on sentman: Cp %.6e (closed form 1/(2s^2) + vr/(2s) = %.6e)\n', Cp0, 1 / (2 * sO^2) + vr / (2 * sO));
fprintf('                 Ctau %.6e (closed form 1/(s sqrt(pi)) = %.6e)\n', Ct0, 1 / (sO * sqrt(pi)));
fprintf('head-on: Cp sentman %.4f, maxwell %.4f, constantCd %.4f\n', CpS(1), CpM(1), CpC(1));
fprintf('leeward (theta = 120 deg): Cp sentman %.3e, constantCd %.3e\n', CpS(121), CpC(121));
figure
subplot(2, 1, 1)
plot(rad2deg(th), [CpS CpM CpC]), grid on, ylabel('C_p (normal)')
legend("sentman", "maxwell 10 % specular", "constantCd 2.2"), title("Plate coefficients for atomic oxygen, s = " + sprintf('%.2f', sO))
subplot(2, 1, 2)
plot(rad2deg(th), [CtS CtM CtC], rad2deg(th), iS.CD, 'k--'), grid on
ylabel('C_\tau (tangential)'), xlabel('angle between normal and relative velocity (deg)')
legend("sentman", "maxwell", "constantCd", "sentman drag projection C_D")

%% Step 4 - Why the mixture is summed species by species
% |kernel.dragPanelForce| calls |gsiCoefficients| once per species with that species' own s and
% weights the coefficients by the species' mass density rho_j = n_j m_j (Doornbos 2011
% eq. 3.57). Here that mixture is formed by hand for the same angles and compared with the
% coefficient of ONE gas at the mean molar mass. On the ram face they agree closely; at grazing
% and leeward angles, where only the slow light species still reach the plate, they do not.
rhoJ = (sp .* mj)';                                    % 1x8 mass densities
CpJ = NonGravForces.kernel.gsiCoefficients("sentman", s', th, alpha, Tw, Tg);   % 181x8 by expansion
CpMix = sum(CpJ .* rhoJ, 2) / sum(rhoJ);
CpOne = NonGravForces.kernel.gsiCoefficients("sentman", sMean, th, alpha, Tw, Tg);
fprintf('%10s %12s %12s %8s\n', 'theta', 'per species', 'mean mass', 'ratio');
for a = [0 60 90 100 120 135]
    i = a + 1;
    fprintf('%7d deg %12.4e %12.4e %8.3f\n', a, CpMix(i), CpOne(i), CpOne(i) / CpMix(i));
end

%% Step 5 - Accommodation: one number, or a surface of two states
% The accommodation coefficient alpha says how far the re-emitted molecules have cooled towards
% the wall. |kernel.accommodation(model, species, T_K, ...)| returns either one alpha per epoch
% (|"fixed"|, the |"pilinski"| fits) or a two-state surface: a fraction |coverage| covered by
% adsorbed oxygen that accommodates fully, and clean substrate with |alphaSubstrate|
% (|"walker2014"|, |"sesam"|). |dragPanelForce| reads either form. Every model reports
% |info.insideFit|: false means the formula is evaluated outside the conditions it was fitted
% to. At GRACE-FO's altitude in solar minimum that is the normal case for the computed models,
% and the toolbox WARNS (|NonGravForces:accommodation:outsideFit|, HANDOVER O-42); the warning
% is silenced here only so that the printout stays readable - the flag says the same thing.
% sesam needs the relative SPEED because its oxygen pressure is a ram pressure.
ws = warning('off', 'NonGravForces:accommodation:outsideFit');
[aF, iF] = NonGravForces.kernel.accommodation("fixed", sp, Tg, 'Alpha', 0.85);
[aP, iP] = NonGravForces.kernel.accommodation("pilinski2010", sp, Tg);
[aW, iW] = NonGravForces.kernel.accommodation("walker2014", sp, Tg, 'SurfaceMass_amu', 263.3223);
[aS, iSe] = NonGravForces.kernel.accommodation("sesam", sp, Tg, 'RelativeSpeed_ms', U);
warning(ws);
fprintf('%-13s %8s %9s %10s %10s\n', 'model', 'alpha', 'coverage', 'alpha_sub', 'insideFit');
fprintf('%-13s %8.4f %9s %10s %10d\n', 'fixed', aF.alpha, '-', '-', iF.insideFit);
fprintf('%-13s %8.4f %9s %10s %10d\n', 'pilinski2010', aP.alpha, '-', '-', iP.insideFit);
fprintf('%-13s %8.4f %9.4f %10.4f %10d\n', 'walker2014', aW.alpha, aW.coverage, aW.alphaSubstrate, iW.insideFit);
fprintf('%-13s %8.4f %9.4f %10.4f %10d\n', 'sesam', aS.alpha, aS.coverage, aS.alphaSubstrate, iSe.insideFit);
fprintf('sesam: sticking coefficient %.3e, oxygen ram pressure %.3e Pa\n', iSe.stickingCoefficient, iSe.pressureO_Pa);

%% Step 6 - The gas force on every panel
% |kernel.dragPanelForce(macro, rho, species, T_K, vrel_body, 'Model', ..., ...)| returns the
% force on EACH panel, 3 x N x K in newtons, body frame. |info.gamma| is cos(theta) per panel,
% |info.speedRatio| the eight speed ratios (the same numbers as Step 2). |kernel.assemble|
% reduces the per-panel array to the total force; dividing by the mass gives the acceleration.
% The effective drag coefficient, referred to the projected area, shows what Sentman does on a
% real shape: the Doornbos (2011) band for GRACE is 3.2-3.8 (tests/tNgfKernel.m), against the
% 2.2 of the fast preset. (|Model| also passes a self-shadowing struct from |geom.flowShadow|
% as |'Shadow'|; it is left out here.)
[Fs, di] = NonGravForces.kernel.dragPanelForce(macro, rho, sp, Tg, vb, 'Model', "sentman", 'Alpha', alpha, 'Tw_K', Tw);
Fm = NonGravForces.kernel.dragPanelForce(macro, rho, sp, Tg, vb, 'Model', "maxwell", 'Alpha', alpha, 'Tw_K', Tw, 'SpecularFraction', 0.1);
Fc = NonGravForces.kernel.dragPanelForce(macro, rho, [], [], vb, 'Model', "constantCd", 'Cd', 2.2);
Fw = NonGravForces.kernel.dragPanelForce(macro, rho, sp, Tg, vb, 'Model', "sentman", 'Alpha', aW, 'Tw_K', Tw);
fprintf('max |speed ratio from dragPanelForce - Step 2| / s = %.1e\n', max(abs(di.speedRatio - s) ./ s));
fprintf('%-16s %8s %12s\n', 'panel', 'gamma', '|F| (N)');
for p = 1:macro.n
    fprintf('%-16s %8.4f %12.4e\n', macro.panels(p).name, di.gamma(p), norm(Fs(:, p)));
end
Aproj = sum(macro.A(:) .* max(macro.N' * (vb / U), 0));
names = ["sentman" "maxwell 10 %" "constantCd 2.2" "sentman, walker2014"];
pages = {Fs, Fm, Fc, Fw};
fprintf('projected area %.3f m^2\n', Aproj);
for i = 1:4
    Ft = NonGravForces.kernel.assemble(pages{i}, macro);
    fprintf('%-20s a = [%9.3e %9.3e %9.3e] m/s^2, Cd on projected area %.3f\n', names(i), ...
        Ft / macro.mass_kg, norm(Ft) / (0.5 * rho * U^2 * Aproj));
end

%% Step 7 - Radiation beams and thermal recoil
% |env.radiationSources(epoch, r, q, options)| turns the Sun (with the eclipse factor inside
% its pressure), the Earth's albedo tiles (band "vis") and its longwave tiles (band "ir") into
% ONE list of beams: |src.dir| 3 x M unit vectors from the satellite to each source (body
% frame), |src.P| the pressure in Pa. The kernel does not know which beam is which effect
% (principle P3). |kernel.radiationPanelForce(macro, src)| applies the flat-plate optical model
% (absorbed, specular, diffuse) to every panel and beam; its |info.absorbed_W_m2| is the
% absorbed flux per panel, from which |geom.thermalState| forms the panel temperatures (here
% the steady state; the standard preset integrates them in time) and |kernel.thermalRecoil|
% the recoil of their Lambertian emission, F_k = -(2/3) eps sigma T^4 A / c n_k.
opt = NonGravForces.Options("standard");
src = NonGravForces.env.radiationSources(e1, r1, q1, opt);
bands = string(src.band);
fprintf('%d beams: %d "vis", %d "ir"; Sun pressure %.3e Pa\n', numel(bands), nnz(bands == "vis"), nnz(bands == "ir"), src.P(1));
[Fr, ri] = NonGravForces.kernel.radiationPanelForce(macro, src);
Tp = NonGravForces.geom.thermalState(macro, src, 'Model', "steadyState");
[Ft, ti] = NonGravForces.kernel.thermalRecoil(macro, Tp);
fprintf('%-16s %10s %8s\n', 'panel', 'absorbed', 'T (K)');
for p = 1:macro.n
    fprintf('%-16s %8.1f W/m^2 %8.1f\n', macro.panels(p).name, ri.absorbed_W_m2(p), Tp(p));
end
aRad = NonGravForces.kernel.assemble(Fr, macro) / macro.mass_kg;
aThm = NonGravForces.kernel.assemble(Ft, macro) / macro.mass_kg;
fprintf('radiation a = [%9.3e %9.3e %9.3e] m/s^2\n', aRad);
fprintf('thermal   a = [%9.3e %9.3e %9.3e] m/s^2 (emitted power %.1f W)\n', aThm, ti.emittedPower_W);

%% Step 8 - One array, two reductions: force and torque
% |[F, tau] = kernel.assemble(Fk, macro)| forms the total force AND the torque
% sum_k (c_k - r_cm) x F_k about the centre of mass |macro.centreOfMass|, from the same
% per-panel array (principle P2), so the torque can never miss a term the force has. Torques
% that act on the body as a whole and never reach an accelerometer come from
% |kernel.bodyTorques|: the gravity gradient 3 GM/|r|^5 r x (I r) needs the position IN THE BODY
% FRAME (R' r) and the inertia tensor of the macro model (JPL D-56935 Table 4); the magnetic
% torque m x B needs a field and a dipole and is left out here.
[~, tauD] = NonGravForces.kernel.assemble(Fs, macro);
[~, tauR] = NonGravForces.kernel.assemble(Fr, macro);
[~, tauT] = NonGravForces.kernel.assemble(Ft, macro);
[tauG, gi] = NonGravForces.kernel.bodyTorques(macro, R' * r1, 'GravityGradient', true);
fprintf('torque about the CoM (N m):\n');
fprintf('  drag             [%10.3e %10.3e %10.3e]\n', tauD);
fprintf('  radiation        [%10.3e %10.3e %10.3e]\n', tauR);
fprintf('  thermal          [%10.3e %10.3e %10.3e]\n', tauT);
fprintf('  gravity gradient [%10.3e %10.3e %10.3e] (enabled %d)\n', tauG, gi.enabled.gravityGradient);

%% Step 9 - The compiled path: kernel.singleEpoch against Model.accel
% |kernel.singleEpoch(macro, gas, src, Tpanel_K, ...)| is the same panel loop written in C, for
% ONE epoch: |gas| holds |rho|, |species|, |T_K|, |vrel_body| and a per-panel |shadow| ([] for
% none); |src| is the beam list; the columns of |F| (N) are drag, radiation and thermal recoil.
% MATLAB is normative and the C is checked against it. Here first against the four MATLAB
% kernel units of Steps 6-7, then end to end against |Model.accel| in the fast preset, composing
% the environment exactly as |tests/tNgfMex.m| does (the fast preset has no self-shadowing and
% no thermal recoil, so only the composition itself is under test).
rep = NonGravForces.buildMex('Check', true, 'Target', "kernel");
if rep.present && rep.upToDate
    gas = struct('rho', rho, 'species', sp, 'T_K', Tg, 'vrel_body', vb, 'shadow', []);
    Fse = NonGravForces.kernel.singleEpoch(macro, gas, src, Tp, 'Model', "sentman", 'Alpha', alpha, 'Tw_K', Tw);
    Fml = [NonGravForces.kernel.assemble(Fs, macro), NonGravForces.kernel.assemble(Fr, macro), NonGravForces.kernel.assemble(Ft, macro)];
    fprintf('singleEpoch vs MATLAB kernel units, per effect: max relative difference %.1e\n', ...
        max(abs(Fse - Fml) ./ max(abs(Fml), [], 1), [], 'all'));

    expo = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
    Mfast = NonGravForces.Model(macro, NonGravForces.Options("fast"), 'Attitude', att, 'Eop', eop, ...
        'SpaceWeather', sw, 'Exponential', expo);
    o = Mfast.options;
    [aRef, ai] = Mfast.accel(e1, r1, v1, q1);
    vrf = NonGravForces.env.wind(o.drag.wind, e1, r1, v1, 'Eop', eop);
    [rf, spf, Tf] = NonGravForces.env.density("exponential", e1, lla, [], 'Exponential', expo);
    srcF = NonGravForces.env.radiationSources(e1, r1, q1, o);
    gasF = struct('rho', rf, 'species', spf, 'T_K', Tf(2), 'vrel_body', R' * vrf, 'shadow', []);
    Ff = NonGravForces.kernel.singleEpoch(macro, gasF, srcF, [], 'Model', o.drag.gsi, 'Cd', o.drag.cd, ...
        'CdScale', o.drag.cdScale, 'CrScale', o.srp.crScale);
    got = Ff / ai.mass_kg;
    want = [ai.drag, ai.radiation, ai.thermal];
    fprintf('singleEpoch vs Model.accel (fast preset): max difference %.1e of the largest effect\n', ...
        max(abs(got - want), [], 'all') / max(abs(want), [], 'all'));
    fprintf('total: compiled [%9.3e %9.3e %9.3e], Model.accel [%9.3e %9.3e %9.3e] m/s^2\n', sum(got, 2), aRef);
else
    fprintf('the compiled kernel is not built here -> run NonGravForces.buildMex(''Target'', "kernel") and re-run this step\n');
end

%% Step 10 - Analytic partials against finite differences
% |kernel.partials(what, wrt, macro, ...)| returns d(F_k)/d(parameter) for every panel, in the
% layout of the force, so |assemble| turns it into the derivative of the total force. It takes
% the same inputs as the force unit, by name. Two of the derivatives are trivial because the
% force is LINEAR in the parameter - a density scale (that is why |app.densityInversion| is a
% division, Tutorial 11) and a radiation pressure scale - and two are not: the accommodation
% alpha and the speed |v_rel|. Each is compared with a central difference of the force itself.
args = {'Rho', rho, 'Species', sp, 'T_K', Tg, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', alpha, 'Tw_K', Tw};
drag = @(da, fU, fR) NonGravForces.kernel.assemble(NonGravForces.kernel.dragPanelForce(macro, rho * fR, sp * fR, Tg, ...
    vb * fU, 'Model', "sentman", 'Alpha', alpha + da, 'Tw_K', Tw), macro);
rel = @(x, y) norm(x - y) / norm(y);
hA = 1e-4; hU = 1e-2; hR = 1e-5;
dA = NonGravForces.kernel.assemble(NonGravForces.kernel.partials("drag", "alpha", macro, args{:}), macro);
dU = NonGravForces.kernel.assemble(NonGravForces.kernel.partials("drag", "speed", macro, args{:}), macro);
dR = NonGravForces.kernel.assemble(NonGravForces.kernel.partials("drag", "densityScale", macro, args{:}), macro);
dP = NonGravForces.kernel.assemble(NonGravForces.kernel.partials("radiation", "pressureScale", macro, 'Sources', src), macro);
nA = (drag(hA, 1, 1) - drag(-hA, 1, 1)) / (2 * hA);
nU = (drag(0, 1 + hU / U, 1) - drag(0, 1 - hU / U, 1)) / (2 * hU);
nR = (drag(0, 1, 1 + hR) - drag(0, 1, 1 - hR)) / (2 * hR);
fprintf('%-28s %12s %12s\n', 'derivative of the force', '|analytic|', 'rel. diff');
fprintf('%-28s %12.4e %12.1e\n', 'drag / alpha (N)', norm(dA), rel(dA, nA));
fprintf('%-28s %12.4e %12.1e\n', 'drag / |v_rel| (N per m/s)', norm(dU), rel(dU, nU));
fprintf('%-28s %12.4e %12.1e\n', 'drag / density scale (N)', norm(dR), rel(dR, nR));
fprintf('%-28s %12.4e %12.1e\n', 'radiation / pressure scale', norm(dP), rel(dP, NonGravForces.kernel.assemble(Fr, macro)));
