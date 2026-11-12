%% Tutorial 5 - The spacecraft: panels, shadows and wall temperatures
% NonGravForces sees a satellite as a set of flat panels ("macro model"): each with an area, an
% outward normal, optical properties in the visible and the infrared, a heat capacity and a
% polygon. Every force is a sum over these panels. This tutorial builds the GRACE-FO and GRACE
% models from their handbooks and a generic box-wing, checks whether the panels close, asks which
% panel shadows which - for a beam of sunlight and for a stream of gas molecules - and computes
% the wall temperatures that drive the thermal recoil, in steady state and with thermal inertia,
% across a real shadow entry.
%
% *What you will learn*
%
% * what a |geom.MacroModel| holds, and how to draw one;
% * what the silhouette defect sum(A n) says about a panel table;
% * how |geom.graceFO|, |geom.grace| and |geom.boxWing| differ;
% * how |geom.selfShadow| computes occlusion by ray casting or from a cube-map table, and what
%   the table's sample count is;
% * why the gas needs its own shadow (|geom.flowShadow|), and the exact and tabulated methods;
% * how |geom.thermalState| turns the radiation sources into wall temperatures, steady state
%   against the transient walls the standard preset uses.
%
% *Requirements:* none - every file this tutorial reads is committed in |data/|. Run time: about
% a minute (the occlusion tables are built once).
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - The panel table of GRACE-C
% |geom.graceFO("C")| transcribes the GRACE-FO handbook (JPL D-56935, Tables 4 and 5): eight
% plates in the Science Reference Frame (+X towards the K-band horn, i.e. towards the other
% satellite, +Z to nadir). The side panels slope at 50 deg and carry the solar arrays on their
% outer faces; their inner faces overhang the nadir plate. The S-band boom is left out: a flat
% plate cannot represent a cylinder (0.27 % of the area). Visible absorptivity absVis and infrared
% emissivity emisIR are per panel (each band's absorbed, specular and diffuse parts sum to 1), as
% is the heat capacity per area of the wall behind it (Woeske 2019, Table 6.1). The MASS is the
% handbook spacecraft mass, the inertia tensor is Table 4's.
macro = NonGravForces.geom.graceFO("C");
P = table([macro.panels.name]', [macro.panels.group]', macro.A', round(macro.N', 4), macro.optVis(1, :)', ...
    macro.optIR(1, :)', macro.heatCapacity', 'VariableNames', {'panel', 'group', 'area_m2', 'normal', 'absVis', 'emisIR', 'C_J_m2K'});
disp(P)
fprintf('%s: %d panels, %.3f m^2 in all, mass %.3f kg, Ixx Iyy Izz = %s kg m^2\n', macro.name, macro.n, sum(macro.A), ...
    macro.mass_kg, mat2str(diag(macro.inertia)', 6));

%% Step 2 - Drawing the panels
% The polygons |macro.V{k}| (3 x V corners in the body frame, wound counter-clockwise about the
% outward normal) exist for occlusion only; every force uses the printed area. Below, each plate
% and its outward normal from its centroid |macro.C(:, k)|. Note the view: +Z points DOWN to the
% Earth.
figure, hold on
col = lines(macro.n);
for k = 1:macro.n
    vk = macro.V{k};
    patch('XData', vk(1, :), 'YData', vk(2, :), 'ZData', vk(3, :), 'FaceColor', col(k, :), 'FaceAlpha', 0.6);
end
quiver3(macro.C(1, :), macro.C(2, :), macro.C(3, :), 0.5 * macro.N(1, :), 0.5 * macro.N(2, :), 0.5 * macro.N(3, :), 0, 'k');
hold off, axis('equal'), grid on, view(35, -25)
set(gca, 'ZDir', 'reverse')
xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m, nadir)')
legend([macro.panels.name], 'Location', 'eastoutside'), title('GRACE-C macro model')

%% Step 3 - Does the panel set close? The silhouette defect
% For ANY closed body, convex or not, the area-weighted sum of the outward normals is zero (the
% divergence theorem). |macro.silhouetteDefect()| returns sum(A_k n_k), so it is a one-line check
% of a panel table. The handbook's printed areas do not close on each other: the defect is not
% zero along z. The polygons, built from the exterior dimensions so that they do close, differ
% from the printed areas by up to 2 % - computed below with the cross-product formula for a
% planar polygon.
fprintf('silhouette defect sum(A n) = %s m^2\n', mat2str(macro.silhouetteDefect()', 4));
fprintf('%-16s %10s %10s %8s\n', 'panel', 'printed', 'polygon', 'diff');
for k = 1:macro.n
    vk = macro.V{k}; s = [0; 0; 0];
    for j = 1:size(vk, 2)
        s = s + cross(vk(:, j), vk(:, mod(j, size(vk, 2)) + 1));
    end
    fprintf('%-16s %10.4f %10.4f %7.2f%%\n', macro.panels(k).name, macro.A(k), norm(s) / 2, 100 * (norm(s) / 2 / macro.A(k) - 1));
end

%% Step 4 - Classic GRACE and the arrays switched off
% |geom.grace("A")| builds GRACE-A/B from the GRACE Product Specification (Bettadpur 2012); its
% surface table is the GRACE-FO handbook's, digit for digit, but the exterior dimensions differ by
% millimetres, the mass is the LAUNCH mass and no inertia is printed (NaN). Both builders take
% |'SolarArrayOperating', false|: a non-operating array absorbs more of the visible light (the
% handbook's note 2).
ga = NonGravForces.geom.grace("A");
fprintf('%-16s %10s   %-16s %10s\n', 'GRACE-C panel', 'area', 'GRACE-A panel', 'area');
for k = 1:macro.n
    fprintf('%-16s %10.4f   %-16s %10.4f\n', macro.panels(k).name, macro.A(k), ga.panels(k).name, ga.A(k));
end
fprintf('GRACE-A: mass %.1f kg (launch), inertia known: %d; silhouette defect %s m^2\n', ga.mass_kg, all(isfinite(ga.inertia(:))), ...
    mat2str(ga.silhouetteDefect()', 4));
off = NonGravForces.geom.graceFO("C", 'SolarArrayOperating', false);
isArr = [macro.panels.group] == "array";
fprintf('array absVis operating %s, not operating %s\n', mat2str(unique(macro.optVis(1, isArr))), mat2str(unique(off.optVis(1, isArr))));

%% Step 5 - A generic box-wing
% |geom.boxWing([Lx Ly Lz], optics)| builds a box of six plates centred at the origin, with the
% caller's optics (|.bus.vis| = [abs spec diff], |.bus.ir| = [emis spec diff]) and, optionally,
% solar-array plates given by area, normal and centroid. Nothing in it is a spacecraft: every
% number is yours. A closed box has zero defect; a ONE-sided array plate opens it by its own area,
% and a two-sided wing (two entries with opposite normals) closes it again. Array plates carry no
% polygon, so a box-wing with arrays cannot be ray cast: |selfShadow| refuses it.
o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]), 'array', struct('vis', [0.65 0.05 0.30], 'ir', [0.81 0.03 0.16]));
bw = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Name', "box", 'Mass_kg', 500);
wing = struct('area_m2', {2, 2}, 'normal', {[0; 0; -1], [0; 0; 1]}, 'centroid_m', {[0; 2; 0], [0; 2; 0]}, 'name', {"wingTop", "wingBottom"});
bw1 = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Array', wing(1));
bw2 = NonGravForces.geom.boxWing([3 1.9 0.7], o, 'Array', wing);
fprintf('box: %d plates, defect %s; one-sided wing: %d plates, defect %s; two-sided: %d plates, defect %s\n', bw.n, ...
    mat2str(bw.silhouetteDefect()'), bw1.n, mat2str(bw1.silhouetteDefect()'), bw2.n, mat2str(bw2.silhouetteDefect()'));
try
    NonGravForces.geom.selfShadow(bw2, [0; 0; 1]);
catch err
    fprintf('selfShadow on the box-wing with arrays: %s (polygons on every plate: %d)\n', err.identifier, bw2.hasVertices);
end

%% Step 6 - Who shadows whom: geom.selfShadow
% For each plate and each direction d (a unit vector FROM the body TOWARD the source), f is the
% fraction of the plate's area from which a ray along d leaves the body without hitting another
% plate. It is occlusion and nothing else: whether the plate FACES d is the consumer's factor
% (n . d)+. On GRACE-FO only the nadir plate is ever shadowed - by the overhanging inner faces of
% the side panels, for directions low over its horizon. Two methods:
%
% * "raycast": an m x m grid of points per plate (|'Samples'|, default 40), exact to about 2/m at
%   a shadow edge;
% * "lut": the ray cast evaluated once at the nodes of a cube map on the sphere (|'LutResolution'|
%   n cells per face edge, the standard preset uses 16) and interpolated - O(1) per direction.
%
% A table is built once and passed back as |'Lut'|; a reused table reports ITS OWN sample count in
% |info.samples|, whatever |'Samples'| says. Below, the nadir plate for directions tilted from
% straight down (+Z) towards +Y.
kN = find([macro.panels.name] == "Nadir");
lut = NonGravForces.geom.selfShadow(macro, [], 'Method', "lut");
el = 0:5:85;
dirs = [zeros(1, numel(el)); sind(el); cosd(el)];
fR = NonGravForces.geom.selfShadow(macro, dirs, 'Plate', kN);
[fL, infoL] = NonGravForces.geom.selfShadow(macro, dirs, 'Method', "lut", 'Lut', lut, 'Samples', 7);
fprintf('table: resolution %d, %d directions, built with %d samples (info.samples %d although Samples 7 was passed)\n', ...
    lut.resolution, size(lut.nodes, 2), lut.samples, infoL.samples);
fprintf('plates the table ever shades: %s\n', strjoin([macro.panels(~all(lut.f == 1, 2)).name], ", "));
fprintf('nadir plate, 85 deg from nadir: ray cast %.3f, table %.3f; largest difference %.3f\n', fR(end), fL(kN, end), max(abs(fR - fL(kN, :))));
figure
plot(el, fR, 'o-', el, fL(kN, :), 's-'), grid on
xlabel('angle of the direction from +Z towards +Y (deg)'), ylabel('unobstructed fraction')
legend("ray cast, 40 samples", "cube-map table", 'Location', 'southwest'), title('Self-shadowing of the GRACE-C nadir plate')

%% Step 7 - The gas has its own shadow: geom.flowShadow
% Molecules do not arrive as a parallel beam. Relative to the satellite the gas is a Maxwellian
% drifting at the relative velocity, so molecules come from a spread of directions of angular
% width about 1/s, s the speed ratio (about 8 for atomic oxygen, 3.5 for helium at 500 km).
% Sentman's coefficients are integrals over that spread, which is why an edge-on or even leeward
% plate feels a force; an obstruction must act INSIDE those integrals, weighted by the flux each
% direction delivers. |geom.flowShadow(macro, vrel_body, T_K, ...)| returns per plate, epoch and
% species the ratios (with obstruction) / (without) of the number, normal-momentum and
% tangential-momentum fluxes, and |S.beam|, the parallel-beam occlusion along v_rel. Method "lut"
% integrates them directly on the occlusion table; Method "table" (the standard preset) reads them
% from a grid over the stream direction and s_O that is built once per session (a few seconds)
% and falls back to "lut" outside it. The inputs here are real: the relative velocity of Tutorial
% 3 (co-rotation removed) rotated into the body frame, and the NRLMSISE-00 gas temperature.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
gnv = f('samples', 'GNV1B_2020-01-01_C_04.txt');
[ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
[~, rI] = NonGravForces.io.readGNV1B(gnv, 'Frame', "ITRS");
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
idx = 1:10:numel(ep.day);
e = referenceSystem.daySecondToEpoch([double(ep.day(idx)) double(ep.sec(idx))], string(ep.scale));
lla = referenceSystem.cartesianToGeodetic(rI(:, idx)', 'Ellipsoid', "WGS84")';
[~, ~, Tgas] = NonGravForces.env.density("msis00", e, lla, sw);
vrel = NonGravForces.env.wind("corotation", e, r(:, idx), v(:, idx), 'Eop', eop);
Rb = referenceSystem.rotationFromQuaternion(att.at(e));                       % BODY -> GCRS
vb = reshape(pagemtimes(Rb, 'transpose', reshape(vrel, 3, 1, []), 'none'), 3, []);
fprintf('v_rel in the body frame: mean [%.1f %.1f %.1f] m/s; gas temperature %.0f K\n', mean(vb, 2), mean(Tgas(2, :)));
tic; [SL, iL] = NonGravForces.geom.flowShadow(macro, vb, Tgas(2, :), 'Method', "lut", 'Lut', lut); tL = toc;
tic; [ST, iT] = NonGravForces.geom.flowShadow(macro, vb, Tgas(2, :), 'Method', "table", 'Lut', lut); tT = toc;
tic; NonGravForces.geom.flowShadow(macro, vb, Tgas(2, :), 'Method', "table", 'Lut', lut); tT2 = toc;
iO = find(NonGravForces.const().speciesName == "O");
fprintf('%d epochs: "lut" %.2f s (%s); "table" %.2f s the first time (builds the grid), %.3f s again; %d fallback epochs\n', ...
    numel(idx), tL, iL.method, tT, tT2, iT.fallbackEpochs);
fprintf('nadir plate, atomic oxygen (speed ratio %.1f): beam %.4f, number %.4f, normal %.4f, tangential %.4f\n', ...
    mean(SL.speedRatio(iO, :)), mean(SL.beam(kN, :)), mean(SL.number(kN, :, iO)), mean(SL.normal(kN, :, iO)), mean(SL.tangential(kN, :, iO)));
fprintf('table vs lut, largest difference over every plate and species: number %.1e, normal %.1e, tangential %.1e\n', ...
    max(abs(ST.number - SL.number), [], 'all'), max(abs(ST.normal - SL.normal), [], 'all'), max(abs(ST.tangential - SL.tangential), [], 'all'));
figure
names = cellstr([macro.panels.name]);
bar(categorical(names, names), squeeze(mean(SL.number(:, :, [1 2 3]), 2))), grid on
ylabel('number-flux ratio (obstructed / free)'), legend(NonGravForces.const().speciesName([1 2 3]), 'Location', 'southwest')
title('Flow shadow per plate along the arc')

%% Step 8 - Wall temperatures: steady state and transient walls
% The panels absorb sunlight, albedo and Earth infrared and re-emit it; the asymmetry of that
% emission is the thermal recoil force. |geom.thermalState(macro, src)| forms the absorbed flux
% per panel from the radiation source list of Tutorial 4 (visible absorptivity for the Sun and the
% albedo, infrared emissivity for the Earth's infrared) and solves for the temperature:
%
% * "steadyState": emitted = absorbed at every epoch, no thermal inertia - a panel that absorbs
%   nothing is at 0 K the moment the light goes;
% * "transient" (the standard and precise presets since 0.18.0): C dT/dt = absorbed - emitted
%   with the heat capacity C of each wall, integrated from |'T0_K'| across the epochs.
%
% Across the shadow entry of the 30-minute arc the difference is the point: the steady walls
% follow the light instantly, the transient ones cool over minutes - which is what ACT1B sees
% (moving the standard preset to transient walls took the radial residual from 3.25 to 2.28
% nm/s^2, CLAUDE.md, PR #126). Here the first epoch's steady state is the initial temperature;
% |Model.accel| takes |'T0_K'| the same way.
md = f('samples', 'multiday');
[ep30, r30] = NonGravForces.io.readGNV1B(fullfile(md, 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq30, q30] = NonGravForces.io.readSCA1B(fullfile(md, 'SCA1B_2020-01-01_C_04.txt'));
att30 = NonGravForces.geom.AttitudeProvider(eq30, q30);
i30 = 1:10:numel(ep30.day);
e30 = referenceSystem.daySecondToEpoch([double(ep30.day(i30)) double(ep30.sec(i30))], string(ep30.scale));
src = NonGravForces.env.radiationSources(e30, r30(:, i30), att30.at(e30), NonGravForces.Options("standard"));
Ts = NonGravForces.geom.thermalState(macro, src);
[Tt, it] = NonGravForces.geom.thermalState(macro, src, 'Model', "transient", 'Epoch', e30, 'T0_K', Ts(:, 1)');
t30 = double(ep30.sec(i30)' - ep30.sec(1));
fprintf('%-16s %12s %12s %16s\n', 'panel', 'steady (K)', 'transient', 'time const. (min)');
for k = 1:macro.n
    fprintf('%-16s %5.0f..%-5.0f %5.0f..%-5.0f %16.0f\n', macro.panels(k).name, min(Ts(k, :)), max(Ts(k, :)), min(Tt(k, :)), max(Tt(k, :)), ...
        median(it.timeConstant_s(k, :)) / 60);
end
show = find(ismember([macro.panels.name], ["Zenith" "Nadir" "PortOuter" "Front"]));
figure
plot(t30, Ts(show, :), '--', t30, Tt(show, :), '-'), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('K')
legend([[macro.panels(show).name] + " steady", [macro.panels(show).name] + " transient"], 'Location', 'eastoutside')
title('Wall temperatures across the shadow entry')
