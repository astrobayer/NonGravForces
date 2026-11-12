%% Tutorial 4 - Sunlight, shadow and Earth radiation
% Radiation pressure is the second largest non-gravitational force on GRACE-FO, and the one with
% the sharpest features: it switches off within seconds when the satellite enters the Earth's
% shadow. This tutorial follows the first 30 minutes of GRACE-C's orbit on 2020-01-01, which
% contain one shadow ENTRY, and builds what the radiation force is made of: the positions of the
% Sun and Moon, the shadow function of three Earth-shadow models across the penumbra, the source
% list (the Sun plus hundreds of Earth tiles that reflect sunlight and emit infrared), the Earth's
% albedo from Knocke's model and from the CERES measurements, and - because GRACE-FO measures it
% too - the Earth's magnetic field from IGRF-14.
%
% *What you will learn*
%
% * where the Sun and Moon come from (|referenceSystem.geocentricBodyState|) and what the beta
%   angle of the orbit is;
% * how |env.eclipseFactor| models the shadow (cylindrical, conical, SOLAARS-CF) and the lunar
%   factor, and how wide the penumbra is on a real shadow entry;
% * what |env.radiationSources| returns: one Sun source and an albedo and an infrared source
%   per Earth tile, all in the body frame;
% * how |env.earthRadiance| tiles the visible cap, and how CERES (|io.importCERES|) replaces
%   Knocke's zonal model with measured monthly fluxes;
% * how |env.magneticField| evaluates IGRF-14 (|io.igrfCoefficients|) and compares with MAG1B.
%
% *Requirements:* none - every file this tutorial reads is committed in |data/|, the CERES
% excerpt and the IGRF-14 coefficients included. Run time: well under a minute.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Thirty minutes of orbit with a shadow entry
% The multiday excerpts hold the first 30 minutes (1800 epochs at 1 s) of GNV1B and SCA1B of
% 2020-01-01; the ten-minute excerpt of Tutorials 1-3 ends before the satellite reaches the
% shadow. The orbit is read into the GCRS (m, m/s) and the attitude into an
% |AttitudeProvider| (BODY -> GCRS), as in Tutorial 2.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
md = f('samples', 'multiday');
[ep, r, v] = NonGravForces.io.readGNV1B(fullfile(md, 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(fullfile(md, 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
t = double(ep.sec(:)' - ep.sec(1));
fprintf('%d epochs from %s GPS, %.0f s\n', numel(t), string(datetime(ep.day(1) + ep.sec(1) / 86400, 'ConvertFrom', 'modifiedjuliandate')), t(end));

%% Step 2 - The Sun and the Moon
% |referenceSystem.geocentricBodyState("sun", epoch)| returns the Earth -> Sun position and
% velocity in the GCRS (m, m/s), from an analytic theory; against JPL's DE440 it differs by
% 4e-9 of the solar radiation pressure (CLAUDE.md, section 7), which is why the standard preset
% uses it and only the precise preset reads DE440. The ecliptic normal, needed by SOLAARS-CF, is
% r_sun x v_sun, normalised. The BETA ANGLE - the elevation of the Sun above the orbit plane -
% decides how long the shadows are: near 0 the satellite crosses the full shadow each orbit, above
% about 70 deg it stays in sunlight.
[rs, vs] = referenceSystem.geocentricBodyState("sun", ep);
rm = referenceSystem.geocentricBodyState("moon", ep);
en = cross(rs, vs, 1); en = en ./ vecnorm(en);
h = cross(r, v, 1); h = h ./ vecnorm(h);
beta = asind(dot(h, rs ./ vecnorm(rs)));
c = NonGravForces.const();
fprintf('Sun at %.5f au, Moon at %.0f km; beta angle %.2f .. %.2f deg\n', mean(vecnorm(rs)) / c.au, mean(vecnorm(rm)) / 1e3, min(beta), max(beta));

%% Step 3 - The shadow function across the penumbra
% |env.eclipseFactor(r, r_sun)| returns nu in [0, 1], the fraction of the solar disc seen from the
% satellite. Three Earth-shadow models:
%
% * "cylindrical": 0 or 1, the shadow of a cylinder - no penumbra at all;
% * "conical" (the standard preset): two discs, the Sun's and the Earth's, overlapping
%   (Montenbruck & Gill 2000, 3.4.2) - a geometric penumbra of a few seconds;
% * "solaarsCF" (the precise preset): Robertson's (2015) fit to a ray-traced shadow with the
%   Earth's oblateness and the ATMOSPHERE (refraction, absorption), which makes the real shadow
%   wider and its edge softer. It needs the ecliptic normal.
%
% Below, for each model, the times at which nu falls below 0.99 and below 0.01 on this entry, and
% the penumbra between them. Against ACT1B, the corrected SOLAARS-CF puts both shadow edges within
% about 1 s, the conical model's entry about 4 s late (CLAUDE.md, PR #129).
models = ["cylindrical" "conical" "solaarsCF"];
nu = zeros(numel(models), numel(t));
for k = 1:numel(models)
    nu(k, :) = NonGravForces.env.eclipseFactor(r, rs, 'Model', models(k), 'EclipticNormal', en);
end
fprintf('%-12s %12s %12s %10s\n', 'model', 'nu < 0.99', 'nu < 0.01', 'penumbra');
for k = 1:numel(models)
    k1 = find(nu(k, :) < 0.99, 1); k0 = find(nu(k, :) < 0.01, 1);
    if isempty(k1) || isempty(k0)
        fprintf('%-12s   no shadow entry in this arc\n', models(k));
    else
        fprintf('%-12s %10.0f s %10.0f s %8.0f s\n', models(k), t(k1), t(k0), t(k0) - t(k1));
    end
end
kE = find(nu(2, :) < 1, 1);
if isempty(kE), win = 1:numel(t); else, win = max(1, kE - 40):min(numel(t), kE + 40); end
figure
plot(t(win), nu(:, win)), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('\nu'), ylim([-0.05 1.05])
legend(models, 'Location', 'southwest'), title('Shadow entry of GRACE-C: three models')

%% Step 4 - The Moon's shadow
% With |'Moon'| the factor also includes a lunar eclipse (always geometric-conical) and multiplies
% the Earth's and the Moon's factors; |info.moon| is the lunar factor alone. Lunar eclipses of a
% LEO satellite are rare and short - the option exists (|srp.lunarEclipse|), the presets leave it
% off. On this arc:
[nuL, infoL] = NonGravForces.env.eclipseFactor(r, rs, 'Moon', rm);
fprintf('lunar factor along the arc: min %.6f; with the Moon the conical nu changes by at most %.1e\n', ...
    min(infoL.moon), max(abs(nuL - nu(2, :))));

%% Step 5 - The source list: Sun, albedo and infrared in the body frame
% |env.radiationSources(epoch, r, q, options)| is the one place where solar and Earth radiation
% become SOURCES (principle P3): the kernel later computes the force on each panel from every
% source without knowing which effect it is. For K epochs it returns M = 1 + 2 x tiles sources:
% the Sun (with its shadow factor), one visible-band (albedo) source per Earth tile and one
% infrared source per tile. |src.dir| (3 x M x K) are unit directions FROM the satellite TO the
% source in the BODY frame; |src.P| (M x K) the radiation pressures in Pa - an irradiance divided
% by the speed of light, so c P is in W/m^2. Unlit or invisible tiles carry 0. Every 10th epoch
% is enough here.
idx = 1:10:numel(t);
e10 = referenceSystem.daySecondToEpoch([double(ep.day(idx)) double(ep.sec(idx))], string(ep.scale));
q10 = att.at(e10);
opt = NonGravForces.Options("standard");
[src, info] = NonGravForces.env.radiationSources(e10, r(:, idx), q10, opt);
cl = info.speedOfLight_m_s; Mt = (src.M - 1) / 2;
Esun = cl * src.P(1, :); Ealb = cl * sum(src.P(2:1 + Mt, :), 1); Eir = cl * sum(src.P(2 + Mt:end, :), 1);
fprintf('options: srp.eclipse %s, erp.backend %s, erp.grid %s; M = %d sources (%d vis, %d ir), %d epochs\n', ...
    info.eclipseModel, info.erpBackend, mat2str(info.grid), src.M, nnz(src.band == "vis"), nnz(src.band == "ir"), src.K);
fprintf('sunlit: Sun %.1f W/m^2; summed over the tiles: albedo up to %.1f W/m^2, infrared %.1f .. %.1f W/m^2\n', ...
    max(Esun), max(Ealb), min(Eir), max(Eir));
figure
plot(t(idx), Esun, t(idx), Ealb, t(idx), Eir), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('W/m^2 (c P summed over the sources)')
legend("Sun", "Earth albedo", "Earth infrared", 'Location', 'east'), title('Radiation reaching GRACE-C')

%% Step 6 - The Earth's cap, tile by tile (Knocke)
% |env.earthRadiance("knocke", epoch, grid, r)| tiles the cap of the Earth visible from the
% satellite in RINGS of equal central angle around the sub-satellite point and SECTORS of equal
% azimuth (|erp.grid| = [rings sectors]). Rings are needed because the sub-satellite tiles
% dominate (the flux falls with 1/distance^2): a global latitude-longitude grid of the same tile
% count misses the flux by half at 490 km. Each tile carries Knocke, Ries & Tapley's (1988) zonal
% albedo and emissivity, with an annual term phased from 22 December. The tile areas sum to the
% cap area exactly.
e1 = referenceSystem.daySecondToEpoch([double(ep.day(1)) double(ep.sec(1))], string(ep.scale));
fprintf('%-8s %6s %16s %16s\n', 'grid', 'tiles', 'sum of areas', 'cap area (km^2)');
for g = {[6 12], [12 24], [24 48]}
    tl = NonGravForces.env.earthRadiance("knocke", e1, g{1}, r(:, 1));
    capArea = 2 * pi * tl.radius_m^2 * (1 - cos(tl.capAngle_rad));
    fprintf('%-8s %6d %16.6e %16.6e\n', mat2str(g{1}), tl.M, sum(tl.area_m2) / 1e6, capArea / 1e6);
end
tK = NonGravForces.env.earthRadiance("knocke", e1, [12 24], r(:, 1));
fprintf('Knocke [12 24]: area-weighted albedo %.3f, emissivity %.3f over the cap (half angle %.1f deg)\n', ...
    sum(tK.albedo .* tK.area_m2) / sum(tK.area_m2), sum(tK.emissivity .* tK.area_m2) / sum(tK.area_m2), rad2deg(tK.capAngle_rad));

%% Step 7 - Measured Earth radiation: CERES EBAF
% NASA's CERES EBAF 4.2.1 gives monthly 1-degree maps of the reflected shortwave (rsut), the
% emitted longwave (rlut) and the incoming solar (rsdt) flux at the top of the atmosphere, in
% W/m^2. |io.importCERES| reads the three files (December 2019 and January 2020 are committed).
% With erp.backend = "ceres" (the precise preset) the same cap tiles take the measured albedo
% rsut/rsdt and the measured longwave exitance instead of Knocke's zonal fit, interpolated in
% space and between the two monthly nodes in time. The tiles must be found in the Earth-fixed
% frame, so the ceres backend needs the ITRS -> GCRS rotation as well.
cer = NonGravForces.io.importCERES(string(f('ceres', ["rsut" "rlut" "rsdt"] + "_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc")));
w = sind(cer.lat_deg(:) + 0.5) - sind(cer.lat_deg(:) - 0.5);          % area weight of a 1-degree row
gm = @(x) sum(double(x) .* w, 'all') / (360 * sum(w));
fprintf('%s %s, %d months; January 2020 global means: SW up %.1f, LW up %.1f, SW down %.1f W/m^2 (albedo %.3f)\n', ...
    cer.product, cer.version, numel(cer.mjd), gm(cer.rsut(:, :, 2)), gm(cer.rlut(:, :, 2)), gm(cer.rsdt(:, :, 2)), ...
    gm(cer.rsut(:, :, 2)) / gm(cer.rsdt(:, :, 2)));
R10 = referenceSystem.rotationItrsToGcrs(e10, 'Eop', eop);
tC = NonGravForces.env.earthRadiance("ceres", e1, [12 24], r(:, 1), 'Ceres', cer, 'RotationItrsToGcrs', R10(:, :, 1));
fprintf('the cap at the first epoch: area-weighted albedo Knocke %.3f, CERES %.3f\n', ...
    sum(tK.albedo .* tK.area_m2) / sum(tK.area_m2), sum(tC.albedo .* tC.area_m2) / sum(tC.area_m2));
optC = NonGravForces.Options("standard", "erp.backend", "ceres");
srcC = NonGravForces.env.radiationSources(e10, r(:, idx), q10, optC, 'Ceres', cer, 'RotationItrsToGcrs', R10);
EalbC = cl * sum(srcC.P(2:1 + Mt, :), 1); EirC = cl * sum(srcC.P(2 + Mt:end, :), 1);
fprintf('along the arc, CERES / Knocke: albedo irradiance %.3f, infrared %.3f (sums over the arc)\n', sum(EalbC) / sum(Ealb), sum(EirC) / sum(Eir));
figure
subplot(2, 1, 1)
imagesc(cer.lon_deg, cer.lat_deg, double(cer.rsut(:, :, 2)) ./ double(cer.rsdt(:, :, 2))), set(gca, 'YDir', 'normal')
caxis([0 1]), colorbar, ylabel('latitude (deg)'), title('CERES EBAF albedo, January 2020')
subplot(2, 1, 2)
plot(t(idx), Ealb, t(idx), EalbC, t(idx), Eir, t(idx), EirC), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('W/m^2')
legend("albedo, Knocke", "albedo, CERES", "infrared, Knocke", "infrared, CERES", 'Location', 'eastoutside')

%% Step 8 - The Earth's magnetic field: IGRF-14 against MAG1B
% Not a force on the centre of mass, but a torque: a magnetorquer current in the Earth's field
% turns the satellite (Tutorial 2, Step 5). |io.igrfCoefficients(file, 'Epoch', year)| reads the
% IGRF-14 table and evaluates the coefficients at one epoch (nanotesla); |env.magneticField|
% synthesises the main field at Earth-fixed positions in TESLA and, with |'Rotation'|, rotates it
% into any frame - here ITRS -> GCRS -> body, to compare with the magnetometer in the SRF. The
% magnetometer is raw and the model is the main field only (no crust above degree 13, no
% ionosphere or magnetosphere), so the intensities agree to a few parts per thousand, not exactly.
% MAG1B samples in pairs every second; each is matched to the GNV1B second it rounds to.
cf = NonGravForces.io.igrfCoefficients(f('igrf', 'igrf14coeffs.txt'), 'Epoch', 2020.0);
[eG, rI] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Frame', "ITRS");
[~, qG] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
[eM, fld] = NonGravForces.io.readMAG1B(f('samples', 'MAG1B_2020-01-01_C_04.txt'));
[tf, loc] = ismember(round(eM.sec(:)'), eG.sec(:)');
k = loc(tf);
eK = referenceSystem.daySecondToEpoch([double(eG.day(k)) double(eG.sec(k))], string(eG.scale));
Rig = referenceSystem.rotationItrsToGcrs(eK, 'Eop', eop);                 % ITRS -> GCRS
Rbg = referenceSystem.rotationFromQuaternion(qG(:, k));                   % BODY -> GCRS
[Bb, infoB] = NonGravForces.env.magneticField(cf, rI(:, k), 'Rotation', pagemtimes(Rbg, 'transpose', Rig, 'none'));
B1 = NonGravForces.env.magneticField(cf, rI(:, k), 'MaxDegree', 1);
ratio = fld.magnitude_T(tf) ./ infoB.intensity_T;
fprintf('IGRF-%d at %.1f, degree %d: |B| %.2f .. %.2f microtesla; MAG1B / IGRF intensity %.4f .. %.4f\n', ...
    14, cf.year, infoB.maxDegree, 1e6 * min(infoB.intensity_T), 1e6 * max(infoB.intensity_T), min(ratio), max(ratio));
fprintf('the dipole alone differs from degree 13 by up to %.1f %% of the intensity\n', 100 * max(vecnorm(B1 - infoB.itrs_T) ./ infoB.intensity_T));
tm = double(eM.sec(tf)' - eM.sec(1));
figure
plot(tm, 1e6 * fld.B_srf_T(:, tf), '.', tm, 1e6 * Bb, '-'), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('\muT')
legend("MAG1B x", "MAG1B y", "MAG1B z", "IGRF x", "IGRF y", "IGRF z", 'Location', 'eastoutside')
title('Magnetic field in the body frame: magnetometer and IGRF-14')
