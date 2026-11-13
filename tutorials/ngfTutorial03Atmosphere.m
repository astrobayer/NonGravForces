%% Tutorial 3 - The atmosphere: density, its observed scale and the wind
% Drag is the largest non-gravitational force on GRACE-FO and the least certain one, because it
% is proportional to the thermosphere density, which no model predicts to better than tens of
% percent. This tutorial takes the atmosphere apart along ten minutes of GRACE-C's real orbit of
% 2020-01-01: the density and composition NRLMSISE-00 gives, the deliberately simple exponential
% atmosphere, the OBSERVED daily scale the standard preset multiplies the model with (from Swarm
% where Swarm flew, predicted from the F30 solar flux elsewhere), an independent check of that
% scale against GRACE-FO's own accelerometer densities, the co-rotating air, and the
% accommodation coefficient that describes how the gas leaves the surface.
%
% *What you will learn*
%
% * how to obtain the GEODETIC position a density model needs, and why the ellipsoid matters;
% * what |env.density| returns: total density, eight species, two temperatures, the drivers;
% * why the exponential atmosphere returns NaN species on purpose;
% * how |env.densityScale| turns Swarm densities into a daily factor, interpolated in height;
% * how |env.densityDriverScale| predicts the factor for any date, and how "combined" chooses
%   between the two per UTC day;
% * how |env.wind| removes the co-rotation of the atmosphere from the velocity;
% * what |kernel.accommodation| computes, and why the computed models warn at GRACE-FO.
%
% *Requirements:* none - every file this tutorial reads is committed in |data/|. NRLMSISE-00
% runs through the MEX that |setupNonGravForces| builds. Run time: well under a minute.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Where the satellite is: GCRS state and geodetic position
% A density model wants latitude, longitude and height above the ELLIPSOID (WGS-84): a spherical
% Earth is 21 km wrong at the poles, which is a factor 1.4 in density. The geodetic position comes
% from the Earth-fixed (ITRS) position through |referenceSystem.cartesianToGeodetic| with the
% ellipsoid named explicitly (its default is GRS80). The drag needs the inertial (GCRS) state as
% well. Both come from |io.readGNV1B| (Tutorial 2). |lla| is 3 x K: latitude and longitude in
% rad, height in m.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
gnv = f('samples', 'GNV1B_2020-01-01_C_04.txt');
[ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
[~, rI] = NonGravForces.io.readGNV1B(gnv, 'Frame', "ITRS");
lla = referenceSystem.cartesianToGeodetic(rI', 'Ellipsoid', "WGS84")';
t = double(ep.sec(:)' - ep.sec(1));
fprintf('%d epochs: latitude %.1f .. %.1f deg, geodetic height %.1f .. %.1f km\n', numel(t), ...
    rad2deg(min(lla(1, :))), rad2deg(max(lla(1, :))), min(lla(3, :)) / 1e3, max(lla(3, :)) / 1e3);
hSphere = vecnorm(rI) - referenceSystem.ellipsoidParameters("WGS84").a;
fprintf('height above a SPHERE of the equatorial radius instead: %.1f .. %.1f km\n', min(hSphere) / 1e3, max(hSphere) / 1e3);

%% Step 2 - NRLMSISE-00 along the arc
% |env.density("msis00", epoch, lla, sw)| evaluates NRLMSISE-00 (Picone et al. 2002) at every
% epoch. The drivers are read from the space-weather struct: the observed F10.7 of the previous
% day, its 81-day centred mean and the 7-element ap history; |info.drivers| shows what was used.
% Outputs: |rho| (1 x K, kg/m^3, total mass density, anomalous oxygen included), |species|
% (8 x K, number densities in 1/m^3, in the order of |NonGravForces.const().speciesName|) and
% |T| (2 x K, K: exospheric temperature and temperature at the satellite). The other backends,
% NRLMSIS 2.1 ("msis21") and CNES's DTM2020 ("dtm2020", academic licence), are compiled from
% sources that |setupNonGravForces| fetches; the loop below uses them when they are built and
% prints the refusal otherwise.
[rho, sp, T, info] = NonGravForces.env.density("msis00", ep, lla, sw);
fprintf('msis00: rho %.3e .. %.3e kg/m^3, exospheric T %.0f K, T at h %.0f K\n', min(rho), max(rho), mean(T(1, :)), mean(T(2, :)));
fprintf('drivers (%s): F10.7 %.1f sfu, F10.7A %.1f sfu, daily Ap %g\n', info.drivers.source, ...
    info.drivers.f107(1), info.drivers.f107a(1), info.drivers.ap(1, 1));
for backend = ["msis21" "dtm2020"]
    try
        rhoB = NonGravForces.env.density(backend, ep, lla, sw);
        fprintf('%s / msis00 along the arc: mean %.3f\n', backend, mean(rhoB ./ rho));
    catch err
        fprintf('%s is not available here: %s\n', backend, err.identifier);
    end
end
figure
yyaxis left, semilogy(t, rho), ylabel('\rho (kg/m^3)')
yyaxis right, plot(t, lla(3, :) / 1e3), ylabel('geodetic height (km)')
grid on, xlabel('seconds since 2020-01-01 00:00:00 GPS'), title('NRLMSISE-00 along GRACE-C')

%% Step 3 - Composition: which gas the satellite flies through
% At GRACE-FO's height in solar minimum the air is mostly atomic oxygen and helium, not "air".
% The drag kernel works per species (each has its own thermal speed, so its own speed ratio),
% which is why the species matter, not just the total. Below: the mean number density of each
% species, its share of the mass (number density times molar mass over Avogadro's number) and
% the mean molar mass of the mixture.
c = NonGravForces.const();
massDensity = sp .* (c.speciesMass(:) / c.NA);                       % 8 x K, kg/m^3
share = sum(massDensity, 2) / sum(massDensity, 'all');
fprintf('%-12s %14s %12s\n', 'species', 'n (1/m^3)', 'mass share');
for k = 1:numel(c.speciesName)
    fprintf('%-12s %14.3e %12.4f\n', c.speciesName(k), mean(sp(k, :)), share(k));
end
fprintf('mean molar mass %.2f g/mol; the species re-form the total density to %.1e (molar masses differ)\n', ...
    1e3 * sum(massDensity, 'all') / sum(sp, 'all') * c.NA, max(abs(sum(massDensity, 1) ./ rho - 1)));
figure
semilogy(t, sp(1:7, :)), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('number density (1/m^3)')
legend(c.speciesName(1:7), 'Location', 'eastoutside'), title('NRLMSISE-00 composition')

%% Step 4 - The exponential atmosphere, and why its species are NaN
% rho = rho0 exp(-(h - h0)/H) is the fast preset's atmosphere: a sizing model with the caller's
% three parameters and no drivers. It models neither composition nor temperature, and says so by
% returning NaN for both - a unit that silently invented them would make every downstream
% gas-surface model look valid. Here its parameters are taken from NRLMSISE-00 itself: rho0 the
% arc's mean density at the arc's mean height h0, and H a local scale height from two MSIS
% evaluations 25 km above and below h0 at the first epoch. The difference that remains along the
% arc is everything an exponential cannot know: latitude, local time, the day-night bulge.
h0 = mean(lla(3, :));
e2 = referenceSystem.daySecondToEpoch([double(ep.day([1 1])) double(ep.sec([1 1]))], string(ep.scale));
rho2 = NonGravForces.env.density("msis00", e2, [lla([1 1], 1)'; lla([2 2], 1)'; h0 + [-25e3 25e3]], sw);
p = struct('rho0_kg_m3', mean(rho), 'h0_m', h0, 'H_m', 50e3 / log(rho2(1) / rho2(2)));
[rhoE, spE, TE] = NonGravForces.env.density("exponential", ep, lla, [], 'Exponential', p);
fprintf('exponential: H = %.1f km; species all NaN: %d, temperatures all NaN: %d\n', p.H_m / 1e3, all(isnan(spE(:))), all(isnan(TE(:))));
fprintf('exponential / msis00 along the arc: %.3f .. %.3f\n', min(rhoE ./ rho), max(rhoE ./ rho));

%% Step 5 - The observed daily scale from Swarm
% NRLMSISE-00 overestimates the density in solar minimum - in 2020-01 by about a factor of two at
% GRACE-FO's height (HANDOVER O-19). The standard preset therefore multiplies it by a factor
% OBSERVED per UTC day. |env.densityScale| reads ESA's Swarm A, B and C POD densities and forms,
% per satellite and day, f = sum(rho_Swarm) / sum(rho_model) over the satellite's valid records
% (every 10th by default), with the model evaluated at Swarm's own positions and the same space
% weather. Swarm A and C fly side by side below GRACE-FO, B above it; |S.at(epoch, h)| interpolates
% log f linearly in height between the lower pair (their geometric mean) and B, and refuses a
% height further outside that bracket than its width. A daily factor removes only the part of the
% model error common to all local times of a day.
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
S = NonGravForces.env.densityScale(swarm, sw);
fprintf('kind %s, model %s, stride %d\n', S.kind, S.model, S.stride);
fprintf('%-6s', 'day'); fprintf('   Swarm %s (f, km)', S.satellite); fprintf('\n');
for dd = 1:numel(S.day)
    fprintf('%-6d', S.day(dd)); fprintf('   %6.3f  %6.1f    ', [S.factor(:, dd)'; S.height_m(:, dd)' / 1e3]); fprintf('\n');
end
hq = 420e3:5e3:560e3;
eQ = referenceSystem.mjdToEpoch(repmat(58849.5, numel(hq), 1), "UTC");
figure
plot(hq / 1e3, S.at(eQ, hq), [S.hLow_m(1) S.hHigh_m(1)] / 1e3, [S.low(1) S.high(1)], 'o'), grid on
xline(mean(lla(3, :)) / 1e3, '--', 'GRACE-C');
xlabel('geodetic height (km)'), ylabel('factor on NRLMSISE-00'), title('Swarm density scale, 2020-01-01')
legend("S.at (log-linear in height)", "the bracket: A+C and B", 'Location', 'best')

%% Step 6 - The predicted scale from the drivers
% Swarm flies since late 2013 and has gaps. |env.densityDriverScale(sw, cls)| PREDICTS the factor
% for any date from drivers that exist since 1951: the 81-day mean of CLS's F30 ending the day
% before, the previous day's deviation from it, the day's Ap, annual and semiannual harmonics and
% a square term - a regression fitted to Swarm, CHAMP, GRACE and GRACE-FO densities over 2000-2026
% (its held-out error is about 18 % per day, its help says). One factor per UTC day, independent
% of height. |'Offset'| adds a satellite's own level in log f (-0.041 for GRACE-FO 1). Below, the
% predicted factor over four months beside Swarm's measured factor at GRACE-C's height.
D = NonGravForces.env.densityDriverScale(sw, cls);
DG = NonGravForces.env.densityDriverScale(sw, cls, 'Offset', -0.041);
days = 58818:58938;                                                   % 2019-12-01 .. 2020-03-30
eD = referenceSystem.mjdToEpoch(days' + 0.5, "UTC");
hD = 490e3 * ones(1, numel(days));
fD = D.at(eD, hD); fG = DG.at(eD, hD);
eS = referenceSystem.mjdToEpoch(S.day' + 0.5, "UTC");
fS = S.at(eS, 490e3 * ones(1, numel(S.day)));
fprintf('driver scale (%s proxy): %.3f .. %.3f over 2019-12..2020-03; Swarm 2020-01-01..03 at 490 km: %s\n', ...
    D.proxy, min(fD), max(fD), mat2str(fS, 3));
fprintf('on 2020-01-01: driver %.3f, with the GRACE-FO offset %.3f, Swarm %.3f\n', fD(days == 58849), fG(days == 58849), fS(1));
tD = datetime(days, 'ConvertFrom', 'modifiedjuliandate');
figure
plot(tD, fD, tD, fG, datetime(S.day, 'ConvertFrom', 'modifiedjuliandate'), fS, 'ko'), grid on
ylabel('factor on NRLMSISE-00'), title('Predicted (driver) and measured (Swarm) density scale')
legend("driver", "driver, GRACE-FO offset", "Swarm at 490 km", 'Location', 'best')

%% Step 7 - "combined": the standard preset's scale
% The standard preset uses drag.densityScale = "combined": Swarm's factor on every UTC day that
% has one, the driver prediction on every other day, never mixed within a day. It is
% |env.densityScale| with the driver scale as |'Fallback'|; |C.source(epoch)| says which one
% serves each epoch. The committed Swarm excerpt ends on 2020-01-03, so from 2020-01-04 on the
% driver takes over. A standard |Model| refuses to evaluate without this object as its
% |'DensityScale'| (Tutorial 1, Step 6) rather than run unscaled. Applied to the arc, the scaled
% density is the model density times |C.at(epoch, h)|.
C = NonGravForces.env.densityScale(swarm, sw, 'Fallback', D);
e5 = referenceSystem.mjdToEpoch((58849:58853)' + 0.5, "UTC");
fprintf('kind %s; 2020-01-01..05 at 490 km: %s from %s\n', C.kind, mat2str(C.at(e5, 490e3 * ones(1, 5)), 3), strjoin(C.source(e5), ", "));
fArc = C.at(ep, lla(3, :));
srcArc = C.source(ep);                % the day of an epoch is its UTC day: GPS 00:00:00 is 2019-12-31 23:59:42 UTC
fprintf('along the arc the factor comes from: %s (%d epochs from the driver)\n', strjoin(unique(srcArc), ", "), nnz(srcArc == "driver"));
fprintf('along the arc: factor %.3f .. %.3f, scaled density %.3e kg/m^3 on average (unscaled %.3e)\n', ...
    min(fArc), max(fArc), mean(rho .* fArc), mean(rho));
opt = NonGravForces.Options("standard");
fprintf('the standard preset: drag.density %s, drag.densityScale %s\n', opt.drag.density, opt.drag.densityScale);

%% Step 8 - An independent check: GRACE-FO's own densities
% TU Delft derives thermosphere densities from the GRACE-FO 1 (GRACE-C) accelerometer itself
% (ESA TOLEOS, 10 s); six hours of 2020-01-01 are committed. They are not the input of the scale
% above, so they test it: if Swarm's factor is right at GRACE-FO, the TOLEOS density over the
% NRLMSISE-00 density at the same positions should be close to |C.at| there. They are not
% expected to agree exactly - TOLEOS densities rest on their own aerodynamic model of the
% satellite, and the GRACE-FO offset of Step 6 is this very difference measured over years.
[eT, llaT, rhoT, auxT] = NonGravForces.io.readSwarmDensity(f('samples', 'toleos', 'GF_OPER_DNS1ACC_2_2020-01-01_00-06.csv'));
use = find(auxT.valid); use = use(1:6:end);                           % one record a minute
eU = referenceSystem.daySecondToEpoch([double(eT.day(use)) double(eT.sec(use))], "UTC");
rhoM = NonGravForces.env.density("msis00", eU, llaT(:, use), sw);
fU = C.at(eU, llaT(3, use));
fprintf('%s (%s), %d records: TOLEOS / msis00 = %.3f; the combined scale there %.3f; TOLEOS / scaled msis00 = %.3f\n', ...
    auxT.satellite, auxT.dataset, numel(use), sum(rhoT(use)) / sum(rhoM), mean(fU), sum(rhoT(use)) / sum(rhoM .* fU));
tU = double(eT.sec(use)) / 3600;
figure
semilogy(tU, rhoT(use), tU, rhoM, tU, rhoM .* fU), grid on
xlabel('hours since 2020-01-01 00:00 UTC'), ylabel('kg/m^3'), title('GRACE-FO 1: observed and modelled density')
legend("TOLEOS (accelerometer)", "NRLMSISE-00", "NRLMSISE-00 x combined scale", 'Location', 'best')

%% Step 9 - The wind: the air co-rotates with the Earth
% The drag acts along the velocity of the satellite RELATIVE TO THE GAS, not along its inertial
% velocity. The atmosphere co-rotates with the Earth, up to about 470 m/s at the equator, which
% turns the relative velocity by a few degrees cross-track. |env.wind("corotation", ...)| removes
% it as v - Rdot (R' r), with R and Rdot from the same |referenceSystem| rotation that made the
% GCRS state - never a hand-written omega x r, which would miss polar motion and the length of
% day. Inputs and outputs are GCRS, m and m/s. "hwm14" adds NRL's Horizontal Wind Model 14
% (quiet-time winds of 40-80 m/s); it is compiled from sources |setupNonGravForces| fetches and
% is the precise preset's choice, so it is not used in this offline tutorial.
[vrel, winfo] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
turn = atan2d(vecnorm(cross(v, vrel)), dot(v, vrel));
fprintf('co-rotation removed: %.1f .. %.1f m/s; |v_rel| %.1f m/s against |v| %.1f m/s; direction turned by up to %.2f deg\n', ...
    min(vecnorm(winfo.corotation_m_s)), max(vecnorm(winfo.corotation_m_s)), mean(winfo.speed_m_s), mean(vecnorm(v)), max(turn));
figure
plot(rad2deg(lla(1, :)), vecnorm(winfo.corotation_m_s), '.'), grid on
xlabel('geodetic latitude (deg)'), ylabel('m/s'), title('Speed of the co-rotating air at GRACE-C')

%% Step 10 - Accommodation: how the gas leaves the surface
% A molecule that hits a panel is re-emitted with part of the wall's temperature: the energy
% accommodation coefficient alpha (0 = it keeps its energy, 1 = it leaves at the wall
% temperature). |kernel.accommodation| returns it per epoch from the gas the density model gives.
% "fixed" is the caller's value - all three presets take March et al. (2021)'s fitted 0.85.
% "pilinski2010" is a fitted Langmuir law in n_O T; "walker2014" and "sesam" describe a surface
% partly covered by adsorbed oxygen (coverage theta), whose effective coefficient is
% theta + (1 - theta) alpha_substrate. The computed models were fitted at lower heights; at
% GRACE-FO they extrapolate, |info.insideFit| says so per epoch, and they WARN
% (|NonGravForces:accommodation:outsideFit|). The warning is silenced here because it is the
% point of the step, then restored. SESAM needs the relative speed: its pressure is a ram pressure.
ws = warning('off', 'NonGravForces:accommodation:outsideFit');
[aF, iF] = NonGravForces.kernel.accommodation("fixed", sp, T(2, :), 'Alpha', 0.85);
[aP, iP] = NonGravForces.kernel.accommodation("pilinski2010", sp, T(2, :));
[aW, iW] = NonGravForces.kernel.accommodation("walker2014", sp, T(2, :), 'SurfaceMass_amu', 263.3223);
[aS, iS] = NonGravForces.kernel.accommodation("sesam", sp, T(2, :), 'RelativeSpeed_ms', winfo.speed_m_s);
warning(ws);
effW = aW.coverage + (1 - aW.coverage) .* aW.alphaSubstrate;
effS = aS.coverage + (1 - aS.coverage) .* aS.alphaSubstrate;
fprintf('%-13s %10s %12s\n', 'model', 'mean alpha', 'inside fit');
row = @(name, alpha, inside) fprintf('%-13s %10.3f %11.0f%%\n', name, mean(alpha), 100 * mean(inside));
row('fixed', aF.alpha, iF.insideFit)
row('pilinski2010', aP.alpha, iP.insideFit)
row('walker2014', effW, iW.insideFit)
row('sesam', effS, iS.insideFit)
fprintf('walker2014 coverage %.3f, sesam coverage %.3f (mean)\n', mean(aW.coverage), mean(aS.coverage));
