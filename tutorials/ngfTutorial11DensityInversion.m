%% Tutorial 11 - Density from the accelerometer
% An accelerometer on a low satellite measures, after calibration, the sum of drag, radiation
% pressure and thermal recoil. Subtract the modelled radiation and recoil, and what is left is the
% drag - and the drag is proportional to the thermosphere's density. |app.densityInversion| turns a
% calibrated acceleration into a density along the orbit. It is a DIVISION, not a search: the gas
% force is exactly linear in every species density (Tutorial 6, Step 10), so one evaluation of the
% model at its own density gives the whole sensitivity. This tutorial shows what the unit does and
% what it reports, proves it on a planted case, applies it to the committed ten minutes of ACT1B,
% and compares the result with the model density with and without the observed Swarm scale.
%
% *What you will learn*
%
% * the inversion s = <a_meas - a_rad - a_thermal, u> / |a_drag|, rho = s rho_model, with u the
%   direction of the modelled drag, and the outputs of |app.densityInversion|;
% * why the inverted density does not depend on the a priori density of the model;
% * the quality flag |residualPerp_m_s2|: what the density cannot explain;
% * what a real accelerometer arc needs first (bias, scale), and why a ten-minute arc cannot fix
%   the level of the density, only its shape;
% * how NRLMSISE-00 and the Swarm-scaled density compare on 2020-01-01.
%
% *Requirements:* none - every file is committed in |data/|. Run time: about twenty seconds.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Two models: NRLMSISE-00 with and without the observed scale
% The standard preset multiplies NRLMSISE-00 by an OBSERVED daily density scale - Swarm's density
% over the model's, interpolated in height to the satellite (Tutorial 3). A second model is the
% same preset with |drag.densityScale = "none"|: plain NRLMSISE-00. |info.rho| (kg/m^3) is the
% density each model flew through, |info.densityScale| the factor applied (1 for "none").
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
[ep, r, vel] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
macro = NonGravForces.geom.graceFO("C");
Ms = NonGravForces.Model(macro, NonGravForces.Options("standard"), 'Attitude', att, 'Eop', eop, ...
    'SpaceWeather', sw, 'DensityScale', S);
Mu = NonGravForces.Model(macro, NonGravForces.Options("standard").set("drag.densityScale", "none"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw);
K = 600; idx = 1:K;
e = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
rr = r(:, idx); vv = vel(:, idx); qq = att.at(e);
t = (ep.sec(idx) - ep.sec(1))';
[aS, iS] = Ms.accel(e, rr, vv, qq);
[aU, iU] = Mu.accel(e, rr, vv, qq);
fprintf('NRLMSISE-00 density, mean          %.4e kg/m^3\n', mean(iU.rho));
fprintf('Swarm-scaled density, mean         %.4e kg/m^3 (scale %.3f)\n', mean(iS.rho), mean(iS.densityScale));
fprintf('modelled drag, mean |a_drag|       %.3e (scaled) / %.3e (unscaled) m/s^2\n', mean(vecnorm(iS.drag)), mean(vecnorm(iU.drag)));
fprintf('radiation + thermal, mean |a|      %.3e m/s^2\n', mean(vecnorm(iS.radiation + iS.thermal)));

%% Step 2 - The call, on a case whose answer is known
% |T = app.densityInversion(model, epoch, r, v, q, aMeasured)|: |epoch| a referenceSystem epoch
% struct of K epochs, |r| and |v| GCRS (m, m/s, 3 x K), |q| the BODY -> GCRS quaternions (or [] for
% the model's attitude provider), |aMeasured| the CALIBRATED acceleration (3 x K, m/s^2, body frame
% unless |'Frame', "gcrs"|). It evaluates the model once per chunk and returns a timetable:
%
% * |rho_kg_m3| the inverted density, |rho_model_kg_m3| the model's own (with its scale), |scale|
%   their ratio;
% * |dragProjected_m_s2| the measurement minus radiation and recoil, projected on the modelled drag
%   direction u, and |dragModel_m_s2| the modelled drag magnitude;
% * |residualPerp_m_s2| the norm of what is left ACROSS u - the quality flag;
% * |usable|: false where the modelled drag is below |'MinDrag'| (default 1e-12 m/s^2); rho is NaN
%   there and the scale is never clipped, so a negative value stays visible.
%
% The only data whose truth is known are the model's own: a "measurement" whose gas part is
% scaled by 0.7 must return a scale of 0.7 and nothing across u.
planted = 0.7;
aPl = iS.radiation + iS.thermal + planted * iS.drag;
[T, ti] = NonGravForces.app.densityInversion(Ms, e, rr, vv, qq, aPl);
disp(head(T, 3))
fprintf('%d of %d epochs usable, %d chunk(s)\n', ti.usable, ti.K, ti.chunks);
fprintf('recovered scale: max |s - %.1f| = %.1e; max residual across u %.1e m/s^2\n', planted, ...
    max(abs(T.scale - planted)), max(T.residualPerp_m_s2));

%% Step 3 - The answer does not depend on the a priori density
% Because the drag is linear in the density, a model with a different density predicts a drag
% proportionally different along the SAME direction, so the ratio s changes by exactly the inverse
% factor and rho = s rho_model does not change at all. The same planted measurement, inverted with
% plain NRLMSISE-00 instead of the Swarm-scaled model, gives the same density - only its |scale|
% column differs, by the observed factor.
Tu = NonGravForces.app.densityInversion(Mu, e, rr, vv, qq, aPl);
fprintf('inverted rho, scaled vs unscaled a priori: max relative difference %.1e\n', ...
    max(abs(Tu.rho_kg_m3 - T.rho_kg_m3) ./ T.rho_kg_m3));
fprintf('scale column: %.4f (scaled a priori) vs %.4f (unscaled a priori), ratio %.4f\n', ...
    mean(T.scale), mean(Tu.scale), mean(Tu.scale ./ T.scale));

%% Step 4 - What the quality flag sees
% Only the component of the measurement along u carries density. An error ACROSS u - here a
% 5 nm/s^2 cross-track offset, the size of a radiation-model error - barely moves rho (u has a small
% cross-track component) but shows in |residualPerp_m_s2|. An error ALONG the drag - a 1 nm/s^2
% along-track offset, e.g. an uncalibrated bias - goes straight into rho, and nothing can flag it:
% that is the warning in the unit's help, made numerical.
dCross = repmat([0; 5e-9; 0], 1, K);
dAlong = repmat([1e-9; 0; 0], 1, K);
Tc = NonGravForces.app.densityInversion(Ms, e, rr, vv, qq, aPl + dCross);
Ta = NonGravForces.app.densityInversion(Ms, e, rr, vv, qq, aPl + dAlong);
fprintf('5 nm/s^2 cross-track: rho moves %.2f %% (mean), residualPerp %.2e m/s^2 (mean)\n', ...
    100 * mean(Tc.rho_kg_m3 ./ T.rho_kg_m3 - 1), mean(Tc.residualPerp_m_s2));
fprintf('1 nm/s^2 along-track: rho moves %.2f %% (mean), residualPerp %.2e m/s^2 (mean)\n', ...
    100 * mean(Ta.rho_kg_m3 ./ T.rho_kg_m3 - 1), mean(Ta.residualPerp_m_s2));

%% Step 5 - The real measurement: ACT1B, with a bias that must come from somewhere
% ACT1B carries an instrument bias of micro-m/s^2 (Tutorial 8), a thousand times the drag. It must
% be removed before inverting, and over ten minutes it cannot be ESTIMATED independently of the
% density: a constant along-track bias and a constant density offset are the same signal (Step 4).
% Here the bias of each axis is therefore CHOSEN, as the mean of ACT1B minus the Swarm-scaled
% model; the scale is taken as 1 (the month baseline of 0.19.0 found S_xx 0.999 for GRACE-C,
% CLAUDE.md section 5). With that choice the arc mean of the inverted density is the Swarm-scaled
% density BY CONSTRUCTION, and the information the measurement adds is the SHAPE along the arc.
% (ACT1B is CRN-filtered and the model inside the inversion is not; Tutorial 8, Step 5 prints how
% small that difference is over this arc.) A proper bias comes from a calibration over many days,
% where a bias per day and the density's orbital variation separate: Tutorial 9.
[epa, act] = NonGravForces.io.readACC1B(f('samples', 'ACT1B_2020-01-01_C_04.txt'));
[tf, loc] = ismember(ep.sec(idx), epa.sec);
k = find(tf);
eK = referenceSystem.daySecondToEpoch([ep.day(k) ep.sec(k)], ep.scale);
am = act(:, loc(tf));
biasS = mean(am - aS(:, k), 2);
aCal = am - biasS;
Tr = NonGravForces.app.densityInversion(Ms, eK, rr(:, k), vv(:, k), qq(:, k), aCal);
fprintf('ACT1B epochs matched: %d; bias removed [%.4e %.4e %.4e] m/s^2\n', numel(k), biasS);
fprintf('inverted density: mean %.4e kg/m^3, scale %.3f +- %.3f (mean +- std over the arc)\n', ...
    mean(Tr.rho_kg_m3, 'omitnan'), mean(Tr.scale, 'omitnan'), std(Tr.scale, 'omitnan'));
fprintf('quality flag: residualPerp RMS %.2e m/s^2 against a projected drag of %.2e m/s^2\n', ...
    sqrt(mean(Tr.residualPerp_m_s2 .^ 2)), mean(abs(Tr.dragProjected_m_s2)));

%% Step 6 - The bias choice sets the level
% The same arc with the bias chosen against PLAIN NRLMSISE-00 instead: the shape is the same, the
% level follows the model the bias was taken from. This is the circularity of Step 5 made visible -
% and the reason the observed density scale of the standard preset is formed from Swarm's
% independent orbit-derived densities, not from GRACE-FO's accelerometer over an arc this short.
biasU = mean(am - aU(:, k), 2);
TrU = NonGravForces.app.densityInversion(Ms, eK, rr(:, k), vv(:, k), qq(:, k), am - biasU);
fprintf('bias against plain NRLMSISE-00: [%.4e %.4e %.4e] m/s^2\n', biasU);
fprintf('inverted density, mean: %.4e kg/m^3 (bias vs Swarm-scaled) and %.4e kg/m^3 (bias vs NRLMSISE-00)\n', ...
    mean(Tr.rho_kg_m3, 'omitnan'), mean(TrU.rho_kg_m3, 'omitnan'));
fprintf('NRLMSISE-00 %.4e, Swarm-scaled %.4e kg/m^3\n', mean(iU.rho(k)), mean(iS.rho(k)));
fprintf('one nm/s^2 of along-track bias is %.1f %% of the density here\n', 100 * 1e-9 / mean(Tr.dragModel_m_s2));

%% Step 7 - Along the arc: three densities
% The inverted density beside the two model densities. The figure also marks the two thruster
% firings of this excerpt (THR1B; Tutorial 8): their impulses, smeared over +-70 s by the CRN
% filter, enter the measurement and therefore the density, mostly through the cross-track and
% radial axes and so through |residualPerp_m_s2|.
[~, thr] = NonGravForces.io.readTHR1B(f('samples', 'THR1B_2020-01-01_C_04.txt'));
tFire = reshape(thr.epoch.sec(any(thr.fired, 1)), 1, []) - ep.sec(1);
tFire = tFire(tFire >= 0 & tFire <= t(end));
tk = t(k);
figure
subplot(2, 1, 1)
plot(tk, Tr.rho_kg_m3, tk, iS.rho(k), tk, iU.rho(k)), grid on, hold on
xline(tFire, 'r:')
ylabel('kg/m^3'), title('density along the arc, 2020-01-01, GRACE-C')
legend("inverted from ACT1B (bias vs Swarm-scaled)", "NRLMSISE-00 x Swarm scale", "NRLMSISE-00", 'Location', 'best')
subplot(2, 1, 2)
plot(tk, 1e9 * Tr.residualPerp_m_s2, tk, 1e9 * Tr.dragModel_m_s2), grid on, hold on
xline(tFire, 'r:')
ylabel('nm/s^2'), legend("residual across the drag (quality flag)", "modelled drag magnitude")
xlabel('seconds since 2020-01-01 00:00:00 GPS')
