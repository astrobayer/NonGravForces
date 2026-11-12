%% Tutorial 2 - Reading the data: Level-1B, space weather, densities and the data pool
% Every number NonGravForces computes starts from a file somebody else wrote: the GRACE-FO
% Level-1B products of JPL (orbit, attitude, accelerometer, thrusters, tank masses, magnetometer,
% gyros), TU Graz's operational attitude, the space-weather indices of CelesTrak, CLS and Kyoto,
% and the thermosphere densities ESA and TU Delft derive from Swarm, CHAMP, GRACE and GRACE-FO.
% This tutorial reads each of them once, with the reader the toolbox provides, and shows what
% comes out: the frame, the time scale, the units, and the pitfalls each reader guards against.
% Nothing is modelled yet; Tutorials 3 to 5 use what is read here.
%
% *What you will learn*
%
% * how |io.readGNV1B| turns the Earth-fixed GPS orbit into the inertial GCRS, and why it
%   needs Earth orientation parameters;
% * why the SCA1B quaternion is conjugated exactly once, and how |geom.checkAttitudeSense|
%   proves the sense against physics;
% * what ACT1B, THR1B, MAS1B, MAG1B and IMU1B contain, and the TANK-mass pitfall;
% * how the star cameras are compared with the gyros (|geom.attitudeQuality|) and how TU
%   Graz's attitude product is read (|AttitudeProvider.fromItsg|);
% * the four space-weather families (SW-All, SOLFSMY/DTCFILE, Kyoto Dst, CLS F30);
% * how the Swarm and TOLEOS density files are read;
% * where full data sets live (the data pool), what can be fetched (|fetchCatalog|), and how
%   every result records where it came from (|provenance|).
%
% *Requirements:* none - every file this tutorial reads is committed in |data/|. Nothing is
% downloaded. Run time: about half a minute.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - The orbit: GNV1B, Earth-fixed in the file, inertial in the toolbox
% GNV1B is the GPS-derived orbit of one satellite at 1 Hz. The file gives position and velocity
% EARTH-FIXED (ITRS, |coord_ref = E|), not inertial - one of the defects of the predecessor
% toolboxes. |io.readGNV1B| rotates the state into the GCRS with |referenceSystem| and the Earth
% orientation parameters (EOP) you pass as |'Eop'|; there is no default series, because the
% rotation cannot be guessed. The velocity then carries the frame-rate term (the rotation of the
% Earth, up to about 470 m/s), which nobody forms by hand. With |'Frame', "ITRS"| you get the
% numbers as printed and no EOP are needed. Outputs: |ep| an epoch struct (|.day| MJD, |.sec|
% seconds of day, |.scale| "GPS"), |r| in m and |v| in m/s as 3 x K arrays, one column per epoch.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
[ep, r, v, auxG] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[~, rI, vI] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Frame', "ITRS");
fprintf('GNV1B GRACE-%s: %d epochs, time scale %s, sampling %g s, frame %s (EOP series %s)\n', ...
    auxG.satellite, numel(ep.day), string(ep.scale), auxG.sampling_s, auxG.frame, auxG.eopSeries);
fprintf('mean |r|, which a rotation keeps: GCRS %.3f km, ITRS %.3f km\n', mean(vecnorm(r)) / 1e3, mean(vecnorm(rI)) / 1e3);
fprintf('mean |v|: GCRS %.4f km/s, ITRS %.4f km/s -- the difference is the rotation of the Earth\n', ...
    mean(vecnorm(v)) / 1e3, mean(vecnorm(vI)) / 1e3);
t = double(ep.sec(:)' - ep.sec(1));
figure
plot(t, vecnorm(v) - vecnorm(vI)), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('m/s')
title('|v| in the GCRS minus |v| in the ITRS')

%% Step 2 - The attitude: SCA1B, conjugated exactly once
% SCA1B stores the star-camera quaternion of the rotation ICRF -> SRF (scalar first). The toolbox
% stores attitude as BODY -> GCRS, so |io.readSCA1B| conjugates the quaternion exactly once; no
% later unit conjugates again. A wrong sense passes every algebraic check (unit norm, round
% trips) and fails only against physics: an Earth-pointing satellite has its body +Z axis
% (nadir plate) towards the Earth's centre. |geom.checkAttitudeSense| measures that angle for
% every epoch, and the same angle with the conjugate, so you see that the other sense would have
% failed. It also tells you whether body +X points along or against the velocity: GRACE-C and
% GRACE-D fly facing each other. |geom.AttitudeProvider| interpolates the quaternions (slerp) to
% any epoch inside the nodes and refuses to extrapolate; its static |nominal| builds the ideal
% Earth-pointing frame from the orbit.
[eq, q, auxS] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
fprintf('SCA1B: %d quaternions, frame %s, largest |norm - 1| in the file %.1e\n', size(q, 2), auxS.frame, auxS.normDeviationMax);
[ok, rep] = NonGravForces.geom.checkAttitudeSense(q, r, v);
fprintf('sense check passed: %d; body +Z vs nadir: max %.2f deg (with the conjugate: mean %.1f deg); body +X %s\n', ...
    ok, rep.maxNadirAngle_deg, mean(rep.nadirAngleConjugate_deg), rep.xSense);
att = NonGravForces.geom.AttitudeProvider(eq, q);
qMid = att.at(referenceSystem.daySecondToEpoch([58849 100.5], "GPS"));
fprintf('provider: %d nodes, method %s; q at 00:01:40.5 GPS = [%.6f %.6f %.6f %.6f]\n', att.nodeCount, att.method, qMid);
% the real attitude against the nominal Earth-pointing frame (GRACE-C: +X against the velocity)
qNom = NonGravForces.geom.AttitudeProvider.nominal(r, v, false);
d = pagemtimes(referenceSystem.rotationFromQuaternion(qNom), 'transpose', referenceSystem.rotationFromQuaternion(q), 'none');
angNom = acosd(min(1, (squeeze(d(1, 1, :) + d(2, 2, :) + d(3, 3, :))' - 1) / 2));
fprintf('SCA1B vs the nominal frame: mean %.2f deg, max %.2f deg\n', mean(angNom), max(angNom));
figure
plot(t, rep.nadirAngle_deg, t, angNom), grid on
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('deg')
legend("body +Z vs nadir", "SCA1B vs nominal frame"), title('Attitude sense and pointing of GRACE-C')

%% Step 3 - The accelerometer and the thrusters: ACT1B and THR1B
% ACT1B is the calibrated accelerometer product of GRACE-FO (for GRACE-D since mid-2018 a
% transplant from GRACE-C's data), linear accelerations in the SRF in m/s^2 at 1 Hz. Nothing is
% calibrated by the reader: bias and scale are the product's (Tutorials 7 and 8 estimate them).
% |io.readACC1B| reads ACC1B, ACT1B and GRACE-D's hybrid ACH1B alike and takes the product name
% from the file header, never from the file name. THR1B lists every thruster activation: per
% record the on-time of each of the 14 thrusters (milliseconds in the file, seconds here). The
% first output is an event table (the format |timeSeriesPreparation.qc| uses for impulse priors),
% |aux| holds the on-times per thruster. THR1B must be tagged in GPS time (|time_ref = G|); a file
% that is not is refused (|readTHR1B:dataTimeRef|) rather than shifted by a guessed offset - on
% GRACE-D, 2020-01-22, such records exist.
[ea, aAcc, auxA] = NonGravForces.io.readACC1B(f('samples', 'ACT1B_2020-01-01_C_04.txt'));
fprintf('%s GRACE-%s: %d epochs at %g s; mean lin_accl [%.3e %.3e %.3e] m/s^2\n', auxA.product, auxA.satellite, ...
    numel(ea.day), auxA.sampling_s, mean(aAcc, 2));
[events, auxT] = NonGravForces.io.readTHR1B(f('samples', 'THR1B_2020-01-01_C_04.txt'));
fprintf('THR1B: %d records over the day, %d with a firing, epochs in %s\n', auxT.numRecords, height(events), string(auxT.epoch.scale));
fprintf('%-16s %8s %12s\n', 'thruster', 'firings', 'on-time (s)');
for k = 1:numel(auxT.thruster)
    fprintf('%-16s %8d %12.3f\n', auxT.thruster(k), nnz(auxT.fired(k, :)), sum(auxT.onTime_s(k, :)));
end
% the firings inside the ten minutes of ACT1B
tFire = double(auxT.epoch.sec(:)') + 86400 * (double(auxT.epoch.day(:)') - 58849);
tFire = tFire(any(auxT.fired, 1) & tFire <= t(end));
fprintf('firings within the ACT1B excerpt at t = %s s\n', mat2str(tFire));
figure
ta = double(ea.sec(:)' - ea.sec(1));
plot(ta, 1e9 * (aAcc(2, :) - mean(aAcc(2, :)))), grid on, hold on
for k = 1:numel(tFire)
    xline(tFire(k), 'r--');
end
hold off
xlabel('seconds since 2020-01-01 00:00:00 GPS'), ylabel('nm/s^2')
title('ACT1B cross-track (mean removed) and THR1B firings (red)')

%% Step 4 - Masses: MAS1B is the TANK, not the spacecraft
% MAS1B reports the COLD-GAS mass in the two tanks, about 15 kg each, hourly. It is NOT the mass
% of the spacecraft (about 601 kg, the handbook value in |NonGravForces.const|). A predecessor
% divided forces by the tank mass and got accelerations 20 to 40 times too large. Therefore
% |io.readMAS1B| refuses any total above 100 kg, and |io.massModel| is the only unit that turns
% tank masses into a spacecraft mass: mode "constant" (the handbook value, standard practice) or
% "mas1b" (the handbook value at the first sample plus the change of the gas since then). It
% refuses any mass at or below 500 kg - try it with the tank mass.
[em, tank] = NonGravForces.io.readMAS1B(f('samples', 'MAS1B_2020-01-01_C_04.txt'));
c = NonGravForces.const();
fprintf('MAS1B: %d samples, tank 1 %.3f kg, tank 2 %.3f kg, total %.3f kg\n', numel(em.day), tank.tank_kg(:, 1), tank.total_kg(1));
fprintf('handbook spacecraft mass GRACE-C %.3f kg, GRACE-D %.3f kg: ratio to the tank total %.1f\n', ...
    c.graceFO.massC, c.graceFO.massD, c.graceFO.massC / tank.total_kg(1));
mC = NonGravForces.io.massModel("C");
mV = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', em, 'GasMass_kg', tank.total_kg);
secDay = [(0:3600:82800)'; 86399];                    % every hour of the day, and its last second
eDay = referenceSystem.daySecondToEpoch([58849 * ones(25, 1) secDay], "GPS");
mt = mV.at(eDay);
fprintf('mass model "%s": %.3f kg; "%s": %.4f .. %.4f kg over the day (source: %s)\n', mC.kind, mC.at(eDay(1)), ...
    mV.kind, min(mt), max(mt), mV.source);
try
    NonGravForces.io.massModel("C", 'Anchor_kg', tank.total_kg(1));
catch err
    fprintf('massModel with the tank mass as spacecraft mass is refused: %s\n', err.identifier);
end
figure
plot((0:24), mt - mt(1), 'o-'), grid on
xlabel('hour of 2020-01-01 (GPS)'), ylabel('kg'), title('Spacecraft mass change from MAS1B (mode "mas1b")')

%% Step 5 - The magnetometer: MAG1B
% MAG1B carries two things: the magnetometer's measurement of the Earth's field in the SRF and
% the commanded currents of the three magnetorquers. A current in a coil in the Earth's field is
% a torque on the satellite - not gravitational, which is why the toolbox reads it. The reader
% converts to SI (microtesla to tesla, milliampere to ampere) and sums each axis's A and B
% branches into the net current; it does NOT turn currents into a torque (the coils' turns-area
% product is in no source the project holds) and does NOT apply the calibration columns, whose
% meaning is not established. Tutorial 4 compares the field with the IGRF model.
[eM, fld, auxM] = NonGravForces.io.readMAG1B(f('samples', 'MAG1B_2020-01-01_C_04.txt'));
fprintf('MAG1B: %d records at %g s, |B| %.2f .. %.2f microtesla, net torquer current up to %.1f mA\n', numel(eM.day), ...
    auxM.sampling_s, 1e6 * min(fld.magnitude_T), 1e6 * max(fld.magnitude_T), 1e3 * max(abs(fld.torquerNet_A(:))));
tm = double(eM.sec(:)' - eM.sec(1));
figure
subplot(2, 1, 1), plot(tm, 1e6 * fld.B_srf_T), grid on, ylabel('\muT'), legend("x", "y", "z"), title('MAG1B field in the SRF')
subplot(2, 1, 2), plot(tm, 1e3 * fld.torquerNet_A), grid on, ylabel('mA'), legend("x", "y", "z"), title('Net magnetorquer current')
xlabel('seconds since 2020-01-01 00:00:00 GPS')

%% Step 6 - The gyros: IMU1B against the star cameras
% IMU1B is NOT a rate. Each record is ONE gyro's filtered, accumulated angle about its own axis
% (|FiltAng|, degrees in the file, radians here); the three gyros are sampled in turn at 8 Hz, so
% |io.readIMU1B| returns one record per sample with the gyro id in |aux.gyroId| and does not invent
% common epochs. |geom.attitudeQuality| reduces both sensors to mean rates over the same 1 s
% intervals: the star cameras' rotation between consecutive SCA1B epochs, and each gyro's angle
% difference. It fits rate_g = a . omega + b per gyro and returns the part the fit cannot absorb,
% the RESIDUAL - attitude noise of either sensor. On an Earth-pointing arc omega is nearly
% constant (the satellite turns about body y once per orbit), so the fitted axis and scale are not
% determined; the |condition| number says so, and only the residual should be read. With
% |'Orbit'| it also compares the attitude with the nominal frame. The multiday excerpts hold
% 30 minutes of 2020-01-01 (1800 SCA1B epochs, 43200 gyro samples).
d1 = f('samples', 'multiday');
[eS, qS] = NonGravForces.io.readSCA1B(fullfile(d1, 'SCA1B_2020-01-01_C_04.txt'));
[eI, aI, auxI] = NonGravForces.io.readIMU1B(fullfile(d1, 'IMU1B_2020-01-01_C_04.txt'));
[~, r30, v30] = NonGravForces.io.readGNV1B(fullfile(d1, 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
fprintf('IMU1B: %d records, gyros %s, %s samples each\n', numel(aI), mat2str(auxI.gyros), ...
    mat2str(arrayfun(@(g) nnz(auxI.gyroId == g), auxI.gyros)));
Q = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId, 'Orbit', struct('r', r30, 'v', v30), 'AlongVelocity', false);
w = median(Q.omega_rad_s, 2);
fprintf('SCA1B median body rate [%.2e %.2e %.2e] rad/s; orbital period from the orbit %.0f s\n', w, Q.period_s);
fprintf('%-5s %14s %14s %12s\n', 'gyro', 'rms rate', 'rms residual', 'condition');
for g = 1:numel(Q.gyro)
    fprintf('%-5d %14.3e %14.3e %12.0f\n', Q.gyro(g).id, Q.gyro(g).rmsRate_rad_s, Q.gyro(g).rms_rad_s, Q.gyro(g).condition);
end
fprintf('attitude vs nominal frame (roll pitch yaw): mean %s mrad, rms %s mrad\n', ...
    mat2str(1e3 * Q.nominal.mean_rad(:)', 3), mat2str(1e3 * Q.nominal.rms_rad(:)', 3));
figure
plot(Q.t_s, 1e6 * vertcat(Q.gyro.residual_rad_s)'), grid on
xlabel('seconds since the first SCA1B epoch'), ylabel('\murad/s')
legend("gyro " + string(Q.gyros)), title('Gyro rates minus their fit to the star-camera rates')

%% Step 7 - TU Graz's attitude product: AttitudeProvider.fromItsg
% TU Graz publishes operational attitude, orbits and modelled non-conservative forces in GROOPS
% instrument files (GPS time as MJD, 10 s). TU Graz names GRACE-A/B "GRACE-1/-2" and GRACE-C/D
% "GRACEFO-1/-2". Its attitude quaternion "SRF -> CRF" is this toolbox's ICRF -> SRF in the other
% quaternion convention, so |fromItsg| conjugates it once, like |readSCA1B|. For GRACE-FO the
% numbers are bit-identical to SCA1B's: it is a second carrier of the same attitude, not a second
% sensor (the gyros of Step 6 are). For classic GRACE (here GRACE-A, 2008-01-01) it is TU Graz's
% own product and the attitude source this toolbox reads.
apFO = NonGravForces.geom.AttitudeProvider.fromItsg(f('samples', 'itsg', 'GRACEFO-1_attitude_2020-01-01.txt'));
k = 1:10:600;
eK = referenceSystem.daySecondToEpoch([double(eq.day(k)) double(eq.sec(k))], string(eq.scale));
dd = pagemtimes(referenceSystem.rotationFromQuaternion(apFO.at(eK)), 'transpose', referenceSystem.rotationFromQuaternion(q(:, k)), 'none');
ang = acos(min(1, (squeeze(dd(1, 1, :) + dd(2, 2, :) + dd(3, 3, :)) - 1) / 2));
fprintf('GRACEFO-1 (TU Graz): %d nodes, scale %s; vs SCA1B at %d shared epochs: max %.1e rad\n', ...
    apFO.nodeCount, apFO.scale, numel(k), max(ang));
apA = NonGravForces.geom.AttitudeProvider.fromItsg(f('samples', 'itsg', 'GRACE-1_attitude_2008-01-01.txt'));
fprintf('GRACE-1 = GRACE-A (TU Graz): %d nodes from MJD %.5f to %.5f (%s)\n', apA.nodeCount, apA.span, apA.scale);

%% Step 8 - Space weather: four families, one reader
% A thermosphere model needs the Sun's activity and the geomagnetic activity. |io.spaceWeather|
% recognises the family from the file CONTENT and returns a struct with a |query| function that
% takes epochs (any scale) and returns the drivers at them:
%
% * "celestrak" (SW-All.csv): F10.7 of the previous day, its 81-day means, Ap, Kp and the
%   NRLMSISE-00 ap history;
% * "jb2008" (SOLFSMY + DTCFILE): the JB2008 indices, never mixed with F10.7;
% * "dst" (Kyoto WDC months): the hourly Dst index;
% * "f30" (CLS): the 30 cm radio flux F30 and its 81-day mean, the driver of the density scale
%   of Tutorial 3.
%
% A query outside the table, or of a missing value, raises an error: a NaN driver would
% silently corrupt a density.
W = f('spaceweather');
sw = NonGravForces.io.spaceWeather(fullfile(W, 'SW-All_2019-09_2020-03.csv'));
sj = NonGravForces.io.spaceWeather(string(fullfile(W, {'SOLFSMY_2019-09_2020-03.TXT', 'DTCFILE_2019-09_2020-03.TXT'})));
sd = NonGravForces.io.spaceWeather(string(fullfile(W, ["dst1912" "dst2001"] + ".for.request")));
cls = NonGravForces.io.spaceWeather(fullfile(W, 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
for s = {sw, sj, sd, cls}
    fprintf('family %-9s %4d days, MJD %d..%d\n', s{1}.format, numel(s{1}.mjd), s{1}.span);
end
fprintf('CLS file adjusted to 1 AU: %d\n', cls.adjusted);
eH = referenceSystem.mjdToEpoch(58849 + (0:71)' / 24 + 1 / 48, "UTC");      % hourly, 2020-01-01..03
dC = sw.query(eH); dD = sd.query(eH); dF = cls.query(eH); dJ = sj.query(eH(1));
fprintf('2020-01-01: F10.7 (previous day) %.1f sfu, 81-day mean %.1f sfu, Ap %g; F30 %.1f sfu (81-day mean %.1f); S10 %.1f sfu\n', ...
    dC.f107(1), dC.f107a(1), dC.ap(1), dF.f30(1), dF.f30Mean81(1), dJ.s10(1));
th = (0:71) + 0.5;
figure
subplot(2, 1, 1), stairs(th, dC.kp), grid on, ylabel('Kp'), title('Geomagnetic activity, 2020-01-01..03 (UTC)')
subplot(2, 1, 2), stairs(th, dD.dst), grid on, ylabel('Dst (nT)'), xlabel('hours since 2020-01-01 00:00 UTC')

%% Step 9 - Thermosphere densities: Swarm and TOLEOS
% ESA derives thermosphere densities from the orbits of Swarm A, B and C (POD-only, 30 s) and,
% in the TOLEOS project with TU Delft, from the accelerometers of CHAMP, GRACE and GRACE-FO
% (10 s). Both come as HAPI CSV files with the same columns, and |io.readSwarmDensity| reads both:
% epochs in UTC, geodetic latitude and longitude in rad and height in m, the density in kg/m^3,
% and in |aux| a validity flag and the satellite. They are OBSERVATIONS of the density - Tutorial 3
% turns Swarm's into a daily correction of the model density and checks it against GRACE-FO's.
sat = ["A" "B" "C"];
figure, hold on
for k = 1:3
    [eD, lla, rho, auxD] = NonGravForces.io.readSwarmDensity(f('samples', 'swarm', "SW_OPER_DNS" + sat(k) + "POD_2_2020-01-01_2020-01-03.csv"));
    fprintf('Swarm %s: %d records (%d valid), mean height %.1f km, mean density %.3e kg/m^3\n', auxD.satellite, ...
        auxD.numRecords, nnz(auxD.valid), mean(lla(3, :)) / 1e3, mean(rho(auxD.valid)));
    td = (double(eD.day(:)') - 58849) * 24 + double(eD.sec(:)') / 3600;
    plot(td(auxD.valid), rho(auxD.valid), '.', 'MarkerSize', 3)
end
for name = ["CH_OPER_DNS_ACC_2_2008-01-01_00-06" "GR_OPER_DNS1ACC_2_2008-01-01_00-06" "GF_OPER_DNS1ACC_2_2020-01-01_00-06"]
    [eD, lla, rho, auxD] = NonGravForces.io.readSwarmDensity(f('samples', 'toleos', name + ".csv"));
    fprintf('TOLEOS %-8s: %d records at %g s, mean height %.1f km, mean density %.3e kg/m^3\n', auxD.satellite, ...
        auxD.numRecords, double(eD.sec(2) - eD.sec(1)), mean(lla(3, :)) / 1e3, mean(rho(auxD.valid)));
    if auxD.satellite == "GRACE-C"
        td = double(eD.sec(:)') / 3600;
        plot(td(auxD.valid), rho(auxD.valid), 'k.', 'MarkerSize', 3)
    end
end
hold off, grid on, set(gca, 'YScale', 'log')
xlabel('hours since 2020-01-01 00:00 UTC'), ylabel('kg/m^3')
legend(["Swarm A" "Swarm B" "Swarm C" "GRACE-FO 1 (TOLEOS)"]), title('Observed thermosphere densities')

%% Step 10 - Provenance, the data pool and the fetch catalogue
% Every reader returns a |provenance| record in its |aux|: the producing unit, the toolbox and
% dependency versions, a hash of the options, the files read and the time. A result without one
% cannot be audited. |NonGravForces.provenance| builds such a record for your own results too.
%
% The committed excerpts are short on purpose. Full days and months live in the DATA POOL, one
% folder shared by every astrobayer toolbox and laid out by product (EOP, SpaceWeather,
% GRACEFO/L1B/RL04/<yyyy>/<yyyy-mm-dd>, Swarm, TOLEOS, ...). |NonGravForces.dataPool| resolves its
% root - an explicit argument, else the environment variable DATAPOOL, else READING_DATA_POOL -
% and refuses when none is set; with |'Required', false| it returns "" instead. Nothing is ever
% written into the toolbox folder as a fallback. |NonGravForces.fetchCatalog| lists everything
% |setupNonGravForces| CAN fetch, with spans, sizes and the licence class this repository has
% recorded, and expands a selection into a plan. It reaches no network and fetches nothing.
p = auxG.provenance;
fprintf('provenance of readGNV1B: unit %s, toolbox %s, option hash %s, %d file(s), created %s\n', p.unit, ...
    p.toolboxVersion, p.optionHash, numel(p.files), string(p.created));
mine = NonGravForces.provenance("tutorial02", 'Files', string(p.files), 'Note', "orbit of Step 1");
fprintf('my own record: unit %s, note "%s", dependencies %s\n', mine.unit, mine.note, strjoin(string({mine.dependencies.name}), ", "));
pool = NonGravForces.dataPool("", 'Required', false);
if pool == ""
    fprintf('no data pool is named on this machine (DATAPOOL and READING_DATA_POOL are unset)\n');
else
    fprintf('the data pool on this machine: %s\n', pool);
end
try
    NonGravForces.dataPool();
    fprintf('dataPool() with Required true returns it too\n');
catch err
    fprintf('dataPool() refuses without a pool: %s\n', err.identifier);
end
T = NonGravForces.fetchCatalog();
disp(T(:, ["id" "group" "destination" "licence" "first" "last" "disk_MB"]))
[Tp, plan] = NonGravForces.fetchCatalog('PoolDays', "2020-01-01..2020-01-03", 'Items', ["itsg" "dst" "f30" "swarm"]);
disp(Tp(Tp.selected, ["id" "count" "totalDisk_MB" "licence"]))
fprintf('the plan expands to %d days: %s\n', numel(plan.PoolDays), strjoin(plan.PoolDays, ", "));
fprintf('constants: S0 = %g W/m^2, au = %.0f m, species %s\n', c.S0, c.au, strjoin(c.speciesName, " "));
