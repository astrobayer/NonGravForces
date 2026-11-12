%% Tutorial 1 - First steps: the toolbox, its data and one modelled acceleration
% This is the first of the first-time tutorials of NonGravForces. It checks the installation,
% looks at the data the toolbox ships, builds the model of GRACE-C and evaluates the
% non-gravitational acceleration -- drag, radiation pressure and thermal recoil -- along ten
% minutes of the real orbit of 2020-01-01, and puts it beside what the accelerometer measured.
% Nothing is calibrated yet; the later tutorials take every step apart, one layer each.
%
% *What you will learn*
%
% * how to check that NonGravForces and its five sibling toolboxes are installed;
% * where the committed data live and what the reference day is;
% * what a preset is, and what the three presets fast, standard and precise choose;
% * how to build a |NonGravForces.Model| and call |accel| on an arc;
% * the units and frames of every input and output (SI, the body frame, GPS epochs).
%
% *Requirements:* none - every file this tutorial reads is committed in |data/|. Run time: about
% ten seconds.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Check the installation
% |setupNonGravForces| (run once per MATLAB session, after the setups of the five sibling
% toolboxes; see the user guide's chapter Installation) puts the toolbox, its tests, tools,
% demos and these tutorials on the path. |NonGravForces.toolboxVersion| tells you which
% release you run - quote it when you report a result - and |NonGravForces.dependencyReport|
% which siblings it found: referenceSystem (time scales, frames, the Sun) and gracefoReading
% (the Level-1B reader) are required; gjToolbox (orbit integration), timeSeriesPreparation
% (the accelerometer filter, quality control) and lsqvce (the least-squares estimator of the
% calibration) are needed by the tutorials that use them.
v = NonGravForces.toolboxVersion();
fprintf('NonGravForces %s of %s\n', v.semantic, v.date)
deps = NonGravForces.dependencyReport('Verbose', false);
disp(struct2table(rmfield(deps, {'setup', 'role'})))

%% Step 2 - The committed data: one real day of GRACE-C
% The toolbox ships EXCERPTS of real data so that every test and tutorial runs without a
% download: the first ten minutes (600 records at 1 Hz) of GRACE-C's Level-1B products of
% 2020-01-01 - GNV1B (orbit), SCA1B (attitude), ACT1B (the calibrated accelerometer) - the whole
% day of THR1B (thruster firings) and MAS1B (tank masses), and 30 minutes of 2020-01-01..03 in
% |data/samples/multiday|. Beside them: the Earth orientation parameters (IERS finals2000A) and
% the space weather (CelesTrak SW-All, CLS F30) of 2019-09..2020-03, and three days of Swarm
% POD densities. Full months live in the DATA POOL (the folder DATAPOOL names), never in the
% toolbox. |f| below is a helper that builds a path into |data/|.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
d = dir(f('samples', '*.txt'));
fprintf('%-28s %8s\n', 'file', 'bytes');
for k = 1:numel(d)
    fprintf('%-28s %8d\n', d(k).name, d(k).bytes);
end

%% Step 3 - Presets: the fidelity ladder
% A |NonGravForces.Options| object fixes every modelling choice. It starts from one of three
% PRESETS and is changed one dotted field at a time (|opt.set("srp.tsi", 1367)|); an unknown
% field or value is refused, so a typo can never select a default silently. "fast" is an
% exponential atmosphere and a constant drag coefficient for quick looks; "standard" is the
% reference for GRACE-FO work (NRLMSISE-00 times an OBSERVED daily density scale from Swarm,
% Sentman's gas-surface model, Earth radiation, transient wall temperatures); "precise" takes
% the slower models (NRLMSIS 2.1, CERES, JPL's DE440 Sun) and needs data the pool provides.
presets = ["fast" "standard" "precise"];
fields = ["drag.density" "drag.densityScale" "drag.gsi" "srp.eclipse" "erp.backend" "thermal.model"];
T = table('Size', [numel(fields) 3], 'VariableTypes', repmat("string", 1, 3), ...
    'VariableNames', presets, 'RowNames', fields);
for p = presets
    s = NonGravForces.Options(p).toStruct();
    for k = 1:numel(fields)
        part = split(fields(k), ".");
        T{fields(k), p} = string(s.(part(1)).(part(2)));
    end
end
disp(T)

%% Step 4 - The spacecraft: a panel model of GRACE-C
% Every force is a sum over flat panels, each with an area, an outward normal, optical
% properties in the visible and the infrared, a heat capacity and a polygon (for shadowing).
% |NonGravForces.geom.graceFO("C")| builds GRACE-C from the GRACE-FO handbook's table. The
% printed handbook areas do not quite enclose a volume - the "silhouette defect", the
% area-weighted sum of the normals, is not zero along z - which Tutorial 5 explains. The mass is
% the handbook's spacecraft mass, NOT MAS1B's tank mass (Tutorial 2).
macro = NonGravForces.geom.graceFO("C");
fprintf('%s: %d panels, total area %.3f m^2, mass %.3f kg\n', macro.name, macro.n, sum(macro.A), macro.mass_kg);
fprintf('silhouette defect sum(A_k n_k) = [%.3f %.3f %.3f] m^2\n', macro.silhouetteDefect());
figure
barh(macro.A), set(gca, 'YTick', 1:macro.n, 'YTickLabel', cellstr([macro.panels.name]), 'TickLabelInterpreter', 'none')
xlabel('area (m^2)'), title('The panels of GRACE-C')

%% Step 5 - Orbit and attitude of the reference day
% |io.readGNV1B| reads the GPS orbit, which GNV1B gives Earth-fixed (ITRS), and rotates it into
% the inertial GCRS with the Earth orientation parameters - positions |r| in m and velocities
% |v| in m/s, 3 x K. Its epochs |ep| are referenceSystem epoch structs (day and second in a
% named time scale, here GPS). |io.readSCA1B| reads the star-camera quaternions and returns them
% as BODY -> GCRS rotations (SCA1B stores the inverse; the reader conjugates exactly once), and
% |geom.AttitudeProvider| interpolates them to any epoch.
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
[ep, r, vel] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
fprintf('GNV1B: %d epochs (%s), mean |r| %.1f km, mean |v| %.4f km/s\n', numel(ep.day), string(ep.scale), ...
    mean(vecnorm(r)) / 1e3, mean(vecnorm(vel)) / 1e3)
fprintf('SCA1B: %d quaternions, interpolated by %s\n', size(q, 2), att.method)

%% Step 6 - The model: configure once, evaluate many
% |NonGravForces.Model| binds the spacecraft, the options and the environment's inputs (attitude,
% Earth orientation, space weather). The standard preset multiplies the NRLMSISE-00 density by
% an observed daily scale, so it needs a |DensityScale| object: Swarm's density ratio where Swarm
% flew, a scale predicted from the F30 solar flux elsewhere (Tutorial 3). A standard Model
% without one is refused at evaluation rather than run unscaled.
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
M = NonGravForces.Model(macro, NonGravForces.Options("standard"), 'Attitude', att, 'Eop', eop, ...
    'SpaceWeather', sw, 'DensityScale', S);
disp(M)

%% Step 7 - Evaluate the acceleration over ten minutes
% |M.accel(e, r, v)| returns the non-gravitational acceleration in the BODY frame (the
% spacecraft's science reference frame: x along track towards the other satellite, z towards
% nadir), in m/s^2, 3 x K, one column per epoch. The second output |info| carries every part:
% |info.drag|, |info.radiation| (Sun and Earth), |info.thermal| and the environment the model
% saw (density, its scale, the relative velocity). One call over the whole arc is much cheaper
% per epoch than one call per epoch (the user guide's section Cost).
idx = 1:600;
e = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
[a, info] = M.accel(e, r(:, idx), vel(:, idx));
rmsOf = @(x) sqrt(mean(vecnorm(x) .^ 2));
fprintf('|a| RMS %.3e m/s^2: drag %.3e, radiation %.3e, thermal %.3e\n', rmsOf(a), ...
    rmsOf(info.drag), rmsOf(info.radiation), rmsOf(info.thermal));
fprintf('density %.3e kg/m^3 (observed scale %.3f), |v_rel| %.1f m/s\n', mean(info.rho), ...
    mean(info.densityScale), mean(vecnorm(info.vrel_body)));
t = (ep.sec(idx) - ep.sec(1))';
figure
ax = ["x (along track)" "y (cross track)" "z (radial)"];
for j = 1:3
    subplot(3, 1, j)
    plot(t, 1e9 * info.drag(j, :), t, 1e9 * info.radiation(j, :), t, 1e9 * info.thermal(j, :), t, 1e9 * a(j, :), 'k')
    ylabel('nm/s^2'), title(ax(j)), grid on
end
xlabel('seconds since 2020-01-01 00:00:00 GPS')
legend(["drag" "radiation" "thermal" "total"], 'Location', 'eastoutside')

%% Step 8 - Beside the accelerometer
% ACT1B is the measured acceleration in the same body frame, after JPL's calibration. It carries
% an unknown bias and scale and is filtered (Tutorials 7 and 8 do this properly). Removing the
% means and comparing the shapes already shows whether the model is right: the cross-track and
% radial correlations should be close to one.
[epa, am] = NonGravForces.io.readACC1B(f('samples', 'ACT1B_2020-01-01_C_04.txt'));
[tf, loc] = ismember(ep.sec(idx), epa.sec);
x = a(:, tf) - mean(a(:, tf), 2);
y = am(:, loc(tf)) - mean(am(:, loc(tf)), 2);
c = sum(x .* y, 2) ./ sqrt(sum(x .^ 2, 2) .* sum(y .^ 2, 2));
fprintf('correlation model vs ACT1B (means removed): x %.3f, y %.3f, z %.3f\n', c)
figure
for j = 2:3
    subplot(2, 1, j - 1)
    plot(t(tf), 1e9 * x(j, :), t(tf), 1e9 * y(j, :)), grid on
    ylabel('nm/s^2'), title(ax(j) + ", means removed"), legend("model", "ACT1B")
end
xlabel('seconds since 2020-01-01 00:00:00 GPS')
