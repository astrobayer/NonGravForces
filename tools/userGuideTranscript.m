function userGuideTranscript()
%USERGUIDETRANSCRIPT Run the user guide's eight-step chain and print what it prints, for the guide's transcript.
%
%   userGuideTranscript
%
%   The user guide's section "What this chain actually prints" quotes the output of its own Steps 1-8
%   on the committed GRACE-C excerpt of 2020-01-01. Until 2026-10-31 that transcript dated from 0.12.0
%   and nothing regenerated it. This tool runs the same calls as the guide's listings -- the standard
%   preset as shipped, with the combined density scale from the committed Swarm and CLS excerpts -- and
%   prints each step in the transcript's format. release-checks.yml runs it ("what: transcript"); the
%   printed block is pasted into doc/src/NonGravForces_UserGuide.tex by hand, with the run it came from.
%
%   Lives in tools/: a documentation instrument, not part of the scientific surface, not gated.
%
% Created 2026-10-31 by Claude (Anthropic) for Matthias Weigelt.
root = fileparts(fileparts(mfilename('fullpath')));
f = @(varargin) fullfile(root, 'data', varargin{:});
v = NonGravForces.toolboxVersion();
fprintf('NonGravForces %s, %s\n\n', v.semantic, string(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm ''UTC''')));

head("Step 1: the auxiliary data");
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw  = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarmFiles = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
fprintf('space weather family %s, CLS family %s, Swarm files %d\n', string(sw.format), string(cls.format), numel(swarmFiles));

head("Step 2: orbit and attitude");
[ep, r, v] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q]    = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att        = NonGravForces.geom.AttitudeProvider(eq, q);
fprintf('GNV1B epochs %d (%s), |r| = %.1f km, |v| = %.4f km/s\n', numel(ep.day), string(ep.scale), mean(vecnorm(r)) / 1e3, mean(vecnorm(v)) / 1e3);
fprintf('SCA1B quaternions %d, provider method %s\n', size(q, 2), att.method);

head("Step 3: the spacecraft");
macro = NonGravForces.geom.graceFO("C");
fprintf('panels %d, area sum %.4f m^2, mass %.3f kg\n', macro.n, sum(macro.A), macro.mass_kg);
fprintf('silhouette defect |sum A_k n_k| = %.3e m^2\n', norm(macro.silhouetteDefect()));
fprintf('   (zero for any CLOSED body; the printed handbook areas do not close on each other)\n');

head("Step 4: the options and the model");
opt = NonGravForces.Options("standard");
S   = NonGravForces.env.densityScale(swarmFiles, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
M   = NonGravForces.Model(macro, opt, 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
fprintf('preset standard: density %s x scale %s, gsi %s, accommodation %s (alpha %.2f), erp %s, thermal %s\n', ...
    opt.drag.density, opt.drag.densityScale, opt.drag.gsi, opt.drag.accommodation, opt.drag.alpha, opt.erp.backend, opt.thermal.model);

head("Step 5: evaluate over an arc");
idx = 1:200;
e   = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
[a, info] = M.accel(e, r(:, idx), v(:, idx));
rmsOf = @(x) sqrt(mean(vecnorm(x) .^ 2));
fprintf('epochs %d, |a| RMS %.4e m/s^2\n', numel(idx), rmsOf(a));
fprintf('   drag      %.4e   (%.1f %% of the total RMS)\n', rmsOf(info.drag), 100 * rmsOf(info.drag) / rmsOf(a));
fprintf('   radiation %.4e\n', rmsOf(info.radiation));
fprintf('   thermal   %.4e\n', rmsOf(info.thermal));
fprintf('info fields: %s\n', strjoin(string(fieldnames(info)), ", "));
fprintf('density %.4e kg/m^3 (scale %.3f), |v_rel| %.1f m/s\n', mean(info.rho), mean(info.densityScale), mean(vecnorm(info.vrel_body)));

head("Step 6: the per-effect budget");
B = M.budget(e, r(:, idx), v(:, idx));
disp(B(:, {'effect', 'rms_norm'}));

head("Step 7: compare with the measured accelerations");
[epa, aMeas] = NonGravForces.io.readACC1B(f('samples', 'ACT1B_2020-01-01_C_04.txt'));
[tf, loc]    = ismember(ep.sec(idx), epa.sec);
am           = aMeas(:, loc(tf));
amod         = a(:, tf);
names = ["cross-track" "radial     "];
for j = 2:3
    x = amod(j, :) - mean(amod(j, :)); y = am(j, :) - mean(am(j, :));
    fprintf('%s correlation with ACT1B  %.4f\n', names(j - 1), sum(x .* y) / sqrt(sum(x .^ 2) * sum(y .^ 2)));
end

head("Step 8: calibrate");
colNames = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
units = struct('lin_accl_x', 'm/s2', 'lin_accl_y', 'm/s2', 'lin_accl_z', 'm/s2');
auxACT1B = struct("obsNames", {colNames}, "variableNames", {colNames}, "variableFieldNames", {colNames}, ...
    "variableUnits", units, "fileName", "ACT1B_2020-01-01_C_04.txt", "product", "ACT1B");
acc = timeSeriesPreparation.Arc(referenceSystem.epochTo(epa, "gpsNanoseconds"), aMeas.', auxACT1B);
try
    res = NonGravForces.app.calibrate(acc, M, 'Orbit', struct('epoch', ep, 'r', r, 'v', v), ...
        'Attitude', att, 'Fit', {'FixCrossTrack', false});
    disp(res.calibration.summaryTable());
catch err
    fprintf('calibrate refused: %s\n  %s\n', err.identifier, err.message);
end
end

function head(title)
fprintf('\n--- %s %s\n', title, repmat('-', 1, max(3, 72 - strlength(title))));
end
