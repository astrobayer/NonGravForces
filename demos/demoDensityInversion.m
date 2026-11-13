%DEMODENSITYINVERSION Recover the thermosphere density that produced an acceleration.
%
%   The gas force is exactly linear in the species densities, so the inversion
%   is a division and not a search. This demo proves that on the only data where
%   the truth is known: the model's own output, fed back in.
%
%   Run:  demoDensityInversion
%
% See also NonGravForces.app.densityInversion.
%
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
root = fileparts(fileparts(mfilename('fullpath')));
f = @(varargin) fullfile(root, varargin{:});

eop = referenceSystem.readEop(f('data', 'eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
[ep, r, v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[epAtt, q] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(epAtt, q);
M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("fast"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, ...
    'Exponential', struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3));

n = 12;
epArc = referenceSystem.daySecondToEpoch([ep.day(1:n) ep.sec(1:n)], ep.scale);
qArc = att.at(epArc);
[~, ai] = M.accel(epArc, r(:, 1:n), v(:, 1:n), qArc);

% A planted factor of 1.7 on the gas part alone: the answer is known exactly.
planted = 1.7;
measured = ai.radiation + ai.thermal + planted * ai.drag;
T = NonGravForces.app.densityInversion(M, epArc, r(:, 1:n), v(:, 1:n), qArc, measured);

fprintf('\nplanted scale        %.4f\n', planted);
fprintf('recovered, worst     %.6f  (relative error %.2e)\n', ...
    max(T.scale), max(abs(T.scale - planted)) / planted);
fprintf('model density        %.4e kg/m^3\n', T.rho_model_kg_m3(1));
fprintf('inverted density     %.4e kg/m^3\n', T.rho_kg_m3(1));
fprintf('residual across drag %.3e m/s^2  (zero: the scaling was along the drag direction)\n', ...
    max(T.residualPerp_m_s2));
fprintf(['\nAgainst real data the perpendicular residual is the honest quality flag: it is what\n' ...
         'the density cannot explain, and it grows when the attitude, the radiation model or\n' ...
         'the calibration is wrong.\n']);
