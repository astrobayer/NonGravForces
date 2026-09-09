%DEMOGRACEFOCALIBRATION Recover an accelerometer's bias and scale from a modelled arc.
%
%   NO ACT1B FILE IS IN THIS REPOSITORY, and that is not an oversight: none is
%   on the owner's host either (HANDOVER obligation O-8), so there is no measured
%   series to calibrate against. This demo therefore does the only honest thing
%   available: it PLANTS a known bias and scale on the model's own output and
%   recovers them, which is how the estimator is tested and is what a real
%   calibration looks like with the truth replaced by a measurement.
%
%   Swap `measured` for a real ACT1B series and nothing else changes.
%
%   Run:  demoGraceFOCalibration
%
% See also NonGravForces.app.calibrate, NonGravForces.sensor.Electrostatic.
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

n = 120;
epArc = referenceSystem.daySecondToEpoch([ep.day(1:n) ep.sec(1:n)], ep.scale);
aModel = M.accel(epArc, r(:, 1:n), v(:, 1:n), att.at(epArc));

trueBias = [1.0e-6; -3.0e-6; 2.0e-7];
trueScale = diag([1.05 0.96 1.01]);
rng(20260909);                                  % a demo that changes with the seed is not a demo
measured = trueBias + trueScale * aModel + 2e-9 * randn(3, n);

cal = NonGravForces.app.calibrate(measured, aModel, epArc, 'Bias', "constant", 'Scale', "diag");

fprintf('\n            planted      recovered     formal sigma\n');
for i = 1:3
    fprintf('bias %c   %11.4e  %11.4e  %11.2e\n', 'w' + i, trueBias(i), cal.bias(i, 1), cal.sigma(i, 1));
end
for i = 1:3
    fprintf('scale %c  %11.4f  %11.4f  %11.2e\n', 'w' + i, trueScale(i, i), cal.scale(i, i), cal.sigma(i, end));
end
fprintf('\nsigma0 (a-posteriori) %.3e m/s^2 against a planted 2.0e-09\n', cal.sigma0(1));
fprintf(['\nThose formal sigmas are OPTIMISTIC. The weights are diagonal and accelerometer noise\n' ...
         'is not white; it rises as 1/f below a few mHz. They compare axes and arcs honestly,\n' ...
         'they are not statements of absolute accuracy. Whiten with qc.residuals first if that\n' ...
         'matters.\n']);
