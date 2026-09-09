%DEMOGRACEFOCALIBRATION Recover an accelerometer's bias and scale from a modelled arc.
%
%   THIS DEMO PLANTS ITS MEASUREMENT, and the reason changed on 2026-09-09. It
%   used to be that no ACT1B file existed anywhere. That is no longer true: the
%   GFZ ISDC archive setupNonGravForces fetches carries ACT1B for the same day,
%   under the same naming as every other product, and NonGravForces.io.readACC1B
%   already reads it.
%
%   It still plants, for a reason that is about DEMONSTRATION rather than data: a
%   planted bias and scale are KNOWN, so this page can show the recovered numbers
%   beside the truth and you can see the estimator work. Against real ACT1B the
%   truth is what you are trying to find, and a demo that prints an answer nobody
%   can check demonstrates nothing.
%
%   To run it on the real series: add ACT1B to the members setupNonGravForces
%   extracts, read it with NonGravForces.io.readACC1B, and put it in `measured`.
%   Nothing else changes.
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
