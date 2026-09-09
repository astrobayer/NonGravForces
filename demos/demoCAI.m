%DEMOCAI A single-axis cold-atom accelerometer on a GRACE-like orbit, and why it needs the attitude.
%
%   The headline number of this demo is not the drag phase. It is the Coriolis
%   phase: at an orbital rate across an orbital speed the rotation term is
%   6.7e7 radians, an apparent 16.7 m/s^2. An interferometer that does not
%   counter-rotate its Raman mirror does not read a small bias on top of the
%   signal, it has no fringe left to read.
%
%   Run:  demoCAI
%
% See also NonGravForces.sensor.AtomInterferometer, NonGravForces.sensor.Electrostatic.
%
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
kEff = 1.610e7;          % rad/m, two-photon Raman on the rubidium D2 line
T = 0.5;                 % s, interrogation time
cycle = 2;               % s, so 1 s of dead time per cycle
fs = 10;                 % Hz
K = 600;
t = (0:K-1) / fs;

s = NonGravForces.sensor.AtomInterferometer(kEff, T, cycle, [1; 0; 0], 'SampleRate', fs);

% A drag of 3e-7 m/s^2 with a modulation at the orbital period's harmonics.
aDrag = [3e-7 + 5e-8 * sin(2 * pi * 0.02 * t); zeros(2, K)];
y = s.observe(struct('a', aDrag));

fprintf('\nreadings per cycle   %d valid of %d samples (dead time %.1f s)\n', ...
    nnz(y.mask), K, y.info.deadTime_s);
fprintf('drag phase           %.4f rad  = %.4e m/s^2 equivalent\n', ...
    mean(y.value(y.mask)), mean(y.value(y.mask)) / (kEff * T^2));
fprintf('wrapped phase        %.4f rad  (what the instrument returns)\n', y.info.wrapped(find(y.mask, 1)));
fprintf('fringe integer       %d  (what the readout cannot tell you)\n', y.info.fringe(find(y.mask, 1)));

% Now switch the rotation on, at a GRACE-like attitude rate and orbital speed.
sr = NonGravForces.sensor.AtomInterferometer(kEff, T, cycle, [1; 0; 0], 'SampleRate', fs, 'Rotation', true);
omega = repmat([0; 0; 1.1e-3], 1, K);
vel = repmat([0; 7600; 0], 1, K);
yr = sr.observe(struct('a', aDrag, 'omega', omega, 'velocity', vel));
coriolis = mean(yr.value(yr.mask)) - mean(y.value(y.mask));

fprintf('\nCoriolis phase       %.4e rad = %.3e fringes\n', coriolis, coriolis / (2 * pi));
fprintf('as an acceleration   %.2f m/s^2, against a drag of %.1e\n', ...
    abs(coriolis) / (kEff * T^2), 3e-7);
fprintf(['\nThat is larger than gravity and eight orders above the signal. The rotation input is\n' ...
         'not a refinement to this instrument; without compensation there is nothing to measure.\n']);
