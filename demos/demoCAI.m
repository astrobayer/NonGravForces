%DEMOCAI A single-axis cold-atom accelerometer on a GRACE-like orbit, and why it needs the attitude.
%
%   The headline number of this demo is not the drag phase. It is the Coriolis
%   phase: at an orbital attitude rate, one millimetre per second of atom
%   velocity RELATIVE TO THE INSTRUMENT is 8.9 radians, an apparent 2.2e-6 m/s^2 --
%   seven times the drag the instrument is there to measure. (The orbital velocity
%   does not enter: atoms and satellite share it. Until 2026-09-15 this demo fed
%   it in and printed 16.7 m/s^2, HANDOVER LP-C finding F6.)
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

% Now switch the rotation on, at a GRACE-like attitude rate, with the atom cloud drifting at
% 1 mm/s across the beam relative to the instrument (a release residual, stated here as an assumption).
sr = NonGravForces.sensor.AtomInterferometer(kEff, T, cycle, [1; 0; 0], 'SampleRate', fs, 'Rotation', true);
omega = repmat([0; 0; 1.1e-3], 1, K);
atomVel = repmat([0; 1e-3; 0], 1, K);
yr = sr.observe(struct('a', aDrag, 'omega', omega, 'atomVelocity', atomVel));
coriolis = mean(yr.value(yr.mask)) - mean(y.value(y.mask));

fprintf('\nCoriolis phase       %.4e rad = %.3e fringes\n', coriolis, coriolis / (2 * pi));
fprintf('as an acceleration   %.2e m/s^2, against a drag of %.1e\n', ...
    abs(coriolis) / (kEff * T^2), 3e-7);
fprintf(['\nOne millimetre per second of atom velocity outweighs the signal several times over. The\n' ...
         'rotation input is not a refinement to this instrument; it has to be measured or compensated.\n']);
