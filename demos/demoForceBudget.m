%DEMO01_FORCEBUDGET Which effect carries the arc, and what does a modelling choice cost?
%
%   A budget is not an acceleration: it answers "which effect dominates, on
%   which axis, and what changes if I turn something off". This runs the same
%   arc through two presets and puts them side by side, which is the comparison
%   the unit exists for.
%
%   Run:  demo01_forceBudget
%
% See also NonGravForces.app.forceBudget, NonGravForces.Model.
%
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
root = fileparts(fileparts(mfilename('fullpath')));
f = @(varargin) fullfile(root, varargin{:});

eop = referenceSystem.readEop(f('data', 'eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
[ep, r, v] = NonGravForces.io.readGNV1B(f('data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[epAtt, q] = NonGravForces.io.readSCA1B(f('data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(epAtt, q);
macro = NonGravForces.geom.graceFO("C");

% Two fidelity tiers on the same arc. The fast preset uses an exponential density,
% a constant drag coefficient, a conical shadow, and neither Earth radiation nor
% thermal recoil nor self-shadowing; the standard one brings all of them. What the
% table shows is what that costs and what it buys.
expo = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
fastModel = NonGravForces.Model(macro, NonGravForces.Options("fast"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'Exponential', expo);
standardModel = NonGravForces.Model(macro, NonGravForces.Options("standard"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw);

% Ten epochs, because a budget is a statement about SHARES and ten of them make the point. (The
% reason this comment gave until 2026-09-18 -- "the standard preset costs 88 ms each" -- was stale
% by a wide margin: WP-C3 fixed that on 2026-09-16 and HANDOVER O-29 was not closed until the figure
% had propagated into four more files. It costs 0.798 ms per epoch at K=100.) A
% demo should finish while you are still looking at it.
n = 10;
epArc = referenceSystem.daySecondToEpoch([ep.day(1:n) ep.sec(1:n)], ep.scale);
B = NonGravForces.app.forceBudget({fastModel, standardModel}, epArc, r(:, 1:n), v(:, 1:n), []);

fprintf('\nPer-effect budget over %d epochs, body frame:\n\n', n);
disp(B(:, {'model', 'effect', 'rms_norm', 'share'}));
fprintf(['The share column is formed from rms_norm, so it can sum to more than one when two\n' ...
         'effects partly cancel in the total. That is a property of the arc, not an error,\n' ...
         'which is why the total row is reported beside them.\n']);
