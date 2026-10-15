%DEMOGRACEFOCALIBRATION Recover an accelerometer's bias and scale from a known truth, through app.calibrate.
%
%   THIS DEMO PLANTS ITS MEASUREMENT, for a reason that is about DEMONSTRATION
%   rather than data: a planted bias and scale are KNOWN, so this page can
%   show the recovered numbers beside the truth and you can see the
%   estimator work. Against real ACT1B the truth is what you are trying to
%   find, and a demo that prints an answer nobody can check demonstrates
%   nothing. (The real day runs in tests/tNgfCalibrate.m, calibrateRealDay.)
%
%   Since N-02 (0.13.0) the truth comes from app.syntheticPair, which passes a
%   planted scale, a drifting bias, coloured noise, spikes and a gap through
%   the REAL Level-1B CRN chain, and the estimate from app.calibrate: the
%   instrument stages of timeSeriesPreparation.run, then lsqvce.fit through
%   app.fitCalibration, with block-bootstrap sigmas.
%
%   Run:  demoGraceFOCalibration
%
% See also NonGravForces.app.calibrate, NonGravForces.app.syntheticPair, NonGravForces.app.Calibration.
%
% Created 2026-09-09 by Claude Opus 5 (Anthropic); rewritten 2026-10-15 by Claude (Anthropic) for N-02c.
[obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 3);
res = NonGravForces.app.calibrate(obs, ref, 'Fit', {'FixCrossTrack', false, 'NBoot', 100});
cal = res.calibration;

fprintf('\n            planted      recovered    bootstrap sigma   formal sigma\n');
for i = 1:3
    fprintf('bias %c   %11.4e  %11.4e\n', 'w' + i, truth.bias0(i), cal.bias.coef(i, 1));
end
for i = 1:3
    fprintf('scale %c  %11.5f  %11.5f      %9.2e       %9.2e\n', 'w' + i, truth.scale(i, i), cal.scale(i, i), ...
        cal.uncertainty.sigmaScaleBootstrap(i, i), cal.uncertainty.sigmaScaleFormal(i, i));
end
fprintf('\nEIV attenuation bound per model axis (1 - lambda): %s\n', mat2str(cal.uncertainty.eiv.scaleBiasBound, 3));
fprintf('stages run: %s\n', strjoin(res.diagnostics.stages, ', '));
fprintf(['\nThe bootstrap sigma is the one the Calibration carries: the CRN-filtered residual is\n' ...
         'correlated, so the formal (white-noise) sigma is optimistic, and their ratio is\n' ...
         'cal.uncertainty.inflation.\n']);
