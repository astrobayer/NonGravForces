% QC  NonGravForces.qc - L7 layer.
% Version 0.13.0  15-Oct-2026
%
%   residuals            - DEPRECATED - forwards to timeSeriesPreparation.qc.residuals (N-01): Whitened residuals measured - modelled with preliminary calibration
%   hampel               - DEPRECATED - forwards to timeSeriesPreparation.qc.hampel (N-01): Robust point-outlier detector (moving median + MAD)
%   detectImpulses       - DEPRECATED - forwards to timeSeriesPreparation.qc.detectImpulses (N-01): Matched filter on the CRN impulse response with THR1B priors
%   detectJumps          - DEPRECATED - forwards to timeSeriesPreparation.qc.detectJumps (N-01): Piecewise-constant/linear change-point detection (PELT + BIC)
%   eventTable           - DEPRECATED - forwards to timeSeriesPreparation.qc.eventTable (N-01): Canonical event table constructor and schema
%   screen               - DEPRECATED - forwards to timeSeriesPreparation.qc.screen (N-01): Orchestrator: residuals -> outliers -> impulses -> jumps -> mask
