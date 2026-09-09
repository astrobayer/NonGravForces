% QC  NonGravForces.qc - L7 layer.
% Version 0.8.4  08-Sep-2026
%
%   residuals            - Whitened residuals measured - modelled with preliminary calibration
%   hampel               - Robust point-outlier detector (moving median + MAD)
%   detectImpulses       - Matched filter on the CRN impulse response with THR1B priors
%   detectJumps          - Piecewise-constant/linear change-point detection (PELT + BIC)
%   eventTable           - Canonical event table constructor and schema
%   screen               - Orchestrator: residuals -> outliers -> impulses -> jumps -> mask
