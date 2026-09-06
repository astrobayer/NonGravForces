% QC  NonGravForces.qc - L7 layer.
% Version 0.1.0  27-Aug-2026
%
%   residuals            - Whitened residuals measured - modelled with preliminary calibration [stub]
%   hampel               - Robust point-outlier detector (moving median + MAD)
%   detectImpulses       - Matched filter on the CRN impulse response with THR1B priors [stub]
%   detectJumps          - Piecewise-constant/linear change-point detection (PELT + BIC) [stub]
%   eventTable           - Canonical event table constructor and schema
%   screen               - Orchestrator: residuals -> outliers -> impulses -> jumps -> mask [stub]
