% APP  NonGravForces.app - L8 layer.
% Version 0.10.0  16-Sep-2026
%
%   DynamicsSource       - The acceleration as gjToolbox's integrator interface (native C entry)
%   calibrate            - Accelerometer calibration: bias basis + 3x3 scale + jumps by separable weighted least squares
%   forceBudget          - Per-effect per-axis acceleration budget over an arc; chunked and across models
%   densityInversion     - Thermosphere density from calibrated accelerations; the inversion is a division
