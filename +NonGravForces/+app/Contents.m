% APP  NonGravForces.app - L8 layer.
% Version 0.12.0  17-Sep-2026
%
%   DynamicsSource       - The acceleration as gjToolbox's integrator interface (native C entry)
%   calibrate            - Accelerometer calibration: bias basis + 3x3 scale + jumps by separable weighted least squares
%   forceBudget          - Per-effect per-axis acceleration budget over an arc; chunked and across models
%   densityInversion     - Thermosphere density from calibrated accelerations; the inversion is a division
%   syntheticPair        - Truth generator: a 10 Hz model arc and its CRN-filtered 1 Hz observation with every planted value
%   Calibration          - An accelerometer calibration: scale + time-varying bias + uncertainty + gates; apply/invert/save/load
%   fitCalibration       - The calibration estimator: CRN-filtered bias basis and model columns per axis through lsqvce.fit; VIF gate + bootstrap + EIV bound
