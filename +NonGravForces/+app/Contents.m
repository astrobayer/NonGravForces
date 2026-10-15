% APP  NonGravForces.app - L8 layer.
% Version 0.15.0  15-Oct-2026
%
%   DynamicsSource       - The acceleration as gjToolbox's integrator interface (native C entry)
%   calibrate            - Accelerometer calibration of an Arc against a Model or a reference Arc: timeSeriesPreparation.run then lsqvce.fit; attribution
%   forceBudget          - Per-effect per-axis acceleration budget over an arc; chunked and across models
%   densityInversion     - Thermosphere density from calibrated accelerations; the inversion is a division
%   syntheticPair        - Truth generator: a 10 Hz model arc and its CRN-filtered 1 Hz observation with every planted value
%   Calibration          - An accelerometer calibration: scale + time-varying bias + uncertainty + gates; apply/invert/save/load
%   fitCalibration       - The calibration estimator: CRN-filtered bias basis and model columns per axis through lsqvce.fit; VIF gate + bootstrap + EIV bound
%   compareCalibrations  - Rank candidate calibration models per axis by blocked cross-validation; AIC/BIC as IC_unreliable
%   perRevolution        - Bias and diagonal scale per orbital revolution; a gated window is NaN with its reason
%   attributeRevolutions - Per-revolution scales regressed on drivers; verdict density | geometry | sensor | mixed | unresolved
