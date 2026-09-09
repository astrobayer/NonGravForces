% SENSOR  NonGravForces.sensor - L6 layer.
% Version 0.8.7  08-Sep-2026
%
%   Sensor               - Abstract measurement model: observe / partials / noiseModel [interface]
%   Electrostatic        - Electrostatic accelerometer: bias scale matrix proof-mass offset coupling FIR filter
%   AtomInterferometer   - Cold-atom interferometer: phase = k_eff a T^2 with dead time, rotation and the 2 pi ambiguity
%   DragFree             - Drag-free residual: sensitivity high-pass cascade with actuator saturation
