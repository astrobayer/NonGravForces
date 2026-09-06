% GEOM  NonGravForces.geom - L3 layer.
% Version 0.1.0  27-Aug-2026
%
%   MacroModel           - Panel geometry frozen and validated once
%   graceFO              - GRACE-FO macro model from JPL D-56935
%   boxWing              - Generic box-wing macro model factory [stub]
%   selfShadow           - Per-panel exposed fraction: exact ray-cast or icosphere LUT [stub]
%   thermalState         - Per-panel temperatures: steady state or transient (Wöske) [stub]
%   AttitudeProvider     - BODY->GCRS quaternions at query epochs (SCA1B or nominal) via referenceSystem.interpolateQuaternions [stub]
%   checkAttitudeSense   - Physics test of quaternion sense (nadir and line-of-sight) [stub]
